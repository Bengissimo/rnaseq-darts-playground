import pandas as pd
from datetime import datetime, timedelta


class DataPreparer:
    def __init__(self, matrix_file, platform_file, sample_file):
        """
        Initialize the DataPreparer with file paths.
        """
        self.matrix_file = matrix_file
        self.platform_file = platform_file
        self.sample_file = sample_file

    def prepare_data(self, gene_symbols):
        """
        Prepare data by loading, filtering, transposing, and adding datetime information.
        Returns a list of DataFrames with datetime, gene_symbol and probe_id.
        """
        # Load expression matrix, platform data and sample metadata
        df = pd.read_csv(self.matrix_file, sep='\t', comment='!', header=0)
        platform_df = pd.read_csv(self.platform_file, sep="\t", comment='#', dtype=str)
        
        if 'GENE_SYMBOL' not in platform_df.columns:
            raise ValueError("The platform DataFrame does not contain 'GENE_SYMBOL' column.")
        filtered = platform_df[platform_df['GENE_SYMBOL'].isin(gene_symbols)]
        probe_ids = filtered['SPOT_ID'].tolist()
        gene_map = dict(zip(filtered['SPOT_ID'], filtered['GENE_SYMBOL']))
        filtered_df = df[df['ID_REF'].isin(probe_ids)]

        samples = pd.read_csv(self.sample_file)
        sample_map = {
            row['Accession']: row['Title'].rsplit('_', 1)
            for _, row in samples.iterrows()
        }

        # Transpose filtered data and add metadata
        df_t = filtered_df.set_index('ID_REF').transpose()
        df_t.index.name = 'Accession'
        df_t = df_t.reset_index()
        df_t['group_and_sample'] = df_t['Accession'].map(lambda x: sample_map.get(x, [None, None])[0])
        df_t['time'] = df_t['Accession'].map(lambda x: sample_map.get(x, [None, None])[1])
        
        # Sort by group and time
        df_t = df_t.sort_values(by=['group_and_sample', 'time'], ascending=[True, True]).reset_index(drop=True)

        # Extract group and sample information
        df_t[['group', 'sample']] = df_t['group_and_sample'].str.extract(r'(\w+)_S#([A-Z]\d)')
        df_t['sample_pair'] = df_t['sample'].str[0]

        # Define base dates for each group
        base_dates = {
            'BDC1': datetime(2024, 1, 1),
            'BDC2': datetime(2024, 1, 11),
            'HDT1': datetime(2024, 1, 21),
            'HDT2': datetime(2024, 1, 31),
            'HDT3': datetime(2024, 2, 10),
            'R': datetime(2024, 2, 20),
        }

        # Create a list to store DataFrames and its corresponding gene symbols
        df_list = []
        for probe_id in probe_ids:
            if probe_id in df_t.columns:
                stats_df = df_t.groupby(['group', 'sample_pair', 'time'])[probe_id].agg(['median', 'std']).reset_index()
                #TODO: handle outliers
                stats_df['sample_index'] = stats_df['sample_pair'].apply(lambda x: ord(x) - ord('A'))
                stats_df['base_date'] = stats_df['group'].map(base_dates)

                # compute datetime
                def compute_datetime(row):
                    if row['base_date'] is not None:
                        date_with_offset = row['base_date'] + timedelta(days=row['sample_index'])
                        time_part = datetime.strptime(row['time'], "%H:%M").time()
                        return datetime.combine(date_with_offset, time_part)
                    return None
                stats_df['datetime'] = stats_df.apply(compute_datetime, axis=1)

                # Add the DataFrame and gene_symbol to the dictionary
                gene_symbol = gene_map[probe_id]
                df_list.append({'df': stats_df, 'gene_symbol': gene_symbol, 'probe_id': probe_id})

        return df_list
