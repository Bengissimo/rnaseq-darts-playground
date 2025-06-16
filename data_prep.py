import pandas as pd
from datetime import datetime, timedelta


class DataPreparer:
    def __init__(self, matrix_file, platform_file):
        """
        Initialize the DataPreparer with file paths.
        """
        self.matrix_file = matrix_file
        self.platform_file = platform_file
        self.df = None
        self.platform_df = None
        self.filtered_df = None
        self.sample_map = None

    def load_expression_and_platform_data(self):
        """
        Load the expression matrix and platform information.
        """
        self.df = pd.read_csv(self.matrix_file, sep='\t', comment='!', header=0)
        #print(self.df.head())  # Debugging line to check the loaded data
        self.platform_df = pd.read_csv(self.platform_file, sep="\t", comment='#', dtype=str)

    def filter_by_gene_symbol(self, gene_symbols):
        """
        Filter the DataFrame by gene symbols and store a mapping of probe_id to gene_symbol.
        """
        if 'GENE_SYMBOL' not in self.platform_df.columns:
            raise ValueError("The platform DataFrame does not contain 'GENE_SYMBOL' column.")
        
        # Filter platform_df by the given gene symbols
        filtered = self.platform_df[self.platform_df['GENE_SYMBOL'].isin(gene_symbols)]
        
        # Store probe IDs and their corresponding gene symbols
        self.probe_ids = filtered['SPOT_ID'].tolist()
        
        # Create a map where the key is probe_id and the value gene_symbol
        self.gene_map = dict(zip(filtered['SPOT_ID'], filtered['GENE_SYMBOL']))
        
        # Filter the expression matrix by the probe IDs
        self.filtered_df = self.df[self.df['ID_REF'].isin(self.probe_ids)]


    def load_sample_metadata(self, sample_file):
        """
        Load the sample map from a file.
        """
        samples = pd.read_csv(sample_file)
        self.sample_map = {
            row['Accession']: row['Title'].rsplit('_', 1)
            for _, row in samples.iterrows()
        }

    def transpose_filtered_data(self):
        """
        Transpose the DataFrame to have samples as rows and genes as columns,
        and add metadata columns.
        """
        if self.filtered_df is None:
            raise ValueError("Data has not been filtered by gene symbols yet.")
        if self.sample_map is None:
            raise ValueError("Sample map has not been loaded yet.")

        # Transpose and add metadata
        df_t = self.filtered_df.set_index('ID_REF').transpose()
        df_t.index.name = 'Accession'
        df_t = df_t.reset_index()

        df_t['group_and_sample'] = df_t['Accession'].map(lambda x: self.sample_map.get(x, [None, None])[0])
        df_t['time'] = df_t['Accession'].map(lambda x: self.sample_map.get(x, [None, None])[1])

        # Sort by metadata
        df_t = df_t.sort_values(by=['group_and_sample', 'time'], ascending=[True, True]).reset_index(drop=True)
        return df_t

    def compute_mean(self):
        """
        Compute mean and standard deviation for a specific gene across groups and times.
        """
        df_t = self.transpose_filtered_data()

        # Extract group and sample information
        df_t[['group', 'sample']] = df_t['group_and_sample'].str.extract(r'(\w+)_S#([A-Z]\d)')
        df_t['sample_pair'] = df_t['sample'].str[0]

        # Group and compute statistics for each probe_id
        values = []
        for probe_id in self.probe_ids:
            if probe_id in df_t.columns:
                values.append(df_t.groupby(['group', 'sample_pair', 'time'])[probe_id]
                              .agg(['mean', 'std']).reset_index().rename(columns={'mean': probe_id}))
        return values

    def add_datetime(self):
        """
        Add datetime information based on group and sample index.
        Returns a list of DataFrames with datetime and gene symbol and probe_id.
        """

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

        # Compute sample index and base date
        for df in self.compute_mean():
            df['sample_index'] = df['sample_pair'].apply(lambda x: ord(x) - ord('A'))
            df['base_date'] = df['group'].map(base_dates)

            # Compute datetime
            def compute_datetime(row):
                if row['base_date'] is not None:
                    date_with_offset = row['base_date'] + timedelta(days=row['sample_index'])
                    time_part = datetime.strptime(row['time'], "%H:%M").time()
                    return datetime.combine(date_with_offset, time_part)
                return None

            df['datetime'] = df.apply(compute_datetime, axis=1)

            probe_id = None
            for column_name in df.columns:
            # if the column name exists in the gene_map as one of the probe_ids 
                if column_name in self.gene_map:
                    probe_id = column_name
                else:
                    continue
        
            # Add the DataFrame and gene_symbol to the dictionary
            df_list.append({'df': df, 'gene_symbol': self.gene_map[probe_id], 'probe_id': probe_id})

        return df_list
