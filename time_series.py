import data_prep as dp

# Example usage
if __name__ == "__main__":
    matrix_file = "data/GSE253864_series_matrix.txt.gz"
    platform_file = "data/GPL15331-30377.txt"
    sample_file = "data/sample.csv"

    d = dp.DataPreparer(matrix_file, platform_file)
    d.load_expression_and_platform_data()

   # TODO give gene symbol as argument
    d.filter_by_gene_symbol(['PER1'])
    d.load_sample_metadata(sample_file)

    df = d.add_datetime()

    print(df.head(30))