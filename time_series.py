"""
Module for time series forecasting using classical and deep learning models.
"""
import numpy as np
import matplotlib.pyplot as plt
from darts import TimeSeries
from darts.models import (
    NaiveSeasonal,
    LinearRegressionModel,
    ExponentialSmoothing,
    Theta,
    ARIMA,
    RandomForest,
    NBEATSModel
)
from darts.metrics import mape, mse
from darts import utils as du
import data_prep as dp


class TimeSeriesForecasting:
    """
    Class for time series forecasting using classical and deep learning models.
    """

    def __init__(self, matrix_file, platform_file, sample_file, gene_symbols):
        """
        Initialize the forecasting class with file paths and gene symbols.

        Args:
            matrix_file (str): Path to the matrix file.
            platform_file (str): Path to the platform file.
            sample_file (str): Path to the sample file.
            gene_symbols (list): List of gene symbols to process.
        """
        self.data_preparer = dp.DataPreparer(matrix_file, platform_file, sample_file)
        self.df_list = self.data_preparer.prepare_data(gene_symbols)

    def run_classical_models(self):
        """
        Run classical forecasting models and plot results.
        """
        for item in self.df_list:
            # Create a TimeSeries object
            ts = TimeSeries.from_dataframe(
                item['df'], time_col='datetime', value_cols=['median'], freq='4h', fill_missing_dates=True
            )
            filled_ts = du.missing_values.fill_missing_values(ts, method='linear')
            filled_ts = filled_ts.astype(np.float32)

            train, test = filled_ts['median'].split_before(0.8)  # 80% train, 20% test

            # Classical models
            models = {
                "LinearRegression": LinearRegressionModel(lags=24),
                "Theta": Theta(seasonality_period=24),
                "ExponentialSmoothing": ExponentialSmoothing(seasonal_periods=24),
                "ARIMA": ARIMA(p=12, d=1, q=24),
                "RandomForest": RandomForest(lags=24, n_estimators=300),
                "NaiveSeasonal": NaiveSeasonal(K=24)
            }

            fig = plt.figure(figsize=(20, 10))
            for i, (name, model) in enumerate(models.items()):
                print(f"Processing {name} for {item['probe_id']} {item['gene_symbol']} ...")
                model.fit(series=train)
                pred = model.predict(n=len(test))

                # Calculate metrics
                mape_value = round(mape(test, pred), 2)
                mse_value = round(mse(test, pred), 2)

                # Plot results
                ax = fig.add_subplot(2, 3, i + 1)
                train.plot(ax=ax, label="Train", lw=1)
                test.plot(ax=ax, label="Test", lw=1)
                pred.plot(ax=ax, label=f"{name} Prediction", lw=1.5)
                ax.set_title(f"{item['probe_id']}_{item['gene_symbol']} - {name}\nMAPE: {mape_value}, MSE: {mse_value}")
                ax.legend()

            fig.tight_layout()
        plt.show()


# Example usage
if __name__ == "__main__":
    MATRIX_FILE = "data/GSE253864_series_matrix.txt.gz"
    PLATFORM_FILE = "data/GPL15331-30377.txt"
    SAMPLE_FILE = "data/sample.csv"
    GENE_SYMBOLS = ['PER1', 'NR1D1']

    ts_forecasting = TimeSeriesForecasting(MATRIX_FILE, PLATFORM_FILE, SAMPLE_FILE, GENE_SYMBOLS)
    ts_forecasting.run_classical_models()
    # ts_forecasting.run_deep_learning_models()