"""
Module for time series forecasting using classical and deep learning models.
"""

import argparse
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
    NBEATSModel,
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

    @staticmethod
    def format_metrics_label(item, name, test, pred):
        """
        Format the label for plots and printouts with metrics.

        Args:
            item (dict): Data info.
            name (str): Model name.
            test (TimeSeries): Test data.
            pred (TimeSeries): Predictions.

        Returns:
            str: Formatted label.
        """
        mape_value = mape(test, pred)
        mse_value = mse(test, pred)
        label = (
            f"{item['probe_id']}_{item['gene_symbol']} - {name}\n"
            f"MAPE: {mape_value:.2f}, MSE: {mse_value:.2f}"
        )
        return label

    def prepare_time_series(self, item):
        """
        Prepare the TimeSeries object and split it into train and test sets.

        Args:
            item (dict): Dictionary containing the DataFrame, probe_id, and gene_symbol.

        Returns:
            tuple: Train and test TimeSeries objects.
        """
        ts = TimeSeries.from_dataframe(
            item['df'],
            time_col='datetime',
            value_cols=['median'],
            freq='4h',
            fill_missing_dates=True,
        )
        filled_ts = du.missing_values.fill_missing_values(ts, method='linear')
        filled_ts = filled_ts.astype(np.float32)
        train, test = filled_ts['median'].split_before(0.8)  # 80% train, 20% test
        return train, test

    def run_models(self, models, plots=True):
        """
        Run forecasting models and plot results.

        Args:
            models (dict): Dictionary of model name to model instance.
            plots (bool): Whether to show plots interactively.
        """
        for item in self.df_list:
            train, test = self.prepare_time_series(item)
            preds = {}
            for i, (name, model) in enumerate(models.items()):
                model.fit(series=train)
                pred = model.predict(n=len(test))
                preds[i] = pred

                # Print metrics with two decimals
                mape_value = mape(test, pred)
                mse_value = mse(test, pred)
                print(
                    f"{name} - Probe: {item['probe_id']}, Gene: {item['gene_symbol']} - "
                    f"MAPE: {mape_value:.2f}, MSE: {mse_value:.2f}"
                )

            if plots:
                fig = plt.figure(figsize=(20, 10))
                for i, (name, model) in enumerate(models.items()):
                    pred = preds[i]
                    label = self.format_metrics_label(item, name, test, pred)
                    ax = fig.add_subplot(2, 3, i + 1)
                    train.plot(ax=ax, label="Train", lw=1)
                    test.plot(ax=ax, label="Test", lw=1)
                    pred.plot(ax=ax, label=f"{name} Prediction", lw=1.5)
                    ax.set_title(label)
                    ax.legend()
                    
                fig.tight_layout()
        if plots:
            plt.show()
                

    def run_classical_models(self, plots=True):
        """
        Run classical forecasting models.
        """
        models = {
            "LinearRegression": LinearRegressionModel(lags=[-1, -2, -3, -6]),
            "Theta": Theta(seasonality_period=6),
            "ExponentialSmoothing": ExponentialSmoothing(seasonal_periods=6),
            "ARIMA": ARIMA(p=1, d=1, q=1, seasonal_order=(1, 1, 1, 6)),
            "RandomForest": RandomForest(lags=[-1, -2, -3, -6], n_estimators=300),
            "NaiveSeasonal": NaiveSeasonal(K=6),
        }
        self.run_models(models, plots)

    def run_deep_learning_models(self, plots=True):
        """
        Run deep learning forecasting models.
        """
        models = {
            "NBEATS": NBEATSModel(
                input_chunk_length=72, output_chunk_length=6, n_epochs=50, random_state=42
            )
        }
        self.run_models(models, plots)


def main():
    """
    Main function to parse arguments and run forecasting.
    """
    parser = argparse.ArgumentParser(
        description="Time series forecasting for gene expression data."
    )
    parser.add_argument(
        "--genes", nargs="+", required=True, help="List of gene symbols to process."
    )
    parser.add_argument(
        "--mode",
        choices=["classical", "deep_learning"],
        default="classical",
        help="Forecasting mode to use.",
    )
    parser.add_argument(
        "--plots", action="store_true", help="Show plots interactively."
    )

    args = parser.parse_args()

    matrix_file = "data/GSE253864_series_matrix.txt.gz"
    platform_file = "data/GPL15331-30377.txt"
    sample_file = "data/sample.csv"

    ts_forecasting = TimeSeriesForecasting(
        matrix_file, platform_file, sample_file, args.genes
    )

    if args.mode == "classical":
        ts_forecasting.run_classical_models(plots=args.plots)
    elif args.mode == "deep_learning":
        ts_forecasting.run_deep_learning_models(plots=args.plots)


if __name__ == "__main__":
    main()