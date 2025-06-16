import data_prep as dp
import numpy as np
from darts import TimeSeries
import matplotlib.pyplot as plt
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


import pandas as pd


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

    df_list = d.add_datetime()

    for df in df_list:
        series = TimeSeries.from_dataframe(df, time_col='datetime', value_cols='mean', freq='4h', fill_missing_dates=True)
        filled = du.missing_values.fill_missing_values(series, method='linear')
        filled = filled.astype(np.float32)


        train, test = filled['mean'].split_before(0.8)  # 80% train, 20% test

        # classical models
        model0 = LinearRegressionModel(lags=24)
        model1 = Theta(seasonality_period=24)
        model2 = ExponentialSmoothing(seasonal_periods=24)
        model3 = ARIMA(p=12, d=1, q=24)
        model4 = RandomForest(lags=24, n_estimators=300)
        model5 = NaiveSeasonal(K=24)

        models = {
            "LinearRegression": model0,
            "Theta": model1,
            "ExponentialSmoothing": model2,
            "ARIMA": model3,
            "RandomForest": model4,
            "NaiveSeasonal": model5
        }

        fig = plt.figure(figsize=(20, 10))
        for i, (name, model) in enumerate(models.items()):
            print(name)
            model.fit(series=train)
            pred = model.predict(n=len(test))

            # Calculate metrics
            mape_value = round(mape(test, pred), 2)
            mse_value = round(mse(test, pred), 2)

            fig.add_subplot(2, 3, i+1)
            train.plot(label="Train", lw=1)
            test.plot(label="Test", lw = 1)
            pred.plot(label=f"{name} Prediction", lw=1.5)
            plt.title(f"Method: {name}, MAPE: {mape_value}, MSE: {mse_value}"); plt.legend()

        fig.tight_layout()
    plt.show()
