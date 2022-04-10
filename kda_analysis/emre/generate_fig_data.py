"""
Module for generating the figure data sets for the 8-state model.
"""

import pandas as pd

from fig_7a import get_7a_dataframe
from fig_7b import get_7b_dataframe
from fig_7d import get_7d_dataframe
from fig_9 import get_9_dataframe


def get_Hussey_datasets():
    # get the dataframes containing the data sets for each figure
    df_7a = get_7a_dataframe()
    df_7b = get_7b_dataframe()
    df_7d = get_7d_dataframe()
    df_9 = get_9_dataframe()
    # combine the dataframes into a single dataframe
    df = pd.concat([df_7a, df_7b, df_7d, df_9], ignore_index=True)
    # return the transpose of the data frame so we can
    # iterate over the rows instead of the columns
    return df


if __name__ == "__main__":

    df = get_Hussey_datasets()
    print("Successfully generated all figure data.")
    fig_data_path = "data/all_figure_data.csv"
    df.to_csv(fig_data_path, index=False)
    print(f"Saving data at location: {fig_data_path}")
