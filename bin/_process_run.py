import os
import pandas as pd
import glob
import numpy as np
from skbio import diversity
 

def add_group_statistics(df, column_to_group, columns):
    """
    Add mean and standard error columns for a given column grouped by another column.

    Parameters:
    - df: pandas DataFrame, the input data.
    - column_to_group: str, the column to group by (e.g., "run").
    - columns: a list of strings, the columns to calculate statistics for (e.g., "Species_standardised").

    Returns:
    - DataFrame with added new columns: `<column>_run_mean` and `<column>_run_std`.
    """
    for column in columns:
        # Group by the specified column
        grouped = df.groupby(column_to_group)[column]

        # Calculate mean and standard error
        means = grouped.mean().rename(f"{column}_run_mean")
        std_errors = grouped.sem().rename(f"{column}_run_std")

        # Merge the results back into the original DataFrame
        df = df.merge(means, on=column_to_group)
        df = df.merge(std_errors, on=column_to_group)

    return df


