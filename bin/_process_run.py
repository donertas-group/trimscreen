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


def filter_runs(df, columns):
    # Initialize an empty dictionary to store the results
    results = {}
    N=50
    P=0.6

    # Loop through each column provided
    for column in columns:
        # Group by 'run' and calculate median and standard deviation
        grouped = df.groupby('run')[column].agg(['median', 'std']).reset_index()

        # Sort by highest median and lowest standard deviation
        grouped_sorted = grouped.sort_values(by=['median', 'std'], ascending=[False, True])

        # Calculate the number of top runs to keep (20% or at least 15)
        top_n = max(N, int(np.ceil(len(grouped_sorted) * P)))

        # Get the top runs based on the criteria
        top_runs = grouped_sorted.head(top_n)['run'].values

        # Store the result in the dictionary
        results[column] = set(top_runs)  # Use a set to handle intersections easily

    # Find the intersection of all top runs across columns
    common_runs = set.intersection(*results.values())

    # Filter the dataframe to keep only the rows corresponding to the common 'run's
    filtered_df = df[df['run'].isin(common_runs)]
    
    print("Evaluation criteria: highest median and lowest standard deviation.\n")
    print(f"Keeping the top {P:.0%} or at least the top {N} runs for each of the following columns:") 
    for c in columns:
        print(c)

    print("\nRuns that satisfy filtering criteria for all evaluated columns:")
    for run in filtered_df['run'].unique():
        print(run)

    return filtered_df

