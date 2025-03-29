import seaborn as sns
import matplotlib.pyplot as plt
import os
import pandas as pd
import math
import numpy as np

def create_scatter_plot(df, X, Y, Z, comb_to_color, outdir):
    # Select only the necessary columns
    relevant_columns = [X, Y, Z, 'trunclenf', 'trunclenr']
    df_relevant = df[relevant_columns]

    # Extract unique combinations of 'trunclenf' and 'trunclenr'
    unique_combs = df_relevant[['trunclenf', 'trunclenr']].drop_duplicates()

    # Map colors for the entire dataset
    colors = df_relevant.apply(lambda row: comb_to_color[(row['trunclenf'], row['trunclenr'])], axis=1)

    # Plot all points at once
    plt.figure(figsize=(20, 6))

    # Normalize the sizes to a range
    min_size, max_size = 10, 300  # Adjust these values as needed
    normalized_sizes = (df_relevant[Z] - df_relevant[Z].min()) / (df_relevant[Z].max() - df_relevant[Z].min())
    sizes = normalized_sizes * (max_size - min_size) + min_size

    plt.scatter(df_relevant[X], df_relevant[Y], c=colors, s=sizes, alpha=0.3)

    # Group by 'trunclenf' and 'trunclenr' and find the top point in each group
    df_relevant['sum_XY'] = df_relevant[X] + df_relevant[Y]
    top_points_by_group = (
        df_relevant.groupby(['trunclenf', 'trunclenr'])
        .apply(lambda group: group.nlargest(1, 'sum_XY'))  # Top point in each group
        .reset_index(drop=True)
    )

    # Select the top 3 points across all groups
    top_3_points = top_points_by_group.nlargest(3, 'sum_XY')

    # Print the top 3 points on the upper-right corner
    annotation_text = "\n".join([
        f"Point {i+1}: X={row[X]:.2f}, Y={row[Y]:.2f}, Z={row[Z]:.2f}, trunclenf={row['trunclenf']}, trunclenr={row['trunclenr']}"
        for i, row in top_3_points.iterrows()
    ])
    plt.text(
        0.99, 0.99, annotation_text,
        transform=plt.gca().transAxes,
        fontsize=10,
        va='bottom', ha='right',
        bbox=dict(boxstyle="round,pad=0.3", edgecolor='gray', facecolor='white', alpha=0.8)
    )

    # Create a legend for unique combinations
    handles = []
    labels = []
    for _, row in unique_combs.iterrows():
        comb = (row['trunclenf'], row['trunclenr'])
        handles.append(plt.Line2D([], [], marker='o', color=comb_to_color[comb], linestyle='', markersize=10))
        labels.append(f"({comb[0]}, {comb[1]})")

    ncol = math.ceil(len(labels) / 20)
    color_legend = plt.legend(
        handles,
        labels,
        title="(trunclenf, trunclenr)",
        loc='upper left',
        bbox_to_anchor=(1, 1),
        ncol=ncol
    )
    plt.gca().add_artist(color_legend)  # Add the color legend explicitly

    # Create a size legend
    size_legend_values = np.linspace(df_relevant[Z].min(), df_relevant[Z].max(), num=5)
    size_legend_sizes = (size_legend_values - df_relevant[Z].min()) / (df_relevant[Z].max() - df_relevant[Z].min())
    size_legend_sizes = size_legend_sizes * (max_size - min_size) + min_size

    size_handles = [
        plt.scatter([], [], c='gray', alpha=0.3, s=size, label=f"{value:.4f}")
        for value, size in zip(size_legend_values, size_legend_sizes)
    ]
    plt.legend(
        handles=size_handles,
        title=Z,
        loc='upper left',
        bbox_to_anchor=(1, 0.45)  # Adjust vertical position for below the color legend
    )

    # Axis labels and title
    plt.xlabel(X)
    plt.ylabel(Y)
    plt.title(f'Scatter Plot of {X} vs {Y} with Colored Trunclen Combinations')

    # Adjust subplots to ensure enough space for the legend
    plt.subplots_adjust(right=0.5)

    # Save the plot as a JPG
    jpg_path = os.path.join(outdir, f'Scatter_Plot_{X}_vs_{Y}.jpg')
    plt.savefig(jpg_path, dpi=400)
    plt.close()  # Close the figure to avoid display
    print(f"Plot saved as {jpg_path}")

def create_scatter_plot_working(df, X, Y, Z, comb_to_color, outdir):
    # Select only the necessary columns
    relevant_columns = [X, Y, Z, 'trunclenf', 'trunclenr']
    df_relevant = df[relevant_columns]

    # Extract unique combinations of 'trunclenf' and 'trunclenr'
    unique_combs = df_relevant[['trunclenf', 'trunclenr']].drop_duplicates()

    # Map colors for the entire dataset
    colors = df_relevant.apply(lambda row: comb_to_color[(row['trunclenf'], row['trunclenr'])], axis=1)

    # Plot all points at once
    plt.figure(figsize=(20, 6))

    # Normalize the sizes to a range
    min_size, max_size = 10, 300  # Adjust these values as needed
    normalized_sizes = (df_relevant[Z] - df_relevant[Z].min()) / (df_relevant[Z].max() - df_relevant[Z].min())
    sizes = normalized_sizes * (max_size - min_size) + min_size

    plt.scatter(df_relevant[X], df_relevant[Y], c=colors, s=sizes, alpha=0.3)

    # Create a legend for unique combinations
    handles = []
    labels = []
    for _, row in unique_combs.iterrows():
        comb = (row['trunclenf'], row['trunclenr'])
        handles.append(plt.Line2D([], [], marker='o', color=comb_to_color[comb], linestyle='', markersize=10))
        labels.append(f"({comb[0]}, {comb[1]})")

    ncol = math.ceil(len(labels) / 20)
    color_legend = plt.legend(
        handles, 
        labels, 
        title="(trunclenf, trunclenr)", 
        loc='upper left', 
        bbox_to_anchor=(1, 1), 
        ncol=ncol
    )
    plt.gca().add_artist(color_legend)  # Add the color legend explicitly

    # Create a size legend
    size_legend_values = np.linspace(df_relevant[Z].min(), df_relevant[Z].max(), num=5)
    size_legend_sizes = (size_legend_values - df_relevant[Z].min()) / (df_relevant[Z].max() - df_relevant[Z].min())
    size_legend_sizes = size_legend_sizes * (max_size - min_size) + min_size

    size_handles = [
        plt.scatter([], [], c='gray', alpha=0.3, s=size, label=f"{value:.4f}") 
        for value, size in zip(size_legend_values, size_legend_sizes)
    ]
    plt.legend(
        handles=size_handles, 
        title=Z, 
        loc='upper left', 
        bbox_to_anchor=(1, 0.45)  # Adjust vertical position for below the color legend
    )

    # Axis labels and title
    plt.xlabel(X)
    plt.ylabel(Y)
    plt.title(f'Scatter Plot of {X} vs {Y} with Colored Trunclen Combinations')

    # Adjust subplots to ensure enough space for the legend
    plt.subplots_adjust(right=0.5)

    # Save the plot as a JPG
    jpg_path = os.path.join(outdir, f'Scatter_Plot_{X}_vs_{Y}.jpg')
    plt.savefig(jpg_path, dpi=400)
    plt.close()  # Close the figure to avoid display
    print(f"Plot saved as {jpg_path}")


def create_violin_plot(df, Y, outdir):
    # Create a new column combining 'trunclenf' and 'trunclenr' as a unique identifier
    df['trunc_combination'] = df['trunclenf'].astype(str) + "_" + df['trunclenr'].astype(str)

    # Plotting
    plt.figure(figsize=(20, 6))

    sns.violinplot(x='trunc_combination', y=Y, data=df,linewidth=0.5)
    plt.xticks(rotation=45)
    plt.xlabel("truncLenf_truncLenf combination")
    plt.ylabel(f"{Y}")
    plt.title(f"Violin Plot of {Y} for Each Trunclen Combination")

    # Adjust subplots to ensure enough space for the legend
    plt.subplots_adjust(right=0.9)

    # Save the plot as a PDF
    jpg_path = os.path.join(outdir, f'Violin_Plot_runs_vs_{Y}.jpg')
    plt.savefig(jpg_path)
    plt.close()  # Close the figure to avoid display
    print(f"Plot saved as {jpg_path}")
