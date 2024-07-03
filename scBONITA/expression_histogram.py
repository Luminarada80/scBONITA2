import pandas as pd
import matplotlib.pyplot as plt

def plot_single_histogram(data):
    array = data.values.flatten()

    filtered_array = array[array != 0]

    # Set up the figure
    plt.figure(figsize=(5, 6))

    plt.hist(filtered_array, bins=250, log=False)
    plt.title('Histogram of scRNAseq Gene Expression')
    plt.xlabel('log2 Normalized Counts')
    plt.ylabel('log10 Frequency')
    plt.xlim((-0.5, 2.5))
    plt.tight_layout()
    plt.show()

def plot_individual_gene_histograms(data):
    # Extract all expression values from the dataset
    expression_values = data.iloc[:, 1:51].values

    # Set up the figure
    plt.figure(figsize=(10, 6))

    # Plot the distribution of each gene with low opacity
    for i in range(expression_values.shape[1]):
        gene_expression_values = expression_values[:, i]
        non_zero_gene_expression_values = gene_expression_values[gene_expression_values > 0]
        plt.hist(non_zero_gene_expression_values, bins=100, alpha=0.1, edgecolor='none', log=False)

    plt.title('Distribution of Expression Values')
    plt.xlabel('Expression Value')
    plt.ylabel('Frequency')
    plt.xlim((-0.5, 5))
    plt.ylim((0, 5000))
    plt.grid(True)
    plt.show()

def plot_total_gene_expression_histogram(data):
    array = data.values.flatten()

    filtered_array = array[array != 0]

    # Create a figure and a set of subplots
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))  # 2 rows, 1 column of subplots

    # Plotting the histogram of the filtered array (non-zero values) with log scale
    axs[0].hist(filtered_array, bins=250, log=False)
    axs[0].set_title('Frequency of Expression Values')
    axs[0].set_xlabel('Expression Value')
    axs[0].set_ylabel('Frequency')
    axs[0].set_xlim((-0.5, 5))


    #Plotting the histogram of all values (including zeros) with log scale
    axs[1].hist(filtered_array, bins=250, log=False)
    axs[1].set_title('Zoomed in View of Expression Frequencies')
    axs[1].set_xlabel('Expression Value')
    axs[1].set_ylabel('Frequency')
    axs[1].set_xlim((-0.5, 5))
    axs[1].set_ylim((0, 0.30e6))

    #Show the plots
    plt.tight_layout()  # Adjusts subplot params so that subplots fit into the figure area.
    plt.show()

if __name__ == '__main__':
    datafile = "../../george_data/hiv_dataset/HIV_dataset_normalized_integrated_counts.csv"
    data = pd.read_csv(datafile, index_col=0)
    plot_single_histogram(data)


