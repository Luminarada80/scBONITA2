import pandas as pd
from io import StringIO
from fastdtw import fastdtw
from scipy.spatial.distance import euclidean
import os
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import statistics
from alive_progress import alive_bar

from file_paths import file_paths

def read_data_from_directory(directory):
    print(f'Reading data from the files')
    dataframes = {}
    for filename in os.listdir(directory):
        print(f'\t{filename}')
        if filename.endswith("_trajectory.csv"):
            filepath = os.path.join(directory, filename)
            df = pd.read_csv(filepath, header=None)
            df.columns = ['Gene'] + [f'Time{i}' for i in range(1, df.shape[1])]
            df.set_index('Gene', inplace=True)
            dataframes[filename] = df
    return dataframes

def extract_time_series(dataframes):
    print('Extracting time series data')
    return {filename: {gene: df.loc[gene].values for gene in df.index} for filename, df in dataframes.items()}

def compute_dtw_distance_pair(file1, file2, time_series_data):
    distances = {}
    for gene in time_series_data[file1].keys():
        if gene in time_series_data[file2]:
            ts1 = time_series_data[file1][gene]
            ts2 = time_series_data[file2][gene]
            distance, _ = fastdtw(ts1, ts2, dist=euclidean)
            distances[gene] = distance
    total_distance = sum(distances.values()) if distances else float('inf')
    return (file1, file2, total_distance)

def compute_dtw_distances(time_series_data, output_directory):
    dtw_distances = {}
    file_names = list(time_series_data.keys())
    total_combinations = len(file_names) * (len(file_names) - 1) // 2
    
    tasks = []
    with ProcessPoolExecutor() as executor, alive_bar(total_combinations) as bar:
        for i in range(len(file_names)):
            for j in range(i + 1, len(file_names)):
                file1 = file_names[i]
                file2 = file_names[j]
                tasks.append(executor.submit(compute_dtw_distance_pair, file1, file2, time_series_data))
                
        for future in as_completed(tasks):
            file1, file2, total_distance = future.result()
            dtw_distances[(file1, file2)] = total_distance
            bar()

    with open(f'{output_directory}/distances.csv', 'w') as outfile:
        for (file1, file2), total_distance in dtw_distances.items():
            file1_cell = file1.split('_trajectory')[0]
            file2_cell = file2.split('_trajectory')[0]
            outfile.write(f'{file1_cell},{file2_cell},{total_distance}\n')

    return dtw_distances

def create_distance_matrix(dtw_distances, file_names):
    # Create a distance matrix
    distance_matrix = np.zeros((len(file_names), len(file_names)))
    
    for (file1, file2), total_distance in dtw_distances.items():
        i = file_names.index(file1)
        j = file_names.index(file2)
        distance_matrix[i, j] = total_distance
        distance_matrix[j, i] = total_distance
    
    return distance_matrix

def find_similar_files(dtw_distances):

    # Extract unique cell names
    cells = set()
    for (cell1, cell2), _ in dtw_distances.items():
        cells.add(cell1.split('_trajectory')[0])
        cells.add(cell2.split('_trajectory')[0])
    cells = sorted(cells)

    # Create a distance matrix
    distance_matrix = pd.DataFrame(np.inf, index=cells, columns=cells)
    for (cell1, cell2), distance in dtw_distances.items():
        distance_matrix.at[cell1.split('_trajectory')[0], cell2.split('_trajectory')[0]] = distance
        distance_matrix.at[cell2.split('_trajectory')[0], cell1.split('_trajectory')[0]] = distance

    # Replace infinity values with a large number
    distance_matrix.replace(np.inf, 1e6, inplace=True)

    # Perform hierarchical clustering
    distance_array = distance_matrix.values[np.triu_indices_from(distance_matrix, k=1)]
    Z = linkage(distance_array, method='ward')


    # Plot the dendrogram
    plt.figure(figsize=(12, 10))
    dendrogram(Z, labels=distance_matrix.index, orientation='top')
    plt.title('Hierarchical Clustering Dendrogram')
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.xlabel('Distance')
    plt.ylabel('Cells')
    plt.show()

    # Set a threshold and get the clusters
    num_clusters = int(input('How many clusters?: '))
    clusters = fcluster(Z, num_clusters, criterion='maxclust')

    # Organize cells by clusters
    cluster_dict = {}
    for cell, cluster_id in zip(cells, clusters):
        if cluster_id not in cluster_dict:
            cluster_dict[cluster_id] = []
        cluster_dict[cluster_id].append(cell)

    return cluster_dict


def summarize_clusters(directory, cell_names, cluster):
    gene_expr_dict = {}

    trajectory_files = []
    for filename in os.listdir(directory):
        if filename.endswith("_trajectory.csv"):
            trajectory_files.append(os.path.join(directory, filename))

    files_to_open = []
    for cell in cell_names:
        for file_name in trajectory_files:
            if cell in file_name:
                files_to_open.append(file_name)

    for file_path in files_to_open:
        with open(file_path, 'r') as sim_file:
            for line in sim_file:
                line = line.strip().split(',')
                gene_name = line[0]
                gene_expression = [int(i) for i in line[1:]]
                
                if gene_name not in gene_expr_dict:
                    gene_expr_dict[gene_name] = []
                
                gene_expr_dict[gene_name].append(gene_expression)

    gene_avg_expr = {}

    for gene, simulation_results in gene_expr_dict.items():

        if gene not in gene_avg_expr:
            gene_avg_expr[gene] = []

        transposed_data = list(map(list, zip(*simulation_results)))

        for i in transposed_data:
            gene_avg_expr[gene].append(statistics.mean(i))
    
    # Convert the dictionary to a DataFrame
    df = pd.DataFrame(gene_avg_expr)

    # Transpose the DataFrame
    df = df.transpose()

    # Create the heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(df, cmap='Greys')
    plt.title(f'Average Gene Expression Heatmap for Cluster {cluster}')
    plt.xlabel('Simulation Time Steps')
    plt.ylabel('Gene')
    plt.show()

def plot_heatmap(distance_matrix, file_names):
    plt.figure(figsize=(8, 8))
    file_names = [i.split('_trajectory')[0] for i in file_names]
    sns.heatmap(distance_matrix, xticklabels=file_names, yticklabels=file_names, cmap='Greys', annot=False)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.title("DTW Distance Heatmap")
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    parser = ArgumentParser()

    parser.add_argument(
        '--dataset_name',
        type=str,
        required=True,
        help='Number of genes to generate'
    )

    parser.add_argument(
        '--network_name',
        type=str,
        required=True,
        help='Number of cells to generate'
    )

    results = parser.parse_args()

    dataset_name = getattr(results, 'dataset_name')
    network_name = getattr(results, 'network_name')

    # Directory containing the files
    directory = f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/text_files'

    dataframes = read_data_from_directory(directory)
    time_series_data = extract_time_series(dataframes)
    dtw_distances = compute_dtw_distances(time_series_data, directory)
    
    file_names = list(dataframes.keys())
    distance_matrix = create_distance_matrix(dtw_distances, file_names)

    cluster_dict = find_similar_files(dtw_distances)

    for cluster, cell_list in cluster_dict.items():
        print(f'Summarizing cluster {cluster}')
        summarize_clusters(directory, cell_list, cluster)
    
    plot_heatmap(distance_matrix, file_names)

    

