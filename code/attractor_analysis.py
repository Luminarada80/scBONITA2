from importance_scores import CalculateImportanceScore
import logging
import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
from scipy.sparse import spmatrix, issparse
import matplotlib.pyplot as plt
import pickle
import os
import glob
from simulate_attractor import *
import argparse
import pickle
import random

from fastdtw import fastdtw
from scipy.spatial.distance import euclidean
import os
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from alive_progress import alive_bar
import statistics

from user_input_prompts import attractor_analysis_arguments
from file_paths import file_paths

def run_attractor_analysis(network, cells):
    """
    Runs the attractor analysis for the dataset and all of the networks.
    """
    # Set up the attractor analysis with the dataset and the parsed networks

    logging.info(f'\nNetwork: {network.name}')
    # Generate attractor
    logging.info(f'\tGenerating attractors...')

        # Specifies the path to the correct network pickle file
    network_pickle_file = f'{file_paths["pickle_files"]}/{dataset_name}_pickle_files/network_pickle_files/{dataset_name}_{network_name}.network.pickle'

    network_attractors, simulated_dataset = generate_attractors(network.nodes, network.dataset)
    attractors_start = [attractor[:, 0] for _, attractor in network_attractors.items()]

    # Convert the sparse dataset to a dense one
    if isinstance(simulated_dataset, spmatrix):
        dense_dataset = simulated_dataset.todense()
    else:
        dense_dataset = simulated_dataset
        
    # Generate Hamming distance matrix between cells and attractors
    logging.info(f'\tCalculating hamming distance between cells and attractors')
    logging.info(f'\tTransposed dataset shape: {dense_dataset.shape}')
    transposed_dataset = transpose_dataset(dense_dataset)
    num_nodes = len(network.nodes)
    num_cells = range(transposed_dataset.shape[0])
    num_attractors = range(len(attractors_start))

    logging.info(f'\t\tNodes: {num_nodes}')
    logging.info(f'\t\tCells: {transposed_dataset.shape[0]}')
    logging.info(f'\t\tAttractors: {len(attractors_start)}')
    cell_attractor_hamming_df = calculate_hamming_distance(num_attractors, num_cells, transposed_dataset, attractors_start)

    # Filter the attractors to only keep those that most closely match the expression of at least one cell
    filtered_attractor_indices = filter_attractors(cell_attractor_hamming_df)
    filtered_attractors = [attractors_start[i] for i in filtered_attractor_indices]
    num_filtered_attractors = range(len(filtered_attractors))

    # Create a distance matrix between each of the attractors
    logging.info(f'\tGenerating attractor distance matrix...')
    attractor_distance_matrix = calculate_hamming_distance(num_filtered_attractors, num_filtered_attractors, filtered_attractors, filtered_attractors)

    # Cluster the attractors using hierarchical agglomerative clustering
    logging.info(f'\tClustering the attractors...')
    clusters, cluster_fig = hierarchical_clustering(attractor_distance_matrix, len(network.nodes), show_plot=False)

    clustered_attractors = {}
    for i, cluster_num in enumerate(clusters):
        if cluster_num not in clustered_attractors:
            clustered_attractors[cluster_num] = {}
        clustered_attractors[cluster_num][filtered_attractor_indices[i]] = filtered_attractors[i]
    
    # Find the Hamming Distance between the cells and the attractors within the clusters
    logging.info(f'\tCalculating Hamming distance between cells and clustered attractors')

    # For each of the cells, find the hamming distance to the best cell
    logging.info(f'Number of cells in the full dataset: {network.dataset.shape[1]}')
    if issparse(network.dataset):
        full_dataset = network.dataset.toarray().T
    else:
        full_dataset = network.dataset.T

    cell_map = {}
    for cell_num, cell in enumerate(full_dataset):
        cell_map[cell_num] = {}

    for cluster, attractors in clustered_attractors.items():
        logging.debug(f'\nCluster {cluster}')

        num_filtered_attractors = range(len(attractors.values()))

        mapped_cluster_attractors = calculate_hamming_distance(num_filtered_attractors, num_cells, attractors.values(), transposed_dataset)

        # transposed_mapped_cluster_attractors = transpose_dataset(mapped_cluster_attractors)
        logging.debug(f'Mapped cluster attractors:\n {mapped_cluster_attractors}')

        # Map the cells to the attractors within the clusters
        filtered_attractor_indices = list(clustered_attractors[cluster].keys())
        logging.debug(f'Filtered attractor indices: {filtered_attractor_indices}')
        mapped_cluster_attractors.index = filtered_attractor_indices # type: ignore
        logging.debug(f'Cluster {cluster} mapped attractors')
        logging.debug(mapped_cluster_attractors)

        # Find the total Hamming distance between each attractor and all of the cells
        total_hamming_distance = mapped_cluster_attractors.sum(axis=1)
        min_sum_index = total_hamming_distance.idxmin()

        representative_attractor = attractors[min_sum_index]

        # Find the hamming distance for each cell for each cluster
        for cell_num, cell in enumerate(full_dataset):
            hamming_distance = 0
            for gene_num, gene in enumerate(cell):
                if cell[gene_num] != representative_attractor[gene_num]:
                   hamming_distance += 1
            
            cell_map[cell_num][cluster] = hamming_distance
    

        logging.debug(f'Representative attractor for cluster {cluster} (index {min_sum_index}):')
        logging.debug(f'\n{attractors[min_sum_index]}')

        network.representative_attractors[cluster] = representative_attractor
    
    # Calculate which cluster each cell should map to
    min_hamming_cluster = {}
    for cell_num, clusters in cell_map.items():

        # Find the cluster with the minimum Hamming distance
        min_cluster = min(clusters, key=clusters.get)
        min_hamming_cluster[cell_num] = min_cluster
        # cell_object = cells[cell_num]

    network.cell_map = min_hamming_cluster

    for cell_object in cells:
        # Add the best attractor for each cell to that cell object
        cell_object.attractor_dict[network.name] = min_hamming_cluster[cell_object.index]

    return cluster_fig

def generate_attractors(nodes, dataset):
    """
    Uses the vectorized simulation function from the importance score calculations to find
    the attractors for the network.
    """

    calculate_importance_score = CalculateImportanceScore(nodes, dataset)

    simulated_dataset = calculate_importance_score.dataset        
    network_attractors = calculate_importance_score.vectorized_run_simulation()

    return network_attractors, simulated_dataset


def extract_node_rows(dataset, nodes):
    """
    Create a dense dataset only containing the nodes that are present in the network from a 
    sparse dataset
    """

    # Convert the sparse dataset to a dense one
    if isinstance(dataset, spmatrix):
        dense_dataset = dataset.todense()
    else:
        dense_dataset = dataset

    # Find the row indices in the dataset for the nodes in the network
    node_indices = []
    for node in nodes:
        if node.dataset_index not in node_indices:
            node_indices.append(node.dataset_index)
        
    return dense_dataset[node_indices]

def transpose_dataset(dataset):
    """
    Converts a dataset to a numpy array and transposes it. Used to order the nodes as columns and 
    the cells as rows.
    """

    try:
        numpy_dataset = np.array(dataset).squeeze() # Get rid of extra dimensions
        logging.info(f'\t\tExtra dimension in dataset, squeezing...')
    except:
        numpy_dataset = np.array(dataset)

    transposed_dataset = np.transpose(numpy_dataset)

    logging.info(f'\t\tTransposed dataset shape: {transposed_dataset.shape}')

    return transposed_dataset

def calculate_hamming_distance(df_index, df_columns, list_1, list_2):
    """
    Iterate through two lists and calculate the Hamming distance between the items, returns
    a matrix of the Hamming distances with list_1 as the columns and list_2 as the rows.
    """

    # Convert lists to numpy arrays for vectorized operations    
    array1 = np.array([list(item) for item in list_1])
    array2 = np.array([list(item) for item in list_2])
    
    # Calculate the hamming distances using broadcasting and vectorization
    # The distances array will be a 2D array where the element at [i, j] is the
    # Hamming distance between array1[i] and array2[j]
    distances = np.sum(array1[:, None, :] != array2[None, :, :], axis=2)
    
    # Create the DataFrame from the distances array
    passed = False

    reduced_dimension = -1
    try:
        df = pd.DataFrame(distances, index=df_index, columns=df_columns)

    # This fixed and issue where the index and columns were not the same dimension, I'm not sure how but I'm leaving it
    except ValueError:
        while passed == False:
            try:
                df_columns_adjusted = df_columns[:reduced_dimension]
                logging.info(f'df_columns_adjusted = {df_columns_adjusted}')
                df = pd.DataFrame(distances, index=df_index, columns=df_columns_adjusted)
                passed = True
            except:
                reduced_dimension -= 1

    return df

def filter_attractors(dataframe):
    """
    Filters a Hamming distance dataframe to extract only the attractors with a minimum
    value for at least one cell, excludes attractors that do not best explain any cells.
    """

    # Ensure that each fo the columns contain numeric values
    for column in dataframe.columns:
        dataframe[column] = pd.to_numeric(dataframe[column], errors='coerce')

    # Exclude the first row and first column
    dataframe_sliced = dataframe.iloc[1:, 1:]

    # Find the row indices with the minimum Hammin distance for each column
    min_indices = dataframe_sliced.idxmin()

    # Keep only unique row indices
    unique_min_indices = min_indices.unique()

    # Filter out rows that do not contain a minimum Hamming distance for any cell
    filtered_df = dataframe.loc[unique_min_indices]

    # Find the attractor indices for these rows
    filtered_attractor_indices = filtered_df.index.tolist()

    return filtered_attractor_indices

def hierarchical_clustering(distance_matrix, num_nodes, show_plot):
    """
    Clusters the distance matrix of the attractors using Hierarchical Agglomerative Clustering.
    The cutoff for the clusters is set as <=25% the number of genes in the network. Returns the
    clusters and allows for a dendrogram to be displayed.
    """

    # Cluster the attractors using Hierarchical Agglomerative Clustering
    # Convert the distance matrix to a numpy array
    distance_matrix = distance_matrix.to_numpy()

    # Convert the square form distance matrix to a condensed form
    condensed_distance_matrix = ssd.squareform(distance_matrix)

    # Perform Hierarchical Agglomerative Clustering
    # 'ward' is one of the methods for calculating the distance between the newly formed cluster
    Z = sch.linkage(condensed_distance_matrix, method='ward')

    cutoff = 0.25 * num_nodes
    clusters = sch.fcluster(Z, t=cutoff, criterion='distance')

    logging.info(f'\t\tClustering cutoff value = {round(cutoff, 2)} (<=20% * number of genes {num_nodes})')
    
    # Plotting the dendrogram
    fig = plt.figure(figsize=(10, 7))
    plt.title("Attractors Hierarchical Clustering Dendrogram")
    
    dendrogram = sch.dendrogram(Z)
    plt.axhline(y=cutoff, color = 'r', linestyle='--')

    plt.ylabel('Distance')
    plt.xlabel('Attractors')

    if show_plot:
        plt.show()
    
    return clusters, fig

def get_starting_state(file):
    starting_state = []

    with open(file, 'r') as network_attractor:
        for line in network_attractor:
            node_starting_state = int([random.choice([0,1]) for i in line])
            starting_state.append(node_starting_state)
    
    return starting_state

def use_cell_starting_state(file):
    with open(file, 'rb') as f:
        network_object = pickle.load(f)
        return network_object

def simulate_network(nodes, filename):
    steps = 20

    starting_state = []
    for i in filename:
        starting_state.append(i[0])

    total_simulation_states = vectorized_run_simulation(nodes, starting_state, steps)
    
    simulation_results = [[int(item) for item in sublist] for sublist in total_simulation_states]

    return simulation_results

def evaluate_expression(data, expression):
    expression = expression.replace('and', '&').replace('or', '|').replace('not', '~')
    if any(op in expression for op in ['&', '|', '~']):
        local_vars = {key: np.array(value).astype(bool) for key, value in data.items()}
    else:
        local_vars = {key: np.array(value) for key, value in data.items()}
    return ne.evaluate(expression, local_dict=local_vars)

def vectorized_run_simulation(nodes, starting_state, steps):
    total_simulation_states = []

    # Run the simulation
    for step in range(steps):
        step_expression = []

        # Iterate through each node in the network
        for node in nodes:

            # Initialize A, B, C to False by default (adjust according to what makes sense in context)
            A, B, C = (False,) * 3
            
            data = {}
            incoming_node_indices = [predecessor_index for predecessor_index in node.predecessors]

            # Get the rows in the dataset for the incoming nodes
            if step == 0:
                if len(incoming_node_indices) > 0:
                    data['A'] = starting_state[incoming_node_indices[0]]
                if len(incoming_node_indices) > 1:
                    data['B'] = starting_state[incoming_node_indices[1]]
                if len(incoming_node_indices) > 2:
                    data['C'] = starting_state[incoming_node_indices[2]]
                if len(incoming_node_indices) > 3:
                    data['D'] = starting_state[incoming_node_indices[3]]
            else:
                if len(incoming_node_indices) > 0:
                    data['A'] = total_simulation_states[step-1][incoming_node_indices[0]]
                if len(incoming_node_indices) > 1:
                    data['B'] = total_simulation_states[step-1][incoming_node_indices[1]]
                if len(incoming_node_indices) > 2:
                    data['C'] = total_simulation_states[step-1][incoming_node_indices[2]]
                if len(incoming_node_indices) > 3:
                    data['D'] = total_simulation_states[step-1][incoming_node_indices[3]]

            next_step_node_expression = evaluate_expression(data, node.calculation_function)

            # Save the expression for the node for this step
            step_expression.append(next_step_node_expression)

        # Save the expression 
        total_simulation_states.append(step_expression)

    return total_simulation_states
    

def visualize_simulation(net, time_series_data, network, show_simulation):
    node_indices = [node.index for node in network.nodes]
    pos = nx.spring_layout(net, k=1, iterations=100)  # Pre-compute layout
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_title(f"Network: {network.name}")
    
    # Initial drawing of the graph
    nx.draw(net, pos, ax=ax, with_labels=True)
    node_collections = ax.collections[0]  # Get the collection of nodes
    frame_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

    def update_graph(frame):
        # Update node colors without redrawing everything
        colors = []
        for i in node_indices:
            if time_series_data[frame-1][i] != time_series_data[frame][i]:
                colors.append('gold' if time_series_data[frame][i] == 1 else 'red')
            else:
                colors.append('green' if time_series_data[frame][i] == 1 else 'red')

        node_collections.set_color(colors)  # Update colors directly
    
        frame_text.set_text(f'Frame: {frame+1}')

    return fig

def create_heatmap(path, title):
    data = []
    gene_names = []

    with open(path, 'r') as file:
        for line in file:
            line = line.strip().split(',')
            gene_name = line[0]
            time_data = [int(i) for i in line[1:]]
            data.append(time_data)
            gene_names.append(gene_name)
    
    num_genes = len(data)
    num_time_steps = len(data[0])

    # Adjusting the data to fit the provided shape
    data_array = np.array(data).reshape((num_genes, num_time_steps))

    # Create a custom colormap
    cmap = mcolors.ListedColormap(['grey', 'green'])
    bounds = [0, 0.5, 1]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    # Create a heatmap
    plot = plt.figure(figsize=(12, 12))
    sns.heatmap(data_array, cmap='Greys', yticklabels=gene_names, xticklabels=True)
    plt.title(title)
    plt.xlabel('Time Steps')
    plt.ylabel('Genes')
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    # plt.tight_layout()

    legend_elements = [
        Patch(facecolor='grey', edgecolor='grey', label='Gene Inactive'),
        Patch(facecolor='black', edgecolor='black', label='Gene Active')
    ]
    plt.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1), title="Legend")


    plt.subplots_adjust(top=0.958, bottom=0.07, left=0.076, right=0.85, hspace=2, wspace=1)

    return plot

def read_data_from_directory(directory, simulated_cells):
    print(f'Reading data from the files')
    dataframes = {}

    for filename in os.listdir(directory):
        # Extract the cell number from the filename
        cell_num = int(''.join(filename.split('_trajectory.csv')[0].split('cell_')[1]))

        if cell_num in simulated_cells and filename.endswith("_trajectory.csv"):
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
            distance, _ = fastdtw(ts1, ts2, radius=1, dist=euclidean)
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
    plt.figure(figsize=(8, 10))
    dendrogram(Z, labels=distance_matrix.index, orientation='top')
    plt.title('Hierarchical Clustering Dendrogram', fontsize=12)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=8, rotation=90)
    plt.xlabel('Cells', fontsize=8)
    plt.ylabel('Distance', fontsize=8)
    plt.tight_layout()

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
    plt.figure(figsize=(8, 10))
    sns.heatmap(df, cmap='Greys', yticklabels=True)
    plt.title(f'Average Gene Expression Heatmap for Cluster {cluster}', fontsize=12)
    plt.xlabel(xlabel='Simulation Time Steps', fontsize=12)
    plt.ylabel(ylabel='Gene', fontsize=12)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.tight_layout()
    plt.show()


def plot_heatmap(distance_matrix, file_names):
    # Convert the square distance matrix to a condensed distance matrix
    condensed_distance_matrix = squareform(distance_matrix)
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_distance_matrix, method='average')
    
    # Create a dendrogram to get the order of the leaves
    dendro = dendrogram(linkage_matrix, no_plot=True)
    order = dendro['leaves']
    
    # Reorder the distance matrix
    reordered_matrix = distance_matrix[np.ix_(order, order)]
    reordered_file_names = [file_names[i].split('_trajectory')[0] for i in order]
    
    plt.figure(figsize=(8, 9))
    sns.heatmap(reordered_matrix, xticklabels=reordered_file_names, yticklabels=reordered_file_names, cmap='Greys', annot=False)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.title("DTW Distance Heatmap")
    plt.tight_layout()
    plt.show()



def save_attractor_simulation(filename, network, simulated_attractor):
    # Save the attractor simulation to a file
    with open(filename, 'w') as file:
        simulated_attractor = np.array(simulated_attractor).T
        for gene_num, expression in enumerate(simulated_attractor):
            file.write(f'{network.nodes[gene_num].name},{",".join([str(i) for i in list(expression)])}\n')


# If you want to run the attractor analysis by itself
if __name__ == '__main__':

    # Set the logging level for output
    logging.basicConfig(format='%(message)s', level=logging.INFO)

    # Allow the user to either add in the dataset name and network name from the command line or as a prompt
    parser = argparse.ArgumentParser()

    dataset_name, show_simulation = attractor_analysis_arguments(parser)

    # Load the network pickle files
    all_networks = []

    # Load the cell and network objects for the dataset
    cell_population = pickle.load(open(f'{file_paths["pickle_files"]}/{dataset_name}_pickle_files/cells_pickle_file/{dataset_name}.cells.pickle', "rb"))
    
    cells = cell_population.cells

    logging.info(f'\nRunning attractor analysis for all networks...')
    pickle_file_path = f'{file_paths["pickle_files"]}/{dataset_name}_pickle_files/network_pickle_files/'
    for pickle_file in glob.glob(pickle_file_path + str(dataset_name) + "_" + "*" + ".network.pickle"):
        if pickle_file:
            logging.info(f'\t\tLoading data file: {pickle_file}')
            network = pickle.load(open(pickle_file, "rb"))
            all_networks.append(network)
        else:
            assert FileNotFoundError("Network pickle file not found")
    
    try:
        dataset_found = all_networks[0].name
        logging.info(f'\nFound dataset!')

    except IndexError:
        error_message = "No networks loaded, check to make sure network pickle files exist in 'pickle_files/'"
        logging.error(error_message)
        raise Exception(error_message)
    
    for network in all_networks:
        # Run the pathway analysis for each of the networks

        logging.info(f'\n----- ATTRACTOR ANALYSIS -----')

        # Convert the network's sparse dataset to a dense one
        dataset = network.dataset
        dense_dataset = np.array(dataset.todense())
        network_name = network.name

        # Run simulations for 1% of the cells (change to something like 10% in the future)
        logging.info('\tSimulating cell trajectories')
        num_simulations = int(round(len(dense_dataset[1]) / 100, 0)) 
        simulated_cells = []
        with alive_bar(num_simulations) as bar:
            for i in range(num_simulations):
                # Select a random column from the network dataset
                cell_index = np.random.choice(dense_dataset.shape[1])

                # Reads in all of the rows for that columns
                selected_column = np.array([random.choice([0,1]) for _ in dense_dataset[:, cell_index]])

                simulated_cells.append(cell_index)

                # Transposes the list of gene expression into a column
                transposed_random_column = selected_column.reshape(-1,1)

                # Specify outfile path for the simulation results
                outfile_folder = f'{file_paths["trajectories"]}/{dataset_name}_{network_name}'
                png_folder = f'{outfile_folder}/png_files'
                text_folder = f'{outfile_folder}/text_files'

                os.makedirs(outfile_folder, exist_ok=True)
                os.makedirs(png_folder, exist_ok=True)
                os.makedirs(text_folder, exist_ok=True)
                
                # Simulate the network
                simulated_attractor = simulate_network(network.nodes, transposed_random_column)

                # Visualize the network simulation results
                fig = visualize_simulation(network.network, simulated_attractor, network, "False")

                # Save the attractor states to a csv file
                save_attractor_simulation(f'{text_folder}/cell_{cell_index}_trajectory.csv', network, simulated_attractor)
                plt.close(fig)

                # Create a heatmap of the expression for easier attractor visualization
                heatmap = create_heatmap(f'{text_folder}/cell_{cell_index}_trajectory.csv', f'Simulation for {dataset_name} {network_name} cell {cell_index} pathway ')
                # heatmap.show()

                # Saves a png of the results
                heatmap.savefig(f'{png_folder}/cell_{cell_index}_trajectory.png', format='png')
                plt.close(heatmap)
                bar()

        # Directory containing the trajectory files
        directory = f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/text_files'

        dataframes = read_data_from_directory(directory, simulated_cells)

        logging.info('\tExtracting time series data and computing dynamic time warping distances')
        time_series_data = extract_time_series(dataframes)
        dtw_distances = compute_dtw_distances(time_series_data, directory)
        
        file_names = list(dataframes.keys())
        distance_matrix = create_distance_matrix(dtw_distances, file_names)

        cluster_dict = find_similar_files(dtw_distances)

        for cluster, cell_list in cluster_dict.items():
            print(f'Summarizing cluster {cluster}')
            summarize_clusters(directory, cell_list, cluster)
        
        plot_heatmap(distance_matrix, file_names)

        cluster_fig = run_attractor_analysis(network, cells)
        cluster_path = f'{file_paths["attractor_analysis_output"]}/{dataset_name}_attractors'
        os.makedirs(cluster_path, exist_ok=True)

        cluster_filename = f'{cluster_path}/{network.name}_cluster.png'

        cluster_fig.savefig(cluster_filename, bbox_inches='tight')
        plt.close(cluster_fig)  # Close the figure after saving


        logging.info(f'\n-----PATHWAY ANALYSIS RESULTS -----')
        logging.info(f'\nNETWORK {network.name.upper()}')

        attractor_counts = {}
        total_cells = 0
        for cell in network.cell_map.keys():
            attractor = network.cell_map[cell]
            if attractor not in attractor_counts:
                attractor_counts[attractor] = 1
            else:
                attractor_counts[attractor] += 1
            total_cells += 1

        for attractor_num, num_cells in attractor_counts.items():
            logging.info(f'\tAttractor {attractor_num} contains {num_cells} cells ({round(num_cells / total_cells * 100, 3)}%)')

        for attractor_num, representative_attractor in network.representative_attractors.items():
            if attractor_num in attractor_counts:

                # Create a dataset- and network-specific directory for the output
                output_directory = f'{file_paths["attractor_analysis_output"]}/{dataset_name}_attractors/{network.name}_attractors/attractor_{str(attractor_num)}'
                os.makedirs(output_directory, exist_ok=True)

                filename = f'{dataset_name}_{network.name}_attractor_{attractor_num}.txt'
                np.savetxt(output_directory + "/" + filename, representative_attractor.T, fmt='%d')

                logging.info(f'\tSaved representative attractor {attractor_num}')
                logging.debug(f'\t{representative_attractor}')
                
                simulation_results = simulate_network(network.nodes, output_directory + "/" + filename)

                svg_output_path = f'{output_directory}/attractor_{attractor_num}_simulation_results.svg'
                png_output_path = f'{output_directory}/attractor_{attractor_num}_simulation_results.png'

                fig = visualize_simulation(network.network, simulation_results, network, show_simulation)

                plt.savefig(svg_output_path, format='svg')
                plt.savefig(png_output_path, format='png')
                plt.close(fig)  # Close the figure after saving

        logging.info(f'\n\tSaved attractor analysis results to "attractor_analysis_output/{dataset_name}_attractors/{network.name}_attractors')


        logging.info(f'\nAdding representative attractor map to network pickle files:')
        network_directory_path = f'{file_paths["pickle_files"]}/{dataset_name}_pickle_files/network_pickle_files'
        os.makedirs(network_directory_path, exist_ok=True)

        network_pickle_file_path = f'{network_directory_path}/{dataset_name}_{network.name}.network.pickle'
        logging.info(f'\tFile: {dataset_name}_{network.name}.network.pickle')
        pickle.dump(network, open(network_pickle_file_path, "wb"))

        cell_pickle_file_path = f'{file_paths["pickle_files"]}/{dataset_name}_pickle_files/cells_pickle_file'
        os.makedirs(cell_pickle_file_path, exist_ok=True)

        cells_pickle_file = f'{dataset_name}.cells.pickle'
        with open(f'{cell_pickle_file_path}/{cells_pickle_file}', "wb") as file:
            pickle.dump(cell_population, file)