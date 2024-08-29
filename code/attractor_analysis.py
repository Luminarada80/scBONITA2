from importance_scores import CalculateImportanceScore
import logging
import pandas as pd
import numpy as np
import scipy
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

# def run_attractor_analysis(network, cells):
#     """
#     Runs the attractor analysis for the dataset and all of the networks.
#     """
#     # Set up the attractor analysis with the dataset and the parsed networks
#
#     logging.info(f'\nNetwork: {network.name}')
#     # Generate attractor
#     logging.info(f'\tGenerating attractors...')
#
#         # Specifies the path to the correct network pickle file
#     network_pickle_file = f'{file_paths["pickle_files"]}/{dataset_name}_pickle_files/network_pickle_files/{dataset_name}_{network_name}.network.pickle'
#
#     network_attractors, simulated_dataset = generate_attractors(network.nodes, network.dataset)
#     attractors_start = [attractor[:, 0] for _, attractor in network_attractors.items()]
#
#     # Convert the sparse dataset to a dense one
#     if isinstance(simulated_dataset, spmatrix):
#         dense_dataset = simulated_dataset.todense()
#     else:
#         dense_dataset = simulated_dataset
#
#     # Generate Hamming distance matrix between cells and attractors
#     logging.info(f'\tCalculating hamming distance between cells and attractors')
#     logging.info(f'\tTransposed dataset shape: {dense_dataset.shape}')
#     transposed_dataset = transpose_dataset(dense_dataset)
#     num_nodes = len(network.nodes)
#     num_cells = range(transposed_dataset.shape[0])
#     num_attractors = range(len(attractors_start))
#
#     logging.info(f'\t\tNodes: {num_nodes}')
#     logging.info(f'\t\tCells: {transposed_dataset.shape[0]}')
#     logging.info(f'\t\tAttractors: {len(attractors_start)}')
#     cell_attractor_hamming_df = calculate_hamming_distance(num_attractors, num_cells, transposed_dataset, attractors_start)
#
#     # Filter the attractors to only keep those that most closely match the expression of at least one cell
#     filtered_attractor_indices = filter_attractors(cell_attractor_hamming_df)
#     filtered_attractors = [attractors_start[i] for i in filtered_attractor_indices]
#     num_filtered_attractors = range(len(filtered_attractors))
#
#     # Create a distance matrix between each of the attractors
#     logging.info(f'\tGenerating attractor distance matrix...')
#     attractor_distance_matrix = calculate_hamming_distance(num_filtered_attractors, num_filtered_attractors, filtered_attractors, filtered_attractors)
#
#     # Cluster the attractors using hierarchical agglomerative clustering
#     logging.info(f'\tClustering the attractors...')
#     clusters, cluster_fig = hierarchical_clustering(attractor_distance_matrix, len(network.nodes))
#
#     clustered_attractors = {}
#     for i, cluster_num in enumerate(clusters):
#         if cluster_num not in clustered_attractors:
#             clustered_attractors[cluster_num] = {}
#         clustered_attractors[cluster_num][filtered_attractor_indices[i]] = filtered_attractors[i]
#
#     # Find the Hamming Distance between the cells and the attractors within the clusters
#     logging.info(f'\tCalculating Hamming distance between cells and clustered attractors')
#
#     # For each of the cells, find the hamming distance to the best cell
#     logging.info(f'Number of cells in the full dataset: {network.dataset.shape[1]}')
#     if issparse(network.dataset):
#         full_dataset = network.dataset.toarray().T
#     else:
#         full_dataset = network.dataset.T
#
#     cell_map = {}
#     for cell_num, cell in enumerate(full_dataset):
#         cell_map[cell_num] = {}
#
#     for cluster, attractors in clustered_attractors.items():
#         logging.debug(f'\nCluster {cluster}')
#
#         num_filtered_attractors = range(len(attractors.values()))
#
#         mapped_cluster_attractors = calculate_hamming_distance(num_filtered_attractors, num_cells, attractors.values(), transposed_dataset)
#
#         # transposed_mapped_cluster_attractors = transpose_dataset(mapped_cluster_attractors)
#         logging.debug(f'Mapped cluster attractors:\n {mapped_cluster_attractors}')
#
#         # Map the cells to the attractors within the clusters
#         filtered_attractor_indices = list(clustered_attractors[cluster].keys())
#         logging.debug(f'Filtered attractor indices: {filtered_attractor_indices}')
#         mapped_cluster_attractors.index = filtered_attractor_indices # type: ignore
#         logging.debug(f'Cluster {cluster} mapped attractors')
#         logging.debug(mapped_cluster_attractors)
#
#         # Find the total Hamming distance between each attractor and all of the cells
#         total_hamming_distance = mapped_cluster_attractors.sum(axis=1)
#         min_sum_index = total_hamming_distance.idxmin()
#
#         representative_attractor = attractors[min_sum_index]
#
#         # Find the hamming distance for each cell for each cluster
#         for cell_num, cell in enumerate(full_dataset):
#             hamming_distance = 0
#             for gene_num, gene in enumerate(cell):
#                 if cell[gene_num] != representative_attractor[gene_num]:
#                    hamming_distance += 1
#
#             cell_map[cell_num][cluster] = hamming_distance
#
#
#         logging.debug(f'Representative attractor for cluster {cluster} (index {min_sum_index}):')
#         logging.debug(f'\n{attractors[min_sum_index]}')
#
#         network.representative_attractors[cluster] = representative_attractor
#
#     # Calculate which cluster each cell should map to
#     min_hamming_cluster = {}
#     for cell_num, clusters in cell_map.items():
#
#         # Find the cluster with the minimum Hamming distance
#         min_cluster = min(clusters, key=clusters.get)
#         min_hamming_cluster[cell_num] = min_cluster
#         # cell_object = cells[cell_num]
#
#     network.cell_map = min_hamming_cluster
#
#     for cell_object in cells:
#         # Add the best attractor for each cell to that cell object
#         cell_object.attractor_dict[network.name] = min_hamming_cluster[cell_object.index]
#
#     return cluster_fig

# def generate_attractors(nodes, dataset):
#     """
#     Uses the vectorized simulation function from the importance score calculations to find
#     the attractors for the network.
#     """
#
#     calculate_importance_score = CalculateImportanceScore(nodes, dataset)
#
#     simulated_dataset = calculate_importance_score.dataset
#     network_attractors = calculate_importance_score.vectorized_run_simulation()
#
#     return network_attractors, simulated_dataset


# def extract_node_rows(dataset, nodes):
#     """
#     Create a dense dataset only containing the nodes that are present in the network from a
#     sparse dataset
#     """
#
#     # Convert the sparse dataset to a dense one
#     if isinstance(dataset, spmatrix):
#         dense_dataset = dataset.todense()
#     else:
#         dense_dataset = dataset
#
#     # Find the row indices in the dataset for the nodes in the network
#     node_indices = []
#     for node in nodes:
#         if node.dataset_index not in node_indices:
#             node_indices.append(node.dataset_index)
#
#     return dense_dataset[node_indices]

# def transpose_dataset(dataset):
#     """
#     Converts a dataset to a numpy array and transposes it. Used to order the nodes as columns and
#     the cells as rows.
#     """
#
#     try:
#         numpy_dataset = np.array(dataset).squeeze() # Get rid of extra dimensions
#         logging.info(f'\t\tExtra dimension in dataset, squeezing...')
#     except:
#         numpy_dataset = np.array(dataset)
#
#     transposed_dataset = np.transpose(numpy_dataset)
#
#     logging.info(f'\t\tTransposed dataset shape: {transposed_dataset.shape}')
#
#     return transposed_dataset

# def calculate_hamming_distance(df_index, df_columns, list_1, list_2):
#     """
#     Iterate through two lists and calculate the Hamming distance between the items, returns
#     a matrix of the Hamming distances with list_1 as the columns and list_2 as the rows.
#     """
#
#     # Convert lists to numpy arrays for vectorized operations
#     array1 = np.array([list(item) for item in list_1])
#     array2 = np.array([list(item) for item in list_2])
#
#     # Calculate the hamming distances using broadcasting and vectorization
#     # The distances array will be a 2D array where the element at [i, j] is the
#     # Hamming distance between array1[i] and array2[j]
#     distances = np.sum(array1[:, None, :] != array2[None, :, :], axis=2)
#
#     # Create the DataFrame from the distances array
#     passed = False
#
#     reduced_dimension = -1
#     try:
#         df = pd.DataFrame(distances, index=df_index, columns=df_columns)
#
#     # This fixed and issue where the index and columns were not the same dimension, I'm not sure how but I'm leaving it
#     except ValueError:
#         while passed == False:
#             try:
#                 df_columns_adjusted = df_columns[:reduced_dimension]
#                 logging.info(f'df_columns_adjusted = {df_columns_adjusted}')
#                 df = pd.DataFrame(distances, index=df_index, columns=df_columns_adjusted)
#                 passed = True
#             except:
#                 reduced_dimension -= 1
#
#     return df
#
# def filter_attractors(dataframe):
#     """
#     Filters a Hamming distance dataframe to extract only the attractors with a minimum
#     value for at least one cell, excludes attractors that do not best explain any cells.
#     """
#
#     # Ensure that each fo the columns contain numeric values
#     for column in dataframe.columns:
#         dataframe[column] = pd.to_numeric(dataframe[column], errors='coerce')
#
#     # Exclude the first row and first column
#     dataframe_sliced = dataframe.iloc[1:, 1:]
#
#     # Find the row indices with the minimum Hammin distance for each column
#     min_indices = dataframe_sliced.idxmin()
#
#     # Keep only unique row indices
#     unique_min_indices = min_indices.unique()
#
#     # Filter out rows that do not contain a minimum Hamming distance for any cell
#     filtered_df = dataframe.loc[unique_min_indices]
#
#     # Find the attractor indices for these rows
#     filtered_attractor_indices = filtered_df.index.tolist()
#
#     return filtered_attractor_indices
#
# def get_starting_state(file):
#     starting_state = []
#
#     with open(file, 'r') as network_attractor:
#         for line in network_attractor:
#             node_starting_state = int([random.choice([0,1]) for i in line])
#             starting_state.append(node_starting_state)
#
#     return starting_state
#
# def use_cell_starting_state(file):
#     with open(file, 'rb') as f:
#         network_object = pickle.load(f)
#         return network_object

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
    

def visualize_simulation(net, cell_trajectory_dict, network, show_simulation):
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
            if cell_trajectory_dict[frame-1][i] != cell_trajectory_dict[frame][i]:
                colors.append('gold' if cell_trajectory_dict[frame][i] == 1 else 'red')
            else:
                colors.append('green' if cell_trajectory_dict[frame][i] == 1 else 'red')

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

def read_data_from_directory(directory, trajectory_files_parsed, num_cells_to_parse):
    """
    Reads in the cell trajectory files and creates a dictionary of dataframes with the file
    name as the key and the dataframe as the falue
    """
    dataframes = {}

    num_cells_parsed = 0
    for filename in os.listdir(directory):
        
        # Keep track of how many files are parsed in this batch and which cells trajectories were analyzed
        if num_cells_parsed <= num_cells_to_parse-1 and filename.endswith("_trajectory.csv") and filename not in trajectory_files_parsed:
            trajectory_files_parsed.append(filename)
            
            # Extract the cell number from the filename
            filepath = os.path.join(directory, filename)
            df = pd.read_csv(filepath, header=None)
            df.columns = ['Gene'] + [f'Time{i}' for i in range(1, df.shape[1])]
            df.set_index('Gene', inplace=True)
            
            # Extract the gene values from the dataset and append them to a dictionary with the gene name as key
            dataframes[filename] = {gene: df.loc[gene].values for gene in df.index}
            num_cells_parsed += 1

    return dataframes, trajectory_files_parsed

def compute_dtw_distance_pair(cell1, cell2, cell_trajectory_dict):
    """
    Computes the dynamic time warping distance between two cell trajectories for each gene.
    
    Parameters:
        cell1: The name of the first cell to be compared
        cell2: The name of the second cell to be compared
        cell_trajectory_dict: A dictionary with cell names as keys and a gene name: trajectory dict as values

    Returns: 
        cell1, cell2: Each of the cells names
        total_distance: The summed DTW distances between the trajectory of each gene
    """
    distances = {}

    for gene in cell_trajectory_dict[cell1].keys():
        if gene in cell_trajectory_dict[cell2].keys():
            ts1 = cell_trajectory_dict[cell1][gene]
            ts2 = cell_trajectory_dict[cell2][gene]
            distance, _ = fastdtw(ts1, ts2, radius=1, dist=2)
            distances[gene] = distance
    total_distance = sum(distances.values()) if distances else float('inf')
    return (cell1, cell2, total_distance)

def compute_dtw_distances(cell_trajectory_dict, output_directory):
    dtw_distances = {}
    cell_names = list(cell_trajectory_dict.keys())
    total_combinations = len(cell_names) * (len(cell_names) - 1) // 2
    
    tasks = []
    # Compute the DTW distance between each pair of cells
    with ProcessPoolExecutor() as executor, alive_bar(total_combinations) as bar:
        for i in range(len(cell_names)):
            for j in range(i + 1, len(cell_names)):
                cell1 = cell_names[i]
                cell2 = cell_names[j]
                tasks.append(executor.submit(compute_dtw_distance_pair, cell1, cell2, cell_trajectory_dict))

        # Once completed, add the distances for each cell to a dictionary
        for future in as_completed(tasks):
            cell1, cell2, total_distance = future.result()
            dtw_distances[(cell1, cell2)] = total_distance
            bar()

    # Write out each cell-cell distance to an outfile
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

def hierarchical_clustering(dtw_distances, num_clusters):

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

    
    plt.savefig(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/dendrogram.png')
    
    # Set a threshold and get the clusters
    if num_clusters == 0:
        logging.info(f'Please open "scBONITA_output/trajectories/{dataset_name}_{network_name}/dendrogram.png"')
        num_clusters = int(input('How many clusters?: '))
    clusters = fcluster(Z, num_clusters, criterion='maxclust')

    # Organize cells by clusters
    cluster_dict = {}
    for cell, cluster_id in zip(cells, clusters):
        if cluster_id not in cluster_dict:
            cluster_dict[cluster_id] = []
        cluster_dict[cluster_id].append(cell)

    plt.close()

    return cluster_dict, num_clusters


def summarize_clusters(directory, cell_names):
    gene_expr_dict = {}

    # Finds the path to all of the trajectory csv files
    trajectory_files = []
    for filename in os.listdir(directory):
        if filename.endswith("_trajectory.csv"):
            trajectory_files.append(os.path.join(directory, filename))

    # finds the files of the cells in the cluster
    files_to_open = []
    for cell in cell_names:
        for file_name in trajectory_files:
            if cell in file_name:
                files_to_open.append(file_name)

    # Reads in the gene expression values from the simulation file
    for file_path in files_to_open:
        with open(file_path, 'r') as sim_file:
            for line in sim_file:
                line = line.strip().split(',')
                gene_name = line[0]
                gene_expression = [int(i) for i in line[1:]]
                
                if gene_name not in gene_expr_dict:
                    gene_expr_dict[gene_name] = []
                
                gene_expr_dict[gene_name].append(gene_expression)
    
    # print(f'Opening {len(files_to_open)} trajectory files')

    gene_avg_expr = {}
    for gene, simulation_results in gene_expr_dict.items():

        if gene not in gene_avg_expr:
            gene_avg_expr[gene] = []

        transposed_data = list(map(list, zip(*simulation_results)))

        # Finds the average gene expression at each point along the simulation
        for i in transposed_data:
            gene_avg_expr[gene].append(statistics.mean(i))
    
    # Convert the average gene expression dictionary to a DataFrame
    df = pd.DataFrame(gene_avg_expr)

    # Transpose the DataFrame
    df = df.transpose()

    return df


def plot_distance_matrix(distance_matrix, file_names):
    # Convert the square distance matrix to a condensed distance matrix
    condensed_distance_matrix = squareform(distance_matrix)
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_distance_matrix, method='average')
    
    # Create a dendrogram to get the order of the leaves
    dendro = dendrogram(linkage_matrix, no_plot=True)
    order = dendro['leaves']
    
    # Reorder the distance matrix
    reordered_matrix = distance_matrix[np.ix_(order, order)]
    reordered_file_names = [i for i in order]
    
    plt.figure(figsize=(8, 9))
    sns.heatmap(reordered_matrix, xticklabels=reordered_file_names, yticklabels=reordered_file_names, cmap='Greys', annot=False)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.title("DTW Distance Heatmap")
    plt.tight_layout()
    plt.show()
    plt.savefig(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/distance_heatmap')
    plt.close()



def save_attractor_simulation(filename, network, simulated_attractor):
    # Save the attractor simulation to a file
    with open(filename, 'w') as file:
        simulated_attractor = np.array(simulated_attractor).T
        for gene_num, expression in enumerate(simulated_attractor):
            file.write(f'{network.nodes[gene_num].name},{",".join([str(i) for i in list(expression)])}\n')


def simulate_cells(dense_dataset, num_simulations):
    # Simulate cell trajectories
    logging.info(f'\tSimulating {num_simulations} cell trajectories')
    simulated_cells = []
    with alive_bar(num_simulations) as bar:
        for i in range(num_simulations):

            # Get the list of existing files in the cell simulation directory
            existing_files = os.listdir(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/text_files')
            existing_indices = set()
            for filename in existing_files:
                if filename.startswith("cell_") and filename.endswith("_trajectory.csv"):
                    index = int(filename.split('_')[1])
                    existing_indices.add(index)

            # Randomly select a cell index that is not in existing_indices
            cell_index = np.random.choice(dense_dataset.shape[1])
            while cell_index in existing_indices:
                cell_index = np.random.choice(dense_dataset.shape[1])

            # Reads in all rows for that columns
            selected_column = np.array([random.choice([0,1]) for _ in dense_dataset[:, cell_index]])

            # Record which cells were simulated
            simulated_cells.append(cell_index)

            # Transposes the list of gene expression into a column
            transposed_random_column = selected_column.reshape(-1,1)
            
            # Simulate the network
            trajectory = simulate_network(network.nodes, transposed_random_column)

            # Visualize the network simulation results
            fig = visualize_simulation(network.network, trajectory, network, "False")

            # Save the attractor states to a csv file
            save_attractor_simulation(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/text_files/cell_{cell_index}_trajectory.csv',
                                        network,
                                        trajectory)
            plt.close(fig)

            # Create a heatmap to visualize the trajectories
            heatmap = create_heatmap(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/text_files/cell_{cell_index}_trajectory.csv',
                                     f'Simulation for {dataset_name} {network_name} cell {cell_index} pathway ')

            # Saves a png of the results
            heatmap.savefig(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/png_files/cell_{cell_index}_trajectory.png', format='png')
            plt.close(heatmap)
            bar()


def plot_average_trajectory(df, cluster):        
    # Create average trajectory directory
    os.makedirs(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/avg_chunks', exist_ok=True)

    # Create the heatmap
    plt.figure(figsize=(8, 10))
    sns.heatmap(df, cmap='Greys', yticklabels=True, vmin=0, vmax=1)
    plt.title(f'Average Gene Expression Heatmap for Cluster {cluster}', fontsize=12)
    plt.xlabel(xlabel='Simulation Time Steps', fontsize=12)
    plt.ylabel(ylabel='Gene', fontsize=12)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.tight_layout()
    plt.savefig(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/avg_chunks/cluster_{cluster}_summary')
    plt.close()


def create_trajectory_chunks(num_chunks, num_clusters, output_directory):
    trajectory_files_parsed = []

    # Chunk the data down to create smaller averaged trajectory clusters
    cluster_chunks = {}
    cells_in_chunks = {}
    for chunk in range(num_chunks):

        logging.info(f'Creating chunk {chunk+1}:')
        # Reads in the cell trajectories from the simulations and creates a dataframe from them
        cell_trajectory_dict, trajectory_files_parsed = read_data_from_directory(output_directory, trajectory_files_parsed,
                                                                             num_cells_per_chunk)

        # Computes the pairwise DTW distances between cells
        dtw_distances = compute_dtw_distances(cell_trajectory_dict, output_directory)

        # Calculates pairwise distance matrix between the cells
        file_names = list(cell_trajectory_dict.keys())
        distance_matrix = create_distance_matrix(dtw_distances, file_names)

        # Performs hierarchical clustering to find the number of clusters
        cluster_dict, num_clusters = hierarchical_clustering(dtw_distances, num_clusters)

        for cluster, cell_list in cluster_dict.items():
            df = summarize_clusters(output_directory, cell_list)

            # Binarize the DataFrame
            df_binarized = pd.DataFrame(np.where(df < 0.5, 0, 1), index=df.index, columns=df.columns)

            # Save the cluster chunk to a dictionary
            cluster_chunks[f'{chunk}:{cluster}'] = df  # df_binarized

            # Keep track of which cells are in which chunks for later
            if chunk not in cells_in_chunks:
                cells_in_chunks[chunk] = {}

            cells_in_chunks[chunk][cluster] = cell_list

    return cluster_chunks, cells_in_chunks, num_clusters

def calculate_dtw(num_files, output_directory, num_cells_per_chunk):
    trajectory_files_parsed = []
    num_clusters = 0
    num_chunks = round(num_files / (num_cells_per_chunk))

    logging.info(f'Creating {num_chunks} chunks ({num_files} cells / {num_cells_per_chunk} cells per chunk)')

    # Ensures there is always one chunk
    if num_chunks == 0:
        num_chunks = 1

    cluster_chunks, cells_in_chunks, num_clusters = create_trajectory_chunks(num_chunks, num_clusters, output_directory)

    logging.info(f'{len(trajectory_files_parsed)} cell trajectories compared, split into {num_chunks} chunks and {num_clusters} clusters')

    logging.info(f'\nComparing chunks')
    chunk_dtw_distances = compute_dtw_distances(cluster_chunks, output_directory)

    # Calculates pairwise distance matrix between the cells
    chunk_names = list(cluster_chunks.keys())
    distance_matrix = create_distance_matrix(chunk_dtw_distances, chunk_names)

    # Performs hierarchical clustering to find the number of clusters
    group_cluster_dict, num_clusters = hierarchical_clustering(chunk_dtw_distances, num_clusters)

    # for chunk, cluster in cells_in_chunks.items():
    #     print(f'Chunk {chunk}')
    #     for clust_num, cells in cluster.items():
    #         print(f'\tCluster {clust_num}')
    #         print(f'\t\tFirst 5 cells: {cells[0:5]}')

    plot_distance_matrix(distance_matrix, group_cluster_dict.keys())

    cells_in_cluster = {}
    for cluster, cell_list in group_cluster_dict.items():
        num_cells_in_cluster = 0
        for chunk_cluster in cell_list:
            
            chunk = chunk_cluster.split(':')[0]
            cluster = chunk_cluster.split(':')[1]
            
            # print(f'\tChunk {chunk}\n\t\t cluster {cluster}\n\t\t\t cells = {cells_in_chunks[int(chunk)][int(cluster)][0:5]}')
            # print(f'Chunk {chunk} cluster {cluster}')
            num_cells_in_cluster += len(cells_in_chunks[int(chunk)][int(cluster)])

            if cluster not in cells_in_cluster:
                cells_in_cluster[cluster] = []
            
            cells_in_cluster[cluster].extend(cells_in_chunks[int(chunk)][int(cluster)])
        

    #     print(f'Total number of cells in Cluster {cluster} = {num_cells_in_cluster}')
    # print(f'Cells in cluster\n{cells_in_cluster}')

    for cluster, cluster_group in group_cluster_dict.items():
        logging.info(f'Summarizing cluster {cluster}')

        df = summarize_clusters(output_directory, cells_in_cluster[str(cluster)])

        # Create average trajectory directory
        os.makedirs(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/avg_chunks', exist_ok=True)
        
        plot_average_trajectory(df, cluster)




# If you want to run the attractor analysis by itself
if __name__ == '__main__':

    # Set the logging level for output
    logging.basicConfig(format='%(message)s', level=logging.INFO)

    # Allow the user to either add in the dataset name and network name from the command line or as a prompt
    parser = argparse.ArgumentParser()

    # Read in the user arguments from the command line or bash script
    dataset_name, show_simulation = attractor_analysis_arguments(parser)

    # Load the cell and network objects for the dataset
    cell_population = pickle.load(open(f'{file_paths["pickle_files"]}/{dataset_name}_pickle_files/cells_pickle_file/{dataset_name}.cells.pickle', "rb"))
    cells = cell_population.cells

    # Load the network pickle files
    all_networks = []
    logging.info(f'\nRunning attractor analysis for all networks...')
    pickle_file_path = f'{file_paths["pickle_files"]}/{dataset_name}_pickle_files/network_pickle_files/'
    for pickle_file in glob.glob(pickle_file_path + str(dataset_name) + "_" + "*" + ".network.pickle"):
        if pickle_file:
            logging.info(f'\t\tLoading data file: {pickle_file}')
            network = pickle.load(open(pickle_file, "rb"))
            all_networks.append(network)
        else:
            assert FileNotFoundError("Network pickle file not found")
    
    # Catch exceptions where no networks are loaded
    try:
        successfully_loaded = all_networks[0].name
        logging.info(f'\nFound dataset!')
    except IndexError:
        error_message = "No networks loaded, check to make sure network pickle files exist in 'pickle_files/'"
        logging.error(error_message)
        raise Exception(error_message)
    
    # Run the pathway analysis for each of the networks
    for network in all_networks:

        logging.info(f'\n----- ATTRACTOR ANALYSIS -----')

        num_cells_per_chunk = 25
        num_cells_to_analyze = 100

        # Convert the network's sparse dataset to a dense one
        dataset = network.dataset
        dense_dataset = np.array(dataset.todense())
        network_name = network.name

        # Specify outfile path for the simulation results
        outfile_dir = f'{file_paths["trajectories"]}/{dataset_name}_{network_name}'
        png_dir = f'{outfile_dir}/png_files'
        text_dir = f'{outfile_dir}/text_files'

        # Make all of the outfile folders
        os.makedirs(outfile_dir, exist_ok=True)
        os.makedirs(png_dir, exist_ok=True)
        os.makedirs(text_dir, exist_ok=True)

        # Finds the number of trajectory files
        num_existing_files = len([file for file in os.listdir(text_dir) if file.endswith('_trajectory.csv')])
        logging.info(f'Found {num_existing_files} trajectory files')

        # Finds the number of cells to simulate based on the number of existing trajectory files
        num_simulations = num_cells_to_analyze - num_existing_files

        # Simulates cells to create the cell trajectory files
        simulate_cells(dense_dataset, num_simulations)

        # Recalculates the number of trajectory files after simulating to ensure there are enough
        num_files = len([file for file in os.listdir(text_dir) if file.endswith('_trajectory.csv')])
        logging.info(f'{num_files} trajectory files after simulation')
        
        # Calculate the dynamic time warping
        logging.info(f'Calculating cell trajectory clusters and averaging the cluster trajectories')
        calculate_dtw(num_files, text_dir, num_cells_per_chunk)