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
import argparse
import pickle
import numexpr as ne
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from fastdtw import fastdtw
import os
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from alive_progress import alive_bar
import statistics
import multiprocessing as mp
from itertools import combinations
import time

from user_input_prompts import attractor_analysis_arguments
from file_paths import file_paths

def vectorized_run_simulation(nodes: list, cell_column: list):
    steps = 20

    # Convert cell_column to NumPy array for faster processing
    starting_state = np.array(cell_column)

    # Preallocate simulation state matrix (steps, nodes)
    total_simulation_states = np.zeros((steps, len(nodes)), dtype=int)

    def evaluate_expression(data, expression):
        expression = expression.replace('and', '&').replace('or', '|').replace('not', '~')
        if any(op in expression for op in ['&', '|', '~']):
            local_vars = {key: np.array(value).astype(bool) for key, value in data.items()}
        else:
            local_vars = {key: np.array(value) for key, value in data.items()}
        return ne.evaluate(expression, local_dict=local_vars)

    # Run the simulation
    for step in range(steps):
        step_expression = []

        # Iterate through each node in the network
        for node_idx, node in enumerate(nodes):
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
                prev_state = total_simulation_states[step - 1]
                if len(incoming_node_indices) > 0:
                    data['A'] = prev_state[incoming_node_indices[0]]
                if len(incoming_node_indices) > 1:
                    data['B'] = prev_state[incoming_node_indices[1]]
                if len(incoming_node_indices) > 2:
                    data['C'] = prev_state[incoming_node_indices[2]]
                if len(incoming_node_indices) > 3:
                    data['D'] = prev_state[incoming_node_indices[3]]

            next_step_node_expression = evaluate_expression(data, node.calculation_function)

            # Save the expression for the node for this step
            step_expression.append(next_step_node_expression)

        # Convert step expression to integer and save the result
        total_simulation_states[step] = np.array(step_expression).astype(int)

    # Return the final matrix of simulation steps where rows are steps and columns are nodes
    return total_simulation_states.tolist()


def simulate_single_cell(cell_index: int, dataset_array: np.ndarray, network: object, dataset_name: str, network_name: str, lock):
    """
    Simulates a single cell trajectory.

    Parameters
    ---------
    cell_index : int
        Index of the cell to simulate
    dataset_array : np.ndarray
        Dense dataset for the cells
    network : object
        The network object containing the nodes
    dataset_name : str
        The name of the dataset
    network_name : str
        The name of the network
    """
    # Reads in the network gene expression for the chosen cell column
    cell_starting_state = np.array([int(gene_expr) for gene_expr in dataset_array[:, cell_index]])

    # Simulate the network using the expression in the selected column as the starting state
    trajectory = vectorized_run_simulation(network.nodes, cell_starting_state)

    # Save the attractor simulation to a csv file
    attractor_sim_path = f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/text_files/cell_trajectories/cell_{cell_index}_trajectory.csv'

    # Makes sure that the cell is not being simulated by multiple threads at once
    with lock:
        if os.path.exists(attractor_sim_path):
            logging.warning(f"File {attractor_sim_path} already exists, skipping write for cell {cell_index}.")
        else:
            with open(attractor_sim_path, 'w') as file:
                trajectory = np.array(trajectory).T
                for gene_num, expression in enumerate(trajectory):
                    file.write(f'{network.nodes[gene_num].name},{",".join([str(i) for i in list(expression)])}\n')

    # Only saves 100 cell simulation trajectory graphs for each network to save space
    if len(os.listdir(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/png_files/cell_trajectories')) <= 100:
        # Creates a heatmap to visualize the trajectories
        heatmap = create_heatmap(
            f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/text_files/cell_trajectories/cell_{cell_index}_trajectory.csv',
            f'Simulation for {dataset_name} {network_name} cell {cell_index} pathway ')

        # Saves a png of the results
        heatmap.savefig(
            f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/png_files/cell_trajectories/cell_{cell_index}_trajectory.png',
            format='png')
        plt.close(heatmap)


def simulate_single_cell_wrapper(args):
    """
    Wrapper to call simulate_single_cell and return 1 for progress tracking.
    This function needs to be in the global scope to be pickled by multiprocessing.
    """
    cell_index, dataset_array, network, dataset_name, network_name, lock = args
    simulate_single_cell(cell_index, dataset_array, network, dataset_name, network_name, lock)
    return 1  # Return 1 to indicate progress for the progress bar


def simulate_cells(dataset_array: np.ndarray, num_simulations: int, network: object, dataset_name: str, network_name: str):
    """
    Simulates the cell trajectories using a network model with parallel processing.

    Parameters
    ---------
    dataset_array : np.ndarray
        The dense dataset for the cells
    num_simulations : int
        The number of simulation steps to run
    network : object
        The network object containing the nodes
    dataset_name : str
        The name of the dataset
    network_name : str
        The name of the network
    """
    logging.info(f'\tSimulating {num_simulations} cell trajectories')
    # Get the list of existing trajectory files
    trajectory_dir = f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/text_files/cell_trajectories'
    existing_files = os.listdir(trajectory_dir)
    existing_indices = set()

    # Find the list of already simulated cell indices from existing files
    for filename in existing_files:
        if filename.startswith("cell_") and filename.endswith("_trajectory.csv"):
            try:
                index = int(filename.split('_')[1])
                existing_indices.add(index)
            except (IndexError, ValueError):
                logging.warning(f"Skipping file {filename}: could not extract cell index.")

    # Get all possible cell indices
    all_indices = set(range(dataset_array.shape[1]))

    # Remove existing indices from the pool of indices to simulate
    indices_to_simulate = list(all_indices - existing_indices)

    # If num_simulations exceeds remaining indices, adjust
    if num_simulations > len(indices_to_simulate):
        logging.warning(
            f"Requested {num_simulations} simulations, but only {len(indices_to_simulate)} unsimulated cells available.")
        num_simulations = len(indices_to_simulate)

    # Randomly select `num_simulations` indices to simulate
    cells_to_simulate = np.random.choice(indices_to_simulate, num_simulations, replace=False)

    # Prepare multiprocessing pool
    pool = mp.Pool(mp.cpu_count())
    manager = mp.Manager()
    lock = manager.Lock()

    # Prepare arguments for each simulation
    arguments = [(cell_index, dataset_array, network, dataset_name, network_name, lock) for cell_index in cells_to_simulate]

    # Use parallel processing and update progress bar
    with alive_bar(num_simulations) as bar:
        # Run simulations in parallel using starmap (blocking version)
        results = pool.imap_unordered(simulate_single_cell_wrapper, arguments)

        for result in results:
            bar()  # Update progress bar for each completed task

        # Close and join the pool after all tasks have been completed
        pool.close()
        pool.join()


def create_heatmap(trajectory_path: str, title: str):
    """
    Creates a trajectory figure using an sns heatmap
    """
    data = []
    gene_names = []

    with open(trajectory_path, 'r') as trajectory_file:
        for line in trajectory_file:
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


def compute_dtw_distance_pair(cell1: str, cell2: str, cell_trajectory_dict: dict):
    """
    Computes the dynamic time warping distance between two cell trajectories for each gene.

    Parameters
    ----------
    cell1: str
        The name of the first cell to be compared
    cell2: str
        The name of the second cell to be compared
    cell_trajectory_dict: dict
        A dictionary with cell names as keys and a dictionary of gene names with trajectories as values

    Returns
    -------
    cell1, cell2 : Str, Str
        Cell names
    total_distance : int | float
        The summed DTW distances between the trajectory of each gene
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

def compute_dtw_distances(cell_trajectory_dict: dict):
    """
    Handles parallel processing of the dynamic time warping calculations.

    Parameters
    ----------
    cell_trajectory_dict: dict
        A dictionary containing the cell names as keys and a dictionary of genes paired to trajectories
        as the values

    Returns
    -------
    dtw_distances: dict
        A dictionary with cell pair tuples as keys and the DTW distance as values
    """
    dtw_distances = {}
    cell_names = list(cell_trajectory_dict.keys())
    cell_pairs = list(combinations(cell_names, 2))  # Get all unique pairs of cells

    # Use multiprocessing to process pairs of cell trajectories
    with mp.Pool(mp.cpu_count()) as pool:
        # Compute DTW distances for all cell pairs in parallel
        tasks = [(cell1, cell2, cell_trajectory_dict) for cell1, cell2 in cell_pairs]
        results = pool.starmap(compute_dtw_distance_pair, tasks)

        # Collect results and update progress bar
        for cell1, cell2, total_distance in results:
            dtw_distances[(cell1, cell2)] = total_distance

    return dtw_distances


def create_distance_matrix(dtw_distances: dict, file_names: list):
    """
    Creates a pairwise distance matrix from the cell-cell distance files.

    Parameters
    ----------
    dtw_distances: dict
        A dictionary with cell pair tuples as keys and the DTW distance as values
    file_names: list
        A list of file names corresponding to the cell pair tuples

    Returns
    ---------
    distance_matrix: np.array
        A numpy array of distances between each cell pair
    """
    # Create a distance matrix
    distance_matrix = np.zeros((len(file_names), len(file_names)))
    for (file1, file2), total_distance in dtw_distances.items():
        i = file_names.index(file1)
        j = file_names.index(file2)
        distance_matrix[i, j] = total_distance
        distance_matrix[j, i] = total_distance

    return distance_matrix


def hierarchical_clustering(dtw_distances: dict, num_clusters: int):
    """ Performs hierarchical clustering on the pairwise DTW distance matrix.

    Parameters
    ----------
    dtw_distances : dict
        A dictionary where each key is a tuple representing a pair of cells,
        and each value is the DTW distance between those cells.

    num_clusters : int
        The number of clusters to split the cells into.

    Returns
    -------
    cluster_labels : list
        A list of cluster labels corresponding to the cells.
    """

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

    # Perform hierarchical clustering
    distance_array = distance_matrix.values[np.triu_indices_from(distance_matrix, k=1)]
    distance_array = np.sqrt(distance_array) # Added sqrt scaling to better visualize clusters
    Z = linkage(distance_array, method='ward')

    # Plot the dendrogram
    plt.figure(figsize=(8, 10))
    dendrogram(Z, labels=distance_matrix.index, orientation='top')
    plt.title('Hierarchical Clustering Dendrogram', fontsize=12)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=8, rotation=90)
    plt.xlabel('Chunk:Cluster', fontsize=8)
    plt.ylabel('Distance', fontsize=8)
    plt.tight_layout()

    plt.savefig(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/png_files/dendrogram.png')

    # Set a threshold and get the clusters
    if num_clusters == 0:
        logging.info(f'Please open "scBONITA_output/trajectories/{dataset_name}_{network_name}/png_files/dendrogram.png"')
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


def summarize_clusters(directory: str, cell_names: list):
    """
    Opens each trajectory files for all cells in the cluster and finds the average gene state at each simulation step.

    Parameters
    ----------
    directory: str
        The path to the trajectory files.
    cell_names: list
        A list of the cell names in the cluster.

    Returns
    -------
    avg_expr_df: pd.DataFrame
        A dataframe containing the average expression at each simulation step for the cells in the cluster.

    """
    gene_expr_dict = {}

    # Finds the path to all trajectory csv files
    trajectory_files = []
    for filename in os.listdir(directory):
        if filename.endswith("_trajectory.csv"):
            trajectory_files.append(os.path.join(directory, filename))

    # finds the files of the cells in the cluster
    files_to_open = []
    for cell in cell_names:
        for file_name in trajectory_files:

            # Finds the exact cell trajectories in the trajectory files
            if cell.split("_")[1] == file_name.split("_")[-2]:
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

    # Finds the average gene expression for the cluster for each gene at each time point
    gene_avg_expr = {}
    for gene, simulation_results in gene_expr_dict.items():

        if gene not in gene_avg_expr:
            gene_avg_expr[gene] = []

        transposed_data = list(map(list, zip(*simulation_results)))

        for i in transposed_data:
            gene_avg_expr[gene].append(statistics.mean(i))

    # Convert the average gene expression dictionary to a DataFrame
    avg_expr_df = pd.DataFrame(gene_avg_expr)

    # Transpose the DataFrame
    avg_expr_df = avg_expr_df.transpose()

    return avg_expr_df


def plot_average_trajectory(df: pd.DataFrame, title: str, path: str):
    """
    Creates an average trajectory figure for a cluster.

    Creates a trajectory simulation figure where the gene expression darkness at each step
    corresponds to how many of the cells in that cluster were expressing the gene at
    that step.

    Parameters
    ----------
    df: pd.DataFrame
        A dataframe containing the average expression at each simulation step.
    title: str
        The title of the figure.
    path: str
        The path specifying where to save the figure.
    """
    # Create the heatmap
    plt.figure(figsize=(8, 10))
    sns.heatmap(df, cmap='Greys', yticklabels=True, vmin=0, vmax=1)
    plt.title(title, fontsize=12)
    plt.xlabel(xlabel='Simulation Time Steps', fontsize=12)
    plt.ylabel(ylabel='Gene', fontsize=12)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def create_trajectory_chunks(num_chunks: int, num_clusters: int, output_directory: str):
    """
    Sorts chunks of cells into clusters based on their pairwise DTW distances.

    This process reads in the individual cell trajectory files and calculates the distance between each cell. The cells
    are clustered based on their distances using hierarchical clustering. Calculating the pairwise distance between each
    cell in the dataset scales exponentially and wouldn't work for large datasets. Instead, the cells are processed into
    clusters in chunks, then the clusters from these chunks are compared like they were individual cells. Basically, the
    clustering process is done twice, once to summarize groups of cells and again to summarize the group chunks into
    n number of clusters. This reduces the number of pairwise comparisons needed for hierarchical clustering while still
    accounting for each cell.

    Parameters
    ----------
    num_chunks : int
        The number of chunks of cells to be processed.
    num_clusters : int
        The number of clusters to be processed.
    output_directory : str
        The directory to save the pairwise distance file to.

    Returns
    ----------
    cluster_chunks: dict
        A dictionary of each chunk in each cluster with the cells in this group
    cells_in_chunks: dict
        A dictionary of which cells are in each chunk
    num_clusters: int
        The number of clusters the cells are split into
    """

    trajectory_files_parsed = []
    cell_txt_traj_dir = f'{output_directory}/cell_trajectories'

    # Chunk the data down to create smaller averaged trajectory clusters
    cluster_chunks = {}
    cells_in_chunks = {}
    logging.info(f'Creating chunks, this may take a while depending on how many cells and the chunk size')
    for chunk in range(num_chunks):
        start = time.time()
        logging.info(f'\tCreating chunk {chunk+1} / {num_chunks}')
        # Reads in the cell trajectories from the simulations and creates a dataframe from them
        cell_trajectory_dict = {}

        num_cells_parsed = 0
        for traj_filename in os.listdir(cell_txt_traj_dir):
            
            # Keep track of how many files are parsed in this batch and which cells trajectories were analyzed
            if num_cells_parsed <= num_cells_per_chunk-1 and traj_filename.endswith("_trajectory.csv") and traj_filename not in trajectory_files_parsed:
                trajectory_files_parsed.append(traj_filename)
                
                # Extract the cell number from the traj_filename
                filepath = os.path.join(cell_txt_traj_dir, traj_filename)
                df = pd.read_csv(filepath, header=None)
                df.columns = ['Gene'] + [f'Time{i}' for i in range(1, df.shape[1])]
                df.set_index('Gene', inplace=True)
                
                # Extract the gene values from the dataset and append them to a dictionary with the gene name as key
                cell_trajectory_dict[traj_filename] = {gene: df.loc[gene].values for gene in df.index}
                num_cells_parsed += 1

        # Computes the pairwise DTW distances between cells
        dtw_distances = compute_dtw_distances(cell_trajectory_dict, cell_txt_traj_dir)

        # Write out the distances or append them if its a different cluster
        if chunk == 0:
            # Write out each cell-cell distance to an outfile
            with open(f'{output_directory}/distances.csv', 'w') as outfile:
                for (cell1, cell2), total_distance in dtw_distances.items():
                    outfile.write(f'{cell1},{cell2},{total_distance}\n')
        else:
            with open(f'{output_directory}/distances.csv', 'a') as outfile:
                for (cell1, cell2), total_distance in dtw_distances.items():
                    outfile.write(f'{cell1},{cell2},{total_distance}\n')


        # Performs hierarchical clustering to find the number of clusters
        cluster_dict, num_clusters = hierarchical_clustering(dtw_distances, num_clusters)

        # Summarize the clusters and track which cells are in which cluster and chunk
        for cluster, cell_list in cluster_dict.items():
            df = summarize_clusters(cell_txt_traj_dir, cell_list)

            # Binarize the DataFrame
            df_binarized = pd.DataFrame(np.where(df < 0.5, 0, 1), index=df.index, columns=df.columns)

            # Save the cluster chunk to a dictionary
            cluster_chunks[f'{chunk}:{cluster}'] = df_binarized

            # Keep track of which cells are in which chunks for later
            if chunk not in cells_in_chunks:
                cells_in_chunks[chunk] = {}

            cells_in_chunks[chunk][cluster] = cell_list

            # # Plots the average trajectory graphs for each cluster in each chunk
            # os.makedirs(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/raw_chunks', exist_ok=True)
            #
            # title: str = f'Average Gene Expression Heatmap for Cluster {cluster}'
            # path: str = f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/raw_chunks/chunk_{chunk}_cluster{cluster}_summary'
            #
            # plot_average_trajectory(df_binarized, title, path)

        end = time.time()
        length = end - start
        if chunk > 0:
            print(
                f'\t\tCompleted in {round(length, 2)} seconds (Est remaining: {round(length * (num_chunks - chunk))}s)')



    return cluster_chunks, cells_in_chunks, num_clusters


def calculate_dtw(num_files: int, output_directory: str, num_cells_per_chunk: int):
    trajectory_files_parsed = []
    num_clusters: int = 0
    num_chunks: int = round(num_files / (num_cells_per_chunk))

    logging.info(f'Creating {num_chunks} chunks ({num_files} cells / {num_cells_per_chunk} cells per chunk)')

    # Ensures there is always one chunk
    if num_chunks == 0:
        num_chunks = 1
    group_dict = {}
    cluster_chunks, cells_in_chunks, num_clusters = create_trajectory_chunks(num_chunks, num_clusters, output_directory)

    logging.info(f'\nComparing chunks')
    chunk_dtw_distances = compute_dtw_distances(cluster_chunks, output_directory)

    # Calculates pairwise distance matrix between the cells
    chunk_names: list = list(cluster_chunks.keys())
    distance_matrix: np.ndarray = create_distance_matrix(chunk_dtw_distances, chunk_names)

    # Performs hierarchical clustering to find the number of clusters
    group_cluster_dict, num_clusters = hierarchical_clustering(chunk_dtw_distances, num_clusters)

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
    
    # Plot the distance heatmap
    plt.figure(figsize=(8, 9))
    sns.heatmap(
        data=reordered_matrix,
        xticklabels=reordered_file_names,
        yticklabels=reordered_file_names,
        cmap='Greys',
        annot=False)
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.title("DTW Distance Heatmap")
    plt.tight_layout()
    plt.savefig(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/png_files/distance_heatmap')
    # plt.show()
    plt.close()

    # Tracks the identity of the cells and number of cells in each cluster and chunk so they can be compiled later
    cells_in_cluster: dict = {}
    for cluster, cell_list in group_cluster_dict.items():
        num_cells_in_cluster: int = 0

        # Extracts the chunks and clusters belonging to the
        for chunk_cluster in cell_list:
            chunk: str = chunk_cluster.split(':')[0]
            cluster: str = chunk_cluster.split(':')[1]

            # Count the number of cells in each cluster
            num_cells_in_cluster += len(cells_in_chunks[int(chunk)][int(cluster)])

            if cluster not in cells_in_cluster:
                cells_in_cluster[cluster]: list = []
            
            cells_in_cluster[cluster].extend(cells_in_chunks[int(chunk)][int(cluster)])

    # Summarize each cluster and plot the average trajectory to compare clusters

    for cluster, cluster_group in group_cluster_dict.items():
        logging.info(f'Summarizing cluster {cluster}')

        # Finds the average signaling state of each gene at each simulation step
        df: pd.DataFrame = summarize_clusters(f'{output_directory}/cell_trajectories', cells_in_cluster[str(cluster)])

        # Create average trajectory directory
        os.makedirs(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/png_files/cluster_summaries', exist_ok=True)
        os.makedirs(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/text_files/cluster_summaries',
                    exist_ok=True)

        title: str = f'Average Gene Expression Heatmap for Cluster {cluster}'
        path: str = f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/png_files/cluster_summaries/cluster_{cluster}_summary'

        # Save the dataframe
        df.to_csv(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/text_files/cluster_summaries/cluster_{cluster}_summary.csv', header=False)

        # Write out df
        plot_average_trajectory(df, title, path)

        # Load the different experimental group pickle files
        groups = []
        pickle_file_path = f'{file_paths["pickle_files"]}/{dataset_name}_pickle_files/network_pickle_files/{dataset_name}_*_pickle_files/'
        for path in glob.glob(pickle_file_path):
            for file in os.listdir(path):
                groups.append(pickle.load(open(f'{path}{file}', 'rb')))

        # Extract the cell numbers for the cells in the cluster
        cell_nums = [str(cell.split('_')[1]) for cell in cells_in_cluster[str(cluster)]]

        # Add the cluster to the group_dict
        if cluster not in group_dict:
            group_dict[cluster] = {}

        # Find which cells are in which groups and add them to the dictionary
        for group_network in groups:
            for cell in group_network.cells:
                cell_str = str(cell).strip()
                if cell_str in cell_nums:
                    if not group_network.name.split('_')[1] in group_dict[cluster]:
                        group_dict[cluster][group_network.name.split('_')[1]] = 0

                    group_dict[cluster][group_network.name.split('_')[1]] += 1

    # Print which cells are in which groups
    with open(f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/text_files/group_data.csv', 'w') as file:
        file.write(f'cluster,group,num_cells\n')
        for cluster, groups in group_dict.items():
            logging.info(f'\nCluster {cluster}')
            for group, num_cells in groups.items():
                logging.info(f'\t{group}: {num_cells} cells')
                file.write(f'{cluster},{group},{num_cells}\n')




# If you want to run the attractor analysis by itself
if __name__ == '__main__':

    # Set the logging level for output
    logging.basicConfig(format='%(message)s', level=logging.INFO)

    # Allow the user to either add in the dataset name and network name from the command line or as a prompt
    parser: argparse.ArgumentParser = argparse.ArgumentParser()

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
            logging.info(f'\tLoading data file: ...{pickle_file[-50:]}')
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

        num_cells_per_chunk: int = 75
        num_cells_to_analyze: int = 10000

        # Convert the network's sparse dataset to a dense one
        dataset = network.dataset
        dense_dataset: np.ndarray = np.array(dataset.todense())
        num_cells_in_dataset = dense_dataset.shape[1]
        logging.info(f'There are {num_cells_in_dataset} cells in dataset')
        network_name: str = network.name

        # Specify outfile path for the simulation results
        outfile_dir = f'{file_paths["trajectories"]}/{dataset_name}_{network_name}'
        png_dir = f'{outfile_dir}/png_files'
        text_dir = f'{outfile_dir}/text_files/'
        txt_traj_dir = f'{text_dir}/cell_trajectories/'
        png_traj_dir = f'{png_dir}/cell_trajectories/'

        # Make all outfile folders
        os.makedirs(txt_traj_dir, exist_ok=True)
        os.makedirs(png_traj_dir, exist_ok=True)

        # Find the number of trajectory files
        num_existing_files: int = len([file for file in os.listdir(txt_traj_dir) if file.endswith('_trajectory.csv')])
        logging.info(f'Found {num_existing_files} trajectory files')

        # Finds the number of cells to simulate based on the number of existing trajectory files
        num_simulations: int = min((num_cells_to_analyze - num_existing_files), num_cells_in_dataset - num_existing_files)
        if num_simulations < 0:
            num_simulations = 0

        # Simulates cells to create the cell trajectory files
        simulate_cells(dense_dataset, num_simulations, network, dataset_name, network_name)

        # Recalculates the number of trajectory files after simulating to ensure there are enough
        num_files: int = len([file for file in os.listdir(txt_traj_dir) if file.endswith('_trajectory.csv')])

        logging.info(f'{num_files} trajectory files after simulation')

        num_files_to_process = min(num_cells_to_analyze, num_files)
        
        # Calculate the dynamic time warping
        logging.info(f'Calculating cell trajectory clusters and averaging the cluster trajectories')
        calculate_dtw(num_files_to_process, text_dir, num_cells_per_chunk)

