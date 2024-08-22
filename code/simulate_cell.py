import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import seaborn as sns
import networkx as nx
import numpy as np
import pickle
import logging
import numexpr as ne
import random
from argparse import ArgumentParser

from user_input_prompts import *
from file_paths import file_paths

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

def save_attractor_simulation(filename, network, simulated_attractor):
    # Save the attractor simulation to a file
    with open(filename, 'w') as file:
        simulated_attractor = np.array(simulated_attractor).T
        for gene_num, expression in enumerate(simulated_attractor):
            file.write(f'{network.nodes[gene_num].name},{",".join([str(i) for i in list(expression)])}\n')

if __name__ == '__main__':

    logging.basicConfig(format='%(message)s', level=logging.INFO)

    parser = ArgumentParser()

    parser.add_argument(
        '--dataset_name',
        type=str,
        required=True,
        help='Name of the dataset to simulate'
    )

    parser.add_argument(
        '--network_name',
        type=str,
        required=True,
        help='Network name to simulate (e.g. hsa04670)'
    )

    parser.add_argument(
        '--num_cells',
        type=str,
        required=True,
        help='number of cells to simulate'
    )

    results = parser.parse_args()

    dataset_name = getattr(results, 'dataset_name')
    network_name = getattr(results, 'network_name')
    num_cells = getattr(results, 'num_cells')

    # Specifies the path to the correct network pickle file
    network_pickle_file = f'{file_paths["pickle_files"]}/{dataset_name}_pickle_files/network_pickle_files/{dataset_name}_{network_name}.network.pickle'

    # Read in the network pickle file from the path
    network = pickle.load(open(network_pickle_file, 'rb'))

    # Convert the network's sparse dataset to a dense one
    dataset = network.dataset
    dense_dataset = np.array(dataset.todense())

    num_simulations = int(num_cells)
    with alive_bar(num_simulations) as bar:
        for i in range(num_simulations):
            # Select a random column from the network dataset
            cell_index = np.random.choice(dense_dataset.shape[1])

            # Reads in all of the rows for that columns
            selected_column = np.array([random.choice([0,1]) for _ in dense_dataset[:, cell_index]])

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