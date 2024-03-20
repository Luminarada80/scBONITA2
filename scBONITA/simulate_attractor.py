from importance_scores import CalculateImportanceScore
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
import networkx as nx
import numpy as np
import pickle
from argparse import ArgumentParser
from setup.user_input_prompts import *
import logging

def get_starting_state(file):
    starting_state = []

    with open(file, 'r') as network_attractor:
        for line in network_attractor:
            node_starting_state = int(line)
            starting_state.append(node_starting_state)
    
    return starting_state

def simulate_network(nodes, filename):
    steps = 100

    starting_state = get_starting_state(filename)

    total_simulation_states = vectorized_run_simulation(nodes, starting_state, steps)
    
    simulation_results = [[int(item) for item in sublist] for sublist in total_simulation_states]

    return simulation_results

def vectorized_run_simulation(nodes, starting_state, steps):
    total_simulation_states = []

    # Run the simulation
    for step in range(steps):
        step_expression = []

        # Iterate through each node in the network
        for node in nodes:
            incoming_nodes = node.best_rule[1][:]
            logic_function = node.best_rule[2]
            inversion_rules = node.best_rule[3]

            # Get the indices for the incoming nodes
            incoming_node_indices = [index for index, name in node.predecessors.items() if name in incoming_nodes]

            # Get the rows in the dataset for the incoming nodes
            if step == 0:
                incoming_node_data = [starting_state[i] for i in incoming_node_indices]
            else:
                incoming_node_data = [total_simulation_states[step-1][i] for i in incoming_node_indices]

            # Find the logic function that should be used for this predicted rule
            calculation_function = node.find_calculation_function(logic_function)

            # Allow for whole rows of the dataset to be passed into the function rather than one at a time
            vectorized_calculation_function = np.vectorize(calculation_function)

            # Set up the argument to be passed into the logic funciton, allows for different number of input nodes
            function_argument = incoming_node_data + inversion_rules

            # Apply the logic function to the data from the incoming nodes to get the predicted output
            next_step_node_expression = vectorized_calculation_function(*function_argument)

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
    
    if show_simulation == "True" or show_simulation == "1":
        ani = animation.FuncAnimation(fig, update_graph, frames=range(1, len(time_series_data)), interval=200)
        plt.show()
    
    else:
        update_graph(10)  # Update to the frame you want to save
        return fig

def save_attractor_simulation(filename, network, simulated_attractor):
    # Save the attractor simulation to a file
    with open(filename, 'w') as file:
        simulated_attractor = np.array(simulated_attractor).T
        for gene_num, expression in enumerate(simulated_attractor):
            file.write(f'{network.nodes[gene_num].name},{",".join([str(i) for i in list(expression)])}\n')

if __name__ == '__main__':
    logging.basicConfig(format='%(message)s', level=logging.INFO)

    # Add in arguments to find the attractor file
    parser = ArgumentParser()

    add_dataset_name_arg(parser)
    add_network_name(parser)

    args = parser.parse_args()

    dataset_name = check_dataset_name(args.dataset_name)
    network_name = args.network_name
    attractor_num = input("Enter attractor number: ")

    # Specify file paths
    network_pickle_file_path = f'pickle_files/{dataset_name}_pickle_files/network_pickle_files/{dataset_name}_{network_name}.network.pickle'
    attractor_path = f'attractor_analysis_output/{dataset_name}_attractors/{network_name}_attractors/attractor_{attractor_num}/{dataset_name}_{network_name}_attractor_{attractor_num}.txt'
    outfile = f'attractor_analysis_output/{dataset_name}_attractors/{network_name}_attractors/attractor_{attractor_num}/{dataset_name}_{network_name}_simulated_attractor_{attractor_num}.txt'

    # Load the network pickle file
    network = pickle.load(open(network_pickle_file_path, "rb"))
    
    # Simulate the network
    simulated_attractor = simulate_network(network.nodes, attractor_path)

    # Visualize the network
    visualize_simulation(network.network, simulated_attractor, network, "True")

    # Save the attractor states to a csv file
    logging.info(f'Saved file to: "{outfile}"')
    save_attractor_simulation(outfile, network, simulated_attractor)





    