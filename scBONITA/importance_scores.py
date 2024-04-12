from boolean_rule_functions import *
import logging
import random
import numpy as np
from alive_progress import alive_bar
from scipy.sparse import csr_matrix, csc_matrix
from network import Network
import os
import pickle
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import csv
from path_manager import PathManager
from setup.user_input_prompts import *

class CalculateImportanceScore():
    def __init__(self, nodes, binarized_matrix):
        self.binarized_matrix = binarized_matrix
        self.nodes = nodes
        self.STEPS = 50
        self.CELLS = 250

        # Ensure you do not sample more columns than exist in the matrix
        n_cols = min(self.CELLS, binarized_matrix.shape[1])

        import random
        # Sample the dataset for cells to simulate
        self.cell_sample_indices = random.sample(range(binarized_matrix.shape[1]), k=n_cols)

        # Check if the matrix is in CSR or CSC format
        if isinstance(binarized_matrix, (csr_matrix, csc_matrix)):
        # Use fancy indexing for sparse matrices
            self.dataset = binarized_matrix[:, self.cell_sample_indices].todense()
        else:
            # Fallback for dense matrices (NumPy arrays)
            self.dataset = binarized_matrix[:, self.cell_sample_indices]
        
        self.zeros_array = np.array([[0] * self.dataset.shape[1]], dtype=np.bool)
        self.ones_array = np.array([[1] * self.dataset.shape[1]], dtype=np.bool)

    def vectorized_run_simulation(self,knockout_node=None, knockin_node=None):
        total_simulation_states = []

        # Run the simulation
        for step in range(self.STEPS):
            step_expression = []

            # Iterate through each node in the network
            for node in self.nodes:
                
                # Find the next step's expression for each cell in the dataset
                if node.name == knockout_node:
                    next_step_node_expression = self.zeros_array

                elif node.name == knockin_node:
                    next_step_node_expression = self.ones_array

                else: 
                    # Get the rows in the dataset for the incoming nodes
                    if step == 0:
                        incoming_node_data = [self.dataset[i] for i in node.incoming_node_indices]
                    else:
                        incoming_node_data = [total_simulation_states[step-1][i] for i in node.incoming_node_indices]

                    # Allow for whole rows of the dataset to be passed into the function rather than one at a time
                    vectorized_calculation_function = np.vectorize(node.calculation_function)

                    # Set up the argument to be passed into the logic funciton, allows for different number of input nodes
                    function_argument = incoming_node_data + node.inversion_rules

                    # Apply the logic function to the data from the incoming nodes to get the predicted output

                    next_step_node_expression = vectorized_calculation_function(*function_argument)

                # Save the expression for the node for this step
                step_expression.append(next_step_node_expression)
            
            # Save the expression 
            total_simulation_states.append(step_expression)

        total_simulation_states = np.squeeze(np.array(total_simulation_states), axis=2)

        attractors = self.calculate_attractors(total_simulation_states)

        return attractors

    def calculate_attractors(self, total_simulation_states):
        # Transpose the matrix
        total_simulation_states = np.transpose(total_simulation_states, (2, 1, 0)) # num_cells, num_nodes, num_steps

        num_cells = total_simulation_states.shape[0]
        num_nodes = total_simulation_states.shape[1]
        num_steps = total_simulation_states.shape[2]

        # For each cell
        attractors = {}
        for cell in range(num_cells):
            # Find the start and end indices of the attractors for the cell
            attractor_start_index, attractor_end_index = self.find_attractors(total_simulation_states, cell, num_steps)
            
            # Append the value of each node within the attractor loop to the dictionary for the cell
            if attractor_start_index is not None and attractor_end_index is not None:
                attractors[cell] = total_simulation_states[cell, :, attractor_start_index:attractor_end_index+1]

        return attractors

    def find_attractors(self, total_simulation_states, cell, num_steps):        
        for i in range(num_steps):
            for j in range(i):
                # Compare the states of each cell across all nodes at two different time steps
                if np.array_equal(total_simulation_states[cell, :, i], total_simulation_states[cell, :, j]):
                    attractor_start_index = j
                    attractor_end_index = i
                    return (attractor_start_index, attractor_end_index)
        return (None, None)

    def perform_knockouts_knockins(self):
        """
        Runs the simulation for each node, knocking in and out each of the nodes sequentially
        """
        normal_signaling = {}
        knockout_results = {}
        knockin_results = {}

        count = 1

        for node in self.nodes:
            incoming_nodes = node.best_rule[1][:]
            logic_function = node.best_rule[2]
            inversion_rules = node.best_rule[3]
            
            node.calculation_function = node.find_calculation_function(logic_function)
            node.incoming_node_indices = [index for index, name in node.predecessors.items() if name in incoming_nodes]
            node.inversion_rules = inversion_rules
            

        # Calculate knockins and knockouts
        with alive_bar(len(self.nodes)) as bar:
            for node in self.nodes:
                # Calculate normal signaling
                normal_signaling[node.name] = self.vectorized_run_simulation(knockout_node=None)

                # Perform knock-out simulation
                knockout_results[node.name] = self.vectorized_run_simulation(knockout_node=node.name) # A list of the attractors calculated for each cell

                # Perform knock-in simulation
                knockin_results[node.name] = self.vectorized_run_simulation(knockin_node=node.name) # A list of the start and stop of each attractor 

                count += 1
                bar()
        
        return normal_signaling, knockout_results, knockin_results

    def calculate_importance_scores(self):
        logging.info(f'\n-----RUNNING NETWORK SIMULATION-----')

        # Perform the knock-out and knock-in simulations for each node
        normal_signaling, knockout_results, knockin_results = self.perform_knockouts_knockins()

        logging.info(f'\n-----CALCULATING IMPORTANCE SCORES-----')
        with alive_bar(len(self.nodes)) as bar:

            raw_importance_scores = {}
            for node in self.nodes:
                total_difference = 0
                
                # For each cell in the sample cell indices
                for cell, _ in enumerate(self.cell_sample_indices):
                    if cell in knockin_results[node.name] and cell in knockout_results[node.name] and cell in normal_signaling[node.name]:
                        # Set the attractors to a NumPy array so they can be aligned
                        normal_attractors = np.array(normal_signaling[node.name][cell])

                        knockin_attractors = np.array(knockin_results[node.name][cell])

                        knockout_attractors = np.array(knockout_results[node.name][cell])

                        # Align the attractors for the knock-in and knock-out conditions based on the shortest list
                        knockin_attractors_aligned, normal_attractors_aligned = self.align_attractors(knockin_attractors, normal_attractors)

                        # Find the sum of the differences between the knock-in and knock-out attractors
                        knockin_difference = np.sum(np.abs(knockin_attractors_aligned ^ normal_attractors_aligned))

                        knockout_attractors_aligned, normal_attractors_aligned = self.align_attractors(knockout_attractors, knockout_attractors)

                        # Find the sum of the differences between the knock-in and knock-out attractors
                        knockout_difference = np.sum(np.abs(knockout_attractors_aligned ^ normal_attractors_aligned))

                        # Debug information
                        if knockin_difference > 0:
                            total_difference += knockin_difference
                        if knockout_difference > 0:
                            total_difference += knockout_difference

                    else:
                        logging.error(f'\nERROR: No knockout attractor found for cell {cell} in node {node.name}. Increase step size\n')
                        
                # Assign the total difference between the knock-in and knock-out to a dictionary with node names as keys
                raw_importance_scores[node.name] = total_difference

                bar()
        
        # Find the node with the maximum importance score for scaling
        max_importance_score = max(raw_importance_scores.values())

        # Scale the importance scores to [0, 1] by dividing each score by the max score
        scaled_importance_scores = {node: score / max_importance_score
                                    for node, score in raw_importance_scores.items()}

        # Update the nodes importance_score attribute
        for node in self.nodes:
            node.importance_score = scaled_importance_scores[node.name]

        return scaled_importance_scores

    def align_attractors(self, ko_attractor, ki_attractor):
        """
        Align the knock-out and knock-in attractor arrays by the minimum size along both axes.
        """
        # Minimum length along the first axis (number of steps/states)
        min_length_steps = min(ko_attractor.shape[0], ki_attractor.shape[0])

        # Minimum length along the second axis (number of features/variables)
        min_length_features = min(ko_attractor.shape[1], ki_attractor.shape[1])

        # Trim both arrays to the minimum length along both axes
        ko_attractor_aligned = ko_attractor[:min_length_steps, :min_length_features]
        ki_attractor_aligned = ki_attractor[:min_length_steps, :min_length_features]

        return ko_attractor_aligned, ki_attractor_aligned

def run_full_importance_score(dataset_name, network_names): 
    # Path to the ruleset pickle file
    ruleset_pickle_file_path = f'pickle_files/{dataset_name}_pickle_files/ruleset_pickle_files/'

    network_name_check = []
    for file_name in os.listdir(ruleset_pickle_file_path):
        network_name = file_name.split(dataset_name+"_")[1].split(".ruleset.pickle")[0]
        network_name_check.append(network_name)
        if network_name in network_names:

            # Check to make sure the ruleset pickle file exists
            logging.info(f'\nLoading: {dataset_name}_{network_name}.ruleset.pickle')

            # Load the ruleset object for the network
            ruleset = pickle.load(open(f'{ruleset_pickle_file_path}/{file_name}', "rb"))

            # Store the network information in a pickle file
            network = Network(name=f'{network_name}')
            network.nodes = ruleset.nodes
            network.rulesets = ruleset.best_ruleset
            network.network = ruleset.network
            network.dataset = ruleset.binarized_matrix

            # Specify the path to the importance score folder
            importance_score_folder = f'{ruleset.dataset_name}_importance_scores/full_dataset_importance_score'
            importance_score_file_name = f'{ruleset.network_name}_{ruleset.dataset_name}_importance_scores'

            # Run the importance score calculation for that ruleset and network
            logging.info(f'Calculating importance score for network {network_name}')
            importance_score(ruleset, network, importance_score_folder, importance_score_file_name)

            logging.info(f'Saving network object as a pickle file')
            network_folder = f'pickle_files/{dataset_name}_pickle_files/network_pickle_files'
            os.makedirs(network_folder, exist_ok=True)

            network_file_path = f'{network_folder}/{dataset_name}_{network_name}.network.pickle'

            pickle.dump(network, open(network_file_path, 'wb'))
            logging.info(f'\tSaved to {network_file_path}')
        else:
            logging.debug(f'Skipping {network_name}')

    # Check to make sure that the network(s) specified have ruleset pickle files
    common_items = [item for item in network_name_check if item in network_names]
    if len(common_items) == 0:
        raise Exception(f'ERROR: Pathways specific do not exist in the ruleset.pickle folder. Check spelling and try again')

def importance_score(ruleset, network, folder_path, file_name):
    importance_score_calculator = CalculateImportanceScore(network.nodes, network.dataset.astype(np.bool))
    
    importance_score_calculator.calculate_importance_scores()

    # Create the paths to the relative abundance results
    content = "\n".join(f"{node.name} = {round(node.importance_score, 3)}" for node in network.nodes)

    # Save the importance scores to a file
    logging.info(f'Saving importance scores to file: {file_name}.txt')
    path_manager.save_text_file(folder_path, file_name, content)

    # Create and save the figure
    fig = ruleset.plot_graph_from_graphml(network.network)

    logging.info(f'Saving importance score figures')
    path_manager.save_png_file(folder_path, file_name, fig)
    path_manager.save_svg_file(folder_path, file_name, fig)

    plt.close(fig)

if __name__ == '__main__':
    # Set the logging level for output
    logging.basicConfig(format='%(message)s', level=logging.INFO)

    path_manager = PathManager('importance_score_output')

    parser = ArgumentParser()

    add_dataset_name_arg(parser)
    add_list_of_kegg_pathways(parser)
    add_organism_code(parser)

    args = parser.parse_args()
   
    dataset_name = check_dataset_name(args.dataset_name)
    network_names = args.list_of_kegg_pathways
    organism_code = args.organism

    for name_index, name in enumerate(network_names):
        org_network_name = organism_code + name
        network_names[name_index] = org_network_name

    txt = f'running importance scores for {dataset_name}'
    logging.info(f' -----{"-" * len(txt)}----- '.center(20))
    logging.info(f'|     {txt.upper()}     |'.center(20))
    logging.info(f' -----{"-" * len(txt)}----- '.center(20))

    run_full_importance_score(dataset_name, network_names)
