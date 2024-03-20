import pickle
import logging
import os
import argparse
import glob

from setup.user_input_prompts import cell_attractor_mapping_arguments

def load_all_network_pickles(dataset_name):
    """
    Loads in all network pickle files for a given dataset_name
    """
    # Load in the network pickle files
    all_networks = []

    pickle_file_path = f'pickle_files/{dataset_name}_pickle_files/network_pickle_files/'
    for pickle_file in glob.glob(pickle_file_path + str(dataset_name) + "_" + "*" + ".network.pickle"):
        if pickle_file:
            logging.info(f'\tLoading data file: {pickle_file}')
            network = pickle.load(open(pickle_file, "rb"))
            all_networks.append(network)
        else:
            assert FileNotFoundError("Network pickle file not found")

    return all_networks

def create_cell_attractor_dictionary(all_networks):
    """
    Creates a dictionary containing each cell and the attractor that it mapped to. The attractor that the cells map to for
    each network is stored in each network object, so this function iterates through each network object and retrieves the 
    attractor that the cell mapped to.
    """
    # Count the total number of cells in all networks
    num_cells = all_networks[0].dataset.shape[1]
    logging.info(f'\nTotal number of cells: {num_cells}')

    # Create a dictionary of cells
    mapped_cells = {}
    for cell in range(num_cells):
        mapped_cells[cell] = []
        for network in all_networks:
            mapped_cells[cell].append(network.cell_map[cell])
    
    return mapped_cells, num_cells

def map_all_cell_attractors(mapped_cells):
    """
    For each cell, this iterates through each of the networks and makes a tuple containing the network name and the attractor that the cell
    mapped to. Returns a list of all cells and the network attractors they mapped to.
    """
    # Make a list of tuples containing which attractor the cell mapped to for each network (network, attractor number)
    cell_states = []
    for cell in mapped_cells:
        states = []
        for network_num, _ in enumerate(mapped_cells[cell]):
            network = all_networks[network_num].name
            attractor_num = mapped_cells[cell][network_num]
            states.append((network, attractor_num))

        cell_states.append(states)
    
    return cell_states

def find_unique_attractor_combinations(cell_states):
    """
    For each network, a cell will map to an attractor. This finds each unique combination of attractors across all networks and counts
    how many cells mapped to that combination.
    """
    # Get the number of unique cluster combinations and count the number of cells belonging to that group
    attractor_combinations = {}
    for cell in cell_states:
        # Sort and convert to tuple to ensure uniqueness regardless of order
        combination = tuple(sorted(cell))
        
        # Count the number of times the combination is seen and append that to a dictionary of the combinations
        if combination in attractor_combinations:
            attractor_combinations[combination] += 1
        else:
            attractor_combinations[combination] = 1
    
    # Sort the combinations from highest to lowest based on the number of cells with that attractor combination
    sorted_combinations = sorted(attractor_combinations.items(), key=lambda item: item[1], reverse=True)
    
    return sorted_combinations

def write_attractor_map(dataset_name, attractor_combinations):
    """
    Saves each attractor combination count to a csv file. Has a header with each network name, the number of cells, and the percent of cells.
    Each row contains the attractor numbers for each network, the number of cells that mapped to that attractor combination, and the percentage of
    the total cells that map to that combination
    """
    # Format and save the combinations to a csv file in the attractors file
    path = f'attractor_analysis_output/{dataset_name}_attractors/{dataset_name}_cell_attractor_states.csv'
    with open(path, 'w') as combination_file:
        pathways = []
        total_attractors = []
        for combinations, _ in attractor_combinations:
            attractors = []
            for pathway, attractor in combinations:
                if pathway not in pathways:
                    pathways.append(pathway)
                attractors.append(attractor)
            total_attractors.append(attractors)
        
        # Header containing pathways
        combination_file.write('\t'.join(pathways) + '\tNumber_of_cells\tPercent_of_cells\n')

        # Print attractors and cell counts
        for pathway_num, pathway in enumerate(total_attractors):
            percent_of_cells = round(attractor_combinations[pathway_num][1] / num_cells * 100, 3)
            combination_file.write('\t'.join(str(i) for i in pathway) + '\t' + str(attractor_combinations[pathway_num][1]) + '\t' + str(percent_of_cells) + '\n')

def write_active_inactive_states(dataset_name, all_networks):
    """
    For each attractor in each network, write out a file containing the node names and whether they are ACTIVE or INACTIVE for that attractor.
    This helps the user to figure out which nodes are active in the network.
    """
    for network in all_networks:
        logging.info(f'Writing out attractor gene activation states for network {network.name}')
        for attractor_num, attractor in network.representative_attractors.items():
            try:
                logging.debug(f'\tattractor_num {attractor_num}')
                path = f'attractor_analysis_output/{dataset_name}_attractors/{network.name}_attractors/attractor_{attractor_num}/attractor_{attractor_num}_node_state.csv'
                with open(path, 'w') as file:
                    header = f'{network.name}   attractor_{attractor_num}\n'
                    file.write(header)
                    for node_index, node_value in enumerate(attractor):
                        node = network.nodes[node_index]
                        if node_value == 1:
                            line = f'{node.name}    ACTIVE\n'
                            file.write(line)
                        elif node_value == 0:
                            line = f'{node.name}    inactive\n'
                            file.write(line)
            except FileNotFoundError:
                pass

def write_only_active_states(dataset_name, all_networks):
    """
    For each attractor in the network, write out a file containing the nodes that are active at the start of the attractor cycle ONLY.
    This can help the user determine which genes are active in larger networks without having to look at every gene.
    """
    for network in all_networks:
        logging.info(f'\nWriting out active genes for network {network.name}')
        for attractor_num, attractor in network.representative_attractors.items():
            try:
                logging.debug(f'\tattractor_num {attractor_num}')
                path = f'attractor_analysis_output/{dataset_name}_attractors/{network.name}_attractors/attractor_{attractor_num}/attractor_{attractor_num}_active_nodes.csv'
                with open(path, 'w') as file:
                    header = f'{network.name}   attractor_{attractor_num}\n'
                    file.write(header)
                    for node_index, node_value in enumerate(attractor):
                        node = network.nodes[node_index]
                        if node_value == 1:
                            line = f'{node.name}    ACTIVE\n'
                            file.write(line)

            except FileNotFoundError:
                pass

if __name__ == '__main__':
    # Set the logging level for output
    logging.basicConfig(format='%(message)s', level=logging.INFO)

    # Allow the user to either add in the dataset name and network name from the command line or as a prompt
    parser = argparse.ArgumentParser()

    dataset_name = cell_attractor_mapping_arguments(parser)

    # Load in the network pickle files for this dataset
    logging.info(f'\nMapping cell states...')
    all_networks = load_all_network_pickles(dataset_name)

    # Create a dictionary of cells
    cells, num_cells = create_cell_attractor_dictionary(all_networks)
    
    # For each cell, make a list of tuples with each network and the attractor the cell maps to for that attractor
    cell_states = map_all_cell_attractors(cells)
    
    # Find each unique combination of attractors that the cells mapped to and count the number of cells that mapped to each combination
    attractor_combinations = find_unique_attractor_combinations(cell_states)    
    
    # Write out a file containing the attractor combinations and the number of cells that mapped
    write_attractor_map(dataset_name, attractor_combinations)

    # Write out files for each attractor in the dataset containing the nodes and whether they are ACTIVE or INACTIVE
    write_active_inactive_states(dataset_name, all_networks)
    
    # Write out a file to each attractor for each network containing only the active genes in the network
    write_only_active_states(dataset_name, all_networks)

            
