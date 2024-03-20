import pickle
import logging
import os
import argparse
import glob
from itertools import product

from setup.user_input_prompts import cell_attractor_mapping_arguments

def load_all_network_pickles(dataset_name: str):
    """
    Loads in all network pickle files for a given dataset_name
    """
    # Load in the network pickle files
    all_networks = []

    pickle_file_path = f'pickle_files/{dataset_name}_pickle_files/network_pickle_files/'
    for pickle_file in glob.glob(pickle_file_path + str(dataset_name) + "_" + "*" + ".network.pickle"):
        if pickle_file:
            logging.info(f'\tLoading data file: {pickle_file}')
            network: object = pickle.load(open(pickle_file, "rb"))
            all_networks.append(network)
        else:
            assert FileNotFoundError("Network pickle file not found")

    return all_networks

def create_cell_index_to_group_dict(network: object):
    """
    Creates a dictionary with each cell index and the group(s) that it maps to. Also
    returns a dictionary with the total number of cells in each group and a list of
    the groups.

    returns cell_group_dict, group_cell_counts, groups
    """
    cell_group_dict = {}
    group_cell_counts = {}
    groups = []
    for group, cell_indices in network.group_cell_indices.items():
        logging.info(f'\tGroup: {group}, {len(cell_indices)} cells')
        for cell_index in cell_indices:
            cell_group_dict[cell_index] = group
        group_cell_counts[group] = len(cell_indices)
        groups.append(group)
    
    return cell_group_dict, group_cell_counts, groups

def create_attractor_group_cell_count_dict(network: object, cell_group_dict: dict):
    """
    Creates an dictionary with each attractor as a key, contains another dictionary
    with each group as the key and the cell counts for that group in the attractor
    as the value

    {attractor_num: {group: count, group: count}}

    returns cluster_counts dictionary
    """
    attractor_counts = {}
    for cell_index, attractor in network.cell_map.items():
        # logging.info(f'Cell: {cell_index}, attractor: {attractor}, Group {cell_group_dict[cell_index]}')
        group = cell_group_dict[cell_index]

        if attractor not in attractor_counts:
            attractor_counts[attractor] = {}

        if group not in attractor_counts[attractor]:
            attractor_counts[attractor][cell_group_dict[cell_index]] = 1
        else:
            attractor_counts[attractor][cell_group_dict[cell_index]] += 1

    return attractor_counts

def map_all_cell_attractors(mapped_cells: dict):
    """
    For each cell, this iterates through each of the networks and makes a tuple containing the network name and the attractor that the cell
    mapped to. Returns a list of all cells and the network attractors they mapped to.
    """
    # Make a list of tuples containing which attractor the cell mapped to for each network (network, attractor number)
    cell_states = []
    for cell in mapped_cells:
        states = []
        for network_num, _ in enumerate(mapped_cells[cell]):
            network_name = all_networks[network_num].name
            attractor_num = mapped_cells[cell][network_num]
            states.append((network_name, attractor_num))

        cell_states.append(states)
    
    return cell_states

def find_unique_attractor_combinations(cell_states: list):
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

def write_attractor_map(dataset_name: str, attractor_combinations: dict):
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

def write_active_inactive_states(dataset_name: str, all_networks: list):
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

def write_only_active_states(dataset_name: str, all_networks: list):
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

    # Find the number of cells in the dataset
    num_cells = all_networks[0].dataset.shape[1]
    logging.info(f'\nTotal number of cells: {num_cells}')

    cell_group_dict = {}

    group_cell_counts = {}

    groups = []

    # Create a dictionary of cell indices to group
    for group, cell_indices in all_networks[0].group_cell_indices.items():
        logging.info(f'\tGroup: {group}, {len(cell_indices)} cells')
        groups.append(group)

        for cell_index in cell_indices:
            cell_group_dict[cell_index] = group
            group_cell_counts[group] = len(cell_indices)

    cell_info = {}
    network_info = {}
    for network in all_networks:

        # Create a dictionary of cell indices to attractor number
        index_to_attractor_dict = {}
        for cell_index in range(num_cells):
            attractor = network.cell_map[cell_index]
            index_to_attractor_dict[cell_index] = attractor

            if cell_index < 5:
                logging.info(f'Cell {cell_index}, Network {network.name}, Attractor {attractor}, Group {cell_group_dict[cell_index]}')
            
            if not cell_index in cell_info:
                cell_info[cell_index] = {}
                
                if not network.name in cell_info[cell_index]:
                    cell_info[cell_index][network.name] = []
            
            cell_info[cell_index][network.name] = [attractor, cell_group_dict[cell_index]]


        attractor_counts = create_attractor_group_cell_count_dict(network, cell_group_dict)

        # Keep attractors that have a significant proportion of the cells (>5%)
        significant_attractors = {}
        for attractor in attractor_counts:
            # Only look at attractors with a significant proportion of the total cells
            attractor_total_counts = sum([attractor_counts[attractor][group] for group in attractor_counts[attractor]])
            if attractor_total_counts / num_cells > 0.05:
                significant_attractors[attractor] = {}

                # Count the number of cells in each group for the attractor
                for group in groups:
                    
                    # If the group is in the attractor, add the cell count
                    if group in attractor_counts[attractor]:       
                        significant_attractors[attractor][group] = attractor_counts[attractor][group]
                        # logging.info(f'\t\t{group}: {attractor_counts[attractor][group]} / {group_cell_counts[group]} cells')
                    
                    # If there are no cells for a group in the attractor, set the cell count to 0
                    else:
                        significant_attractors[attractor][group] = 0
                        # logging.info(f'\t\t{group}: 0 / {group_cell_counts[group]} cells')

        network_info[network.name] = significant_attractors
    


    """
    {combination: {'count' : num_cells, 'group_counts' : {'Healthy': num_healthy, 'HIV' : num_HIV}}}
    
    """

    network_combos = {}

    attractor_combo_strings = []
    for cell_num, cell_info in cell_info.items():

        attractor_combo = []
        
        for network_name, info in cell_info.items():
            attractor_combo.append(info[0])
            group_name = info[1]

        attractor_combo_string = ",".join([str(i) for i in attractor_combo])
        attractor_combo_strings.append(attractor_combo_string)

        if attractor_combo_string not in network_combos:
            network_combos[attractor_combo_string] = {}

        if 'count' not in network_combos[attractor_combo_string]:
            network_combos[attractor_combo_string]['count'] = 0
        
        if 'group_counts' not in network_combos[attractor_combo_string]:
            network_combos[attractor_combo_string]['group_counts'] = {}

        for group in groups:
            if group not in network_combos[attractor_combo_string]['group_counts']:
                network_combos[attractor_combo_string]['group_counts'][group] = 0
            
            if group_name == group:
                network_combos[attractor_combo_string]['group_counts'][group] += 1


        network_combos[attractor_combo_string]['count'] += 1

    print(network_combos)

    for combination in network_combos:
        logging.info(f'Combination: {combination}')
        logging.info(f'\tNumber of cells: {network_combos[combination]["count"]}')

        for group in network_combos[combination]["group_counts"]:
            logging.info(f'\t\t{group}: {network_combos[combination]["group_counts"][group]}')


    # Iterate through the networks and find the list of attractors for each network
    network_attractors = []
    for network_name, significant_attractors in network_info.items():
        attractors_list = []
        logging.info(f'\nNetwork {network_name}')

        for attractor_num, groups in significant_attractors.items():
            attractor_cell_count = sum(groups.values())
            logging.info(f'\tAttractor {attractor_num}: ({attractor_cell_count} total cells)')
            attractors_list.append((network_name, attractor_num, attractor_cell_count))


            for group, num_cells in groups.items():
                logging.info(f'\t\t{group}: {num_cells} cells')
        network_attractors.append(attractors_list)

    # Step 2: Generate all combinations of attractors between networks
    attractor_combinations = list(product(*network_attractors))

    # Step 3: Calculate total cell counts for each combination and sort
    combination_cell_counts = []
    for combo in attractor_combinations:
        total_cells = sum(attractor[2] for attractor in combo)  # Summing the cell counts
        combination_cell_counts.append((combo, total_cells))

    # Sort by total cell count
    sorted_combinations = sorted(combination_cell_counts, key=lambda x: x[1], reverse=True)

    # Print sorted combinations and their total cell counts
    for combo, total_cells in sorted_combinations:
        attractor_info = ', '.join(f'Network {attractor[0]} Attractor {attractor[1]}' for attractor in combo)
        print(f'Combination: [{attractor_info}], Total Cells: {total_cells}')
    
    network_counts = {}

    # for cell_num, cell_info in cell_info.items():
    #     # print(f'Cell {cell_num}')
    #     for network_name, info in cell_info.items():
    #         if not network_name in network_counts:
    #             network_counts[network_name] = {}
    #         # print(f'\tNetwork: {network_name}, Attractor {info[0]}, Group {info[1]}')

    # Format and save the combinations to a csv file in the attractors file
    path = f'attractor_analysis_output/{dataset_name}_attractors/{dataset_name}_cell_attractor_states.csv'
    with open(path, 'w') as combination_file:
        pathways = []
        total_attractors = []

        for combo, total_cells in sorted_combinations:
            attractors = []

            for info in combo:
                pathway = info[0]
                attractor = info[1]

                if pathway not in pathways:
                    pathways.append(pathway)

                attractors.append(attractor)

            total_attractors.append(attractors)
        
        # Header containing pathways
        combination_file.write('\t'.join(pathways) + '\tNumber_of_cells\tTotal_Cells\tPercent_of_cells\n')

        # Print attractors and cell counts
        for pathway_num, pathway in enumerate(total_attractors):
            percent_of_cells = round(sorted_combinations[pathway_num][1] / (num_cells * len(all_networks)) * 100, 3)
            combination_file.write('\t'.join(str(i) for i in pathway) + '\t' + str(sorted_combinations[pathway_num][1]) + '\t' + str(num_cells) + '\t' + str(percent_of_cells) + '\n')

                

    # # Get the number of unique cluster combinations and count the number of cells belonging to that group
    # attractor_combinations = {}
    # for cell in cell_states:
    #     # Sort and convert to tuple to ensure uniqueness regardless of order
    #     combination = tuple(sorted(cell))
        
    #     # Count the number of times the combination is seen and append that to a dictionary of the combinations
    #     if combination in attractor_combinations:
    #         attractor_combinations[combination] += 1
    #     else:
    #         attractor_combinations[combination] = 1
    
    # # Sort the combinations from highest to lowest based on the number of cells with that attractor combination
    # sorted_combinations = sorted(attractor_combinations.items(), key=lambda item: item[1], reverse=True)
        




    # # Find each unique combination of attractors that the cells mapped to and count the number of cells that mapped to each combination
    # attractor_combinations = find_unique_attractor_combinations(index_to_attractor_dict)    
    
    # # Write out a file containing the attractor combinations and the number of cells that mapped
    # write_attractor_map(dataset_name, attractor_combinations)

    # # Write out files for each attractor in the dataset containing the nodes and whether they are ACTIVE or INACTIVE
    # write_active_inactive_states(dataset_name, all_networks)
    
    # # Write out a file to each attractor for each network containing only the active genes in the network
    # write_only_active_states(dataset_name, all_networks)

            
