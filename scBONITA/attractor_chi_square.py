
import logging
from setup.user_input_prompts import *
import argparse
import scipy.stats as stats
import pickle

from file_paths import file_paths

# Load in the full network pickle file with the cell maps

def load_all_network_pickles(dataset_name):
    """
    Loads in all network pickle files for a given dataset_name
    """
    # Load in the network pickle files
    all_networks = []

    pickle_file_path = f'{file_paths["pickle_files"]}/{dataset_name}_pickle_files/network_pickle_files/'
    for pickle_file in glob.glob(pickle_file_path + str(dataset_name) + "_" + "*" + ".network.pickle"):
        if pickle_file:
            logging.info(f'\tLoading data file: {pickle_file}')
            network = pickle.load(open(pickle_file, "rb"))
            all_networks.append(network)
        else:
            assert FileNotFoundError("Network pickle file not found")

    return all_networks

def create_cell_index_to_group_dict(network):
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

def create_attractor_group_cell_count_dict(network, cell_group_dict):
    attractor_counts = {}
    for cell_index, attractor in network.cell_map.items():
        # logging.info(f'Cell: {key+1}, attractor: {value}')
        logging.info(f'Cell: {cell_index}, attractor: {attractor}, Group {cell_group_dict[cell_index]}')
        group = cell_group_dict[cell_index]

        if attractor not in attractor_counts:
            attractor_counts[attractor] = {}

        if group not in attractor_counts[attractor]:
            attractor_counts[attractor][cell_group_dict[cell_index]] = 1
        else:
            attractor_counts[attractor][cell_group_dict[cell_index]] += 1

    return attractor_counts

if __name__ == '__main__':
    logging.basicConfig(format='%(message)s', level=logging.INFO)
    parser = argparse.ArgumentParser()

    add_dataset_name_arg(parser)

    args = parser.parse_args()

    dataset_name = check_dataset_name(args.dataset_name)

    all_networks = load_all_networks(dataset_name)

    for network in all_networks:

        # Creates a dictionary where each cell index is the key and each value is the group of that cell
        cell_group_dict, group_cell_counts, groups = create_cell_index_to_group_dict(network)
        
        # Creates a count of the number of cells in each group for each attractor. {attractor: {group: count, group: count}}
        cluster_counts = create_attractor_group_cell_count_dict(network, cell_group_dict)

        logging.info(cluster_counts)

        total_cell_counts = sum([group_cell_counts[group] for group in groups])

        contingency_table = []
        for cluster in cluster_counts:
            # Only look at attractors with a significant proportion of the total cells
            cluster_total_counts = sum([cluster_counts[cluster][group] for group in cluster_counts[cluster]])
            if cluster_total_counts / total_cell_counts > 0.05:
                print(f'Attractor {cluster}')
                cluster_observed = []
                cluster_expected = []

                # Count the number of cells in each group for the attractor
                for group in groups:
                    # If the group is in the attractor, add the cell count
                    if group in cluster_counts[cluster]:                
                        print(f'\tGroup {group}, {cluster_counts[cluster][group]} / {group_cell_counts[group]} cells')
                        cluster_observed.append(cluster_counts[cluster][group])
                        cluster_expected.append(group_cell_counts[group])
                    
                    # If there are no cells for a group in the attractor, set the cell count to 0
                    else:
                        print(f'\t{group}: 0 / {group_cell_counts[group]} cells')

                        cluster_observed.append(0)
                        cluster_expected.append(group_cell_counts[group])
                contingency_table.append(cluster_observed)


        # print(f'Observed = {cluster_observed}')
        # print(f'Expected = {cluster_expected}')

        # Perform the chi-squared test
        chi2_stat, p_val, dof, expected = stats.chi2_contingency(contingency_table)

        print(f"Chi-square statistic: {chi2_stat}")
        print(f"P-value: {p_val}")
        print(f"Degrees of freedom: {dof}")
        print(f"Expected frequencies:\n{expected}")

        if p_val < 0.05:
            logging.info(f'p value = {p_val}. There is a significant association between attractor mapping and disease state')
