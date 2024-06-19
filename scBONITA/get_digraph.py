import pickle
import glob
import logging
import numpy as np

from file_paths import file_paths

def load_all_network_pickles(dataset_name: str):
    """
    Loads in all network pickle files for a given dataset_name
    """
    # Load in the network pickle files
    all_networks = []

    pickle_file_path = f'{file_paths["pickle_files"]}/{dataset_name}_pickle_files/network_pickle_files/'
    for pickle_file in glob.glob(pickle_file_path + str(dataset_name) + "_" + "*" + ".network.pickle"):
        if pickle_file:
            logging.info(f'\tLoading data file: {pickle_file}')
            network: object = pickle.load(open(pickle_file, "rb"))
            all_networks.append(network)
        else:
            assert FileNotFoundError("Network pickle file not found")

    return all_networks

all_networks = load_all_network_pickles('george_hiv')

network0 = all_networks[0]

vertices = []
edges = []

for node in network0.nodes:
    vertices.append(node.index)

for node in network0.nodes:
    if len(node.predecessors.keys()) > 1:
        for predecessor in node.predecessors.keys():
            edges.append(tuple([predecessor, node.index]))

print(vertices)

print(edges)
