from rule_inference import RuleInference
from pathways import Pathways
from network import Network
from importance_scores import CalculateImportanceScore
from cell_class import CellPopulation

# Setup
from setup.scbonita_banner import make_banner
from setup.user_input_prompts import rule_inference_arguments

# General packages
import glob
import os
import pickle
import logging
from datetime import timedelta
import pandas as pd
from sklearn.preprocessing import MaxAbsScaler
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser


class Pipeline():
    def __init__(self):
        # Set the logging level for output
        logging.basicConfig(format='%(message)s', level=logging.INFO)

        parser = ArgumentParser()

        # Read in arguments from shell script using argparse
        parser_argument = rule_inference_arguments(parser)
        
        # Store the results from argparse
        self.dataset_name = getattr(parser_argument, 'dataset_name', "")
        self.network_files = getattr(parser_argument, 'network_files', [])
        self.max_nodes = getattr(parser_argument, 'max_nodes', 20000)
        self.data_file = parser_argument.data_file
        self.datafile_sep = getattr(parser_argument, 'datafile_sep', ',')
        self.get_kegg_pathways = self._convert_string_to_boolean(getattr(parser_argument, 'get_kegg_pathways', True))
        self.list_of_kegg_pathways = getattr(parser_argument, 'list_of_kegg_pathways', [])
        self.pathway_list = getattr(parser_argument, 'pathway_list', [])
        self.organism = getattr(parser_argument, 'organism', 'hsa')
        self.cv_threshold = getattr(parser_argument, 'cv_threshold', 0.001)
        self.sample_cells = self._convert_string_to_boolean(getattr(parser_argument, 'sample_cells', True))
        self.binarize_threshold = float(getattr(parser_argument, 'binarize_threshold', 0.001))
        self.display_title = getattr(parser_argument, 'display_title', True)

        # Other variables
        self.write_graphml = True
        self.node_indices = []

        if self.display_title:
            logging.info(make_banner())
        
        self.rule_inference()
    
    def rule_inference(self):
        # Create a Pathways object to store the network information
        pathways = Pathways(
            self.dataset_name, self.cv_threshold, self.data_file, self.datafile_sep, self.write_graphml, self.organism)

        # Use a list of kegg pathways specified by the user
        logging.info(f'-----PARSING NETWORKS-----')
        if self.list_of_kegg_pathways is None:
            logging.info(f'\tNo KEGG Pathways listed...')
            self.list_of_kegg_pathways = []
        else:
            logging.info(f'\tKEGG pathways = {self.list_of_kegg_pathways}')


        # Get all kegg pathways matching the genes in the dataset
        if self.get_kegg_pathways == True:
            logging.info(f'\tFinding and formatting KEGG Pathways...')

            pathways.pathway_graphs = pathways.find_kegg_pathways(
                kegg_pathway_list=self.list_of_kegg_pathways,
                write_graphml=self.write_graphml,
                organism=self.organism,
                minimumOverlap=25
            )

            pathways.add_pathways(pathways.pathway_graphs, minOverlap=1, organism=self.organism)
        
        # Use the pathway(s) specified by pathway_list
        else:
            if isinstance(self.pathway_list, str):
                logging.info(f'\tPathways = {self.pathway_list}')
                pathways.add_pathways([self.pathway_list], minOverlap=1, organism=self.organism)

            elif isinstance(self.pathway_list, list):
                logging.info(f'\tPathways = {self.pathway_list}')
                pathways.add_pathways(self.pathway_list, minOverlap=1, organism=self.organism)

            # Else throw an exception if no pathways are specified
            else:
                msg = f'ERROR: get_kegg_pathways = {self.get_kegg_pathways} and pathway_list = {self.pathway_list}. ' \
                    f'If get_kegg_pathways is False, specify a list or string of pathways to use'
                assert Exception(msg)

        if len(self.network_files) > 0:
            logging.info(f'\tCustom pathway: {self.network_files}')
            pathways.add_pathways(self.network_files, minOverlap=1, organism=self.organism)

        # Get the information from each pathway and pass the network information into a ruleset object
        for pathway in pathways.pathway_graphs:
            node_indices = []
            graph = pathways.pathway_graphs[pathway]
            
            for node in graph.nodes():
                node_indices.append(pathways.gene_list.index(node))
                

            # node_indices = set(node_indices)  # Only keep unique values
            self.node_indices = list(node_indices)  # Convert back to a list

            logging.info(f'\n-----RULE INFERENCE-----')
            logging.info(f'Pathway: {pathway}')
            logging.info(f'Num nodes: {len(node_indices)}')

            # Generate the rule inference object for the pathway
            ruleset = self.generate_ruleset(pathway, node_indices, pathways.gene_list)

            processed_graphml_path = f'graphml_files/{self.dataset_name}/{self.organism}{pathway}_processed.graphml'

            self.infer_rules(pathway, processed_graphml_path, ruleset)
    
    def generate_ruleset(self, network, node_indices, gene_list):
        # Create RuleInference object
        ruleset = RuleInference(
            self.data_file,
            self.dataset_name,
            network,
            self.datafile_sep,
            node_indices,
            gene_list,
            self.max_nodes,
            self.binarize_threshold,
            self.sample_cells)

        logging.debug(f'Pipeline: Created ruleset object')
        
        # Filter the data based on the cutoff value threshold
        ruleset.filterData(threshold=self.cv_threshold)
        
        # Set the parameters in ruleset based on the dataset
        logging.info(f'\tSetting ruleset parameters')
        logging.debug(f'Pipeline: Setting up the ruleset parameters')
        ruleset.max_samples = 15000 #max_samples
        ruleset.gene_list = [ruleset.gene_list[node] for node in node_indices]
        ruleset.node_list = ruleset.gene_list
        ruleset.node_positions = [ruleset.gene_list.index(node) for node in ruleset.node_list]

        return ruleset

    def infer_rules(self, network_name, processed_graphml_path, ruleset):
        logging.info(f'\tRunning rule inference for {network_name}')

        # Specify and create the folder for the dataset pickle files
        data_pickle_folder = f'pickle_files/{self.dataset_name}_pickle_files/ruleset_pickle_files'
        os.makedirs(data_pickle_folder, exist_ok=True)

        # Specify the name of the pickle file for the dataset and network ruleset object
        data_pickle_file_path = f'{data_pickle_folder}/{self.dataset_name}_{self.organism}{network_name}.ruleset.pickle'

        # Run the rule inference
        ruleset.rule_determination(graph=str(processed_graphml_path))

        # Write out the cells objects to a pickle file
        cell_population = CellPopulation(ruleset.cells)

        cell_pickle_file_path = f'pickle_files/{self.dataset_name}_pickle_files/cells_pickle_file'
        os.makedirs(cell_pickle_file_path, exist_ok=True)

        cells_pickle_file = f'{self.dataset_name}.cells.pickle'
        with open(f'{cell_pickle_file_path}/{cells_pickle_file}', "wb") as file:
            pickle.dump(cell_population, file)
        
        logging.info(f'\nRule inference complete, saving to {data_pickle_file_path.split("/")[-1]}')
        pickle.dump(ruleset, open(data_pickle_file_path, "wb"))

    def _convert_string_to_boolean(self, variable):
        if variable == "True" or variable is True:
            return True
        else:
            return False

if __name__ == "__main__":
    Pipeline()

