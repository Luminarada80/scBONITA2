from network_setup import *
from keggParser import *
import scipy.sparse as sparse
import numpy as np
from os import path
import matplotlib.pyplot as plt
import time
from deap_class import CustomDeap
import csv
import networkx as nx
from cell_class import Cell

class RuleInference(NetworkSetup):

    """Class for single-cell experiments"""

    def __init__(
        self,
        data_file,
        dataset_name,
        network_name,
        sep,
        node_indices,
        gene_list,
        max_nodes=15000,
        binarize_threshold=0.001,
        sample_cells=True,
    ):
        # Check is the data file exists
        if path.isfile(data_file):
            pass
        else:
            raise FileNotFoundError(f'File not found: {data_file}')
        
        self.data_file = data_file
        self.sample_cells = sample_cells
        self.node_indices = node_indices
        self.dataset_name = dataset_name
        self.network_name = network_name
        self.binarize_threshold = binarize_threshold
        self.max_samples = 15000
        self.cells = []

        logging.info(f'\n-----EXTRACTING AND FORMATTING DATA-----')
        # Extract the data from the data file based on the separator, sample the cells if over 15,000 cells
        logging.info(f'Extracting cell expression data from "{data_file}"')
        self.cell_names, self.gene_names, self.data = self._extract_data(data_file, sep, sample_cells, node_indices)


        logging.info(f'\tFirst 2 genes: {self.gene_names[:2]}')
        logging.info(f'\tFirst 2 cells: {self.cell_names[:2]}')
        
        logging.info(f'\tNumber of genes: {len(self.gene_names)}')
        logging.info(f'\tNumber of cells: {len(self.cell_names)}')

        self.sparse_matrix = sparse.csr_matrix(self.data)
        logging.info(f'\tCreated sparse matrix')
        logging.debug(f'\tShape: {self.sparse_matrix.shape}')
        
        self.gene_names = list(self.gene_names)
        self.cell_names = list(self.cell_names)
        self.sparse_matrix.eliminate_zeros()

        # Check if there are at least 1 sample selected
        if self.data.shape[0] > 0:
            # Binarize the values in the sparse matrix
            logging.info(f'\tBinarized sparse matrix')
            self.binarized_matrix = preprocessing.binarize(self.sparse_matrix, threshold=binarize_threshold, copy=True)
            logging.debug(f'{self.binarized_matrix[:5,:5]}')
            
        else:
            # Handle the case where no samples were selected
            raise ValueError("No samples selected for binarization")
        
        full_matrix = self.binarized_matrix.todense()

        # Create cell objects
        for cell_index, cell_name in enumerate(self.cell_names):
            cell = Cell(cell_index)
            cell.name = cell_name
            for row_num, row in enumerate(full_matrix):
                row_array = np.array(row).flatten()
                cell.expression[self.gene_names[row_num]] = row_array[cell_index]

            self.cells.append(cell)


        self.max_nodes = max_nodes
        self.pathway_graphs = {}
        self.node_list = []
        self.node_positions = []
        self.gene_list = gene_list

    def _extract_data(self, data_file, sep, sample_cells, node_indices):
        """
        Extract the data from the data file
        Parameters
        ----------
        data_file
        sep
        sample_cells

        Returns
        -------
        cell_names, data
        """
        with open(data_file, 'r') as file:
            reader = csv.reader(file, delimiter=sep)

            # Extract the header (cell_names)
            cell_names = next(reader)[1:]

            cell_count = len(cell_names)

            # Randomly sample the cells in the dataset
            if cell_count >= self.max_samples or sample_cells:
                logging.info(f'\tRandomly sampling {self.max_samples} cells...')
                sampled_cell_indices = np.random.choice(
                    range(cell_count),
                    replace=False,
                    size=min(self.max_samples, cell_count),
                )
                logging.info(f'\t\tNumber of cells: {len(sampled_cell_indices)}')

            else:
                sampled_cell_indices = range(cell_count)
                logging.info(f'\tLoading all {len(sampled_cell_indices)} cells...')

            # Data extraction
            data_shape = (len(node_indices), len(sampled_cell_indices))
            data = np.empty(data_shape, dtype="float")
            gene_names = []
            data_row_index = 0  # Separate index for data array

            for i, row in enumerate(reader):
                if (i + 1) in node_indices:  # Adjust index for skipped header
                    gene_names.append(row[0])
                    # Offset cell indices by 1 to skip the gene name column
                    selected_data = [float(row[cell_index + 1]) for cell_index in sampled_cell_indices]
                    data[data_row_index, :] = selected_data
                    data_row_index += 1

            # Convert the filtered data to a NumPy array
            logging.info("\tConverting filtered data to numpy array...")
            # print(f'Raw data:\n\t{data}')

            return cell_names, gene_names, data

    def filterData(self, threshold):
        """Filters the data to include genes with high variability (genes with a std dev / mean ratio above the cv_cutoff threshold)"""
        self.cv_genes = []
        if threshold is not None:
            for i in range(0, self.sparse_matrix.get_shape()[0]):
                rowData = list(self.sparse_matrix.getrow(i).todense())
                if np.std(rowData) / np.mean(rowData) >= threshold:
                    self.cv_genes.append(self.gene_names[i])
        else:
            self.cv_genes = copy.deepcopy(self.gene_names)
    
    def genetic_algorithm(self, net):
        # Genetic algorithm
        custom_deap = CustomDeap(
            net,
            self.network_name,
            self.dataset_name,
            self.binarized_matrix,
            self.nodeList,
            self.nodes,
            self.deap_individual_length,
            self.nodeDict,
            self.successorNums
            )

        # raw_fitnesses, population, logbook = custom_deap.genetic_algorithm()
        # best_ruleset = custom_deap.find_best_individual(population, raw_fitnesses)

        return custom_deap.chunked_data_numpy

    def rule_determination(self, graph):
        """Main function that performs rule determination and node scoring in preparation for pathway analysis"""

        # Load the processed graphml file as a NetworkX object
        start_time = time.time()
        if path.exists(graph):
            logging.info(f'\t\tLoading: {graph.split("/")[-1]}')
            self.network = nx.read_graphml(graph)
        else:
            msg = f'File "{graph} not found"'
            raise FileNotFoundError(msg)

        # Find the network gene names
        netGenes = [self.gene_names.index(gene) for gene in list(self.network) if gene in self.gene_names]

        logging.debug(f'Network Genes: {netGenes}')

        self.__inherit(
            self.network,
            removeSelfEdges=False,
            restrictIncomingEdges=True,
            maxIncomingEdges=3,
            groundTruth=False,
        )

        chunked_data_numpy = self.genetic_algorithm(self.network)
        
        # Runs the genetic algorithm and rule refinement
        self.best_ruleset = self.fitness_calculation(self.nodes, chunked_data_numpy)

    def fitness_calculation(self, nodes, dataset):
        best_rule_dict = {}

        # Calculate the error for each node in the individual
        print(f'\n Calculate error for each node in the dataset')
        for node in nodes:

            errors = {}

            print(f'{node.name}')
            # print(f'\tName: {node.name}')
            # print(f'\tIndex: {node.index}')
            # print(f'\tPredecessors: {node.predecessors}')
            # print(f'\tIncoming node indices: {node.incoming_node_indices}')
            # print(f'\tInversions: {node.inversions}')
            
            # get the indices in the dataset for the incoming nodes
            incoming_node_inversions = []
            # print(f'Incoming node indices:')
            for incoming_node in node.predecessors.keys():
                node.incoming_node_indices.append(incoming_node)
                incoming_node_inversions.append(node.inversions[incoming_node])
            
            # Get the row in the dataset for the node being evaluated
            node_evaluated = dataset[node.index]

            # Get the dataset row indices for the incoming nodes included in this rule
            incoming_node_indices = [index for index in node.predecessors]
            incoming_node_names = [value for value in node.predecessors.values()]
            incoming_node_inversions = [node.inversions[index] for index in incoming_node_indices]

            # Get the row(s) in the dataset for the incoming node(s)
            incoming_node_data = [dataset[i] for i in incoming_node_indices]
            # print(incoming_node_data)

            # for each node in the dataset, calculate the error for each rule combination
            
            print(f'\tIncoming nodes: {incoming_node_indices}')
            print(f'\tNumber of incoming nodes: {len(incoming_node_indices)}')
            print(f'\tInversions: {incoming_node_inversions}')
            possible_indices = []

            A = None
            B = None
            C = None

            not_a = None
            not_b = None
            not_c = None
            
            # Dont try to calculate 3 node rules if there are only 2 incoming nodes
            # Find each of the possible indices 
            # (for three nodes, its all of the possibilities)
            # print(f'len(incoming_node_data): {len(incoming_node_data)}')

            if len(incoming_node_indices) == 0:
                best_rule_dict[node.name] = node.name
            
            else:

                if len(incoming_node_indices) == 3:
                    possible_indices = range(len(node.possible_rules.keys()))
                    A = incoming_node_data[0]
                    B = incoming_node_data[1]
                    C = incoming_node_data[2]

                    not_a = incoming_node_inversions[0]
                    not_b = incoming_node_inversions[1]
                    not_c = incoming_node_inversions[2]
                    
                # (for two nodes, its only the two and one node possibilities)    
                elif len(incoming_node_indices) == 2:
                    possible_indices = [13, 14, 15, 16]
                    A = incoming_node_data[0]
                    B = incoming_node_data[1]

                    not_a = incoming_node_inversions[0]
                    not_b = incoming_node_inversions[1]
                    
                # (for one node, its only the one node possibilities)      
                elif len(incoming_node_indices) == 1:
                    possible_indices = [16]
                    A = incoming_node_data[0]

                    not_a = incoming_node_inversions[0]
                
                # print(f'len(incoming_node_data[0]): {len(incoming_node_data[0])}')
                for column, _ in enumerate(incoming_node_data[0]):      
                    # print(f'Column: {column}')

                    # If there are possible rules (i.e. there are incoming nodes), calculate error
                    if len(possible_indices) > 0:
                        # Calculate the error for each possible rule
                        for i in node.possible_rules:
                            if i in possible_indices:
                                if A is not None:
                                    # print(f'\tA = {A[column]} not_a = {not_a}')
                                    A_new = (not A[column] if not_a else A[column])
                                if B is not None:
                                    # print(f'\tB = {B[column]} not_b = {not_b}')
                                    B_new = (not B[column] if not_b else B[column])
                                if C is not None:
                                    # print(f'\tC = {C[column]} not_c = {not_c}')
                                    C_new = (not C[column] if not_c else C[column])
                                # print(f'\t{node.possible_rules[i]}')
                                
                                # I want to ignore cases where none of the rules are 1
                                if A_new + B_new + C_new > 0:
                                    if i == 0:
                                        # A and B and C
                                        result = int(A_new and B_new and C_new)

                                    elif i == 1:
                                        # (A and B) or C
                                        result = int((A_new and B_new) or C_new)
                                    
                                    elif i == 2:
                                        result = int(A_new and (B_new or C_new))
                                    
                                    elif i == 3:
                                        result = int((A_new or B_new) and C_new)
                                    
                                    elif i == 4:
                                        result = int(A_new or (B_new and C_new))
                                    
                                    elif i == 5:
                                        result = int((A_new and C_new) or B_new)
                                    
                                    elif i == 6:
                                        result = int((A_new or C_new) and B_new)
                                    
                                    elif i == 7:
                                        result = int(A_new or C_new or B_new)
                                    
                                    elif i == 8:
                                        result = int(A_new and C_new)
                                    
                                    elif i == 9:
                                        result = int(A_new or C_new)
                                    
                                    elif i == 10:
                                        result = int(B_new and C_new)
                                    
                                    elif i == 11:
                                        result = int(B_new or C_new)
                                    
                                    elif i == 12:
                                        result = int(C_new)
                                    
                                    elif i == 13:
                                        result = int(A_new and B_new)
                                    
                                    elif i == 14:
                                        result = int(A_new or B_new)
                                    
                                    elif i == 15:
                                        result = int(B_new)
                                    
                                    elif i == 16:
                                        result = int(C_new)
                                    
                                    # print(f'\t\tPredicted = {result}')
                                    # print(f'\t\tActual =  {node_evaluated[column]}')

                                    for x in possible_indices:
                                        if not x in errors:
                                            errors[x] = 0

                                    if node_evaluated[column] != result:
                                        errors[i] += 1
                                        # print(f'\t\tERROR')

                formatted_errors = {}
                print(f'\tRule Errors:')

                for key, value in errors.items():
                    if key in possible_indices:
                        rule_template = node.possible_rules[key]
                        
                        # Dictionary to hold the actual values to substitute into the rule
                        substitutions = {}
                        
                        # Iterate over each incoming node to determine the substitution
                        for i, inversion in enumerate(incoming_node_inversions):
                            placeholder = chr(65 + i)  # 'A', 'B', 'C' corresponding to 0, 1, 2
                            node_name = incoming_node_names[i]

                            if inversion:
                                substitutions[placeholder] = f'not {node_name}'
                            else:
                                substitutions[placeholder] = node_name

                        # Use str.format_map() to replace placeholders with actual values
                        rule = rule_template.format_map(substitutions)

                        print(f'\t\t{rule}: {value}')
                        
                        formatted_errors[rule] = value
                    
                # Assuming 'errors' is already defined and 'node.possible_rules' maps rule indices to rule descriptions

                # Step 1: Calculate the average error
                total_error = sum(errors.values())
                num_rules = len(errors)
                average_error = total_error / num_rules if num_rules else 0

                # Step 2: Define the threshold
                threshold = 0.75 * average_error

                # Step 3: Filter the rules
                significantly_lower_errors = {key: value for key, value in formatted_errors.items() if value < threshold}

                # Step 4: Output results
                print(f'\t\tAverage Error: {average_error}')

                if len(significantly_lower_errors.items()) > 0:
                    print("\t\tRules with significantly lower errors than average:")
                    min_val = min(significantly_lower_errors.values())
                    print(f'\t\tmin_val = {min_val}')
                    for rule, error in significantly_lower_errors.items():
                        print(f'\t\t\tRule {rule}: {error}')
                        if error == min_val and not node.name in best_rule_dict:
                            best_rule_dict[node.name] = rule

                            for i, (formatted_rule, value) in enumerate(formatted_errors.items()):
                                if formatted_rule == rule:
                                    print(f'\t\t\tRule index = {i}')

                                    node.best_rule = i
                                    print(f'\t\t\t\tBEST RULE: {node.best_rule}')
                            
                    
                            print(f'\t\t\t\tmin_rule = {best_rule_dict[node.name]}')

                else:
                    print(f'\t\tNo low error rules for {node.name}')
                    best_rule_dict[node.name] = node.name

        print(f'\n')
        for key, value in best_rule_dict.items():
            print(f'{key} = {value}')
        
        self.best_rule_dict = best_rule_dict

        return best_rule_dict

    def plot_graph_from_graphml(self, network):
        # Load the graph
        G = network
        logging.debug(f'\n-----IMPORTANCE SCORE GRAPH-----')

        # Assuming each node has an attribute 'value' that is a number between 0 and 1
        values = [node.importance_score for node in self.nodes]
        logging.debug(f'Values: {values}')

        # Choose a colormap (e.g., 'viridis', 'plasma', 'inferno', 'magma', 'cividis')
        cmap = plt.cm.Greys

        # Map 'values' to colors
        node_colors = [cmap(value) for value in values]
        logging.debug(f'node_color: {node_colors}')

        pos = nx.spring_layout(G, k=1)

        # Invert the text color if the value is above a certain value, so the gene name can still be read
        text_colors = ['white' if value > 0.8 else 'black' for value in values]


        # Draw the graph
        fig = plt.figure(figsize=(12, 8))
        # Draw the graph
        fig = plt.figure(figsize=(12, 8))
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=500)
        nx.draw_networkx_edges(G, pos, edge_color='gray')
        nx.draw_networkx_labels(G, pos, font_color=text_colors, font_size=10)
        plt.title("Importance score for each node in the network")
        
        return fig
    
    def __inherit(
        self,
        graph,
        removeSelfEdges=False,
        restrictIncomingEdges=True,
        maxIncomingEdges=3,
        groundTruth=False,
        graphName=""):
        super().__init__(
            graph,
            removeSelfEdges,
            restrictIncomingEdges,
            maxIncomingEdges,
            groundTruth,
            graphName,
        )

