import numpy as np
import random
from boolean_rule_functions import *
import logging

class ErrorCalculation:
    def __init__(self,nodes,individual_length,nodeDict):

        self.nodes = nodes
        self.individual_length = individual_length
        self.nodeDict = nodeDict

        # ----- Gate choice list/dictionary -----
        self.gate_choices_twogenes = [
            A_AND_B,
            A_OR_B
        ]

        self.gate_choices_threegenes = [
            A_AND_B_AND_C,
            A_AND_B_OR_C,
            A_OR_B_OR_C
        ]

        # ----- Create dictionary of the functions, so I can call them using their names -----
        self.gate_dict_twogenes = {func.__name__: func for func in self.gate_choices_twogenes}
        self.gate_dict_threegenes = {func.__name__: func for func in self.gate_choices_threegenes}

    def evaluate_rule(self, A, B, C, predicted_rule):
        return eval(predicted_rule)
    
    def calculate_error(self, A, B, C, node, predicted_rule):
        predicted = np.vectorize(self.evaluate_rule(A, B, C, predicted_rule))



    def rule_calculation_1input(self, A, node_evaluated, invert_node):
        if invert_node:
            A = np.logical_not(A)

        difference = np.sum(A != node_evaluated)

        return difference, len(A)

    # --- New rule calculation ---
    def rule_calculation_2input(self, A, B, not_a, not_b, node_evaluated, predicted_logic_function, gate_dict_twogenes):
        predicted_logic_function = np.vectorize(gate_dict_twogenes.get(predicted_logic_function))
        
        predicted = predicted_logic_function(A, B, not_a, not_b)
        difference = np.sum(((A + B) >= 1) & (predicted != node_evaluated)) + np.sum((A + B == 0) & (node_evaluated == 1))
        return difference, len(A)

    # --- New rule calculation ---
    def rule_calculation_3input(self, A, B, C, not_a, not_b, not_c, node_evaluated, predicted_logic_function, gate_dict_threegenes):
        # Vectorize the function
        predicted_logic_function = np.vectorize(gate_dict_threegenes.get(predicted_logic_function))

        # Calculate predicted values
        predicted = predicted_logic_function(A, B, C, not_a, not_b, not_c)

        # Identify indices where at least one of A, B, C is 1
        one_gene_active = A + B + C >= 1

        # Calculate differences where at least one of A, B, C is 1
        difference_active = np.sum((predicted[one_gene_active] != node_evaluated[one_gene_active]) &
                                   ((predicted[one_gene_active] + node_evaluated[one_gene_active]) % 2 == 1))

        # Identify indices where all of A, B, C are 0 but the node is active
        all_zero_node_active = (A + B + C == 0) & (node_evaluated == 1)

        # Calculate differences where all of A, B, C are 0 but the node is active
        difference_zero = np.sum(all_zero_node_active)

        total_difference = difference_active + difference_zero # Sum the differences

        return total_difference, len(A)

    def calculate_error(self, total_rules, dataset):
        total_difference = 0
        total_count = 0

        # # Timer for diagnosing error
        # timer = 0

        for node_predictions in total_rules:
            for connection in node_predictions: # For each predicted connection in the graph
                node_name = connection[0] # The node being evaluated is the first element of connection
                incoming_nodes = connection[1] # The incoming nodes and predicted rule is the second element
                rule = connection[2]
                inversion = connection[3]

                node_data = dataset[self.nodeDict[node_name]]

                # For nodes with only 1 input
                if len(incoming_nodes) == 1:
                    input_gene1 = incoming_nodes
                    difference, count = self.rule_calculation_1input(
                        dataset[self.nodeDict[input_gene1[0]]], # The data for the input gene
                        incoming_nodes, # The data for the node being looked at
                        inversion[0]
                    )

                    total_difference += difference
                    total_count += count

                # For nodes with 2 inputs
                elif len(incoming_nodes) == 2:
                    
                    input_gene1, input_gene2, = incoming_nodes
                    difference, count = self.rule_calculation_2input(
                        dataset[self.nodeDict[input_gene1]], # A
                        dataset[self.nodeDict[input_gene2]], # B
                        inversion[0],  # not_a
                        inversion[1],  # not_b
                        node_data, # node_evaluated
                        rule,
                        self.gate_dict_twogenes
                    )

                # For nodes with 3 inputs
                elif len(incoming_nodes) == 3:
                    input_gene1, input_gene2, input_gene3 = incoming_nodes
                    difference, count = self.rule_calculation_3input(
                        dataset[self.nodeDict[input_gene1]], # A
                        dataset[self.nodeDict[input_gene2]], # B
                        dataset[self.nodeDict[input_gene3]], # C
                        inversion[0], # not_a
                        inversion[1], # not_b
                        inversion[2], # not_c
                        node_data,
                        rule,
                        self.gate_dict_threegenes
                    )

                total_difference += difference
                total_count += count

        if total_count > 0:
            total_error = total_difference / total_count
            # print('Total Error = total difference / total count: ', str(total_difference), '/', str(total_count), ' = ',  (total_difference / total_count * 100), '%')

        else:
            total_error = 0
            # print('TOTAL COUNT = 0')

        return total_error

    def calculate_node_error(self, predicted_rule, dataset):
        node_name = predicted_rule[0] # The node being evaluated is the first element of connection
        incoming_nodes = predicted_rule[1] # The incoming nodes and predicted rule is the second element
        rule = predicted_rule[2]
        inversion = predicted_rule[3]

        node_data = dataset[self.nodeDict[node_name]]
        
        # Check if an index exists
        def index_exists(lst, index):
            try:
                _ = lst[index]
                return True
            except IndexError:
                return False
        
        # Set A, B, and C to 0
        A = 0
        B = 0
        C = 0
        
        # Update the values of A, B, and C based on the gene index in the dataset
        if index_exists(incoming_nodes, 0):
            A = dataset[incoming_nodes[0]]
        if index_exists(incoming_nodes, 1):
            B = dataset[incoming_nodes[1]]
        if index_exists(incoming_nodes, 2):
            C = dataset[incoming_nodes[2]]
        
        predicted = np.vectorize(eval(predicted_rule))
        expected = dataset[node.dataset_index]

        difference = 0
        count = 0

        

        if count > 0:
            total_error = difference / count
        else:
            total_error = 0

        return total_error