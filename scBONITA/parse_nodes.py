import logging
from boolean_rule_functions import *

class Node:
    def __init__(self,
                 name,
                 index,
                 three_node_combos,
                 two_node_combos,
                 one_node_combos,
                 predecessors):

        self.name = name # Gene name
        self.index = index # Index in the network
        self.dataset_index = None

        # Possible rules for each of the incoming node combinations
        self.three_node_combos = three_node_combos or []
        self.two_node_AND_combos = two_node_combos or []
        self.two_node_OR_combos = two_node_combos or []
        self.one_node_combos = one_node_combos or []

        # Combining all of the combinations
        self.total_combos = three_node_combos + two_node_combos + one_node_combos
        self.total_combo_length = len(self.total_combos)
        self.rvalues = None

        # Inversion rules for each of the incoming node combinations
        self.inversions = []
        self.inversion_rules = []

        self.predecessors = predecessors # incoming nodes for this node
        self.incoming_node_indices = []
        self.rule_length = self.find_rule_length()

        # Indices in the individuals where the rules for this node start and end
        self.rule_start_index = None
        self.rule_end_index = None

        # Rule Predictions
        self.rule_predictions = None

        self.three_node_predictions = None
        self.two_node_AND_predictions = None
        self.two_node_OR_predictions = None
        self.one_node_predictions = None
        self.node_rules = []
        self.best_rule_index = None
        self.best_rule = None
        self.best_rule_error = None

        self.logic_function = None

        # Compressed best rule information
        self.calculation_function = None

        # Importance scores
        self.importance_score = 0
        self.importance_score_stdev = 0.0

        # Relative abundance
        self.relative_abundance = None

    # Find the length of the rule combinations for the current node
    def find_rule_length(self):
        """
        Find the length of the rule combinations for the current node (for determining how many bits should be in
        an individual)

        :Return: rule_length
        """
        rule_length = len(self.three_node_combos)\
                      + len(self.two_node_AND_combos) \
                      + len(self.two_node_OR_combos)\
                      + len(self.one_node_combos)
        return rule_length

    # Return the possible combinations that are predicted for the individual
    def find_rule_predictions(self):
        """
        Looks at the rule predictions for the current node in the individual. Returns the predicted rules for each
        combination of incoming nodes.
        Returns
        -------

        """

        # Finds the rules predicted for this node in the deap individual
        if self.rule_predictions:
            self.three_node_predictions = [rule for i, rule in enumerate(self.three_node_combos) if self.rule_predictions[i] == 1]
            self.two_node_AND_predictions = [rule for i, rule in enumerate(self.two_node_AND_combos) if self.rule_predictions[i + len(self.three_node_combos)] == 1]
            self.two_node_OR_predictions = [rule for i, rule in enumerate(self.two_node_OR_combos) if self.rule_predictions[i + len(self.three_node_combos) + len(self.two_node_AND_combos)] == 1]
            self.one_node_predictions = [rule for i, rule in enumerate(self.one_node_combos) if self.rule_predictions[i + len(self.three_node_combos) + len(self.two_node_AND_combos) + len(self.two_node_OR_combos)] == 1]

    # Format the AND OR and NOT rules
    def make_AND_OR_NOT_rules(self):
        if self.rule_predictions:

            # Make the rules for three-incoming node combinations
            if self.three_node_predictions:
                # print(f'\tThree node rules:')
                three_node_OR_nodes = []

                # Find the incoming nodes that aren't in the AND connection, set these as the OR connections
                for combination in self.three_node_predictions:
                    # Create sets of the incoming nodes and the combination
                    incoming_nodes = set(self.predecessors.keys())
                    rule_combination = set(combination)

                    # Find which nodes are not in the combination
                    incoming_nodes_not_in_combination = incoming_nodes - rule_combination

                    # The incoming nodes not in the combination are set as the OR rules
                    three_node_OR_nodes.append([node for node in incoming_nodes_not_in_combination])

                # Put together the rules for three incoming node combinations
                rules = []
                for i, prediction in enumerate(self.three_node_predictions):
                    # Retrieve the names of the nodes from the predecessors dictionary
                    and_rules = []
                    or_rules = []
                    placeholder = ['A', 'B', 'C']
                    and_nodes = prediction
                    or_nodes = three_node_OR_nodes[i]

                    # Process AND nodes
                    for i, node in enumerate(and_nodes):
                        and_rules.append(placeholder[i])
                    rule_string = '_AND_'.join(and_rules)

                    if or_nodes:
                        rule_string = rule_string + '_OR_'
                        for i, node in enumerate(or_nodes):
                            or_rules.append(placeholder[i + len(and_rules)])
                        rule_string = rule_string + '_OR_'.join(or_rules)

                    # Append the rule string to the list
                    rules.append(rule_string)

                # Create the rules for the node from the predictions
                for i, prediction in enumerate(self.three_node_predictions):
                    # Extract the relevant information for the rules
                    and_node_names = [self.predecessors[node] for node in prediction]
                    and_inversion_rules = [self.inversions[node] for node in prediction]
                    or_node_names = [self.predecessors[node] for node in three_node_OR_nodes[i]]
                    or_inversion_rules = [self.inversions[node] for node in three_node_OR_nodes[i]]

                    # Combine the names for AND rule predictions and OR rule predictions
                    incoming_node_names = and_node_names + or_node_names

                    # Combine the inversion rules for the AND predictions adn OR predictions
                    inversion_rules = and_inversion_rules + or_inversion_rules

                    # print(f'\t\t\tRULES: {node_names}, {rules[i]}')
                    # Combine the info into the correct format
                    self.node_rules.append([self.name, incoming_node_names, rules[i], inversion_rules])

            # Make the rules for two-incoming node AND combinations
            if self.two_node_AND_predictions:
                # print(f'\tTwo node rules:')
                # print(f'\t\tAND_nodes_list: {self.two_node_AND_predictions}')
                rules = f'A_AND_B'
                for i, prediction in enumerate(self.two_node_AND_predictions):
                    # Extract the relevant information for the rules
                    incoming_and_node_names = [self.predecessors[node] for node in prediction]
                    inversion_rules = [self.inversions[node] for node in prediction]

                    # Combine the incoming node names and the rules
                    self.node_rules.append([self.name, incoming_and_node_names, rules, inversion_rules])
                    # print(f'\t\t\tRULES: {and_node_names}, {rules}')

            # Make the rules for two-incoming node OR combinations
            if self.two_node_OR_predictions:
                # print(f'\t\tOR_nodes_list: {self.two_node_OR_predictions}')
                rules = f'A_OR_B'
                for i, prediction in enumerate(self.two_node_OR_predictions):
                    # Extract the relevant information for the rules
                    incoming_or_node_names = [self.predecessors[node] for node in prediction]
                    inversion_rules = [self.inversions[node] for node in prediction]

                    # Combine the incoming node names and the rules
                    self.node_rules.append([self.name, incoming_or_node_names, rules, inversion_rules])
                    # print(f'\t\t\tRULES: {or_node_names}, {rules}')


            # Make the rules for one-incoming node combinations
            if self.one_node_predictions:
                # print(f'\tOne node rules:')
                # print(f'\t\tOR_nodes_list: {self.one_node_predictions}')
                for i, prediction in enumerate(self.one_node_predictions):
                    node_name = [self.predecessors[node] for node in prediction]
                    inversion_rules = [self.inversions[node] for node in prediction]
                    # print(f'\t\t\tRULES: {node_name}')
                    self.node_rules.append([self.name, node_name, 'A', inversion_rules])
        else:
            print(f'\tNo predictions')

    def find_calculation_function(self, rule):
        # ----- Gate choice list/dictionary -----
        possible_rules = [
            A_AND_B_AND_C,
            A_AND_B_OR_C,
            A_OR_B_OR_C,
            A_AND_B,
            A_OR_B,
            A
        ]

        # Find the function based on the name of the function
        rules = {func.__name__: func for func in possible_rules}

        # Find which rule function is in the best rule
        if rule in rules:
            calculation_function = rules[rule]

            return calculation_function

        else:
            msg = f'ERROR: rule "{rule}" is not in the list of possible rules'
            assert Exception(msg)
            

    # Print information about the node
    def print_info(self):
        logging.info(f'\nNode {self.name} Length: {self.rule_length}')
        logging.info(f'\tindex {self.index}')
        logging.info(f'\tpredecessors: {self.predecessors}')
        logging.info(f'\tinversions: {self.inversions}')
        if self.rule_predictions:
            logging.debug(f'\tpredictions {self.rule_predictions} (length {len(self.rule_predictions)})')
            logging.debug(f'\tAll Predictions: {self.three_node_combos + self.two_node_AND_combos + self.two_node_OR_combos + self.one_node_combos}')
            logging.debug(f'\tthree_node_predictions {self.three_node_predictions}')
            logging.debug(f'\ttwo_node_AND_predictions {self.two_node_AND_predictions}')
            logging.debug(f'\ttwo_node_OR_predictions {self.two_node_OR_predictions}')
            logging.debug(f'\tone_node_predictions {self.one_node_predictions}')

            logging.info(f'\tnode_rules: {self.node_rules}')
        else:
            logging.info(f'\tNo predictions')

    # Reset the state of the node
    def reset_state(self):
        """
        Reset the state of the Node instance to clear any data that is specific to an individual's rule set.
        This method should be called before processing a new individual.
        """
        self.rule_predictions = None
        self.three_node_predictions = None
        self.two_node_AND_predictions = None
        self.two_node_OR_predictions = None
        self.one_node_predictions = None
        self.node_rules = []  # Clearing the list of rules

        # Any other attributes that need to be reset can be added here



