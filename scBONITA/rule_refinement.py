import logging
from itertools import product
import numpy as np

def find_min_error_individuals(population, fitnesses):
    """Finds the individuals with the minimum error from the population"""
    # logging.info(f'Finding the minimum error individuals')
    min_error_individuals = [ind for ind in population if ind.fitness.values == min(fitnesses)]

    # Set to keep track of unique value combinations
    seen_values = set()

    # List comprehension to filter out duplicates
    unique_min_error_individuals = []
    for ind in min_error_individuals:
        # Convert the individual's values to a tuple (or another immutable type if necessary)
        values_tuple = tuple(ind[1])

        # Check if we've already seen these values
        if values_tuple not in seen_values:
            unique_min_error_individuals.append(ind)
            seen_values.add(values_tuple)
    
    return unique_min_error_individuals

def get_third_element(item):
    """
    Returns the third item in a list, gets around the restriction of 
    using lambda when pickling objects
    """
    return item[2]

def handle_self_loop(node):
    """
    Handles instances when a node only connects to itself or has
    no incoming nodes
    """
    # If the node has only an incoming connection from itself, skip the rule refinement
    # Add in self signaling for the node
    # logging.info(f'\tNode {node.name} signals only to itself')
    # logging.info(f'\t[{node.name}, [{node.index}]')
    node.node_rules = [node.name, [node.index], 'A']
    node.predecessors[node.index] = node.name
    node.inversions[node.index] = False

    # Set the node's incoming node rules to itself
    best_rule = ([node.name, [node.index], 'A'], 0, 0)

    return best_rule

def calculate_refined_errors(node, chunked_dataset):
    """
    1) Finds all possible rule combinations for the node
    2) Calculates the error for each rule based on chunked_dataset
    3) Finds the rules with the minimum error
    4) Finds the simplest rules with the minimum error (fewest incoming nodes)
    5) Finds the simplest rules with the greatest number of 'or' connections

    returns best_rules
    """

    # Calculate the error for each prediction based on the entire dataset (not chunked)
    # logging.info(f'\tCalculate_refined_errors function:')
    prediction_errors = []
    most_incoming_node_rules = []
    best_rule_errors = []
    best_rule_indices = []
    best_rules = []

    # logging.info(f'\t\tPredicted rule errors:')

    # # Calculates the error for each possible rule for the node
    # logging.info(f'\t\t\tNode rules: {node.node_rules}')
    rules = node.find_all_rule_predictions()

    def generate_not_combinations(rule):

        variables = []

        if 'A' in rule:
            variables.append('A')
        if 'B' in rule:
            variables.append('B')
        if 'C' in rule:
            variables.append('C')
        if 'D' in rule:
            variables.append('D')

        # Remove the existing 'not' statements
        rule = rule.replace('not ', '')

        # Generate all possible combinations of the variables and their negations
        combinations = list(product(*[(var, f'not {var}') for var in variables]))
        # logging.info(f'\t\t\tNode combinations: {combinations}')
        
        results = []
        for combo in combinations:
            # Create a temporary rule with the current combination
            temp_rule = rule
            for var, replacement in zip(variables, combo):
                # Replace variables with the current combination of variable or its negation
                temp_rule = temp_rule.replace(var, replacement)
            results.append(temp_rule)
        
        return results
    
    # Generate all 'not' combinations for each rule
    
    not_combinations = [generate_not_combinations(rule[2]) for rule in rules]
    # print('NOT COMBINATIONS')
    for combination in not_combinations:
        for rule in combination:
            # print(rule)
            incoming_node_indices = [predecessor_index for predecessor_index in node.predecessors]
            formatted_rule = [node.name, incoming_node_indices, rule]
            if formatted_rule not in rules:
                rules.append([node.name, incoming_node_indices, rule])

    # logging.info(f'\n\n\n RULES')
    # for rule in rules:
    #     print(rule)

    # Calculate the error for each rule
    for rule in rules:
        # logging.info(f'\t\t\tPredicted rule: {rule}')
        difference, count = calculate_error(node, rule[2], chunked_dataset)
        prediction_error = difference / count
        # logging.info(f'\t\t\t\tError: {prediction_error}')
        prediction_errors.append(prediction_error)

    # Find all of the rules that have the minimum error
    if prediction_errors: # Excludes nodes with no predictions
        # logging.info(f'\t\t\tPrediction errors: {prediction_errors}')
        minimum_error = min(prediction_errors)
        minimum_error_indices = [index for index, value in enumerate(prediction_errors) if
                                    value == minimum_error]

        # logging.info(f'\n\t\tMinimum error indices: {minimum_error_indices}, minimum error: {minimum_error}')

        # Now that we have found the rules with the minimum error, we want to minimize the number of incoming nodes
        # (Find the simplest rule that explains the data)
        
        # Create a list of the number of incoming nodes for each of the minimum rules
        num_incoming_nodes_list = []
        for index in minimum_error_indices: 
            num_incoming_nodes = 0
            if 'A' in rules[index][2]:
                num_incoming_nodes += 1
            if 'B' in rules[index][2]:
                num_incoming_nodes += 1
            if 'C' in rules[index][2]:
                num_incoming_nodes += 1
            num_incoming_nodes_list.append(num_incoming_nodes)
        # logging.info(f'\t\t\tNum incoming nodes list: {num_incoming_nodes_list}')

        # Find the maximum number of incoming nodes
        max_incoming_nodes = max(num_incoming_nodes_list)

        # Compare to see if the current node has the same number of nodes as the max
        for i, index in enumerate(minimum_error_indices): 

            num_incoming_nodes = num_incoming_nodes_list[i]

            # logging.info(f'\n\t\t\tRule: {rules[index]}')
            # logging.info(f'\t\t\tminimum error index: {index}')
            # logging.info(f'\t\t\t\tmax_incoming_nodes = {max_incoming_nodes}')
            # logging.info(f'\t\t\t\tnum_incoming_nodes = {num_incoming_nodes}')

            def append_best_rule(index):
                # logging.info(f'\t\t\t\tbest_rule: {rules[index]}, Index: {index}, error: {prediction_errors[index]}')
                most_incoming_node_rules.append(rules[index])
                best_rule_indices.append(index)
                best_rule_errors.append(prediction_errors[index])

            # Prioritize nodes with 1 incoming rule
            if max_incoming_nodes == 1 and num_incoming_nodes == 1:
                append_best_rule(index)
            # If there are no minimum error rules with one incoming node, look for rules with two incoming nodes
            elif max_incoming_nodes == 2 and num_incoming_nodes == 2:
                append_best_rule(index)
            # If no two incoming node rules, look for rules with three incoming nodes
            elif max_incoming_nodes == 3 and num_incoming_nodes == 3:
                append_best_rule(index)
        
        # logging.info(f'\t\t\tMost incoming node rules: {most_incoming_node_rules}')

        # Find the rules with the fewest 'and' statements
        min_rules = min([rule.count("and") for rule in most_incoming_node_rules])

        # logging.info(f'\t\tmost_incoming_node_rules: {most_incoming_node_rules}')
        for i, rule in enumerate(most_incoming_node_rules):
            if rule.count("and") == min_rules:
                best_rules.append((most_incoming_node_rules[i], best_rule_indices[i], best_rule_errors[i]))
        
        # logging.info('\t\tbest rules:')
        # for i in best_rules:
        #     logging.info(f'\t\t\t{i[0]}')
        #     logging.info(f'\t\t\tError = {i[2]}')

    return best_rules

def calculate_error(node, predicted_rule, dataset):
    # logging.info(f'CustomDeap calculate_error:')
    # Get the row in the dataset for the node being evaluated
    node_evaluated = dataset[node.index]
    # logging.info(f'\tNode {node.name}, Dataset index {node.index}')
    
    # Get the dataset row indices for the incoming nodes included in this rule
    # logging.info(f'\tPredicted rule: {predicted_rule}')
    incoming_node_indices = [predecessor_index for predecessor_index in node.predecessors]
    # logging.info(f'\tIncoming node indices: {incoming_node_indices}')

    # Initialize A, B, C to False by default (adjust according to what makes sense in context)
    A, B, C = (False,) * 3
    
    data = {}

    # Map dataset values to A, B, C based on their indices
    if len(incoming_node_indices) > 0:
        data['A'] = dataset[incoming_node_indices[0]]
    if len(incoming_node_indices) > 1:
        data['B'] = dataset[incoming_node_indices[1]]
    if len(incoming_node_indices) > 2:
        data['C'] = dataset[incoming_node_indices[2]]
    if len(incoming_node_indices) > 3:
        data['D'] = dataset[incoming_node_indices[3]]

    def evaluate_expression(var_data, expression):
        # logging.info(f'\tEvaluate_expression function:')
        # logging.info(f'\t\tvar_data = {var_data}')
        # logging.info(f'\t\texpression = {expression}')

        def eval_func(*args):
            local_vars = {name: arg for name, arg in zip(var_data.keys(), args)}
            # logging.info(f'\t\texpression: {expression}, local_vars = {local_vars}')
            return eval(expression, {}, local_vars)
        
        vectorized_eval = np.vectorize(eval_func)

        # Prepare the argument list in the correct order, corresponding to the keys in var_data
        arg_list = [var_data[key] for key in var_data.keys()]
        return vectorized_eval(*arg_list)

    # logging.info(f'\tdata = {data}')
    predicted = evaluate_expression(data, predicted_rule)

    # # Evaluate the rule using numpy's vectorize function
    # predicted = np.vectorize(eval(predicted_rule))
    
    # logging.info(f'Predicted: {predicted}')
    difference = np.sum(predicted != node_evaluated)
    # logging.info(f'\tDifferences: {difference}')
    # logging.info(f'\tError = {difference / len(predicted) * 100}%')

    count = len(predicted)

    return difference, count



