import random
import numpy as np
import re
from create_test_network import CreateTestNetwork
import networkx as nx
import time
from datetime import timedelta
from alive_progress import alive_bar

def generate_random_state(length):
    return np.random.randint(2, size=length)

def evaluate_expression(expression, state):
    parsed_expression = re.sub(r'\b\d+\b', lambda match: f'{state[int(match.group(0))]}', expression)
    try:
        return int(eval(parsed_expression))
    except Exception as e:
        raise ValueError(f"Error evaluating expression: {expression}") from e

def evaluate_rule(rules, state):
    return [evaluate_expression(rule, state) for rule in rules]

def generate_states(rules, num_cells, num_genes):
    cells = []
    max_iterations = 5
    max_initial_states_tried = 5
    
    for _ in range(num_cells):
        initial_state = generate_random_state(num_genes)
        iterations = 0
        initial_states_tried = 0

        while iterations <= max_iterations:
            new_genes = evaluate_rule(rules, initial_state)
            if np.array_equal(new_genes, initial_state) or initial_states_tried > max_initial_states_tried:
                break
            if iterations >= max_iterations:
                initial_state = generate_random_state(num_genes)
                iterations = 0
                initial_states_tried += 1
            else:
                initial_state = new_genes
                iterations += 1

        cells.append(initial_state)
    
    return np.array(cells)

def preprocess_rules(rules):
    return [re.sub(r'\b\d+\b', lambda match: f'state[{int(match.group(0))}]', rule) for rule in rules]

def evaluate_preprocessed_expression(parsed_expression, state):
    try:
        return int(eval(parsed_expression))
    except Exception as e:
        raise ValueError(f"Error evaluating expression: {parsed_expression}") from e
    
def vectorized_evaluate_rules(rules, dataset):
    preprocessed_rules = preprocess_rules(rules)
    num_cells = dataset.shape[1]
    num_genes = dataset.shape[0]

    # Initialize a matrix to store the results of rule evaluations
    results = np.zeros((num_genes, num_cells), dtype=int)

    # Evaluate rules for each cell
    for gene_index, parsed_rule in enumerate(preprocessed_rules):
        for cell in range(num_cells):
            state = dataset[:, cell]
            results[gene_index, cell] = evaluate_preprocessed_expression(parsed_rule, state)

    return results

def check_rules(dataset, rules):
    results = vectorized_evaluate_rules(rules, dataset)
    mismatches = np.sum(results != dataset)
    return mismatches

def generate_network(num_genes, network_object):
    network_object.graph = network_object.create_network(num_genes)
    network = network_object.graph
    rule_dict = network_object.generate_network_rules()

    return rule_dict, network

def generate_dataset(chunks, ruleset, num_cells, num_genes):
    matrix_chunks = []
    total_mismatches = 0
    
    with alive_bar(chunks) as bar:
        for i in range(chunks):
            
            cells = generate_states(ruleset, num_cells, num_genes)
            matrix_chunk = cells.T
            total_mismatches += check_rules(matrix_chunk, ruleset)
                
            matrix_chunks.append(matrix_chunk)
            bar()
            
    return np.hstack(matrix_chunks), total_mismatches

def simulate_network():
    num_genes = 50
    num_cells = 20
    chunks = 100

    print('Generating simulated data...')
    print(f'\tCreating network')
    print(f'\tAttempting to generate dataset')

    test_network = CreateTestNetwork(num_genes=num_genes)
    
    rule_dict, network = generate_network(num_genes, test_network)
    ruleset = [rule.replace('Gene', '') for _, rule in rule_dict.items()]
    matrix, total_mismatches = generate_dataset(chunks, ruleset, num_cells, num_genes)
    percent_error = round(total_mismatches / (num_cells * chunks * num_genes) * 100)

    print(f'Number of mismatches = {total_mismatches}')
    print(f'Percent rule mismatch = {percent_error}')

    rules_filename = f"/home/emoeller/github/scBONITA/network_simulation_data/test_network_rules_{num_genes}g_{num_cells*chunks}c.txt"
    network_filename = f"/home/emoeller/github/scBONITA/network_simulation_data/test_network_{num_genes}g_{num_cells*chunks}c.graphml"
    data_filename = f"/home/emoeller/github/scBONITA/network_simulation_data/test_data_file_{num_genes}g_{num_cells*chunks}c.csv"
    percent_error_file_path = f'/home/emoeller/github/scBONITA/network_simulation_data/percent_error_log.txt'

    nx.write_graphml(network, network_filename)

    with open(rules_filename, 'w') as ruleset_file:
        for target, rule in rule_dict.items():
            line = f'{target} = {rule}\n'
            ruleset_file.write(line)
    
    mutation_rate = 0.1 
    num_mutations = int(num_genes * num_cells * mutation_rate)
    rows_to_mutate = np.random.choice(num_genes, num_mutations, replace=True)
    cols_to_mutate = np.random.choice(num_cells * chunks, num_mutations, replace=True)
    matrix[rows_to_mutate, cols_to_mutate] = 1 - matrix[rows_to_mutate, cols_to_mutate]

    row_labels = ["Gene" + str(i) for i in range(num_genes)]
    labeled_matrix = np.column_stack((row_labels, matrix))
    column_labels = [""] + ["Cell" + str(i + 1) for i in range(num_cells * chunks)]
    labeled_matrix_with_columns = np.vstack([column_labels, labeled_matrix])

    np.savetxt(data_filename, labeled_matrix_with_columns, delimiter=',', fmt='%s')

    print(f'Number of cells: {num_cells * chunks}')
    print(f'Number of genes: {num_genes}')
    print(f'Matrix shape = {labeled_matrix_with_columns.shape}')
    print(f'CSV dataset file "{data_filename}" created')
    print(f'Network graphml file "{network_filename}" created')
    print(f'Rules text file "{rules_filename}" created')

    with open(percent_error_file_path, 'a') as percent_error_file:
        percent_error_file.write(f'network_{num_genes}g_{num_cells*chunks}c\t->\t{percent_error}')


if __name__ == '__main__':
    start_time = time.time()
    simulate_network()
    end_time = time.time()
    elapsed_time = end_time - start_time
    formatted_time = str(timedelta(seconds=round(elapsed_time)))

    print(f'Time Elapsed: {formatted_time}')