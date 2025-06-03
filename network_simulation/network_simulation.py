from create_test_network import CreateTestNetwork
import random
import numpy as np
import argparse

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate scBONITA rules vs. simulated rules.")
    parser.add_argument(
        "--dataset_name",
        type=str,
        required=True,
        help="Path to scBONITA output rules file (one gene per line)."
    )
    parser.add_argument(
        "--num_genes",
        type=str,
        required=True,
        help="Path to scBONITA output rules file (one gene per line)."
    )
    return parser.parse_args()

def generate_random_state(length):
    return [random.randint(0, 1) for _ in range(length)]

def evaluate_rule(rules, state):
    """
    Evaluate a single rule based on the current state.

    Args:
    rule (list): The rule to be evaluated.
    state (list): The current state of the values.

    Returns:
    int: The result of the rule evaluation.
    """

    ruleset = []
    for rule in rules:
        result = None
        operator = None
        for item in rule:
            if item == 'and':
                operator = 'and'
            elif item == 'or':
                operator = 'or'
            elif item == 'not':
                operator = 'not'
            else:
                value = state[int(item)]
                if operator == 'not':
                    value = not value
                    operator = None

                if result is None:
                    result = value
                else:
                    if operator == 'and':
                        result = result and value
                    elif operator == 'or':
                        result = result or value
        ruleset.append(1 if result else 0)
    return ruleset

def generate_states(rules, num_cells, num_genes):
    cells = []
    for cell in range(num_cells):
        initial_state = generate_random_state(num_genes)
        # print(f'cell{cell} initial_state: \t{initial_state}')

        iterations = 0
        initial_states_tried = 0

        while True:
            # print(f'Cell {cell} Iteration: {iterations} Initial State: {initial_states_tried}')
            new_genes = evaluate_rule(rules, initial_state)
            # print(f'\tstate: \t\t{new_genes}')
            if iterations > 10:
                initial_state = generate_random_state(num_genes)
                iterations = 0
                initial_states_tried += 1
            if new_genes == initial_state or initial_states_tried > 5:
                break
            initial_state = new_genes

            iterations += 1
        # print(f'\tcell{cell} = \t{initial_state}\n')
        cells.append(initial_state)
    return cells

def check_rules(dataset, rules):
    """
    Check if the dataset follows the given rules.

    Args:
    dataset (numpy array): The dataset containing the states of the genes.
    rules (list of lists): The rules for each gene.

    Returns:
    bool: True if the dataset follows the rules, False otherwise.
    """
    num_cells = dataset.shape[1]
    num_genes = dataset.shape[0]
    mismatches = []

    for cell in range(num_cells):
        for gene_index, rule in enumerate(rules):
            result = evaluate_generated_dataset(rule, dataset[:, cell])
            if result != dataset[gene_index, cell]:
                mismatches.append((gene_index, cell, result, dataset[gene_index, cell]))
    return mismatches

def evaluate_generated_dataset(rule, state):
    """
    Evaluate a single rule based on the current state.

    Args:
    rule (list): The rule to be evaluated.
    state (list): The current state of the values.

    Returns:
    int: The result of the rule evaluation.
    """
    result = None
    operator = None
    for item in rule:
        if item == 'and':
            operator = 'and'
        elif item == 'or':
            operator = 'or'
        elif item == 'not':
            operator = 'not'
        else:
            value = state[int(item)]
            if operator == 'not':
                value = not value
                operator = None

            if result is None:
                result = value
            else:
                if operator == 'and':
                    result = result and value
                elif operator == 'or':
                    result = result or value
    return 1 if result else 0

def parse_rules(rules_filename):
    # Read the rules for the test network and create a dataset
    with open(rules_filename, 'r') as datafile:
        ruleset = []
        for line in datafile:

            # Get the gene names
            gene_name = line.split(' = ')[0]
            gene_name = gene_name.strip('Gene')

            # Get and format the rules
            gene_rules = line.split(' = ')[1] # Specify the rules as being the right side of the '=' sign
            # print(gene_rules)
            gene_rules = gene_rules.split(' ') # Split elements based on spaces
            gene_rules = [i.strip('Gene') for i in gene_rules] # Strip 'Gene' from each element to get gene numbers
            gene_rules = [i.strip() for i in gene_rules] # Get rid of newline characters
            gene_rules = [i.lower() for i in gene_rules] # Set rules to lowercase

            ruleset.append(gene_rules) # Append the ruleset for the file

        return ruleset

def create_network(num_genes, network_filename, rules_filename):
    #Instantiate the network class
    test_network = CreateTestNetwork(num_genes=num_genes)
    
    # Export the graphml file for the network
    test_network.export_network_graphml(network_filename)
    
    # Export the network rules to a text file
    test_network.export_network_rules(filename=rules_filename)

    # Visualize the network
    # test_network.visualize_network()

    # Extract the rules from the rules file
    ruleset = parse_rules(rules_filename)

    return ruleset

def generate_dataset(chunks, ruleset, num_cells, num_genes):
    matrix_chunks = []

    num_tries = 0
    
    for _ in range(chunks):
        while True:
            # Generate Cells
            cells = generate_states(ruleset, num_cells, num_genes)

            # Put the simulated data into a matrix
            matrix_chunk = np.array(cells).T

            # Check if the dataset follows the rules
            mismatches = check_rules(matrix_chunk, ruleset)

            # Break if there are no mismatches, else try another set of data
            if not mismatches:
                num_tries = 0
                break
                
            else:
                num_tries += 1
                if num_tries == 100:
                    return False

        matrix_chunks.append(matrix_chunk)
            
    matrix = np.concatenate(matrix_chunks, axis = 1)

    return matrix
                
def simulate_network():
    args = parse_args()
    dataset_name = args.dataset_name
    
    num_genes = int(args.num_genes)
    rules_filename = f"/home/emoeller/github/scBONITA2/network_simulation/data/{dataset_name}_{num_genes}_rules.txt"
    network_filename = f"/home/emoeller/github/scBONITA2/input/custom_graphml_files/{dataset_name}_{num_genes}.graphml"
    data_filename = f"/home/emoeller/github/scBONITA2/input/{dataset_name}_{num_genes}.csv"

    # Generate a dataset based on the rules
    num_cells = 25
    chunks = 100

    # Main loop
    print(f'Generating simulated data...')
    while True:
        # Create the network and the rules
        ruleset = create_network(num_genes, network_filename, rules_filename)

        # Generate the data
        matrix = generate_dataset(chunks, ruleset, num_cells, num_genes)

        if isinstance(matrix, np.ndarray):
            break
    
    # Randomly choose indices to mutate
    mutation_rate = 0.1 
    num_mutations = int(num_genes * num_cells * mutation_rate)
    rows_to_mutate = np.random.choice(num_genes, num_mutations, replace=True)
    cols_to_mutate = np.random.choice(num_cells, num_mutations, replace=True)

    # Perform the mutations
    for row, col in zip(rows_to_mutate, cols_to_mutate):
        matrix[row, col] = 1 - matrix[row, col]  # Flip the bit

    # Add the column and row labels
    row_labels = ["Gene" + str(i) for i in range(matrix.shape[0])]
    labeled_matrix = np.column_stack((row_labels, matrix))
    column_labels = [""] + ["Cell" + str(i + 1) for i in range(matrix.shape[1])]
    labeled_matrix_with_columns = np.vstack([column_labels, labeled_matrix])
    # print(labeled_matrix_with_columns)

    # Save the numpy array to a CSV file without joining rows
    np.savetxt(data_filename, labeled_matrix_with_columns, delimiter=',', fmt='%s')

    print(f'Number of cells: {num_cells * chunks}')
    print(f'Number of genes: {num_genes}\n')

    print(f'Matrix shape = ({labeled_matrix_with_columns.shape[0]},{labeled_matrix_with_columns.shape[1]})')
    print(f'\tCSV dataset file "{data_filename}" created')
    print(f'\tNetwork graphml file "{network_filename}" created')
    print(f'\tRules text file {rules_filename} created')

if __name__ == '__main__':
   simulate_network()



