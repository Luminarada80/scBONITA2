import numpy as np
import csv
import pandas as pd

predicted_rules_filename = '../scBONITA/output_rules'
generated_rules_filename = '../data/network_rules.txt'

data_file = '../data/test_data_file.csv'


with open(data_file, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    data = list(reader)
    for row_number, row in enumerate(data): 
        for column_number, value in enumerate(row):
            if row_number > 0 and column_number > 0:
                data[row_number][column_number] = int(data[row_number][column_number])
    
    data_matrix = pd.DataFrame(data[1:], columns=data[0])
    data_matrix.set_index(data_matrix.columns[0], inplace=True)

    print(f'\nData\n\t{data_matrix}')

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

    print(f'\nCells: {num_cells}, Genes: {num_genes}\n')
    mismatches = []

    for cell in range(num_cells):
        print('\n')
        print(data_matrix)
        print(dataset.iloc[:, cell].values)
        for gene_index, rule in enumerate(rules):

            result = evaluate_generated_dataset(rule, dataset.iloc[:, cell].values)
            print(f'Cell{cell+1}, Gene{gene_index}, Result {result}, Dataset {dataset.iloc[gene_index, cell]}')
            if result != dataset.iloc[gene_index, cell]:
                print(f'\tCell{cell+1}, Gene{gene_index}: Expected {result}, Found {dataset.iloc[gene_index, cell]}')
                mismatches.append((gene_index, cell, result, dataset.iloc[gene_index, cell]))
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


def check_accuracy(rule_file):
    with open(rule_file, 'r') as datafile:
        dataset_rules = ''
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

        print(f'\nRuleset:\n\t{ruleset}')
    
        mismatches = check_rules(data_matrix, ruleset)

        # if mismatches:
        #     print("The dataset does not follow the rules. Mismatches found:")
        #     for mismatch in mismatches:
        #         gene_index, cell, expected, actual = mismatch
        #         print(f"Cell {cell + 1}, Gene{gene_index}: Expected {expected}, Found {actual}")
            
        # else:
        #     print("The dataset follows the rules.")

check_accuracy(generated_rules_filename)