import numpy as np
import random

def check_rules(dataset, rules, raw_rule_text):
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
    total = 0
    mismatch_count = 0

    rule_error_dict = {}
    

    for cell in range(num_cells):
        for gene_index, rule in enumerate(rules):
            if raw_rule_text[gene_index] not in rule_error_dict:
                rule_error_dict[str(raw_rule_text[gene_index])] = 0

            result = check_logic(rule, dataset[:, cell])

            if result != dataset[gene_index, cell]:
                mismatch_count += 1
                rule_error_dict[str(raw_rule_text[gene_index])] += 1
        
            total += 1
    
    print(f'\nRules with mismatches:')
    for rule, mismatches in rule_error_dict.items():
        line_new = '{:>12}  {:>12}'.format(rule, mismatches)
        if mismatches != 0:
            print(f'\tMismatches: {mismatches}\tRule: {rule}')
    
    if mismatch_count == 0:
        print(f'\tNo mismatches')
    
    error =  mismatch_count / total * 100

    # print(f'\nMismatches = {mismatch_count} (total = {total})')
    # print(f'Error = {error}\n')

    return mismatch_count, total, error

def check_logic(rule, state):
    result = None
    operator = None
    print(f"Evaluating Rule: {rule} with State: {state}")  # Debugging output
    for item in rule:
        if item in ['and', 'or', 'not']:
            operator = item
        else:
            value = state[int(item)]
            if operator == 'not':
                value = not value
                operator = None  # Reset operator after use

            print(f"Item: {item}, Operator: {operator}, Value: {value}")  # Debugging output

            if result is None:
                result = value
            else:
                if operator == 'and':
                    result = result and value
                elif operator == 'or':
                    result = result or value

            print(f"Intermediate Result: {result}")  # Debugging output

    final_result = 1 if result else 0
    print(f"Final Result: {final_result}")  # Debugging output
    return final_result

def get_rules(path):
    raw_ruletext = []
    with open(path, 'r') as rules_file:
        ruleset = []
        for line in rules_file:

            # Get the gene names
            gene_name = line.split(' = ')[0]
            gene_name = gene_name.strip('Gene')

            # Get and format the rules
            try:
                gene_rules = line.split(' = ')[1] # Specify the rules as being the right side of the '=' sign
                # print(gene_rules)
                gene_rules = gene_rules.split(' ') # Split elements based on spaces
                gene_rules = [i.strip('Gene') for i in gene_rules] # Strip 'Gene' from each element to get gene numbers
                gene_rules = [i.strip() for i in gene_rules] # Get rid of newline characters
                gene_rules = [i.lower() for i in gene_rules] # Set rules to lowercase
                raw_ruletext.append(line.strip())
            except IndexError:
                continue

            ruleset.append(gene_rules) # Append the ruleset for the file
    return ruleset, raw_ruletext

def get_matrix(path):
    with open(path, 'r') as data_file:
        matrix = []
        file_lines = []
        for line in data_file:
            line = line.strip()
            line = line.split(',')
            file_lines.append(line)
        for entry in file_lines[1:]:
            integer_data = [int(i) for i in entry[1:]]
            matrix.append(integer_data)
    return matrix

def evaluate_scbonita_output():
    scbonita_rules_file_path = '/home/emoeller/github/scBONITA/scBONITA/rules_output/test_data_rules/test_network.graphml_test_data_ind_1_rules.txt'
    simulated_rules_file_path = './data/network_rules.txt'
    data_file_path = './data/test_data_file.csv'

    scbonita_ruleset, scbonita_rules = get_rules(scbonita_rules_file_path)
    simulated_ruleset, simulated_rules = get_rules(simulated_rules_file_path)

    data = get_matrix(data_file_path)

    matrix = np.array(data)

    print('scBONITA')
    scbonita_mismatch_count, scbonita_total, scbonita_error = check_rules(matrix, scbonita_ruleset, scbonita_rules)

    print('Simulated Rules')
    simulated_mismatch_count, simulated_total, simulated_error = check_rules(matrix, simulated_ruleset, simulated_rules)

    print(f'\nscBONITA Rules')
    print(f'Mismatches = {scbonita_mismatch_count} (total = {scbonita_total})')
    print(f'Error = {round(scbonita_error, 2)}%')

    print(f'\nSimulated Rules')
    print(f'Mismatches = {simulated_mismatch_count} (total = {simulated_total})')
    print(f'Error = {simulated_error}%')


if __name__ == '__main__':
    evaluate_scbonita_output()