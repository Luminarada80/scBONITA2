import re
import statistics
import random
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency
import numpy as np
from alive_progress import alive_bar
import os
import sys

# Get the files from the parent dir so that file_paths can be imported
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)

from file_paths import file_paths, sim_data_file_paths

def evaluate_expression(expression, state):
    # Replace indices with their boolean values from the state
    def replace_indices(match):
        index = int(match.group(0).replace('Gene', ''))
        return str(bool(state[index]))
    
    # Use regex to find all indices in the expression and replace them
    parsed_expression = re.sub(r'\bGene\d+\b', replace_indices, expression)
    
    # Evaluate the parsed expression safely
    try:
        result = eval(parsed_expression)
    except Exception as e:
        raise ValueError(f"Error evaluating expression: {expression}") from e
    
    return result

def evaluate_ruleset(ruleset, state):
    results = []
    for rule in ruleset:
        result = evaluate_expression(rule, state)
        results.append(result)
    return results

def compare_rulesets(ruleset1, ruleset2, state):
    results1 = evaluate_ruleset(ruleset1, state)
    results2 = evaluate_ruleset(ruleset2, state)
    
    match_count = sum(1 for r1, r2 in zip(results1, results2) if r1 == r2)
    match_percentage = (match_count / len(ruleset1)) * 100
    
    return match_percentage, results1, results2

def parse_ruleset(rule_file_path):
    ruleset = []
    with open(rule_file_path, 'r') as rule_file:
        indegree = []
        for line in rule_file:
            rule_indegree = 0
            if ' = ' in line:
                rule = line.strip().split(' = ')[1]

                # Count the number of logic gates in the rule, indegree is one more than the number of gates
                indegree.append(sum(line.count(gate) for gate in ['and', 'or']) + 1)
                ruleset.append(rule)

    return ruleset, indegree


num_genes = 10
num_cells = 2000

while num_genes <= 100:
    true_ruleset_path = f'{sim_data_file_paths["sim_ruleset"]}/sim_rules_{num_genes}g_{num_cells}c.txt'
    test_ruleset_path = f'{file_paths["rules_output"]}/sim_data_{num_genes}g_{num_cells}c_rules/sim_network_{num_genes}g_{num_cells}c.graphml_sim_data_{num_genes}g_{num_cells}c_ind_1_rules.txt'
    results_output_path = f'{sim_data_file_paths["sim_ruleset"]}'
    
    if os.path.exists(true_ruleset_path) and os.path.exists(test_ruleset_path):

        print(f'Checking rules for simulated data with {num_genes} genes and {num_cells} cells')

        true_ruleset, true_indegree = parse_ruleset(true_ruleset_path)
        test_ruleset, test_indegree = parse_ruleset(test_ruleset_path)


        def per_indegree(indegree_list, indegree):
            return len([i for i in indegree_list if i == indegree]) / len(indegree_list) * 100

        num_trials = 1000
        sim_steps = 100
        num_genes = len(true_ruleset)

        # Compare rulesets
        percent_matches = []

        TP = 0
        TN = 0
        FP = 0
        FN = 0

        with alive_bar(num_trials) as bar:
            for i in range(num_trials):
                state = [random.choice([0,1]) for i in true_ruleset]
                # print(f'State {state}')
                state_matches = []
                for j in range(sim_steps):
                    # print(f'\tSimulation step {j}')
                    match_percentage, results1, results2 = compare_rulesets(true_ruleset, test_ruleset, state)

                    # True positive
                    for true, test in zip(results1, results2):
                        if true == 1 and test == 1:
                            TP += 1
                        
                        # True negative
                        elif true == 0 and test == 0:
                            TN += 1
                        
                        # False positive
                        elif true == 0 and test == 1:
                            FP += 1
                        
                        # False negative
                        elif true == 1 and test == 0:
                            FN += 1

                    state_matches.append(match_percentage)
                    state = results1
                    # print(f'\tPercent Match = {match_percentage}')
                percent_matches.append(statistics.mean(state_matches))
                bar()

        precision = TP / (TP + FP)
        recall = TP / (TP + FN)

        percent_true_positives = TP / (num_trials * num_genes * sim_steps) * 100
        percent_true_negatives = TN / (num_trials * num_genes * sim_steps) * 100
        percent_false_positives = FP / (num_trials * num_genes * sim_steps) * 100
        percent_false_negatives = FN / (num_trials * num_genes * sim_steps) * 100

        true_avg_indegree = statistics.mean(true_indegree)
        scbonita_avg_indegree = statistics.mean(test_indegree)

        true_percent_one_indegree = per_indegree(true_indegree, 1)
        true_percent_two_indegree = per_indegree(true_indegree, 2)
        true_percent_three_indegree = per_indegree(true_indegree, 3)

        scbonita_percent_one_indegree = per_indegree(test_indegree, 1)
        scbonita_percent_two_indegree = per_indegree(test_indegree, 2)
        scbonita_percent_three_indegree = per_indegree(test_indegree, 3)

        # print(f'True Positives (TP): {percent_true_positives}%')
        # print(f'True Negatives (TN): {percent_true_negatives}%')
        # print(f'False Positives (FP): {percent_false_positives}%')
        # print(f'False Negatives (FN): {percent_false_negatives}%\n')
        # print(f'Precision = {precision}')
        # print(f'Recall = {recall}\n')

        # print(f'True Ruleset indegree:')
        # print(f'\t1 indegree: {true_percent_one_indegree}')
        # print(f'\t2 indegree: {true_percent_two_indegree}')
        # print(f'\t2 indegree: {true_percent_three_indegree}')

        # print(f'scBONITA Ruleset indegree:')
        # print(f'\t1 indegree: {scbonita_percent_one_indegree}')
        # print(f'\t2 indegree: {scbonita_percent_two_indegree}')
        # print(f'\t2 indegree: {scbonita_percent_three_indegree}')

        contingency_table = np.array([[TP, FN], [FN, TN]])

        chi2, p, dof, ex = chi2_contingency(contingency_table)

        # print(f'Chi-square value = {chi2}, p = {p}')
        # print(f'DoF = {dof}, Expected frequencies = {ex}\n')
                
        avg = statistics.mean(percent_matches)
        stdev = statistics.stdev(percent_matches)
        minimum = min(percent_matches)
        maximum = max(percent_matches)

        # print(f"Average Percent Match: {avg}%")
        # print(f'Stdev = {round(stdev,2)}%')
        # print(f'Min = {minimum}%')
        # print(f'Max = {maximum}%')

        # plt.figure(figsize=(12, 6))
        # plt.boxplot(percent_matches, vert=True, patch_artist=True)
        # plt.title(f'Percent of scBONITA rules matching True rules over {num_trials} trials for a pathway with {num_genes} genes')
        # plt.ylabel('Percentage')
        # plt.xticks([])  # Remove the x-axis tick (1)

        # # Adding legend with statistics
        # legend_text = (
        #     f"Avg = {avg:.2f}%\n"
        #     f"Stdev = {stdev:.1f}%\n"
        #     f"Min = {minimum}%\n"
        #     f"Max = {maximum}%"
        # )
        # plt.figtext(0.85, 0.5, legend_text, verticalalignment='center', fontsize=10, bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})

        # plt.grid(True)
        # plt.ylim((0,100))
        # plt.subplots_adjust(right=0.8)
        # plt.show()

        with open(f'{results_output_path}/results_{num_genes}_{num_cells}.txt', 'w') as results_file:
            results_file.write(f'avg_percent_match\t{avg}\n')
            results_file.write(f'stdev\t{round(stdev,2)}\n')
            results_file.write(f'min\t{minimum}\n')
            results_file.write(f'max\t{maximum}\n')
            
            results_file.write(f'TP\t{percent_true_positives}\n')
            results_file.write(f'TN\t{percent_true_negatives}\n')
            results_file.write(f'FP\t{percent_false_positives}\n')
            results_file.write(f'FN\t{percent_false_positives}\n')

            results_file.write(f'precision\t{precision}\n')
            results_file.write(f'recall\t{recall}\n')

            results_file.write(f'true_avg_indegree\t{true_avg_indegree}\n')
            results_file.write(f'scbonita_avg_indegree\t{scbonita_avg_indegree}\n')

            results_file.write(f'true_percent_one_indegree\t{true_percent_one_indegree}\n')
            results_file.write(f'true_percent_two_indegree\t{true_percent_two_indegree}\n')
            results_file.write(f'true_percent_three_indegree\t{true_percent_three_indegree}\n')

            results_file.write(f'scbonita_percent_one_indegree\t{scbonita_percent_one_indegree}\n')
            results_file.write(f'scbonita_percent_two_indegree\t{scbonita_percent_two_indegree}\n')
            results_file.write(f'scbonita_percent_three_indegree\t{scbonita_percent_three_indegree}\n')

            results_file.write(f'Chi_square_value\t{chi2}\n')
            results_file.write(f'p_value\t{p}\n')

    else:
        print(f'Rule files do not exist for network of size {num_genes} genes x {num_cells} cells')

    num_genes += 10