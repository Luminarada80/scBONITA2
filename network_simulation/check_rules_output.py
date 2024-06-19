import re
import statistics
import random
import matplotlib.pyplot as plt
<<<<<<< HEAD
<<<<<<< HEAD
=======
from scipy.stats import chi2_contingency
import numpy as np
>>>>>>> 8d47ab02881342d43c2a52c7151cbc2cdbe854aa
=======
from scipy.stats import chi2_contingency
import numpy as np
>>>>>>> 8d47ab02881342d43c2a52c7151cbc2cdbe854aa
from alive_progress import alive_bar

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
        for line in rule_file:
            if ' = ' in line:
                rule = line.strip().split(' = ')[1]
                ruleset.append(rule)
    return ruleset

<<<<<<< HEAD
<<<<<<< HEAD
num_genes = 100
num_cells = 5000

true_ruleset_path = f'/home/emoeller/github/scBONITA/scBONITA/network_rules_{num_genes}_genes_{num_cells}_cells.txt'
test_ruleset_path = f'/home/emoeller/github/scBONITA/scBONITA/rules_output/test_data_{num_genes}_genes_rules/test_network_{num_genes}_genes_{num_cells}_cells.graphml_test_data_{num_genes}_genes_ind_1_rules.txt'
=======
=======
>>>>>>> 8d47ab02881342d43c2a52c7151cbc2cdbe854aa
num_genes = 50
num_cells = 2000
per_err = 44

true_ruleset_path = f'/home/emoeller/github/scBONITA2/scBONITA/test_network_rules_{num_genes}g_{num_cells}c_{per_err}e.txt'
test_ruleset_path = f'/home/emoeller/github/scBONITA2/scBONITA/rules_output/test_data_{num_genes}g_{num_cells}c_{per_err}e_rules/test_network_{num_genes}g_{num_cells}c_{per_err}e.graphml_test_data_{num_genes}g_{num_cells}c_{per_err}e_ind_1_rules.txt'
<<<<<<< HEAD
>>>>>>> 8d47ab02881342d43c2a52c7151cbc2cdbe854aa
=======
>>>>>>> 8d47ab02881342d43c2a52c7151cbc2cdbe854aa

true_ruleset = parse_ruleset(true_ruleset_path)
test_ruleset = parse_ruleset(test_ruleset_path)

num_trials = 1000
sim_steps = 100
num_genes = len(true_ruleset)

# Compare rulesets
percent_matches = []

<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
>>>>>>> 8d47ab02881342d43c2a52c7151cbc2cdbe854aa
TP = 0
TN = 0
FP = 0
FN = 0

<<<<<<< HEAD
>>>>>>> 8d47ab02881342d43c2a52c7151cbc2cdbe854aa
=======
>>>>>>> 8d47ab02881342d43c2a52c7151cbc2cdbe854aa
with alive_bar(num_trials) as bar:
    for i in range(num_trials):
        state = [random.choice([0,1]) for i in true_ruleset]
        # print(f'State {state}')
        state_matches = []
        for j in range(sim_steps):
            # print(f'\tSimulation step {j}')
            match_percentage, results1, results2 = compare_rulesets(true_ruleset, test_ruleset, state)
<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
>>>>>>> 8d47ab02881342d43c2a52c7151cbc2cdbe854aa

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

<<<<<<< HEAD
>>>>>>> 8d47ab02881342d43c2a52c7151cbc2cdbe854aa
=======
>>>>>>> 8d47ab02881342d43c2a52c7151cbc2cdbe854aa
            state_matches.append(match_percentage)
            state = results1
            # print(f'\tPercent Match = {match_percentage}')
        percent_matches.append(statistics.mean(state_matches))
        bar()
<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
>>>>>>> 8d47ab02881342d43c2a52c7151cbc2cdbe854aa

precision = TP / (TP + FP)
recall = TP / (TP + FN)

print(f'True Positives (TP): {TP / (num_trials * num_genes * sim_steps) * 100}%')
print(f'True Negatives (TN): {TN / (num_trials * num_genes * sim_steps) * 100}%')
print(f'False Positives (FP): {FP / (num_trials * num_genes * sim_steps) * 100}%')
print(f'False Negatives (FN): {FN / (num_trials * num_genes * sim_steps) * 100}%\n')
print(f'Precision = {precision}')
print(f'Recall = {recall}')

contingency_table = np.array([[TP, FN], [FN, TN]])

chi2, p, dof, ex = chi2_contingency(contingency_table)

print(f'Chi-square value = {chi2}, p = {p}')
print(f'DoF = {dof}, Expected frequencies = {ex}\n')
<<<<<<< HEAD
>>>>>>> 8d47ab02881342d43c2a52c7151cbc2cdbe854aa
=======
>>>>>>> 8d47ab02881342d43c2a52c7151cbc2cdbe854aa
        
avg = statistics.mean(percent_matches)
stdev = statistics.stdev(percent_matches)
minimum = min(percent_matches)
maximum = max(percent_matches)

print(f"Average Percent Match: {avg}%")
print(f'Stdev = {round(stdev,2)}%')
print(f'Min = {minimum}%')
print(f'Max = {maximum}%')

plt.figure(figsize=(12, 6))
plt.boxplot(percent_matches, vert=True, patch_artist=True)
plt.title(f'Percent of scBONITA rules matching True rules over {num_trials} trials for a pathway with {num_genes} genes')
plt.ylabel('Percentage')
plt.xticks([])  # Remove the x-axis tick (1)

# Adding legend with statistics
legend_text = (
    f"Avg = {avg:.2f}%\n"
    f"Stdev = {stdev:.1f}%\n"
    f"Min = {minimum}%\n"
    f"Max = {maximum}%"
)
plt.figtext(0.85, 0.5, legend_text, verticalalignment='center', fontsize=10, bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 10})

plt.grid(True)
plt.ylim((0,100))
plt.subplots_adjust(right=0.8)
plt.show()