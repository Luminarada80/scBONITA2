import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

# Get the files from the parent dir so that file_paths can be imported
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)

from file_paths import sim_data_file_paths

# Data from the provided table
data = {}

results_output_path = f'{sim_data_file_paths["sim_ruleset"]}'

num_genes = 10
num_cells = 2000

while num_genes <= 100:
    
    results_path = f'{results_output_path}/results_{num_genes}_{num_cells}.txt'
    if os.path.exists(results_path):

        print(f'Reading in simulated data results for {num_genes} genes x {num_cells} cells')

        if 'num_genes' not in data:
            data['num_genes'] = []
        data['num_genes'].append(num_genes)

        with open(results_path, 'r') as results_file:
            for line in results_file:
                line = line.strip().split('\t')
                data_name = line[0]
                data_value = float(line[1])
                
                if data_name not in data:
                    data[data_name] = []
                
                data[data_name].append(data_value)
    
    num_genes += 10

# Convert to DataFrame
df = pd.DataFrame(data)

# Prepare the data for box plots
box_plot_data = []
for i in range(len(df)):
    num_genes = df.iloc[i]["num_genes"]
    avg = df.iloc[i]["avg_percent_match"]
    stdev = df.iloc[i]["stdev"]
    min_val = df.iloc[i]["min"]
    max_val = df.iloc[i]["max"]
    true_positive = df.iloc[i]["TP"]
    true_negative = df.iloc[i]["TN"]
    false_positive = df.iloc[i]["FP"]
    false_negative = df.iloc[i]["FN"]
    precision = df.iloc[i]["precision"]
    recall = df.iloc[i]["recall"]
    true_avg_indegree = df.iloc[i]["true_avg_indegree"]
    scbonita_avg_indegree = df.iloc[i]["scbonita_avg_indegree"]
    true_percent_one_indegree = df.iloc[i]["true_percent_one_indegree"]
    true_percent_two_indegree = df.iloc[i]["true_percent_two_indegree"]
    true_percent_three_indegree = df.iloc[i]["true_percent_three_indegree"]
    scbonita_percent_one_indegree = df.iloc[i]["scbonita_percent_one_indegree"]
    scbonita_percent_two_indegree = df.iloc[i]["scbonita_percent_two_indegree"]
    scbonita_percent_three_indegree = df.iloc[i]["scbonita_percent_three_indegree"]
    chi_square = df.iloc[i]["Chi_square_value"]
    p_value = df.iloc[i]["p_value"]

    # Simulate the data for the box plot
    values = [avg - stdev, avg, avg + stdev]
    values = [min(max(val, min_val), max_val) for val in values]  # Ensure values are within min and max
    box_plot_data.append(values)

# Prepare data for bar plot
indegree_data = {
    "num_genes": df["num_genes"],
    "true_avg_indegree": df["true_avg_indegree"],
    "scbonita_avg_indegree": df["scbonita_avg_indegree"],
}

indegree_df = pd.DataFrame(indegree_data)

# Melt the DataFrame to long format for easier plotting
melted_df = indegree_df.melt(id_vars=["num_genes"], 
                             value_vars=["true_avg_indegree",
                                         "scbonita_avg_indegree"],
                             var_name="indegree_type", value_name="percent")

# Create subplots
fig, ax = plt.subplots(2, 1, figsize=(6, 8))

# Create the box plot
ax[0].boxplot(box_plot_data, positions=df["num_genes"], widths=5)
# ax[0].set_xticks(df["num_genes"])
ax[0].set_xticklabels(df["num_genes"], fontsize=10)
# ax[0].set_xlabel('Number of Genes', fontsize=10)
ax[0].set_ylabel('Accuracy', fontsize=10)
ax[0].set_title('scBONITA Rule Recovery Accuracy', fontsize=12)
ax[0].set_ylim(0, 100)
ax[0].grid(True)

# Create the bar plot
bar_width = 0.1
num_bars = len(melted_df["indegree_type"].unique())
indices = range(len(indegree_df["num_genes"].unique()))

for i, indegree_type in enumerate(melted_df["indegree_type"].unique()):
    subset = melted_df[melted_df["indegree_type"] == indegree_type]
    ax[1].bar([index + bar_width * i for index in indices], subset["percent"], bar_width, label=indegree_type)

ax[1].set_xlabel('Number of Genes', fontsize=10)
ax[1].set_ylabel('Percent Indegree', fontsize=10)
ax[1].set_title('Average Indegree Percentages by Number of Genes', fontsize=12)
ax[1].set_xticks([index + bar_width * (num_bars / 2) for index in indices])
ax[1].set_xticklabels(indegree_df["num_genes"].unique(), fontsize=10)
ax[1].set_ylim(0, 3)
ax[1].legend(fontsize=10)
ax[1].grid(True)

plt.tight_layout()
plt.show()