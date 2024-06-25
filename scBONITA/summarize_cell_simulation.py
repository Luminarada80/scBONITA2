import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pickle
from setup.user_input_prompts import *
import logging
from heatmap import create_heatmap
import seaborn as sns
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import numexpr as ne
import random
from argparse import ArgumentParser
import statistics

from file_paths import file_paths

if __name__ == '__main__':

    logging.basicConfig(format='%(message)s', level=logging.INFO)

    parser = ArgumentParser()

    parser.add_argument(
        '--dataset_name',
        type=str,
        required=True,
        help='Number of genes to generate'
    )

    parser.add_argument(
        '--network_name',
        type=str,
        required=True,
        help='Number of cells to generate'
    )

    results = parser.parse_args()

    dataset_name = getattr(results, 'dataset_name')
    network_name = getattr(results, 'network_name')

    simulation_result_dir = f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/text_files/'

    gene_expr_dict = {}

    for file in os.listdir(simulation_result_dir):
        file_path = os.path.join(simulation_result_dir, file)
        with open(file_path, 'r') as sim_file:
            for line in sim_file:
                line = line.strip().split(',')
                gene_name = line[0]
                gene_expression = [int(i) for i in line[1:]]
                
                if gene_name not in gene_expr_dict:
                    gene_expr_dict[gene_name] = []
                
                gene_expr_dict[gene_name].append(gene_expression)

    gene_avg_expr = {}

    for gene, simulation_results in gene_expr_dict.items():

        if gene not in gene_avg_expr:
            gene_avg_expr[gene] = []

        transposed_data = list(map(list, zip(*simulation_results)))

        for i in transposed_data:
            gene_avg_expr[gene].append(statistics.mean(i))
    
    # Convert the dictionary to a DataFrame
    df = pd.DataFrame(gene_avg_expr)

    # Transpose the DataFrame
    df = df.transpose()

    # Create the heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(df, cmap='Greys')
    plt.title('Gene Expression Heatmap')
    plt.xlabel('Sample')
    plt.ylabel('Gene')
    plt.show()