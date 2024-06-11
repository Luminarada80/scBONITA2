import pandas as pd
import numpy as np
import random

class GenerateDataset:
    def __init__(self, rules_file, num_cells):
        self.rules_file = rules_file
        self.num_cells = num_cells
        self.rules = self.read_rules()
        self.gene_dictionary, self.genes_with_values = self.build_gene_dictionary()

    def read_rules(self):
        rules = {}
        with open(self.rules_file, 'r') as file:
            for line in file:
                parts = line.strip().split(' = ')
                gene = parts[0]
                expression = parts[1]
                rules[gene] = expression
        return rules

    def build_gene_dictionary(self):
        # Initialize a DataFrame with random 0s and 1s
        gene_dictionary = {}
        for gene_name, incoming_nodes in self.rules.items():
            incoming_node_name = [i for i in self.rules[gene_name].split(' ') if "Gene" in i]
            rules = [i for i in self.rules[gene_name].split(' ') if not "Gene" in i]
            formatted_rules = self.combine_not(rules)

            gene_dictionary[gene_name] = [incoming_node_name, formatted_rules]

        genes_with_values = []
        for gene in gene_dictionary:
            if len(gene_dictionary[gene][1]) == 0:
                genes_with_values.append(gene)
                gene_dictionary[gene].append(random.choice([0, 1]))

        return gene_dictionary, genes_with_values


    def generate_cell_data(self):
        print(self.gene_dictionary)
        for node in self.genes_with_values:
            print(node)
        for gene in self.gene_dictionary:
            if len(self.gene_dictionary[gene][0]) == 1:
                incoming_nodes = self.gene_dictionary[gene][0]
                for incoming_node in incoming_nodes:
                    if incoming_node in self.genes_with_values:
                        print(f'{gene}: incoming_node {incoming_node}, incoming_nodes = {incoming_nodes}\n')
                        expression = self.gene_dictionary[incoming_node][2]
                        self.gene_dictionary[gene].append(expression)
                        print(self.gene_dictionary[gene])
        print('\n')

    def combine_not(self, operators):
        i = 0
        while i < len(operators):
            if operators[i] == 'not' and i + 1 < len(operators):
                operators[i] = 'not ' + operators[i + 1]
                del operators[i + 1]
            else:
                i += 1
        return operators