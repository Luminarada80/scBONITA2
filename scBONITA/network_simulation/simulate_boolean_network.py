import random as rand
import itertools
import os
import time
import csv
import networkx as nx

def get_truth_table():
    inputs = [(1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
              (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)]
    outputs = rand.choices([0, 1], k=8)
    return dict(zip(inputs, outputs))

def truth_table_to_rule(truth_table):
    expressions = []
    for inputs, output in truth_table.items():
        if output == 1:
            term = []
            for i, val in enumerate(inputs):
                term.append(f"{'NOT ' if val == 0 else ''}A{i+1}")
            expressions.append(f"({' AND '.join(term)})")
    return ' OR '.join(expressions)

def apply_rule(neighbor_states, truth_table):
    # Use the neighbor states directly
    return truth_table[tuple(neighbor_states)]

def print_state(matrix):
    for row in matrix:
        print(' '.join(map(str, row)))

def create_graph(nodes, connections, truth_tables):
    G = nx.DiGraph()
    for node in range(nodes):
        gene_name = f'gene{node+1}'
        G.add_node(gene_name, rule=gene_name)

    for node, conn in enumerate(connections):
        gene_name = f'gene{node+1}'
        for target in conn:
            target_gene = f'gene{target+1}'
            interaction = 'a'
            G.add_edge(gene_name, target_gene, interaction=interaction)    
    
    return G

def output_rules_to_file(node_rules, filename="./boolean_network_conditions/node_rules.txt"):
    with open(filename, "w") as file:
        for node, rule in node_rules.items():
            file.write(f"{node}: {rule}\n")

def output_simulation_to_csv(simulation, nodes, filename="./boolean_network_conditions/network_simulation.csv"):
    with open(filename, "w", newline="") as file:
        writer = csv.writer(file)
        
        # Write header
        header = [''] + [f'cell{i+1}' for i in range(len(simulation))]
        writer.writerow(header)

        # Transpose the simulation matrix to write each gene's state across cells
        transposed_simulation = list(zip(*simulation))
        
        # Write each gene's states
        for gene_idx, gene_states in enumerate(transposed_simulation):
            row = [f'gene{gene_idx+1}'] + list(gene_states)
            writer.writerow(row)


nodes = 25
iterations = 1000
truth_tables = [get_truth_table() for _ in range(nodes)]
initial_state = rand.choices([1, 0], k=nodes)
connections = rand.choices(list(itertools.combinations(range(nodes), 3)), k=nodes)
network_state = [initial_state]

# Output the rules for each node to a text file
node_rules = {f'gene{node+1}': truth_table_to_rule(truth_tables[node]) for node in range(nodes)}
output_rules_to_file(node_rules)

for i in range(iterations - 1):
    # os.system('clear')
    state = network_state[-1]

    # Apply the rule for each node using its corresponding truth table
    new_state = []
    for node in range(nodes):
        neighbor_states = [state[index] for index in connections[node]]
        new_state.append(apply_rule(neighbor_states, truth_tables[node]))

    network_state.append(new_state)
    # print_state(network_state)

# Output the simulation results to a CSV file
output_simulation_to_csv(network_state, nodes)

# Create the graph and write to GraphML
G = create_graph(nodes, connections, truth_tables)
nx.write_graphml(G, "./boolean_network_conditions/random_boolean_network.graphml")
