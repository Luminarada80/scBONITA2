import networkx as nx
import random
import itertools
import matplotlib.pyplot as plt

class CreateTestNetwork():
    def __init__(self, num_genes):
        self.num_genes = num_genes  # number of nodes in the network
        self.graph = self.create_network()

    def create_network(self):
        G_genes = nx.DiGraph()

        # Add gene nodes using the gene name as the node identifier
        for i in range(self.num_genes):
            gene_name = f'Gene{i}'
            G_genes.add_node(gene_name, name=gene_name)

        # Add random edges between genes
        for node_num in range(self.num_genes):
            gene_name = f'Gene{node_num}'
            num_connections = random.randint(1, 3)  # Random number of connections
            possible_targets = [f'Gene{j}' for j in range(self.num_genes) if j != node_num and G_genes.in_degree(f'Gene{j}') < 3]
            if possible_targets:
                connections = random.sample(possible_targets, min(num_connections, len(possible_targets)))
                for target in connections:
                    interaction = 'a' if random.random() < 0.8 else 'i'
                    G_genes.add_edge(gene_name, target, interaction=interaction)

        # Add self-loops to genes without predecessors
        for node_num in range(self.num_genes):
            gene_name = f'Gene{node_num}'
            if len(list(G_genes.predecessors(gene_name))) == 0:
                G_genes.add_edge(gene_name, gene_name, interaction='a')

        return G_genes

    def is_cyclic(self, graph, source, target):
        """
        Check if adding an edge from source to target creates a cycle.
        """
        if source == target:  # Allow self-loops
            return False
        
        visited = set()

        def visit(vertex):
            if vertex in visited:
                return False
            visited.add(vertex)
            for neighbor in graph.successors(vertex):
                if neighbor == target or not visit(neighbor):
                    return True
            visited.remove(vertex)
            return False

        return visit(source)
    
    def export_network_rules(self, filename):
        """Export the network rules to a text file."""
        with open(filename, 'w') as file:
            for node in self.graph.nodes():
                gene_connections = []
                operators = []
                not_operators = []
                for source, target, data in self.graph.in_edges(node, data=True):
                    gene_connections.append(self.graph.nodes[source]['name'])
                    not_operator = 'NOT ' if data['interaction'] == 'i' else ''
                    operator = str(random.choice([' AND ', ' OR ']))
                    not_operators.append(not_operator)
                    operators.append(operator)

                expression = ''

                # Build the expression - Intersperse the operators between gene names
                for i, (not_op, op, gene) in enumerate(zip(not_operators, operators, gene_connections)):
                    if i < len(gene_connections)-1: # If there are more genes being connected after this one
                        expression += ''.join([not_op + gene + op])
                    elif i == len(gene_connections): # If this is the last gene
                        expression += ''.join([not_op + gene])
                    else: # If there is only one connection
                        expression += ''.join([not_op + gene]) 
                
                # Write the formatted expression for the current node
                if expression != '':
                    line = f"{self.graph.nodes[node]['name']} = {expression}\n" 

                # If there are no nodes connected, create a self connection
                else:
                    line = f"{self.graph.nodes[node]['name']} = {self.graph.nodes[node]['name']}\n"

                # print(line)
                file.write(line)

    def export_network_graphml(self, filename):
        """Export the network to a GraphML file."""
        nx.write_graphml(self.graph, filename)
    
    def visualize_network(self):
        plt.figure(figsize=(12, 8))
        pos = nx.spring_layout(self.graph)  # positions for all nodes

        # Draw nodes
        nx.draw_networkx_nodes(self.graph, pos, node_size=700)

        # Draw edges
        nx.draw_networkx_edges(self.graph, pos, edgelist=self.graph.edges(), arrowstyle='->', arrowsize=20)

        # Draw labels
        nx.draw_networkx_labels(self.graph, pos, font_size=12, font_family='sans-serif')

        # Show edge attributes (optional)
        edge_labels = nx.get_edge_attributes(self.graph, 'edge_type')
        nx.draw_networkx_edge_labels(self.graph, pos, edge_labels=edge_labels)

        plt.axis('off')  # Turn off axis
        # plt.show()