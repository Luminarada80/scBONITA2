import networkx as nx
import random
import itertools
import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations, product
import re

class CreateTestNetwork():
    def __init__(self, num_genes):
        self.num_genes = num_genes  # number of nodes in the network
        self.graph = None

    def create_network(self, num_genes):
        network = nx.DiGraph()

        # Add gene nodes using the gene name as the node identifier
        for i in range(num_genes):
            gene_name = f'Gene{i}'
            network.add_node(gene_name, name=gene_name)

        # Add random edges between genes
        for node_num in range(num_genes):
            gene_name = f'Gene{node_num}'
            num_connections = random.randint(1, 3)  # Random number of connections
            possible_targets = [f'Gene{j}' for j in range(num_genes) if j != node_num and network.in_degree(f'Gene{j}') < 3]
            if possible_targets:
                connections = random.sample(possible_targets, min(num_connections, len(possible_targets)))
                for target in connections:
                    interaction = 'a' if random.random() < 0.8 else 'i'
                    network.add_edge(gene_name, target, interaction=interaction)

        # Add self-loops to genes without predecessors
        for node_num in range(num_genes):
            gene_name = f'Gene{node_num}'
            if len(list(network.predecessors(gene_name))) == 0:
                network.add_edge(gene_name, gene_name, interaction='a')
        
        self.num_genes = num_genes
        self.graph = network

        return network

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
    
    def enumerate_possibilities(self, num_incoming_nodes):
        cache = dict()
        or0, or1, and0, and1 = "  or  ", " or ", "  and  ", " and "  

        def cacheResult(keys, result=None):
            if not result:
                return [r.format(*keys) for r in cache.get(len(keys), [])]   
            cache[len(keys)] = resFormat = []
            result = sorted(result, key=lambda x: x.replace("  ", " "))
            for r in result:
                r = r.replace("and", "&").replace("or", "|")
                for i, k in enumerate(keys):
                    r = r.replace(k, "{" + str(i) + "}")
                r = r.replace("&", "and").replace("|", "or")
                resFormat.append(r)
            return result

        def boolCombo(keys):
            if len(keys) == 1: 
                return list(keys)
            
            result = cacheResult(keys) or set()
            if result: 
                return result
            
            def addResult(left, right):
                OR = or0.join(sorted(left.split(or0) + right.split(or0)))
                result.add(OR.replace(and0, and1))
                if or0 in left:  
                    left = f"({left})"
                if or0 in right: 
                    right = f"({right})"
                AND = and0.join(sorted(left.split(and0) + right.split(and0)))
                result.add(AND.replace(or0, or1))
                    
            seenLeft = set()
            for leftSize in range(1, len(keys) // 2 + 1):
                for leftKeys in combinations(keys, leftSize):
                    rightKeys = [k for k in keys if k not in leftKeys]
                    if len(leftKeys) == len(rightKeys):
                        if tuple(rightKeys) in seenLeft: continue
                        seenLeft.add(tuple(leftKeys))
                    for left, right in product(*map(boolCombo, (leftKeys, rightKeys))):
                        addResult(left, right)
            return cacheResult(keys, result)
                
        if num_incoming_nodes == 4:
            possibilities = np.array(boolCombo("ABCD") + boolCombo("ABC") + boolCombo("AB") + boolCombo("AC") + boolCombo("BC") + ["A"] + ["B"] + ["C"])
        if num_incoming_nodes == 3:
            possibilities = np.array(boolCombo("ABC") + boolCombo("AB") + boolCombo("AC") + boolCombo("BC") + ["A"] + ["B"] + ["C"])
        elif num_incoming_nodes == 2:
            possibilities = np.array(boolCombo("AB") + ["A"] + ["B"])
        elif num_incoming_nodes == 1 or num_incoming_nodes == 0:
            possibilities = np.array(["A"])
        else:
            assert IndexError('Num incoming nodes out of range')

        random_rule = random.choice(possibilities).replace("  ", " ")

        return random_rule

    def generate_network_rules(self):
        rule_dict = {}
        for node in self.graph.nodes():
            gene_connections = []
            operators = []
            not_operators = []
            for source, target, data in self.graph.in_edges(node, data=True):
                gene_connections.append(self.graph.nodes[source]['name'])
                not_operator = True if data['interaction'] == 'i' else False
                operator = str(random.choice([' and ', ' or ']))
                not_operators.append(not_operator)
                operators.append(operator)
            
            random_rule = self.enumerate_possibilities(len(gene_connections))
            
            A = ''
            B = ''
            C = ''
            if len(gene_connections) > 0:
                A = f'not {gene_connections[0]}' if not_operators[0] else gene_connections[0]
            if len(gene_connections) > 1:
                B = f'not {gene_connections[1]}' if not_operators[1] else gene_connections[1]
            if len(gene_connections) > 2:
                C = f'not {gene_connections[2]}' if not_operators[2] else gene_connections[2]
            
            formatted_rule = random_rule.replace('A', A).replace('B', B).replace('C', C)
            
            expression = formatted_rule
            line = f"{self.graph.nodes[node]['name']} = {expression}\n" if expression else f"{self.graph.nodes[node]['name']} = {self.graph.nodes[node]['name']}\n"

            if node not in rule_dict:
                rule_dict[node] = None
            rule_dict[node] = expression
        
        return rule_dict
    
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
        plt.show()

