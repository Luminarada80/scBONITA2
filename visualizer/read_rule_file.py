import re
from node import Node
from gate import Gate
import random

def parse_equation(equation):
    """Parse the equation into LHS and RHS components."""
    match = re.match(r"^\s*([^\s=]+)\s*=\s*(.+)\s*$", equation)
    if not match:
        raise ValueError(f"Invalid equation format: {equation}")
    lhs, rhs = match.groups()
    return lhs.strip(), rhs.strip()

def parse_rhs(rhs):
    """Parse the RHS into a list of terms and logical operators."""
    terms = re.split(r'\s*(and|or|not)\s*', rhs)
    return terms

def create_nodes_and_gates(file_path, game):
    """Read the file, create nodes and gates, and set up connections."""
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    node_dict = {}
    gate_list = []

    for line in lines:
        line = line.strip().replace('(', '').replace(')', '')
        if line == "":
            continue
        
        try:
            lhs, rhs = parse_equation(line)
        except ValueError as e:
            print(e)
            continue
        
        if lhs not in node_dict:
            node = Node(lhs, (game.WIDTH / 2 + random.randint(-200, 200), game.HEIGHT / 2 + random.randint(-200, 200)), "light blue")
            game.nodes.append(node)
            node_dict[lhs] = node
        
        rhs_terms = parse_rhs(rhs)
        
        incoming_nodes = []
        gate = None
        
        i = 0
        while i < len(rhs_terms):
            term = rhs_terms[i]
            if term == "and" or term == "or" or term == "not":
                if gate is None:
                    gate_type = term.upper()
                    gate = Gate(gate_type, (game.WIDTH / 2 + random.randint(-200, 200), game.HEIGHT / 2 + random.randint(-200, 200)))
                    game.gates.append(gate)
                    gate_list.append(gate)
                else:
                    gate_type = term.upper()
                    new_gate = Gate(gate_type, (game.WIDTH / 2 + random.randint(-200, 200), game.HEIGHT / 2 + random.randint(-200, 200)))
                    game.gates.append(new_gate)
                    gate_list.append(new_gate)
                    gate.outgoing_connections.append(new_gate)
                    gate = new_gate
            else:
                if term not in node_dict:
                    node = Node(term, (game.WIDTH / 2 + random.randint(-200, 200), game.HEIGHT / 2 + random.randint(-200, 200)), "light blue")
                    game.nodes.append(node)
                    node_dict[term] = node
                
                if gate is not None:
                    node_dict[term].outgoing_connections.append(gate)
                    gate.incoming_connections.append(node_dict[term])
                else:
                    incoming_nodes.append(node_dict[term])
            i += 1

        if gate is not None:
            for node in incoming_nodes:
                node.outgoing_connections.append(gate)
                gate.incoming_connections.append(node)
            node_dict[lhs].incoming_connections.append(gate)
            gate.outgoing_connections.append(node_dict[lhs])
        else:
            for node in incoming_nodes:
                node.outgoing_connections.append(node_dict[lhs])
                node_dict[lhs].incoming_connections.append(node)
    
    return node_dict, gate_list

# Example usage:
# game = Game()  # Ensure you have an instance of your Game class
# create_nodes_and_gates('path_to_your_file.txt', game)
# Now game.nodes and game.gates should be populated based on the file content
