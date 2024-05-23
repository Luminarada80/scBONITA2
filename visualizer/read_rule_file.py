import re
from node import Node
from gate import Gate
import random
import pygame

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
            if term in {"and", "or", "not"}:
                gate_type = term.upper()
                if gate is None:
                    gate = Gate(gate_type, (game.WIDTH / 2 + random.randint(-200, 200), game.HEIGHT / 2 + random.randint(-200, 200)))
                    game.gates.append(gate)
                    gate_list.append(gate)
                else:
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
    
    # Update game attributes
    game.node_ids = [node.id for node in game.nodes]
    game.gate_ids = [gate.id for gate in game.gates]
    game.nodes_group = pygame.sprite.Group(*game.nodes)
    game.gates_group = pygame.sprite.Group(*game.gates)
    game.objects_group = pygame.sprite.Group(*game.nodes, *game.gates)
    
    for obj in game.objects_group:
        game.uuids[obj.id] = obj
        obj.node_ids = game.node_ids
        obj.gate_ids = game.gate_ids

    return node_dict, gate_list