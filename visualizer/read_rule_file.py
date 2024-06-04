import re
from node import Node
from gate import Gate
import random
import pygame
from alive_progress import alive_bar

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

def obj_from_line(target, expression, game, node_dict, gate_list):

    if target not in node_dict:
        node = Node(target, (game.WIDTH / 2 + random.randint(-200, 200), game.HEIGHT / 2 + random.randint(-200, 200)), "light blue")
        game.nodes.append(node)
        node_dict[target] = node
    
    rhs_terms = parse_rhs(expression)
    
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
        node_dict[target].incoming_connections.append(gate)
        gate.outgoing_connections.append(node_dict[target])
    else:
        for node in incoming_nodes:
            node.outgoing_connections.append(node_dict[target])
            node_dict[target].incoming_connections.append(node)
    
    return node_dict, gate_list

def create_nodes_and_gates(file_path, game):
    """Read the file, create nodes and gates, and set up connections."""
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    node_dict = {}
    gate_list = []

    with alive_bar(len(lines)) as bar:
        for line in lines:
            # Handles parentheses, basically builds the parentheses rules then adds on the non-parentheses rules
            if "(" in line:
                target = line.strip().split(' = ')[0].strip()
                rules = line.strip().split(' = ')[1].strip()
                parentheses = line.strip().split("(")[1].split(")")[0].strip().replace('(', '').replace(')', '')
                par_terms = re.split(r'\s*(and|or|not)\s*', parentheses)
                out_parentheses = rules.strip().replace('(', '').replace(')', '').replace(f'{parentheses}', '').replace(f'{target}', '').strip()
                not_par_terms = re.split(r'\s*(and|or|not)\s*', out_parentheses)

                if target not in node_dict:
                    node = Node(target, (game.WIDTH / 2 + random.randint(-200, 200), game.HEIGHT / 2 + random.randint(-200, 200)), "light blue")
                    game.nodes.append(node)
                    node_dict[target] = node

                incoming_nodes = []
                gate = None

                i = 0
                while i < len(par_terms):
                    term = par_terms[i]

                    # If the term is a Gate:
                    if term in {"and", "or", "not"}:
                        gate_type = term.upper()

                        # If no gate exists, create it
                        if gate is None:
                            gate = Gate(gate_type, (game.WIDTH / 2 + random.randint(-200, 200), game.HEIGHT / 2 + random.randint(-200, 200)))
                            game.gates.append(gate)
                            gate_list.append(gate)
                        
                        # If there is already a gate, make a new one and add its connections
                        else:
                            new_gate = Gate(gate_type, (game.WIDTH / 2 + random.randint(-200, 200), game.HEIGHT / 2 + random.randint(-200, 200)))
                            game.gates.append(new_gate)
                            gate_list.append(new_gate)
                            gate.outgoing_connections.append(new_gate)
                            new_gate.incoming_connections.append(gate)
                            gate = new_gate
                    
                    # Else if the term is a Node:
                    else:
                        # If the node does not exist, create it
                        if term not in node_dict:
                            node = Node(term, (game.WIDTH / 2 + random.randint(-200, 200), game.HEIGHT / 2 + random.randint(-200, 200)), "light blue")
                            game.nodes.append(node)
                            node_dict[term] = node
                        
                        # If a gate exists in this expression, add the gate to the nodes outgoing connections and the node to the gates incoming connections
                        if gate is not None:
                            node_dict[term].outgoing_connections.append(gate)
                            gate.incoming_connections.append(node_dict[term])
                        
                        # If the node does exist, add it to the list of incoming nodes for the gate 
                        else:
                            incoming_nodes.append(node_dict[term])
                    i += 1

                # If there is a gate:
                if gate is not None:
                    
                    # For each incoming node in the expression
                    for node in incoming_nodes:

                        # Add an edge from the nodes to the gates
                        node.outgoing_connections.append(gate)
                        gate.incoming_connections.append(node)
                    

                    # Now look at the rest of the expression outside of the parentheses

                    # Set the gate as the target and reset the gate to None
                    incoming_gate = gate
                    gate = None

                    i = 0
                    while i < len(not_par_terms):
                        term = not_par_terms[i]

                        # If the term is a Gate:
                        if term in {"and", "or", "not"}:
                            gate_type = term.upper()

                            # If no gate exists, create it
                            if gate is None:
                                gate = Gate(gate_type, (game.WIDTH / 2 + random.randint(-200, 200), game.HEIGHT / 2 + random.randint(-200, 200)))
                                game.gates.append(gate)
                                gate_list.append(gate)
                            
                            # If there is already a gate, make a new one and add its connections
                            else:
                                new_gate = Gate(gate_type, (game.WIDTH / 2 + random.randint(-200, 200), game.HEIGHT / 2 + random.randint(-200, 200)))
                                game.gates.append(new_gate)
                                gate_list.append(new_gate)
                                gate.outgoing_connections.append(new_gate)
                                new_gate.incoming_connections.append(gate)
                                gate = new_gate
                        
                        # Else if the term is a Node:
                        else:
                            # If the node does not exist, create it
                            if term not in node_dict:
                                node = Node(term, (game.WIDTH / 2 + random.randint(-200, 200), game.HEIGHT / 2 + random.randint(-200, 200)), "light blue")
                                game.nodes.append(node)
                                node_dict[term] = node

                            # If a gate exists in this expression, add the gate to the nodes outgoing connections and the node to the gates incoming connections
                            if gate is not None:
                                node_dict[term].outgoing_connections.append(gate)
                                gate.incoming_connections.append(node_dict[term])
                            
                            # If the node does exist, add it to the list of incoming nodes for the gate 
                            else:
                                incoming_nodes.append(node_dict[term])
                        i += 1

                    # If there is a gate outside of the parentheses:
                    if gate is not None:
                        
                        # Add an edge to the target gate
                        incoming_gate.outgoing_connections.append(gate)
                        gate.incoming_connections.append(incoming_gate)
                    
                    # Connect the gate to the target
                    gate.outgoing_connections.append(node_dict[target])
                    node_dict[target].incoming_connections.append(gate)


                # Else if there is no gate, connect the node straight to the target node
                else:
                    for node in incoming_nodes:
                        node.outgoing_connections.append(node_dict[target])
                        node_dict[target].incoming_connections.append(node)

                
            else:
                line = line.strip().replace('(', '').replace(')', '')
                if line == "":
                    continue
                
                try:
                    lhs, rhs = parse_equation(line)
                except ValueError as e:
                    print(e)
                    continue
                
                if lhs != rhs:
                
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
                                new_gate.incoming_connections.append(gate)
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
            bar()
    
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