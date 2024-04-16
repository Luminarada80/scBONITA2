import pygame, sys
from node import Node
from gate import Gate
from display_box import DisplayBox
import random
import pickle

class Game:
    def __init__(self):

        pygame.init()
        self.WIDTH = 1920
        self.HEIGHT = 1080

        self.screen = pygame.display.set_mode((self.WIDTH,self.HEIGHT))
        pygame.display.set_caption('Graph visualizer')
        self.clock = pygame.time.Clock()

        num_nodes = 7
        num_and_gates = 25
        num_or_gates = 25
        num_not_gates = 10

        relative_abundance = False

        self.nodes = []
        self.gates = []

        # node_list = {
        #     "CDH5",
        #     "PECAM1",
        #     "ITGB2",
        #     "CTNND1",
        #     "BCAR1",
        #     "SIPA1",
        #     "F11R",
        #     "ITGAL",
        #     "ITGB1",
        #     "ITGA4",
        #     "JAM2",
        #     "ITGAM",
        #     "CTNNB1",
        #     "PTK2",
        #     "RHOA",
        #     "JAM3",
        #     "CYBA",
        #     "NCF2",
        #     "PXN",
        #     "NCF1",
        #     "RAC1",
        #     "NCF4",
        #     "RASSF5",
        #     "PTK2B",
        #     "CXCL12",
        #     "CXCR4",
        #     "THY1",
        #     "VCAM1",
        #     "ICAM1",
        #     "OCLN",
        #     "ESAM",
        #     "MMP2",
        #     "MMP9",
        #     "MYL2",
        #     "MYL5",
        #     "MYL7",
        #     "MYL9",
        #     "MYL10",
        #     "MYLPF",
        #     "ROCK1",
        #     "ROCK2",
        #     "PRKCA",
        #     "PRKCB",
        #     "PRKCG",
        #     "ARHGAP35",
        #     "ARHGAP5",
        #     "RAP1A",
        #     "RAP1B",
        #     "MSN",
        #     "EZR",
        #     "RAPGEF4",
        #     "RAPGEF3",
        #     "PLCG1",
        #     "PLCG2",
        #     }

        # # Nodes
        # for i, name in enumerate(node_list):
        #     if i < 10:
        #         self.nodes.append(Node(f'{name}', (self.WIDTH/2+350,self.HEIGHT/2+(25*i)), "light blue"))
        #     if 11 <= i <= 20:
        #         self.nodes.append(Node(f'{name}', (self.WIDTH/2+400,self.HEIGHT/2+(25*(i-10))), "light blue"))
        #     if 21 <= i <= 30:
        #         self.nodes.append(Node(f'{name}', (self.WIDTH/2+450,self.HEIGHT/2+(25*(i-20))), "light blue"))
        #     if 31 <= i <= 40:
        #         self.nodes.append(Node(f'{name}', (self.WIDTH/2+500,self.HEIGHT/2+(25*(i-30))), "light blue"))
        #     if 41 <= i <= 50:
        #         self.nodes.append(Node(f'{name}', (self.WIDTH/2+550,self.HEIGHT/2+(25*(i-40))), "light blue"))
        #     if 51 <= i <= 60:
        #         self.nodes.append(Node(f'{name}', (self.WIDTH/2+600,self.HEIGHT/2+(25*(i-50))), "light blue"))
        
        # counter = 1
        # with open('visualizer/connections.txt', 'r') as connection_file:
        #     for line in connection_file:
        #         line = line.strip()
        #         target_name = line.split(' = ')[0]
        #         connections = line.split(' = ')[1].split(' ')
        #         incoming_gates = [gate for gate in connections if gate in ("OR", "AND", "NOT")]
        #         incoming_nodes = [gate for gate in connections if gate not in ("OR", "AND", "NOT")]

        #         #print(f'Node: {target_name}') 
        #         #print(f'\tincoming gates: {incoming_gates}')
        #         #print(f'\tincoming nodes: {incoming_nodes}')

        #         # Make the target node
        #         target_node = Node(f'{target_name}', (self.WIDTH/2+350,self.HEIGHT/2+(25*counter)), "light blue")

        #         node_objects = []

        #         # Get a list of the nodes that are already in the nodes list
        #         existing_nodes = [node.name for node in self.nodes]

        #         # For each incoming node to the target node:
        #         for inc_node in incoming_nodes:             
        #             # If the incoming node is not in the list of existing node names, create it           
        #             if inc_node not in existing_nodes:
        #                 node = Node(f'{inc_node}', (self.WIDTH/2+350,self.HEIGHT/2+(25*counter)), "light blue")
        #                 self.nodes.append(node)
        #                 node_objects.append(node)

        #             # If the incoming node is in the list, reference it
        #             else:
        #                 node = self.nodes[existing_nodes.index(inc_node)]
        #                 node_objects.append(node)

        #         gates = []
        #         # If there is more than one incoming node
        #         if len(node_objects) > 0:
        #             if len(incoming_gates) > 0:
        #                 if incoming_gates[0] == 'AND':
        #                     and_gate = Gate('AND', (self.WIDTH/2+200, self.HEIGHT/2-150 - (15*counter)))
        #                     self.gates.append(and_gate)
        #                     gates.append(and_gate)

        #                     node_objects[0].outgoing_connections.append(and_gate)
        #                     node_objects[1].outgoing_connections.append(and_gate)

        #                     and_gate.incoming_connections.append(node_objects[0])
        #                     and_gate.incoming_connections.append(node_objects[1])
                        
        #                 elif incoming_gates[0] == 'OR':
        #                     or_gate = Gate('OR', (self.WIDTH/2+300, self.HEIGHT/2-150 - (15*counter)))
        #                     self.gates.append(or_gate)
        #                     gates.append(or_gate)

        #                     node_objects[0].outgoing_connections.append(or_gate)
        #                     node_objects[1].outgoing_connections.append(or_gate)

        #                     or_gate.incoming_connections.append(node_objects[0])
        #                     or_gate.incoming_connections.append(node_objects[1])
                        
        #                 if incoming_gates[0] == 'NOT':
        #                     not_gate = Gate('AND', (self.WIDTH/2+200, self.HEIGHT/2-150 - (15*counter)))
        #                     self.gates.append(not_gate)
        #                     gates.append(not_gate)

        #                     node_objects[0].outgoing_connections.append(not_gate)
        #                     node_objects[1].outgoing_connections.append(not_gate)

        #                     not_gate.incoming_connections.append(node_objects[0])
        #                     not_gate.incoming_connections.append(node_objects[1])
                        
        #                 if len(incoming_gates) == 1:
        #                     target_node.incoming_connections.append(gates[0])
        #                     gates[0].outgoing_connections.append(target_node)

        #             #print(f'Num incoming gates: {len(incoming_gates)}')
        #             if len(incoming_gates) == 2:
        #                 if incoming_gates[1] == 'AND':
        #                     and_gate = Gate('AND', (self.WIDTH/2+200, self.HEIGHT/2-150 - (15*counter)))
        #                     self.gates.append(and_gate)
        #                     gates.append(and_gate)

        #                     gates[1].outgoing_connections.append(and_gate)

        #                     node_objects[2].outgoing_connections.append(and_gate)

        #                     and_gate.incoming_connections.append(gates[1])
        #                     and_gate.incoming_connections.append(node_objects[1])

        #                 if incoming_gates[1] == 'OR':
        #                     or_gate = Gate('OR', (self.WIDTH/2+300, self.HEIGHT/2-150 - (15*counter)))
        #                     self.gates.append(or_gate)
        #                     gates.append(or_gate)

        #                     gates[1].outgoing_connections.append(or_gate)

        #                     node_objects[2].outgoing_connections.append(or_gate)

        #                     or_gate.incoming_connections.append(gates[1])
        #                     or_gate.incoming_connections.append(node_objects[1])
                        
        #                 if incoming_gates[1] == 'NOT':
        #                     not_gate = Gate('AND', (self.WIDTH/2+200, self.HEIGHT/2-150 - (15*counter)))
        #                     self.gates.append(not_gate)
        #                     gates.append(not_gate)

        #                     gates[1].outgoing_connections.append(not_gate)

        #                     node_objects[2].outgoing_connections.append(not_gate)

        #                     not_gate.incoming_connections.append(gates[1])
        #                     not_gate.incoming_connections.append(node_objects[1])
                        
        #                 # If there is one rule
        #                 if len(incoming_gates) == 1:
        #                     # Set the incoming nodes 
        #                     target_node.incoming_connections.append(gates)
        #                     gates[0].outgoing_connections.append(target_node)

        #                 # If there are two rules
        #                 elif len(incoming_gates) == 2:
        #                     # Set the outgoing connection of the first rule to the second rule
        #                     gates[0].outgoing_connections.append(gates[1])

        #                     # Set the outgoing connection of the second rule to the target
        #                     gates[1].outgoing_connections.append(target_node)

        #                     # Set the target nodes incoming connection to the second gate
        #                     target_node.incoming_connections.append(gates)
                            
                    
        #             #print(f'\tIncoming connections: {target_node.incoming_connections}')
        #             counter += 1

        counter = 1
        counter2 = 1
        with open('visualizer/connections.txt', 'r') as connection_file:
            for line in connection_file:
                line = line.strip()
                target_name, connection_expr = line.split(' = ')
                elements = connection_expr.split()

                incoming_nodes = [ele for ele in elements if ele not in ("OR", "AND", "NOT")]
                gates = [ele for ele in elements if ele in ("OR", "AND", "NOT")]

                #print(f'Node: {target_name}')
                #print(f'\tincoming gates: {gates}')
                #print(f'\tincoming nodes: {incoming_nodes}')

                target_node = Node(target_name, (self.WIDTH / 2 + 350 + (50*counter2), self.HEIGHT / 2 + (15 * counter)), "light blue")
                self.nodes.append(target_node)

                existing_nodes = {node.name: node for node in self.nodes}
                gate_objects = []
                last_output = None

                def parse_not_gate(node_position):
                    # NOT gates should only have one input
                    node_name = incoming_nodes[node_position]
                    node = existing_nodes.get(node_name, Node(node_name, gate_position, "light blue"))
                    if node not in self.nodes:
                        self.nodes.append(node)

                    not_gate = Gate('NOT', gate_position)
                    self.gates.append(not_gate)
                    not_gate.incoming_connections.append(node)
                    node.outgoing_connections.append(not_gate)

                    gate_objects.append(not_gate)
                
                count = 0
                gates_no_not = [ele for ele in elements if ele not in ("OR", "AND")]
                skip = False
                for i, gate_type in enumerate(gates):
                    
                    if skip == False:
                        # Determine the correct position for visualization
                        gate_position = (self.WIDTH / 2 + 200 + (50*counter2), self.HEIGHT / 2 - 150 + (15 * counter))
                        
                        if gate_type == 'NOT':
                            parse_not_gate(node_position=0)

                        else:
                            # AND/OR gates
                            gate = Gate(gate_type, gate_position)
                            self.gates.append(gate)
                            #print(f'\n\tGate type: {gate_type}')

                            # Get the incoming node names
                            node1_name = incoming_nodes[0+count]
                            #print(f'\t\tnode1_name = {node1_name}')
                            node1 = existing_nodes.get(node1_name, Node(node1_name, gate_position, "light blue"))
                            # Append the node and it's connections
                            self.nodes.append(node1)
                            gate.incoming_connections.append(node1)
                            node1.outgoing_connections.append(gate)
                            count += 1
                            
                        
                            # Check if the next node has a NOT gate
                            if count < len(gates)-1:
                                if gates[count] == 'NOT':
                                    parse_not_gate(count)
                                    skip = True

                            if i == 0:
                                # Get the second node for the gate
                                # Get the incoming node names
                                node2_name = incoming_nodes[count]
                                #print(f'\t\tnode2_name = {node2_name}')
                                node2 = existing_nodes.get(node2_name, Node(node2_name, gate_position, "light blue"))
                                self.nodes.append(node2)
                                # Append the node and it's connections
                                gate.incoming_connections.append(node2)
                                node2.outgoing_connections.append(gate)
                                count += 1
                                
                            
                            else:
                                #print(f'\t\tincoming gate name = {gate_objects[-1].name}')
                                gate.incoming_connections.append(gate_objects[-1])
                                gate_objects[-1].outgoing_connections.append(gate)
                                

                            # #print(f'\tGate incoming connections:')
                            # for node in gate.incoming_connections:
                            #     #print(f'\t\t{node.name}')
                            


                            gate_objects.append(gate)

                    else:
                        skip = True
                
                if len(gate_objects) > 0:
                    target_node.incoming_connections.append(gate_objects[-1])
                    gate_objects[-1].outgoing_connections.append(target_node)
                    
                if len(incoming_nodes) == 1:
                    node_name = incoming_nodes[0]
                    node = existing_nodes.get(node_name, Node(node_name, gate_position, "light blue"))
                    self.nodes.append(node)
                    target_node.incoming_connections.append(node)
                    node.outgoing_connections.append(target_node)

                # for i, gate in enumerate(gate_objects):
                    #print(f'\n\tGate name: {gate.name}')
                    #print(f'\t\tIncoming connections:')
                    # for node in gate.incoming_connections:
                        #print(f'\t\t\t{node.name}')

                    #print(f'\t\tOutgoing connections:')
                    # for node in gate.outgoing_connections:
                        #print(f'\t\t\t{node.name}')
                
                if counter == 25:
                    counter = 0
                    counter2 += 1
                counter += 1
                
                



            
        
        # Boolean Nodes
        # for i in range(1, num_nodes+1):
        #     self.nodes.append(Node('', (self.WIDTH/2+350,self.HEIGHT/2+(50*i)), "light blue"))

        # Change the size and color of the nodes for relative abundance. Manual entry required
        if relative_abundance:
            colors = [(239, 224, 215), (255, 208, 198), (235, 103, 75), (206, 219, 255), (255, 98, 98), "blue", "red"]
            sizes = [60, 35, 40, 35, 55, 65, 70]
            for node_num, node in enumerate(self.nodes):
                node.color = colors[node_num]
                node.set_size(sizes[node_num])
                node.font = pygame.font.Font('arial.ttf', int(round(node.size / 2, 1)))

        
        # # AND Gates
        # for i in range(num_and_gates):
        #     self.gates.append(Gate('AND', (self.WIDTH/2+200, self.HEIGHT/2-150 - (15*i))))

        # # OR Gates
        # for i in range(num_or_gates):
        #     self.gates.append(Gate('OR', (self.WIDTH/2+300, self.HEIGHT/2-150 - (15*i))))

        # # NOT Gates
        # for i in range(num_not_gates):
        #     self.gates.append(Gate('NOT', (self.WIDTH/2+400, self.HEIGHT/2-150 - (15*i))))

        self.nodes_group = pygame.sprite.Group(*self.nodes)
        self.gates_group = pygame.sprite.Group(*self.gates)

        self.connections = 0

        self.objects_group = pygame.sprite.Group([self.gates_group, self.nodes_group])

        self.display_box = DisplayBox()

        self.update_num = 0

        self.states = {}

        self.node_being_moved = None  # No node is being moved initially

    def run(self):
        while True:
            events = pygame.event.get()
            for event in events:
                if event.type == pygame.QUIT:
                    pygame.quit()
                    sys.exit()

            self.screen.fill("white")      
            for node in self.nodes_group:
                node.draw_connections()
                node.update_object(events, self.connections, self.objects_group, self)
                if node.is_drawing_line == True:
                    self.connections += 1  
            
            for gate in self.gates_group:
                gate.draw_connections()
                gate.update_object(events, self.connections, self.objects_group, self)
                if gate.is_drawing_line == True:
                    self.connections += 1

            # Reset the number of connections if the 1 key is not pressed
            keys = pygame.key.get_pressed()
            if not keys[pygame.K_1]:
                self.connections = 0

            self.nodes_group.draw(self.screen)

            for node in self.nodes_group:
                node.draw_lock((node.position[0] - 37, node.position[1]))

            keys = pygame.key.get_pressed()
            if keys[pygame.K_TAB]:
                self.display_box.display_text('Node State', (self.WIDTH / 4 + 25, self.HEIGHT - 400))
                self.display_box.display_text(f'Update {len(self.states)}', (self.WIDTH / 4 + 25, self.HEIGHT - 375))

                if self.nodes[0].update_num <= 35:
                    self.states[self.nodes[0].update_num] = []
                    for node_num, node in enumerate(self.nodes_group):
                        position_adjustment = 15 * node_num
                        update_adjustment = 10 * self.nodes[0].update_num - 150
                        self.states[self.nodes[0].update_num].append((node.state, (self.WIDTH / 4 + update_adjustment, self.HEIGHT - 350 + position_adjustment)))


                # Displaying states
                for update in self.states:
                    for node_num, node in enumerate(self.nodes_group):
                        position_adjustment = 25 * node_num
                        update_adjustment = 15 * len(self.states) - 100
                        self.display_box.display_text(f'{self.states[update][node_num][0]}', self.states[update][node_num][1])

            else:
                self.states = {}

            pygame.display.update()
            self.clock.tick(60)

if __name__ == "__main__":
    game = Game()
    game.run()
