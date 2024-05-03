import pygame, sys
from node import Node
from gate import Gate
from display_box import DisplayBox
import random
import pickle

class SaveState():
    def __init__(self):
        self.object_types = []
        self.object_names = []
        self.object_positions = []
        self.object_color = []
        self.object_ids = []
        self.outgoing_connections = []
        self.incoming_connections = []


    def save_objects(self, objects):
        # I want the name, position, size, outgoing_connections, and incoming_connections
        for object in objects:
            self.object_types.append(type(object).__name__)
            self.object_names.append(object.name)
            self.object_positions.append(object.position)
            self.object_color.append(object.color)
            self.object_ids.append(object.id)

            object_outgoing_connections = []
            object_incoming_connections = []

            for connected_object in object.outgoing_connections:
                object_outgoing_connections.append(connected_object.id)

            for connected_object in object.incoming_connections:
                object_incoming_connections.append(connected_object.id)

            self.outgoing_connections.append(object_outgoing_connections)
            self.incoming_connections.append(object_incoming_connections)
    
    def load_objects(self):
        nodes = []
        gates = []
        uuid_dict = {}  # Dictionary to map UUIDs to object instances

        # Instantiate all objects and map them
        for i, object_type in enumerate(self.object_types):
            if object_type == "Node":
                node = Node(self.object_names[i], self.object_positions[i], self.object_color[i])
                node.id = self.object_ids[i]
                nodes.append(node)
                uuid_dict[node.id] = node
            elif object_type == "Gate":
                gate = Gate(self.object_names[i], self.object_positions[i])
                gate.id = self.object_ids[i]
                gates.append(gate)
                uuid_dict[gate.id] = gate

        # Restore connections using the UUID dictionary
        for node in nodes:
            node.outgoing_connections = [uuid_dict[uid] for uid in self.outgoing_connections[i]]
            node.incoming_connections = [uuid_dict[uid] for uid in self.incoming_connections[i]]
        for gate in gates:
            gate.outgoing_connections = [uuid_dict[uid] for uid in self.outgoing_connections[i]]
            gate.incoming_connections = [uuid_dict[uid] for uid in self.incoming_connections[i]]


        # Group the objects into sprite groups
        nodes_group = pygame.sprite.Group(*nodes)
        gates_group = pygame.sprite.Group(*gates)
        objects_group = pygame.sprite.Group([gates_group, nodes_group])

        return nodes_group, gates_group, objects_group

def save_game_state(save_state, file_name):
    try:
        with open(file_name, 'wb') as file:
            pickle.dump(save_state, file)
            print("Game state saved successfully!")
    except IOError:
        print("Error: Unable to save game state.")

def load_game_state(file_name):
    try:
        with open(file_name, 'rb') as file:
            game_state = pickle.load(file)
            print("Game state loaded successfully!")
            return game_state
    except (IOError, pickle.UnpicklingError):
        print("Error: Unable to load game state.")
        

class Game:
    def __init__(self):

        pygame.init()
        self.WIDTH = 1920
        self.HEIGHT = 1080

        self.screen = pygame.display.set_mode((self.WIDTH,self.HEIGHT))
        pygame.display.set_caption('Graph visualizer')
        self.clock = pygame.time.Clock()

        self.uuids = {}

        num_nodes = 10
        num_and_gates = 10
        num_or_gates = 10
        num_not_gates = 4

        relative_abundance = False

        self.nodes = []
        self.gates = []

        # counter = 1
        # counter2 = 1
        # with open('visualizer/connections.txt', 'r') as connection_file:
        #     for line in connection_file:
        #         line = line.strip()
        #         target_name, connection_expr = line.split(' = ')
        #         elements = connection_expr.split()

        #         incoming_nodes = [ele for ele in elements if ele not in ("OR", "AND", "NOT")]
        #         gates = [ele for ele in elements if ele in ("OR", "AND", "NOT")]

        #         #print(f'Node: {target_name}')
        #         #print(f'\tincoming gates: {gates}')
        #         #print(f'\tincoming nodes: {incoming_nodes}')

        #         target_node = Node(target_name, (self.WIDTH / 2 + 350 + (50*counter2), self.HEIGHT / 2 + (15 * counter)), "light blue")
        #         self.nodes.append(target_node)

        #         existing_nodes = {node.name: node for node in self.nodes}
        #         gate_objects = []
        #         last_output = None

        #         def parse_not_gate(node_position):
        #             # NOT gates should only have one input
        #             node_name = incoming_nodes[node_position]
        #             node = existing_nodes.get(node_name, Node(node_name, gate_position, "light blue"))
        #             if node not in self.nodes:
        #                 self.nodes.append(node)

        #             not_gate = Gate('NOT', gate_position)
        #             self.gates.append(not_gate)
        #             not_gate.incoming_connections.append(node)
        #             node.outgoing_connections.append(not_gate)

        #             gate_objects.append(not_gate)
                
        #         count = 0
        #         gates_no_not = [ele for ele in elements if ele not in ("OR", "AND")]
        #         skip = False
        #         for i, gate_type in enumerate(gates):
                    
        #             if skip == False:
        #                 # Determine the correct position for visualization
        #                 gate_position = (self.WIDTH / 2 + 200 + (50*counter2), self.HEIGHT / 2 - 150 + (15 * counter))
                        
        #                 if gate_type == 'NOT':
        #                     parse_not_gate(node_position=0)

        #                 else:
        #                     # AND/OR gates
        #                     gate = Gate(gate_type, gate_position)
        #                     self.gates.append(gate)
        #                     #print(f'\n\tGate type: {gate_type}')

        #                     # Get the incoming node names
        #                     node1_name = incoming_nodes[0+count]
        #                     #print(f'\t\tnode1_name = {node1_name}')
        #                     node1 = existing_nodes.get(node1_name, Node(node1_name, gate_position, "light blue"))
        #                     # Append the node and it's connections
        #                     self.nodes.append(node1)
        #                     gate.incoming_connections.append(node1)
        #                     node1.outgoing_connections.append(gate)
        #                     count += 1
                            
                        
        #                     # Check if the next node has a NOT gate
        #                     if count < len(gates)-1:
        #                         if gates[count] == 'NOT':
        #                             parse_not_gate(count)
        #                             skip = True

        #                     if i == 0:
        #                         # Get the second node for the gate
        #                         # Get the incoming node names
        #                         node2_name = incoming_nodes[count]
        #                         #print(f'\t\tnode2_name = {node2_name}')
        #                         node2 = existing_nodes.get(node2_name, Node(node2_name, gate_position, "light blue"))
        #                         self.nodes.append(node2)
        #                         # Append the node and it's connections
        #                         gate.incoming_connections.append(node2)
        #                         node2.outgoing_connections.append(gate)
        #                         count += 1
                                
                            
        #                     else:
        #                         #print(f'\t\tincoming gate name = {gate_objects[-1].name}')
        #                         gate.incoming_connections.append(gate_objects[-1])
        #                         gate_objects[-1].outgoing_connections.append(gate)
                                
        #                     # #print(f'\tGate incoming connections:')
        #                     # for node in gate.incoming_connections:
        #                     #     #print(f'\t\t{node.name}')

        #                     gate_objects.append(gate)

        #             else:
        #                 skip = True
                
        #         if len(gate_objects) > 0:
        #             target_node.incoming_connections.append(gate_objects[-1])
        #             gate_objects[-1].outgoing_connections.append(target_node)
                    
        #         if len(incoming_nodes) == 1:
        #             node_name = incoming_nodes[0]
        #             node = existing_nodes.get(node_name, Node(node_name, gate_position, "light blue"))
        #             self.nodes.append(node)
        #             target_node.incoming_connections.append(node)
        #             node.outgoing_connections.append(target_node)

                
        #         if counter == 25:
        #             counter = 0
        #             counter2 += 1
        #         counter += 1

        # Boolean Nodes
        # for i in range(1, num_nodes+1):
        #     self.nodes.append(Node('', (self.WIDTH/2+350,self.HEIGHT/2+(50*i)), "light blue"))

        for i in range(1, num_nodes+1):
            self.nodes.append(Node(f'Node {i}', (self.WIDTH/2+350,self.HEIGHT/2+(50*i)), "light blue"))

        # Change the size and color of the nodes for relative abundance. Manual entry required
        if relative_abundance:
            colors = [(239, 224, 215), (255, 208, 198), (235, 103, 75), (206, 219, 255), (255, 98, 98), "blue", "red"]
            sizes = [60, 35, 40, 35, 55, 65, 70]
            for node_num, node in enumerate(self.nodes):
                node.color = colors[node_num]
                node.set_size(sizes[node_num])
                node.font = pygame.font.Font('arial.ttf', int(round(node.size / 2, 1)))

        # AND Gates
        for i in range(num_and_gates):
            self.gates.append(Gate('AND', (self.WIDTH/2+200, self.HEIGHT/2-150 - (15*i))))

        # OR Gates
        for i in range(num_or_gates):
            self.gates.append(Gate('OR', (self.WIDTH/2+300, self.HEIGHT/2-150 - (15*i))))

        # NOT Gates
        for i in range(num_not_gates):
            self.gates.append(Gate('NOT', (self.WIDTH/2+400, self.HEIGHT/2-150 - (15*i))))

        # Create a list of the unique IDs for the objects
        self.node_ids = [node.id for node in self.nodes]
        self.gate_ids = [gate.id for gate in self.gates]

        self.nodes_group = pygame.sprite.Group(*self.nodes)
        self.gates_group = pygame.sprite.Group(*self.gates)

        self.objects_group = pygame.sprite.Group(*self.nodes, *self.gates)


        for object in self.objects_group:
            self.uuids[object.id] = object
            object.node_ids = self.node_ids
            object.gate_ids = self.gate_ids

        self.connections = 0        

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

            keys = pygame.key.get_pressed()
            if keys[pygame.K_s]:
                game_state = SaveState()
                game_state.save_objects(self.objects_group)
                save_game_state(game_state, 'save_game.pickle')

            if keys[pygame.K_l]:
                game_state = load_game_state('save_game.pickle')
                self.nodes_group, self.gates_group, self.objects_group = game_state.load_objects()

            self.screen.fill("white")      
            for node in self.nodes_group:
                # node.draw_connections()
                node.update_object(events, self.connections, self.gate_ids, self.node_ids, self.uuids, self)
                if node.is_drawing_line == True:
                    self.connections += 1  
            
            for gate in self.gates_group:
                # gate.draw_connections()
                gate.update_object(events, self.connections, self.gate_ids, self.node_ids, self.uuids, self)
                if gate.is_drawing_line == True:
                    self.connections += 1

            # Reset the number of connections if the 1 key is not pressed
            
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
