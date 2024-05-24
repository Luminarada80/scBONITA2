import pygame, sys
from node import Node
from gate import Gate
from display_box import DisplayBox
import random
import pickle
import read_rule_file
from datetime import datetime

class SaveState():
    def __init__(self):
        self.object_types = []
        self.object_names = []
        self.object_positions = []
        self.object_color = []
        self.object_ids = []
        self.object_image = []
        self.outgoing_connections = []
        self.incoming_connections = []
        self.nodes = []
        self.gates = []
        self.uuid_dict = {}  # Dictionary to map UUIDs to object instances
        
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
        # Instantiate all objects and map them
        for i, object_type in enumerate(self.object_types):
            if object_type == "Node":
                node = Node(self.object_names[i], self.object_positions[i], self.object_color[i])
                node.id = self.object_ids[i]
                node.rect = node.image.get_rect(center=node.position)
                node.outgoing_connections = self.outgoing_connections[i]
                node.incoming_connections = self.incoming_connections[i]

                if node not in self.nodes:
                    self.nodes.append(node)
                
                if node.id not in self.uuid_dict:
                    self.uuid_dict[node.id] = node

            elif object_type == "Gate":
                gate = Gate(self.object_names[i], self.object_positions[i])
                gate.id = self.object_ids[i]
                gate.inactive_image, gate.active_image = gate.choose_gate_image()
                gate.image = gate.inactive_image
                gate.rect = gate.image.get_rect(center=gate.position)
                gate.outgoing_connections = self.outgoing_connections[i]
                gate.incoming_connections = self.incoming_connections[i]

                if gate not in self.gates:
                    self.gates.append(gate)
                
                if gate.id not in self.uuid_dict:
                    self.uuid_dict[gate.id] = gate

        # Restore connections using the UUID dictionary
        for i, object_type in enumerate(self.object_types):
            object = self.uuid_dict[self.object_ids[i]]
            object.outgoing_connections = [self.uuid_dict[uid] for uid in self.outgoing_connections[i]]
            object.incoming_connections = [self.uuid_dict[uid] for uid in self.incoming_connections[i]]


        # Group the objects into sprite groups
        nodes_group = pygame.sprite.Group(*self.nodes)
        gates_group = pygame.sprite.Group(*self.gates)
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

        # num_nodes = 10
        # num_and_gates = 10
        # num_or_gates = 10
        # num_not_gates = 4

        # relative_abundance = False

        self.nodes = []
        self.gates = []

        # read_rule_file.create_nodes_and_gates('visualizer/04670.txt', self)
        read_rule_file.create_nodes_and_gates('visualizer/05171.txt', self)

        # Boolean Nodes
        # for i in range(1, num_nodes+1):
        #     self.nodes.append(Node('', (self.WIDTH/2+350,self.HEIGHT/2+(50*i)), "light blue"))

        # for i in range(1, num_nodes+1):
        #     self.nodes.append(Node(f'Node {i}', (self.WIDTH/2+350,self.HEIGHT/2+(50*i)), "light blue"))

        # Change the size and color of the nodes for relative abundance. Manual entry required
        # if relative_abundance:
        #     colors = [(239, 224, 215), (255, 208, 198), (235, 103, 75), (206, 219, 255), (255, 98, 98), "blue", "red"]
        #     sizes = [60, 35, 40, 35, 55, 65, 70]
        #     for node_num, node in enumerate(self.nodes):
        #         node.color = colors[node_num]
        #         node.set_size(sizes[node_num])
        #         node.font = pygame.font.Font('arial.ttf', int(round(node.size / 2, 1)))

        # # AND Gates
        # for i in range(num_and_gates):
        #     self.gates.append(Gate('AND', (self.WIDTH/2+200, self.HEIGHT/2-150 - (15*i))))

        # # OR Gates
        # for i in range(num_or_gates):
        #     self.gates.append(Gate('OR', (self.WIDTH/2+300, self.HEIGHT/2-150 - (15*i))))

        # # NOT Gates
        # for i in range(num_not_gates):
        #     self.gates.append(Gate('NOT', (self.WIDTH/2+400, self.HEIGHT/2-150 - (15*i))))

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

        self.save_interval = 180000  # Save every 3 minutes (180 seconds)
        self.last_save_time = pygame.time.get_ticks()
    
    def get_current_datetime(self):
        now = datetime.now()
        formatted_datetime = now.strftime("%m/%d/%Y %H:%M")
        return formatted_datetime

    def run(self):
        while True:
            events = pygame.event.get()
            for event in events:
                if event.type == pygame.QUIT:
                    pygame.quit()
                    sys.exit()

            save_file = 'hsa05171_save_game.pickle'

            # Press "S" to save the game
            keys = pygame.key.get_pressed()
            if keys[pygame.K_s]:
                game_state = SaveState()
                game_state.save_objects(self.objects_group)
                print(f'Loaded {save_file}')
                save_game_state(game_state, save_file)
                self.last_save_time = current_time

            # Automatically saves the game every 2 minutes
            current_time = pygame.time.get_ticks()
            if current_time - self.last_save_time >= self.save_interval:
                print(f'Autosaving to {save_file}: {self.get_current_datetime()}')
                game_state = SaveState()
                game_state.save_objects(self.objects_group)
                save_game_state(game_state, save_file)
                self.last_save_time = current_time  # Reset the last save time


            if keys[pygame.K_l]:
                # Clear existing objects
                self.objects_group.empty()
                self.nodes_group.empty()
                self.gates_group.empty()
                
                # Clear lists
                self.objects_group = pygame.sprite.Group()
                self.nodes_group = pygame.sprite.Group()
                self.gates_group = pygame.sprite.Group()

                # Clear nodes and gates
                self.nodes = []
                self.gates = []

                game_state = load_game_state(save_file)
                self.nodes_group, self.gates_group, self.objects_group = game_state.load_objects()
                # Update UUIDs for event handling
                self.uuids = {obj.id: obj for obj in self.objects_group}
                self.node_ids = {obj.id: obj for obj in self.nodes_group}
                self.gate_ids = {obj.id: obj for obj in self.gates_group}


            self.screen.fill("white")      
            for node in self.nodes_group:
                if node not in self.nodes:
                    self.nodes.append(node)
                # node.draw_connections()
                node.update_object(events, self.connections, self.gate_ids, self.node_ids, self.uuids, self)
                if node.is_drawing_line == True:
                    self.connections += 1  
            
            for gate in self.gates_group:
                if gate not in self.gates:
                    self.gates.append(gate)
                # gate.draw_connections()
                gate.update_object(events, self.connections, self.gate_ids, self.node_ids, self.uuids, self)
                if gate.is_drawing_line == True:
                    self.connections += 1

            # Reset the number of connections if the 1 key is not pressed
            if self.connections != 0 and not keys[pygame.K_1]:
                self.connections = 0
            self.nodes_group.draw(self.screen)

            for node in self.nodes_group:
                node.draw_lock((node.position[0] - 37, node.position[1]))

            mouse_pos = pygame.mouse.get_pos()

            if keys[pygame.K_DELETE]:
                for object in self.objects_group:
                    if object.rect.collidepoint(mouse_pos):
                        if object in self.uuids:
                            del self.uuids[object.id]
                        if object in self.gate_ids:
                            del self.gate_ids[object.id]
                        if object in self.node_ids:
                            del self.gate_ids[object.id]
                        object.kill()
                        object.remove_all_connections(self.objects_group, just_killed=True)

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
