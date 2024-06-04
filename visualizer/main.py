import pygame, sys
from node import Node
from gate import Gate
from display_box import DisplayBox
import random
import pickle
import read_rule_file
from datetime import datetime
import math

class SaveState():
    def __init__(self):
        self.object_types = []
        self.object_names = []
        self.object_positions = []
        self.object_color = []
        self.object_ids = []
        self.object_image = []
        self.object_states = []
        self.outgoing_connections = []
        self.incoming_connections = []
        self.active_connections = []
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
            self.object_states.append(object.state)

            object_outgoing_connections = []
            object_incoming_connections = []
            object_active_connections = []

            for connected_object in object.outgoing_connections:
                object_outgoing_connections.append(connected_object.id)
            
            for active_object in object.active_incoming_connections:
                object_active_connections.append(active_object.id)

            for connected_object in object.incoming_connections:
                object_incoming_connections.append(connected_object.id)

            self.outgoing_connections.append(object_outgoing_connections)
            self.incoming_connections.append(object_incoming_connections)
            self.active_connections.append(object_active_connections)
    
    def load_objects(self):
        # Instantiate all objects and map them
        for i, object_type in enumerate(self.object_types):
            if object_type == "Node":
                node = Node(self.object_names[i], self.object_positions[i], self.object_color[i])
                node.id = self.object_ids[i]
                node.rect = node.image.get_rect(center=node.position)
                node.outgoing_connections = self.outgoing_connections[i]
                node.incoming_connections = self.incoming_connections[i]
                node.active_incoming_connections = self.active_connections[i]
                node.state = self.object_states[i]

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
                gate.active_incoming_connections = self.active_connections[i]
                gate.state = self.object_states[i]

                if gate not in self.gates:
                    self.gates.append(gate)
                
                if gate.id not in self.uuid_dict:
                    self.uuid_dict[gate.id] = gate

        # Restore connections using the UUID dictionary
        for i, object_type in enumerate(self.object_types):
            object = self.uuid_dict[self.object_ids[i]]
            object.outgoing_connections = [self.uuid_dict[uid] for uid in self.outgoing_connections[i]]
            object.incoming_connections = [self.uuid_dict[uid] for uid in self.incoming_connections[i]]
            object.active_incoming_connections = [self.uuid_dict[uid] for uid in self.active_connections[i]]

        # Group the objects into sprite groups
        nodes_group = pygame.sprite.Group(*self.nodes)
        gates_group = pygame.sprite.Group(*self.gates)
        objects_group = pygame.sprite.Group([gates_group, nodes_group])

        return nodes_group, gates_group, objects_group



        

class Game:
    def __init__(self):

        pygame.init()
        self.WIDTH = 2560
        self.HEIGHT = 1440

        self.screen = pygame.display.set_mode((self.WIDTH,self.HEIGHT))
        pygame.display.set_caption('Graph visualizer')
        self.clock = pygame.time.Clock()

        self.selection_box_start = None
        self.selection_box_end = None
        self.selected_objects = pygame.sprite.Group()
        self.dragging_selection = False
        self.moving_objects = False

        self.save_file = 'hsa05171_save_game_huge.pickle'
        self.mouse_pos = pygame.mouse.get_pos()
        self.keys = pygame.key.get_pressed()

        self.uuids = {}
        self.last_mouse_pos = None

        self.is_panning = False  # Flag to track panning state
        self.pan_start_pos = (0, 0)  # Store the initial position when panning starts
        self.offset = [0, 0]  # Offset for panning

        # num_nodes = 10
        # num_and_gates = 10
        # num_or_gates = 10
        # num_not_gates = 4

        # relative_abundance = False

        self.nodes = []
        self.gates = []

        # read_rule_file.create_nodes_and_gates('visualizer/04670.txt', self)
        read_rule_file.create_nodes_and_gates('visualizer/05171_huge.txt', self)

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
    
    def save_game_state(self):
        game_state = SaveState()
        game_state.save_objects(self.objects_group)
        try:
            with open(self.save_file, 'wb') as file:
                pickle.dump(game_state, file)
                print("Game state saved successfully!")
        except IOError:
            print("Error: Unable to save game state.")
        current_time = pygame.time.get_ticks()
        self.last_save_time = current_time

    def load_game_state(self):
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

        # Load the saved pickle file
        try:
            with open(self.save_file, 'rb') as file:
                game_state = pickle.load(file)
                print("Game state loaded successfully!")

        except (IOError, pickle.UnpicklingError):
            print("Error: Unable to load game state.")

        # Load the objects
        self.nodes_group, self.gates_group, self.objects_group = game_state.load_objects()

        # Update UUIDs for event handling
        self.uuids = {obj.id: obj for obj in self.objects_group}
        self.node_ids = {obj.id: obj for obj in self.nodes_group}
        self.gate_ids = {obj.id: obj for obj in self.gates_group}
    
    def autosave(self):
        # Automatically saves the game every 2 minutes
        current_time = pygame.time.get_ticks()
        if current_time - self.last_save_time >= self.save_interval:
            print(f'Autosaving to {self.save_file}: {self.get_current_datetime()}')
            self.save_game_state(self.save_file)
    
    def delete_all_connections(self):
        for object in self.objects_group:
            if object.rect.collidepoint(self.mouse_pos):
                if object in self.uuids:
                    del self.uuids[object.id]
                if object in self.gate_ids:
                    del self.gate_ids[object.id]
                if object in self.node_ids:
                    del self.gate_ids[object.id]
                object.kill()
                object.remove_all_connections(self.objects_group, just_killed=True)

    def draw_selection_box(self):
        if self.selection_box_start and self.selection_box_end:
            x1, y1 = self.selection_box_start
            x2, y2 = self.selection_box_end
            rect = pygame.Rect(min(x1, x2), min(y1, y2), abs(x2 - x1), abs(y2 - y1))
            pygame.draw.rect(self.screen, (0, 255, 0), rect, 2)

    def update_selected_objects(self):
        self.selected_objects.empty()
        if self.selection_box_start and self.selection_box_end:
            x1, y1 = self.selection_box_start
            x2, y2 = self.selection_box_end
            selection_rect = pygame.Rect(min(x1, x2), min(y1, y2), abs(x2 - x1), abs(y2 - y1))
            for obj in self.objects_group:
                if selection_rect.colliderect(obj.rect):
                    self.selected_objects.add(obj)

    def check_selected_objects(self, mouse_pos):
        mouse_pos_difference = pygame.Vector2(mouse_pos) - pygame.Vector2(self.last_mouse_pos)
        for obj in self.selected_objects:
            if obj.rect.collidepoint(mouse_pos):
                # obj.moving = True
                # self.node_being_moved = obj
                # obj.drag_offset = obj.position - pygame.Vector2(mouse_pos)
                obj.position = obj.position - mouse_pos_difference
                break

    def highlight_selected_objects(self):
        for obj in self.selected_objects:
            # Draw a golden circle around the object
            pygame.draw.circle(self.screen, (255, 215, 0), obj.rect.center, max(obj.rect.width, obj.rect.height) // 2 + 10, 3)
    
    def auto_organize_nodes(self):
        gravity = 0.01
        repulsive_force_modifier = 50
        attractive_force_modifier = 0.01

        # There should be a repulsive force away from 
        # Pull objects toward the center of the screen
        center_x = self.WIDTH / 2
        center_y = self.HEIGHT / 2
        center = pygame.Vector2(center_x, center_y)

        for obj in self.objects_group:
            
            dir_to_center = center - obj.position
            distance_to_center = dir_to_center.length()

            if distance_to_center != 0:
                dir_to_center.normalize_ip()

            attraction = gravity * ((math.sqrt(distance_to_center)**1.5))
            repulsion = -30/math.sqrt(distance_to_center)

            # Move the object towards the center
            obj.position += dir_to_center * (attraction + repulsion)
            obj.rect.center = obj.position      

            repulsive_force = 1 + len(obj.incoming_connections) * repulsive_force_modifier
            attractive_force = 1 + len(obj.outgoing_connections) * attractive_force_modifier

            # Push incoming connections away from the node
            for conn_obj in obj.incoming_connections + obj.outgoing_connections:
                # Calculate the direction vector between the centers of the two sprites
                direction_vector = pygame.math.Vector2(obj.rect.centerx - conn_obj.rect.centerx + 0.001,
                                                        obj.rect.centery - conn_obj.rect.centery + 0.001)

                # Find the distance between the objects
                distance = direction_vector.length()
                # Calculate the overlap distance
                overlap = distance - (obj.rect.width / 2 + conn_obj.rect.width / 2 + 5)

                # Make sure that nodes don't collide a ton
                if not pygame.sprite.collide_rect(obj, conn_obj) and distance > 20:
                    
                    conn_repulsion = -repulsive_force/(math.sqrt(distance))
                    conn_attraction = attractive_force * (math.sqrt(distance)**1.1)
                    move_vector = direction_vector.normalize() * (conn_attraction + conn_repulsion)
                    # Move the object by the move_vector
                    conn_obj.position += move_vector
                    conn_obj.rect.center = conn_obj.position

    def run(self):
        while True:
            events = pygame.event.get()
            for event in events:
                if event.type == pygame.QUIT:
                    pygame.quit()
                    sys.exit()
                elif event.type == pygame.MOUSEBUTTONDOWN:
                    if event.button == 3:  # Right mouse button for selection box
                        self.selection_box_start = pygame.Vector2(event.pos)
                        self.selection_box_end = pygame.Vector2(event.pos)
                        self.dragging_selection = True
                    elif event.button == 1:  # Left mouse button for moving
                        self.check_selected_objects(event.pos)

                    elif event.button == 2:  # Middle mouse button
                        self.is_panning = True
                        self.pan_start_pos = pygame.mouse.get_pos()

                elif event.type == pygame.MOUSEBUTTONUP:
                    if event.button == 3 and self.dragging_selection:
                        self.dragging_selection = False
                        self.update_selected_objects()
                        self.selection_box_start = None  # Reset the selection box start
                        self.selection_box_end = None    # Reset the selection box end
                    elif event.button == 1:
                        for obj in self.selected_objects:
                            obj.moving = False
                        self.node_being_moved = None

                    elif event.button == 2:  # Middle mouse button
                        self.is_panning = False

                elif event.type == pygame.MOUSEMOTION:
                    if self.dragging_selection:
                        self.selection_box_end = pygame.Vector2(event.pos)

                    elif self.node_being_moved:
                        mouse_pos_difference = pygame.Vector2(pygame.mouse.get_pos()) - pygame.Vector2(self.last_mouse_pos)
                        for obj in self.selected_objects:
                            obj.position = obj.position + mouse_pos_difference
                            object.rect.center = obj.position
                        self.last_mouse_pos = pygame.mouse.get_pos()  # Update the last mouse position

                    elif self.is_panning:
                        mouse_pos_difference = pygame.Vector2(pygame.mouse.get_pos()) - pygame.Vector2(self.last_mouse_pos)
                        for object in self.objects_group:
                            object.position = object.position + mouse_pos_difference
                            object.rect.center = object.position
                        self.last_mouse_pos = pygame.mouse.get_pos()  # Update the last mouse position
                
            self.mouse_pos = pygame.mouse.get_pos()
            self.keys = pygame.key.get_pressed()

            if self.keys[pygame.K_s]:
                self.save_game_state()

            if self.keys[pygame.K_l]:
                self.load_game_state()
            
            if self.keys[pygame.K_DELETE]:
                self.delete_all_connections()
            
            if self.keys[pygame.K_SPACE]:
                self.auto_organize_nodes()

            self.screen.fill("white")      

            

            for object in self.objects_group:

                if object.is_node and object not in self.nodes:
                    self.nodes.append(object)

                elif not object.is_node and object not in self.gates:
                    self.gates.append(object)

                object.update_object(events, self.connections, self.gate_ids, self.node_ids, self.uuids, self, self.keys, self.mouse_pos)
                if object.is_drawing_line:
                    self.connections += 1
                
                # object.position = (object.position[0] + self.offset[0], object.position[1] + self.offset[1])
            
            # if not self.keys[pygame.K_SPACE]:
            # for node in self.nodes_group:
            #     node.move_if_colliding(self.nodes_group)  # Move and resolve collisions


            # Reset the number of connections if the 1 key is not pressed
            if self.connections != 0 and not self.keys[pygame.K_1]:
                self.connections = 0

            self.nodes_group.draw(self.screen)

            # Draw the lock for the object (after drawn to screen so lock is in front)
            for object in self.nodes_group:
                object.draw_lock()

            # if self.keys[pygame.K_TAB]:
            #     self.display_box.display_text('Node State', (self.WIDTH / 4 + 25, self.HEIGHT - 400))
            #     self.display_box.display_text(f'Update {len(self.states)}', (self.WIDTH / 4 + 25, self.HEIGHT - 375))

            #     if self.nodes[0].update_num <= 35:
            #         self.states[self.nodes[0].update_num] = []
            #         for node_num, node in enumerate(self.nodes_group):
            #             position_adjustment = 15 * node_num
            #             update_adjustment = 10 * self.nodes[0].update_num - 150
            #             self.states[self.nodes[0].update_num].append((node.state, (self.WIDTH / 4 + update_adjustment, self.HEIGHT - 350 + position_adjustment)))

            #     # Displaying states
            #     for update in self.states:
            #         for node_num, node in enumerate(self.nodes_group):
            #             position_adjustment = 25 * node_num
            #             update_adjustment = 15 * len(self.states) - 100
            #             self.display_box.display_text(f'{self.states[update][node_num][0]}', self.states[update][node_num][1])

            # else:
            #     self.states = {}

            self.draw_selection_box()  # Draw the selection box
            self.highlight_selected_objects()  # Highlight selected objects

            # Update the last mouse pos
            self.last_mouse_pos = self.mouse_pos
            pygame.display.update()
            self.clock.tick(60)

if __name__ == "__main__":
    game = Game()
    game.run()
