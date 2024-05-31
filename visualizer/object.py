import pygame
import math
import uuid


class Object(pygame.sprite.Sprite):
    def __init__(self, name, position, color):
        super().__init__()

        self.size = 50
        self.name = name
        self.position = position
        self.color = color
        self.outgoing_connections = []
        self.incoming_connections = []
        self.active_incoming_connections = []
        self.num_active_incoming_connections = len(self.active_incoming_connections)
        self.activation_threshold = 0
        self.locked = False

        self.node_ids = []
        self.gate_ids = []
        
        self.id = uuid.uuid4()

        self.simulation_running = False

        self.display_surface = pygame.display.get_surface()
        self.moving = False

        self.cooldown = 200
        self.update_time = None
        self.can_update = True

        self.is_drawing_line = False

        self.state = 0

        self.highlight_color = "white"
        self.highlight_size = self.size - 4

        self.line_color = "black"

        self.update_num = 0
        self.velocity = pygame.math.Vector2(0, 0)

        self.keys = pygame.key.get_pressed()
        self.mouse_pos = pygame.mouse.get_pos()

    def draw_lock(self):
        if self.keys[pygame.K_3] and self.rect.collidepoint(self.mouse_pos):
            self.lock_state(self.rect)
        
        if self.locked:
            if not hasattr(self, 'lock_scaled_image'):
                self.lock_image = pygame.image.load('visualizer/images/lock.png').convert_alpha()
                self.lock_scaled_image = pygame.transform.smoothscale(self.lock_image, (self.size + 10, self.size + 10))
                self.lock_rect = self.lock_scaled_image.get_rect(center=self.position)
            self.display_surface.blit(self.lock_scaled_image, (self.position[0] - 37, self.position[1]))

    def draw_connections(self):
        if self.outgoing_connections:
            for obj in self.outgoing_connections:
                direction = (obj.position[0] - self.position[0], obj.position[1] - self.position[1])
                distance = math.sqrt(direction[0] ** 2 + direction[1] ** 2) + 1e-5
                unit_direction = (direction[0] / distance, direction[1] / distance)

                start_point = (self.position[0] + unit_direction[0] * self.size / 2, self.position[1] + unit_direction[1] * self.size / 2)
                end_point = (obj.position[0] - unit_direction[0] * obj.size / 2, obj.position[1] - unit_direction[1] * obj.size / 2)

                pygame.draw.line(self.display_surface, self.line_color, start_point, end_point, 3)
                
                arrow_length = 6
                arrow_degrees = math.radians(30)

                arrow_body_end = (end_point[0] - unit_direction[0] * arrow_length, end_point[1] - unit_direction[1] * arrow_length)
                angle = math.atan2(direction[1], direction[0])
                arrow_point_1 = (arrow_body_end[0] - arrow_length * math.cos(angle + arrow_degrees), 
                                 arrow_body_end[1] - arrow_length * math.sin(angle + arrow_degrees))
                
                arrow_point_2 = (arrow_body_end[0] - arrow_length * math.cos(angle - arrow_degrees), 
                                 arrow_body_end[1] - arrow_length * math.sin(angle - arrow_degrees))

                pygame.draw.polygon(self.display_surface, self.line_color, [end_point, arrow_point_1, arrow_point_2])

    def move(self, events, rect, game):
        for event in events:
            if event.type == pygame.MOUSEBUTTONDOWN and rect.collidepoint(self.mouse_pos) and game.node_being_moved is None and event.button == 1:
                self.moving = True
                game.node_being_moved = self

            elif event.type == pygame.MOUSEBUTTONUP and self.moving and event.button == 1:
                self.moving = False
                game.node_being_moved = None
                
            if self.moving:
                self.position = self.mouse_pos
                rect.center = self.position

    def line_with_arrow(self):
        direction = (self.mouse_pos[0] - self.position[0], self.mouse_pos[1] - self.position[1])
        distance = math.sqrt(direction[0] ** 2 + direction[1] ** 2) + 1e-5
        unit_direction = (direction[0] / distance, direction[1] / distance)

        node_edge = (self.position[0] + unit_direction[0] * self.size, self.position[1] + unit_direction[1] * self.size)
        mouse = (self.mouse_pos[0] - unit_direction[0], self.mouse_pos[1] - unit_direction[1])

        pygame.draw.line(self.display_surface, self.line_color, node_edge, mouse, 2)

        arrow_length = 7
        arrow_degrees = math.radians(30)

        angle = math.atan2(direction[1], direction[0])
        arrow_point_1 = (mouse[0] - arrow_length * math.cos(angle + arrow_degrees), 
                         mouse[1] - arrow_length * math.sin(angle + arrow_degrees))
        
        arrow_point_2 = (mouse[0] - arrow_length * math.cos(angle - arrow_degrees), 
                         mouse[1] - arrow_length * math.sin(angle - arrow_degrees))

        pygame.draw.polygon(self.display_surface, self.line_color, [mouse, arrow_point_1, arrow_point_2])

    def draw_edge(self, connections, rect):
        if self.keys[pygame.K_1] and rect.collidepoint(self.mouse_pos) and connections < 2:
            self.is_drawing_line = True
        
        if self.is_drawing_line:
            if not self.keys[pygame.K_1]:
                self.is_drawing_line = False
            else:
                self.line_with_arrow()

    def connect_to_object(self, events, objects):
        if self.is_drawing_line and self.keys[pygame.K_1]:
            for obj in objects:
                if obj.rect.collidepoint(self.mouse_pos):
                    for event in events:
                        if event.type == pygame.MOUSEBUTTONDOWN:
                            if obj not in self.outgoing_connections:
                                self.outgoing_connections.append(obj)
                            if self not in obj.incoming_connections:
                                obj.incoming_connections.append(self)
                            self.is_drawing_line = False

    def disconnect_objects(self, objects):
        if self.is_drawing_line and self.keys[pygame.K_t]:
            for obj in objects:
                if obj.rect.collidepoint(self.mouse_pos):
                    self.remove_connection(obj)
    
    def remove_all_connections(self, objects, just_killed):
        if self.keys[pygame.K_y] or just_killed:
            for obj in objects:
                if obj.rect.collidepoint(self.mouse_pos):
                    self.remove_connection(obj)

    def remove_connection(self, obj):
        if obj in self.outgoing_connections:
            self.outgoing_connections.remove(obj)
        
        if obj in self.active_incoming_connections:
            self.active_incoming_connections.remove(obj)
        
        if obj in self.incoming_connections:
            self.incoming_connections.remove(obj)
        
        if self in obj.outgoing_connections:
            obj.outgoing_connections.remove(self)
        
        if self in obj.active_incoming_connections:
            obj.active_incoming_connections.remove(self)
        
        if self in obj.incoming_connections:
            obj.incoming_connections.remove(self)

    def run_simulation(self, draw_object_function):
        if self.keys[pygame.K_TAB] and self.can_update:
            self.update_time = pygame.time.get_ticks()

            if self in self.incoming_connections:
                self.incoming_connections.remove(self)
            if self in self.outgoing_connections:
                self.outgoing_connections.remove(self)
            if self in self.active_incoming_connections:
                self.active_incoming_connections.remove(self)

            for obj in self.outgoing_connections:
                if not obj.locked:
                    if self.state == 1:
                        if self not in obj.active_incoming_connections:
                            obj.active_incoming_connections.append(self)
                    else:
                        if self in obj.active_incoming_connections:
                            obj.active_incoming_connections.remove(self)
                    obj.num_active_incoming_connections = len(obj.active_incoming_connections)
                    obj.state = 1 if obj.num_active_incoming_connections >= obj.activation_threshold else 0
                    draw_object_function()
            self.update_num += 1
            if self.is_node:
                self.can_update = False
        else:
            self.update_num = 0
    
    def auto_organize_nodes(self):
        if self.keys[pygame.K_SPACE]:

            # There should be a repulsive force away from 

            # Calculate the gravity
            # The gravity pulling toward the node will be equal to the number of connections the node has
            # Pull connected objects toward this object
            if len(self.incoming_connections) + len(self.outgoing_connections) > 0:
                gravity = len(self.incoming_connections) + len(self.outgoing_connections)
                for obj in self.incoming_connections + self.outgoing_connections:
                    if obj.is_node:

                        # Calculate the direction vector between the centers of the two sprites
                        direction_vector = pygame.math.Vector2(self.rect.centerx - obj.rect.centerx + 0.001,
                                                            self.rect.centery - obj.rect.centery + 0.001)

                        # Find the distance between the objects
                        distance = direction_vector.length()
                        # Calculate the overlap distance
                        overlap = distance - (self.rect.width / 2 + obj.rect.width / 2 + 5)
                        
                        # Make sure that nodes dont collide a ton
                        if not pygame.sprite.collide_rect(self, obj) and overlap > 1:
                            move_vector = (direction_vector.normalize() * gravity) * math.sqrt(distance)
                            # Move the object by the move_vector
                            obj.position += move_vector
                            obj.rect.center = obj.position
                    
                    

            


    def simulation_step_cooldown(self):
        if not self.can_update and pygame.time.get_ticks() - self.update_time >= self.cooldown:
            self.can_update = True

    def toggle_state(self, rect):
        if self.keys[pygame.K_2] and rect.collidepoint(self.mouse_pos) and self.can_update and type(self).__name__ == "Node":
            self.update_time = pygame.time.get_ticks()
            self.state = 1 - self.state
            self.can_update = False

    def lock_state(self, rect):
        if self.keys[pygame.K_3] and rect.collidepoint(self.mouse_pos) and self.can_update:
            self.update_time = pygame.time.get_ticks()
            self.locked = not self.locked
            self.can_update = False

    def update_activation_highlight(self, draw_object_function):
        self.highlight_color = "gold" if self.state == 1 else "black"
        self.highlight_size = self.size if self.state == 1 else self.size - 4
        draw_object_function()

    def move_if_colliding(self, objects):
        # Check for collisions with other objects
        for sprite in objects:
            if sprite != self and pygame.sprite.collide_rect(self, sprite):
                self.resolve_collision(sprite)

    def resolve_collision(self, other_sprite):
        # Calculate the direction vector between the centers of the two sprites
        direction_vector = pygame.math.Vector2(self.rect.centerx - other_sprite.rect.centerx + 0.001,
                                            self.rect.centery - other_sprite.rect.centery + 0.001)
        # Prevent division by zero
        distance = direction_vector.length() or 1
        
        # Calculate the overlap distance
        overlap = 0.5 * (distance - (self.rect.width / 2 + other_sprite.rect.width / 2 + 5))

        # If there is an overlap, resolve the collision
        if overlap < 0:
            adjustment_factor = 0.5
            move_vector = direction_vector.normalize() * (-overlap)
            self.position += move_vector * adjustment_factor
            other_sprite.position -= move_vector * adjustment_factor

        # Update the rect positions to match the new positions
        self.rect.center = self.position
        other_sprite.rect.center = other_sprite.position



    def request_uuid_objects(self, uuid_dict, object_uuids):
        return [uuid_dict.get(id) for id in object_uuids if id in uuid_dict]
            
    def update(self, events, connections, gate_ids, node_ids, uuid_dict, rect, draw_object_function, game, keys, mouse_pos):
        self.keys = keys
        self.mouse_pos = mouse_pos
        gates = self.request_uuid_objects(uuid_dict, gate_ids)
        nodes = self.request_uuid_objects(uuid_dict, node_ids)
        objects = nodes + gates
        self.auto_organize_nodes()
        self.move(events, rect, game)
        self.connect_to_object(events, objects)
        self.draw_connections()
        self.lock_state(rect)
        self.simulation_step_cooldown()
        self.draw_edge(connections, rect)
        self.disconnect_objects(objects)
        self.remove_all_connections(objects, False)
        self.run_simulation(draw_object_function)
        self.toggle_state(rect)
        