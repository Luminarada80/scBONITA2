import pygame
import math
import uuid


class Object(pygame.sprite.Sprite):
    def __init__(self, name, position, color):
        super().__init__()

        self.size = 100
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
        self.highlight_size = self.size-4

        self.line_color = "black"

        self.update_num = 0
        self.velocity = pygame.math.Vector2(0, 0)  # Assuming each sprite has velocity
    
    def draw_lock(self, position):
        if self.locked:
            self.lock_image = pygame.image.load('visualizer/images/lock.png').convert_alpha()  # Adjusted for gate drawing space

            self.lock_scaled_image = pygame.transform.smoothscale(self.lock_image, (self.size + 10, self.size + 10))
            # Update the lock_rect position before drawing
            self.lock_rect = self.lock_scaled_image.get_rect(center=self.position)
            # Draw the scaled lock image at the updated position
            self.display_surface.blit(self.lock_scaled_image, position)

    
       # Draw connections to nodes in the "connections" list
    def draw_connections(self):
        """Draw a white arrow between this object and all other objects in the self.outgoing_connections list"""

        if self.outgoing_connections is not None:
            for object in self.outgoing_connections:
                # Vector from this node to the connected node
                direction = (object.position[0] - self.position[0], object.position[1] - self.position[1])

                # Distance between the nodes
                distance = math.sqrt(direction[0] ** 2 + direction[1] ** 2) + 0.00001

                # Normalize the direction vector
                unit_direction = (direction[0] / distance, direction[1] / distance)

                # Calculate start and end points on the edges of the nodes, not the centers
                start_point = (self.position[0] + unit_direction[0] * self.size / 2, self.position[1] + unit_direction[1] * self.size / 2)
                end_point = (object.position[0] - unit_direction[0] * object.size / 2, object.position[1] - unit_direction[1] * object.size / 2)

                # Draw the line
                pygame.draw.line(self.display_surface, self.line_color, start_point, end_point, 4)
                
                # Assuming arrow_length and arrow_degrees define the size and angle of the arrowhead
                arrow_length = 12
                arrow_degrees = math.radians(30)

                # Calculate the end point of the arrow body (slightly before the actual end_point)
                arrow_body_end = (end_point[0] - unit_direction[0] * arrow_length, end_point[1] - unit_direction[1] * arrow_length)

                # Calculate points for the arrowhead
                angle = math.atan2(direction[1], direction[0])
                arrow_point_1 = (arrow_body_end[0] - arrow_length * math.cos(angle + arrow_degrees), 
                                arrow_body_end[1] - arrow_length * math.sin(angle + arrow_degrees))
                
                arrow_point_2 = (arrow_body_end[0] - arrow_length  * math.cos(angle - arrow_degrees), 
                                arrow_body_end[1] - arrow_length * math.sin(angle - arrow_degrees))

                # Draw the arrowhead
                pygame.draw.polygon(self.display_surface, self.line_color, [end_point, arrow_point_1, arrow_point_2])
        
    # Move the node by clicking and dragging with the mouse
    def move(self, events, rect, game):
        """
        Move the object by clicking and dragging with the mouse. Ensure only one object can be moved at a time.
        
        Keyword arguments:
        events -- the events variable from the main game loop (pygame.event.get())
        rect -- the rectangle defining the object's boundaries
        game -- reference to the game instance controlling the movement
        """

        mouse_pos = pygame.mouse.get_pos()

        for event in events:
            if event.type == pygame.MOUSEBUTTONDOWN and rect.collidepoint(mouse_pos):
                # Check if no other node is being moved
                if game.node_being_moved is None:
                    self.moving = True
                    game.node_being_moved = self  # Set this node as the one being moved
            elif event.type == pygame.MOUSEBUTTONUP:
                if self.moving:
                    self.moving = False
                    game.node_being_moved = None  # No node is being moved now
                
            if self.moving:
                self.position = pygame.mouse.get_pos()
                rect.center = self.position
    
    def line_with_arrow(self):
        mouse_pos = pygame.mouse.get_pos()
        adjustment = 0.0000001

        # Vector from this node to the connected node
        direction = (mouse_pos[0] - self.position[0] + adjustment, mouse_pos[1] - self.position[1] + adjustment)

        # Distance between the nodes
        distance = math.sqrt(direction[0] ** 2 + direction[1] ** 2)

        # Normalize the direction vector
        unit_direction = (direction[0] / distance, direction[1] / distance)

        node_edge = (self.position[0] + unit_direction[0] * self.size, self.position[1] + unit_direction[1] * self.size)
        mouse = (mouse_pos[0] - unit_direction[0], mouse_pos[1] - unit_direction[1])

        pygame.draw.line(self.display_surface, self.line_color, node_edge, mouse, 2)

        arrow_length = 10
        arrow_degrees = math.radians(30)

        # Calculate points for the arrowhead
        angle = math.atan2(direction[1], direction[0])
        arrow_point_1 = (mouse[0] - arrow_length * math.cos(angle + arrow_degrees), 
                        mouse[1] - arrow_length * math.sin(angle + arrow_degrees))
        
        arrow_point_2 = (mouse[0] - arrow_length * math.cos(angle - arrow_degrees), 
                        mouse[1] - arrow_length * math.sin(angle - arrow_degrees))

        # Draw the arrowhead
        pygame.draw.polygon(self.display_surface, self.line_color, [mouse, arrow_point_1, arrow_point_2])

    # Draw an arrow from this object by hovering over it with the mouse and holding down 1
    def draw_edge(self, connections, rect):
        mouse_pos = pygame.mouse.get_pos()
        keys = pygame.key.get_pressed()

        if keys[pygame.K_1] and rect.collidepoint(mouse_pos) and connections < 2:
            self.is_drawing_line = True
        
        if self.is_drawing_line == True:
            if not keys[pygame.K_1]:
                self.is_drawing_line = False
            
            else:
                self.line_with_arrow()
        
        
    
    # Connect to an object by pressing the '1' key while drawing a line and hovering over the object
    def connect_to_object(self, events, objects):
        """
        Connect the current object to another object by hovering over it with the mouse while
        drawing a line and pressing the '1' key
        """

        mouse_pos = pygame.mouse.get_pos()

        for object in objects:
            keys = pygame.key.get_pressed()

            if self.is_drawing_line and keys[pygame.K_1] and object.rect.collidepoint(mouse_pos):
                for event in events:
                    if event.type == pygame.MOUSEBUTTONDOWN:
                        if object not in self.incoming_connections:
                            self.outgoing_connections.append(object)
                            print(self.outgoing_connections)
                        if self not in object.incoming_connections:
                            object.incoming_connections.append(self)
                            print(object.incoming_connections)
                        self.is_drawing_line = False

    
    # Activate the remove_connection method if 'T' is pressed while drawing a line and hovering over another object
    def disconnect_objects(self, objects):
        """
        Activate the remove_connection method if 'T' is pressed while drawing 
        a line and hovering over another object
        """

        mouse_pos = pygame.mouse.get_pos()

        for object in objects:
            keys = pygame.key.get_pressed()

            if self.is_drawing_line and keys[pygame.K_t] and object.rect.collidepoint(mouse_pos):
                self.remove_connection(object)

    # Remove a connection between two nodes
    def remove_connection(self, object):
        """Remove an object from the self.outgoing_connections variable"""
        
        # Remove the connection to the other object from this objects connections list
        if object in self.outgoing_connections:
            self.outgoing_connections.remove(object)
        
        if object in self.active_incoming_connections:
            self.active_incoming_connections.remove(object)
        
        if object in self.incoming_connections:
            self.incoming_connections.remove(object)
        
        # Remove the connection to this object from the other objects connections list
        if self in object.outgoing_connections:
            object.outgoing_connections.remove(self)
        
        if self in object.active_incoming_connections:
            object.active_incoming_connections.remove(self)
        
        if self in object.incoming_connections:
            object.incoming_connections.remove(self)


    # Activate or deactivate the connections based on the connections between nodes
    def run_simulation(self, draw_object_function):
        """Activate the simulation by holding down the 'Tab' key"""

        keys = pygame.key.get_pressed()

        # I need to implement a timer to that it doesnt all update at once
        if keys[pygame.K_TAB]:
            if self.can_update:
                self.update_time = pygame.time.get_ticks()
                
                # Signal a 1 to the other object if this objects state is 1
                if self.state == 1:
                    for object in self.outgoing_connections:
                        # print(f'\t{self.name} locked: {self.locked}, {object.name} locked: {object.locked}')
                        if not object.locked:
                            
                            if self not in object.active_incoming_connections:
                                object.active_incoming_connections.append(self)
                            object.num_active_incoming_connections = len(object.active_incoming_connections)

                            if object.num_active_incoming_connections >= object.activation_threshold:
                                object.state = 1
                            else:
                                object.state = 0
                            draw_object_function()


                # Signal a 0 to the other object if this objects state is 0
                elif self.state == 0:
                    for object in self.outgoing_connections:
                        # print(f'\t{self.name} locked: {self.locked}, {object.name} locked: {object.locked}')
                        if not object.locked:
                            

                            if self in object.active_incoming_connections:
                                object.active_incoming_connections.remove(self)
                            object.num_active_incoming_connections = len(object.active_incoming_connections)

                            if object.num_active_incoming_connections >= object.activation_threshold:
                                object.state = 1
                            else:
                                object.state = 0

                            draw_object_function()

            
                self.update_num += 1

                # Only the gates can update (allows for logic to update correctly if multiple gates are strung together)
                if self in self.nodes_group:
                    self.can_update = False
        else:
            self.update_num = 0

    def simulation_step_cooldown(self):
        """Create a delay between steps in the simulation so that the simulation doesn't update at every frame"""

        current_time = pygame.time.get_ticks()

        if not self.can_update:
            if current_time - self.update_time >= self.cooldown:
                self.can_update = True

    # Press 2 while hovering over the object to activate it
    def toggle_state(self, rect):
        """Set the state of the object to 1 by hovering over it and pressing the '2' key"""
        mouse_pos = pygame.mouse.get_pos()
        keys = pygame.key.get_pressed()

        if keys[pygame.K_2] and rect.collidepoint(mouse_pos) and self.can_update:
            self.update_time = pygame.time.get_ticks()
            if self.state == 1:
                self.state = 0
                self.can_update = False

            elif self.state == 0:
                self.state = 1
                self.can_update = False
    
    def lock_state(self, rect):
        mouse_pos = pygame.mouse.get_pos()
        keys = pygame.key.get_pressed()

        if keys[pygame.K_3] and rect.collidepoint(mouse_pos) and self.can_update and not self.locked:
            self.update_time = pygame.time.get_ticks()
            self.locked = True
            self.can_update = False
            # print(f'{self.name} locked = {self.locked}')

        
        elif keys[pygame.K_3] and rect.collidepoint(mouse_pos) and self.can_update and self.locked:
            self.update_time = pygame.time.get_ticks()
            self.locked = False
            self.can_update = False
            # print(f'{self.name} locked = {self.locked}')

    def update_activation_highlight(self, draw_object_function):
        """
        Updates the highlight around the node based on self.state. 
        Calls the draw function again to re-draw the object to the surface
        """
        # Update the gate's appearance based on the current state
        if self.state == 0:
            self.highlight_color = "black"
            self.highlight_size = self.size-4
        
        elif self.state == 1:
            self.highlight_color = "gold"
            self.highlight_size = self.size

        # Redraw the gate and text
        draw_object_function()

    def move_if_colliding(self, objects):
        # Move sprite by its velocity
        self.rect.x += self.velocity.x
        self.rect.y += self.velocity.y
        
        # Check for collisions with all other sprites in the group
        for sprite in objects:
            if sprite == self:
                continue  # Skip checking collision with itself
            if pygame.sprite.collide_rect(self, sprite):
                # Handle collision
                self.resolve_collision(sprite)

    def resolve_collision(self, other_sprite):
        # Calculate the vector between the sprite centers
        direction_vector = pygame.math.Vector2(self.rect.centerx - other_sprite.rect.centerx,
                                            self.rect.centery - other_sprite.rect.centery)
        distance = direction_vector.length()
        if distance == 0:  # Prevent division by zero
            direction_vector = pygame.math.Vector2(1, 0)  # Arbitrary direction
            distance = 1

        # Calculate overlap, considering a small buffer for separation
        overlap = 0.5 * (distance - (self.rect.width / 2 + other_sprite.rect.width / 2 + 5))  # 5 pixels buffer

        # Adjust positions based on overlap and direction vector, ensuring they move apart
        if overlap < 0:  # Only adjust positions if sprites are overlapping
            move_vector = direction_vector.normalize() * (-overlap)
            self.rect.x += move_vector.x
            self.rect.y += move_vector.y
            other_sprite.rect.x -= move_vector.x
            other_sprite.rect.y -= move_vector.y

        # Optionally, reverse or adjust velocities for a more dynamic response
        self.velocity = -self.velocity
        other_sprite.velocity = -other_sprite.velocity
    
    def request_uuid_objects(self, uuid_dict, object_uuids):
        """Requests the objects from main.py in a group"""
        objects = []
        for id in object_uuids:
            objects.append(uuid_dict[id])

        return objects
            
    def update(self, events, connections, gate_ids, node_ids, uuid_dict, rect, draw_object_function, game):
        gates = self.request_uuid_objects(uuid_dict, gate_ids)
        nodes = self.request_uuid_objects(uuid_dict, node_ids)
        objects = nodes + gates
        self.draw_connections()
        self.lock_state(rect)
        self.move(events, rect, game)
        self.simulation_step_cooldown()
        self.draw_edge(connections, rect)
        self.connect_to_object(events, objects)
        self.disconnect_objects(objects)
        self.run_simulation(draw_object_function)
        # self.move_if_colliding(objects)
        self.toggle_state(rect)
