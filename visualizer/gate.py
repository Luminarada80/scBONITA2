import pygame
import math
from object import Object
import os

class Gate(Object):
    def __init__(self, gate_type, position):
        super().__init__(gate_type, position, color='grey')
        # Gate specific attributes
        self.gate_type = gate_type
        self.activation_threshold = self.incoming_connections

        # General attributes
        self.size = self.size/1.5 # This represents half the width/height of the gate for drawing purposes
        self.draw_object_function = self.draw_gate
        self.label_offset = 0

        # AND gate images
        self.inactive_and_image = pygame.image.load('visualizer/images/and.png').convert_alpha()  # Adjusted for gate drawing space
        self.active_and_image = pygame.image.load('visualizer/images/and.png').convert_alpha()

        # OR gate images
        self.inactive_or_image = pygame.image.load('visualizer/images/or.png').convert_alpha()  # Adjusted for gate drawing space
        self.active_or_image = pygame.image.load('visualizer/images/or.png').convert_alpha()

        # NOT gate images
        self.inactive_not_image = pygame.image.load('visualizer/images/not.png').convert_alpha()  # Adjusted for gate drawing space
        self.active_not_image = pygame.image.load('visualizer/images/not.png').convert_alpha()

        self.inactive_image, self.active_image = self.choose_gate_image()

        self.inactive_scaled_image = pygame.transform.smoothscale(self.inactive_image, (self.size+5, self.size))
        self.active_scaled_image = pygame.transform.smoothscale(self.active_image, (self.size+5, self.size))
        
        # Initial image and rect setup
        self.image = self.inactive_image
        self.rect = self.image.get_rect(center=self.position)
        self.original_image = self.image  # This will hold the currently displayed image

        self.font_size = 16
        self.font = pygame.freetype.Font('visualizer/arial.ttf', self.font_size)

        self.angle = 0
    
    def choose_gate_image(self):
        if self.gate_type == 'AND':
            inactive_image = self.inactive_and_image
            active_image = self.active_and_image
            self.label_offset = 0

        elif self.gate_type == 'OR':
            inactive_image = self.inactive_or_image
            active_image = self.active_or_image
            self.label_offset = 0

        elif self.gate_type == 'NOT':
            inactive_image = self.inactive_not_image
            active_image = self.active_not_image
            self.label_offset = 15
        
        return inactive_image, active_image

    def draw_gate(self):
        # Decide which image to use based on state
        if self.state == 0:
            self.image = self.inactive_scaled_image
        elif self.state == 1:
            self.image = self.active_scaled_image
        # Update original image for rotation
        self.original_image = self.image.copy()
        # Apply rotation if needed
        self.image = self.rotate(self.angle)


        text_surface, text_rect = self.font.render(self.gate_type, pygame.Color('black'))

        if self.gate_type == 'NOT' and self.angle == 180:
            text_rect.center = (self.size/2 + self.label_offset + 5, self.size/2)
        else:
            text_rect.center = (self.size/2 - self.label_offset, self.size/2)
            
        self.image.blit(text_surface, text_rect)
        


        self.rect = self.image.get_rect(center=self.position)

        # Blit the rotated image
        self.display_surface.blit(self.image, self.rect)

    def adjust_rotation(self):
        # Directly adjust to the desired orientation for simplicity in this example
        # You could implement smoothing or gradual rotation here
        keys = pygame.key.get_pressed()
        mouse_pos = pygame.mouse.get_pos()

        if self.rect.collidepoint(mouse_pos):
            if keys[pygame.K_d]:
                angle = 0
                self.rotate(angle)
            
            elif keys[pygame.K_w]:
                angle = -90
                self.rotate(angle)
            
            elif keys[pygame.K_a]:
                angle = 180
                self.rotate(angle)
            
            elif keys[pygame.K_s]:
                angle = 90
                self.rotate(angle)

            
    def rotate(self, angle):
        """Rotate the sprite image to the given angle and apply flipping for 180 degrees."""
        self.angle = angle
        rotated_image = pygame.transform.rotate(self.original_image, -self.angle)

        if angle == 180:
            # Flip horizontally to mirror the image when rotated 180 degrees
            self.image = pygame.transform.flip(rotated_image, False, True)

        else:
            self.image = rotated_image
        
        return self.image

    def simulate_and_logic(self):
        # Get which objects are connected to this one        
        self.activation_threshold = len(self.incoming_connections)
    
    def simulate_or_logic(self):
        self.activation_threshold = 1
    
    def simulate_not_logic(self):
        if len(self.active_incoming_connections) > 0:
            self.activation_threshold = 999
        elif len(self.active_incoming_connections) == 0:
            self.activation_threshold = 0

    def update_object(self, events, connections, gate_ids, node_ids, uuid_dict, game):
        self.update(events, connections, gate_ids, node_ids, uuid_dict, self.rect, self.draw_gate, game)
        # self.calculate_desired_orientation()

        if self.gate_type == 'AND':
            self.simulate_and_logic()
        elif self.gate_type == 'OR':
            self.simulate_or_logic()
        elif self.gate_type == 'NOT':
            self.simulate_not_logic()

        self.adjust_rotation()
        self.update_activation_highlight(self.draw_gate)


