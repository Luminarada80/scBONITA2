import pygame
import math
from object import Object
import pygame.freetype

class Node(Object):
    def __init__(self, name, position, color):
        super().__init__(name, position, color)

        pygame.freetype.init()
        
        # Create a surface for drawing the node
        self.image = pygame.Surface((self.size/2, self.size/2), pygame.SRCALPHA)

        self.inactive_image = pygame.image.load('visualizer/images/node_inactive.png').convert_alpha()
        self.active_image = pygame.image.load('visualizer/images/node_active.png').convert_alpha()

        self.inactive_scaled_image = pygame.transform.scale(self.inactive_image, (self.size, self.size))
        self.active_scaled_image = pygame.transform.scale(self.active_image, (self.size, self.size))

        self.activation_threshold = 1

        self.image = self.inactive_image
        self.rect = self.image.get_rect(center=self.position)
        self.original_image = self.image

        self.font_size = 32
        self.font = pygame.freetype.Font(None, self.font_size)

        # Circle and text will be drawn on self.image
        self.draw_circle()
        self.draw_object_function = self.draw_circle

        self.angle = 0

    def draw_circle(self):
        # Clear the image
        if self.state == 0:
            self.image = self.inactive_scaled_image
        elif self.state == 1:
            self.image = self.active_scaled_image


        # Apply rotation if needed
        self.rect = self.image.get_rect(center=self.position)
        
        text_surface, text_rect = self.font.render(self.name, pygame.Color('black'))
        text_rect.center = (self.size/2, self.size/2)
        self.image.blit(text_surface, text_rect)

        if self.name == '':
            activation_surface, activation_rect = self.font.render(str(self.state), pygame.Color('black'))
            activation_rect.center = (self.size/2, self.size/2)
            self.image.blit(activation_surface, activation_rect)

        self.original_image = self.image.copy()  # Update the original image
    
    def set_size(self, new_size):
        self.size = new_size
        self.image = pygame.Surface((self.size, self.size), pygame.SRCALPHA)
        self.draw_circle()
        self.rect = self.image.get_rect(center=self.position)
        self.original_image = self.image.copy()

    
    def update_object(self, events, connections, nodes, game):
        self.update(events, connections, nodes, self.rect, self.draw_circle, game)
        self.update_activation_highlight(self.draw_circle)

        # Display the updated image to the surface
        # self.display_surface.blit(self.image, self.rect)
        


