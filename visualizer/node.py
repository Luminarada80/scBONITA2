import pygame
import math
from object import Object

class Node(Object):
    def __init__(self, name, position, color):
        super().__init__(name, position, color)
        
        # Create a surface for drawing the node
        self.image = pygame.Surface((self.size * 2, self.size * 2), pygame.SRCALPHA)

        self.rect = self.image.get_rect(center=self.position)
        self.original_image = self.image.copy()

        self.activation_threshold = 1
        
        # Circle and text will be drawn on self.image
        self.draw_circle()
        self.draw_object_function = self.draw_circle

    def draw_circle(self):
        # Clear the image
        self.image.fill((0, 0, 0, 0))

        pygame.draw.circle(self.image, self.highlight_color, (self.size, self.size), self.highlight_size)
        pygame.draw.circle(self.image, self.color, (self.size, self.size), self.size-8)
        
        black_text = self.font.render(self.name, True, "black")

        black_text_rect = black_text.get_rect(center=(self.size, self.size))
        self.image.blit(black_text, black_text_rect)

        self.original_image = self.image.copy()  # Update the original image
    
    def set_size(self, new_size):
        self.size = new_size
        self.image = pygame.Surface((self.size * 2, self.size * 2), pygame.SRCALPHA)
        self.draw_circle()
        self.rect = self.image.get_rect(center=self.position)
        self.original_image = self.image.copy()

    
    def update_object(self, events, connections, nodes, game):
        self.update(events, connections, nodes, self.rect, self.draw_circle, game)
        self.update_activation_highlight(self.draw_circle)

        # Display the updated image to the surface
        # self.display_surface.blit(self.image, self.rect)
        


