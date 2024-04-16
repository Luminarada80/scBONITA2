import pygame

class DisplayBox():
    def __init__(self):

        self.display_surface = pygame.display.get_surface()
    
        self.font = pygame.font.Font(None, 24)

    def display_text(self, text, position):
        if text == "1":
            text = self.font.render(text, True, "black")
        elif text == "0":
            text = self.font.render(text, True, (201, 201, 201))
        else:
            text = self.font.render(text, True, "black")

        text_rect = text.get_rect()

        text_rect.center = position

        self.display_surface.blit(text, text_rect)
