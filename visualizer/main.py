import pygame, sys
from node import Node
from gate import Gate
from display_box import DisplayBox

class Game:
    def __init__(self):

        pygame.init()
        self.WIDTH = 1080
        self.HEIGHT = 720

        self.screen = pygame.display.set_mode((self.WIDTH,self.HEIGHT))
        pygame.display.set_caption('Graph visualizer')
        self.clock = pygame.time.Clock()

        num_nodes = 7
        num_and_gates = 1
        num_or_gates = 2
        num_not_gates = 1

        relative_abundance = False

        self.nodes = []
        # Nodes
        for i in range(1, num_nodes+1):
            self.nodes.append(Node(f'Gene {i}', (self.WIDTH/2+350,self.HEIGHT/2+(50*i)), "light blue"))

        # Change the size and color of the nodes for relative abundance. Manual entry required
        if relative_abundance:
            colors = [(239, 224, 215), (255, 208, 198), (235, 103, 75), (206, 219, 255), (255, 98, 98), "blue", "red"]
            sizes = [60, 35, 40, 35, 55, 65, 70]
            for node_num, node in enumerate(self.nodes):
                node.color = colors[node_num]
                node.set_size(sizes[node_num])
                node.font = pygame.font.Font(None, int(round(node.size / 2, 1)))

        self.gates = []
        # AND Gates
        for i in range(num_and_gates):
            self.gates.append(Gate('AND', (self.WIDTH/2+200, self.HEIGHT/2-150 - (50*i))))

        # OR Gates
        for i in range(num_or_gates):
            self.gates.append(Gate('OR', (self.WIDTH/2+300, self.HEIGHT/2-150 - (50*i))))

        # NOT Gates
        for i in range(num_not_gates):
            self.gates.append(Gate('NOT', (self.WIDTH/2+400, self.HEIGHT/2-150 - (50*i))))

        self.nodes_group = pygame.sprite.Group([self.nodes])
        self.gates_group = pygame.sprite.Group([self.gates])

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
                self.display_box.display_text('Node State', (self.WIDTH / 6 + 25, self.HEIGHT - 200))
                self.display_box.display_text(f'Update {len(self.states)}', (self.WIDTH / 6 + 25, self.HEIGHT - 175))

                if self.nodes[0].update_num <= 35:
                    self.states[self.nodes[0].update_num] = []
                    for node_num, node in enumerate(self.nodes_group):
                        position_adjustment = 15 * node_num
                        update_adjustment = 10 * self.nodes[0].update_num - 150
                        self.states[self.nodes[0].update_num].append((node.state, (self.WIDTH / 6 + update_adjustment, self.HEIGHT - 150 + position_adjustment)))


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