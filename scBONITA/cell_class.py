class Cell:
    def __init__(self, index):
    
        self.index: int = index
        self.name: str = ""
        self.attractor_dict: dict = {}
        self.expression: dict = {}
        self.simulation_results: list= []
        self.attractor_barcode: list = []
        self.transcription_factors: dict = {}
    
    def create_attractor_barcode(self):
        self.attractor_barcode = self.attractor_dict.values()

class CellPopulation:
    def __init__(self, cells):

        self.cells = cells