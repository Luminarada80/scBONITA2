from pathlib import Path
import pickle

class PathManager:
    def __init__(self, base_dir):

        # Set up the base directory
        self.base_dir = Path(base_dir)

    # Saving files
    def save_text_file(self, file_path, filename, content):
        path = self.base_dir / file_path / 'text_files' / f'{filename}.txt'
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('w', encoding='utf-8') as file:
            file.write(content)

    def save_png_file(self, file_path, filename, fig):
        path = self.base_dir / file_path / 'png_files' / f'{filename}.png'
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, bbox_inches='tight', format='png')
    
    def save_svg_file(self, path, filename, fig):
        path = self.base_dir / path / 'svg_files' / f'{filename}.svg'
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, bbox_inches='tight', format='svg')
    
    def save_pickle_file(self, obj, file_path):
        path = self.base_dir / 'pickle_files' / file_path
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('wb') as file:
            pickle.dump(obj, file)

