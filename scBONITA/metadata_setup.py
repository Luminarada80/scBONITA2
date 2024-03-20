
directory_path = 'george_data/hiv_dataset'

data_file = 'HIV_dataset_normalized_integrated_counts.csv'
meta_file = 'hiv_meta.txt'

with open(f'{directory_path}/{data_file}', 'r') as datafile, open(f'{directory_path}/{meta_file}', 'w') as meta_file:
    first_line = datafile.readline()
    first_line = first_line.strip().split(',')
    for cell_num, cell_name in enumerate(first_line[1:]):
        disease_state = cell_name.split('_')[1]
        meta_file.write(f'"{cell_num}" "{cell_name}" "{disease_state}"')
        meta_file.write('\n')

