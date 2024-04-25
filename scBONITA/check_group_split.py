

file1_path = '../../kazer_data/merged_data_AS-.csv'
file2_path = '../../kazer_data/merged_data_AS+.csv'

with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
    header1 = file1.readline().split(',')
    header2 = file2.readline().split(',')

    for names1, names2 in zip(header1, header2):
        print(f'Name1: {names1}, Name2: {names2}')
