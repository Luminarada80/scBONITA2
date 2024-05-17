import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def create_heatmap(path, title):
    data = []
    gene_names = []
    with open(path, 'r') as file:
        for line in file:
            line = line.strip().split(',')
            gene_name = line[0]
            time_data = [int(i) for i in line[1:]]
            data.append(time_data)
            gene_names.append(gene_name)
    
    num_genes = len(data)
    num_time_steps = len(data[0])

    # Adjusting the data to fit the provided shape
    data_array = np.array(data).reshape((num_genes, num_time_steps))

    # Create a custom colormap
    cmap = mcolors.ListedColormap(['grey', 'green'])
    bounds = [0, 0.5, 1]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    # Create a heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(data_array, cmap=cmap, norm=norm, cbar=False, yticklabels=gene_names, xticklabels=True)
    plt.title(title)
    plt.xlabel('Time Steps')
    plt.ylabel('Genes')
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.show()

atherosclerosis_main_path = 'scBONITA/attractor_analysis_output/atherosclerosis_attractors/hsa05166_attractors/'
george_hiv_main_path = 'scBONITA/attractor_analysis_output/george_hiv_attractors/hsa04010_attractors/'


create_heatmap(f'{george_hiv_main_path}/attractor_15/george_hiv_hsa04010_simulated_attractor_15.txt', 'George HIV hsa04010 attractor 15')
create_heatmap(f'{george_hiv_main_path}/attractor_14/george_hiv_hsa04010_simulated_attractor_14.txt', 'George HIV hsa04010 attractor 14')
create_heatmap(f'{george_hiv_main_path}/attractor_5/george_hiv_hsa04010_simulated_attractor_5.txt', 'George HIV hsa04010 attractor 5')




