import pandas as pd
import matplotlib.pyplot as plt

datafile = "../../kazer_data/merged_data.csv"

data = pd.read_csv(datafile, index_col=0)

array = data.values.flatten()

filtered_array = array[array != 0]

# Create a figure and a set of subplots
fig, axs = plt.subplots(1, 2, figsize=(12, 6))  # 2 rows, 1 column of subplots

# Plotting the histogram of the filtered array (non-zero values) with log scale
axs[0].hist(array, bins=250, log=False)
axs[0].set_title('Frequency of Expression Values')
axs[0].set_xlabel('Expression Value')
axs[0].set_ylabel('Frequency')
axs[0].set_xlim((-0.5, 5))


# Plotting the histogram of all values (including zeros) with log scale
axs[1].hist(filtered_array, bins=250, log=False)
axs[1].set_title('Zoomed in View of Expression Frequencies')
axs[1].set_xlabel('Expression Value')
axs[1].set_ylabel('Frequency')
axs[1].set_xlim((-0.5, 5))
axs[1].set_ylim((0, 0.30e6))

# Show the plots
plt.tight_layout()  # Adjusts subplot params so that subplots fit into the figure area.
plt.show()

# import matplotlib.pyplot as plt
# import numpy as np
# import math

# # Data Setup
# pathways = ['Apoptosis', 'Cell adhesion\nmolecules (CAMs)', 'Complement and \ncoagulation cascades', 'Glycolysis / Gluconeogenesis']
# data = {
#     'BON': [1.33, 0.10, 2.56, 0.05],
#     'CLP': [0.64, 0.64, 0.39, 0.51],
#     'CAM': [0.42, 1.46, 0.01, 1.44]
# }

# p_value = 0.05
# negative_log10_p_value = -math.log10(p_value)
# negative_log10_p_value

# # X locations for the groups
# ind = np.arange(len(pathways))
# width = 0.25  # width of the bars

# fig, ax = plt.subplots()

# # Creating bars
# bars1 = ax.bar(ind - width, data['BON'], width, label='BONITA')
# bars2 = ax.bar(ind, data['CLP'], width, label='CLIPPER')
# bars3 = ax.bar(ind + width, data['CAM'], width, label='CAMERA')

# # Adding labels, title, and custom x-axis tick labels
# ax.set_ylabel('-log10 p-values')
# ax.set_title('Differentially regulated pathways between infants with mild and severe RSV infection')
# ax.set_xticks(ind)
# ax.set_xticklabels(pathways)
# ax.legend()

# # Adding a horizontal dashed line at y=1.30
# ax.axhline(y=negative_log10_p_value, color='r', linestyle='--')
# ax.legend()

# # Adding some text for labels, title, and custom x-axis tick labels, etc.
# plt.tight_layout()

# plt.show()


