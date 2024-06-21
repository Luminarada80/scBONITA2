import pandas as pd
import matplotlib.pyplot as plt

# Data from the provided table
data = {
    "Num Genes": [10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
    "Avg": [99.34, 92.41, 88.26, 96.43, 87.96, 94.69, 96.78, 93.09, 95.51, 94.68],
    "Stdev": [0.3, 9.2, 3.4, 2.7, 4.0, 4.2, 1.9, 3.5, 2.4, 4.04],
    "Min": [98.5, 74.85, 76.7, 92.0, 78.22, 84.95, 88.64, 83.85, 86.66, 85.1],
    "Max": [100.0, 100.0, 96.4, 99.825, 97.86, 99.85, 99.91, 98.59, 99.92, 99.75]
}

# Convert to DataFrame
df = pd.DataFrame(data)

# Prepare the data for box plots
box_plot_data = []
for i in range(len(df)):
    num_genes = df.iloc[i]["Num Genes"]
    avg = df.iloc[i]["Avg"]
    stdev = df.iloc[i]["Stdev"]
    min_val = df.iloc[i]["Min"]
    max_val = df.iloc[i]["Max"]
    
    # Simulate the data for the box plot
    values = [avg - stdev, avg, avg + stdev]
    values = [min(max(val, min_val), max_val) for val in values]  # Ensure values are within min and max
    box_plot_data.append(values)

# Create the box plots with the further increased font size
plt.figure(figsize=(12, 8))
plt.boxplot(box_plot_data, positions=df["Num Genes"], widths=5)
plt.xticks(df["Num Genes"], fontsize=14)
plt.xlabel('Number of Genes', fontsize=18)
plt.ylabel('Accuracy', fontsize=18)
plt.title('scBONITA Rule Recovery Accuracy', fontsize=20)
plt.ylim(0, 100)
plt.grid(True)
plt.show()
