import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

data = {"Gene":[], "Expression Change":[]}
with open('/home/emoeller/github/scBONITA/scBONITA/relative_abundance_output/george_hiv/Healthy_vs_HIV/text_files/hsa04670_george_hiv_Healthy_vs_HIV_relative_abundance.txt', 'r') as abundance_file:
    next(abundance_file)
    for line in abundance_file:
        line = line.strip().split(',')
        gene_name = line[0]
        relative_abundance = line[2]

        print(gene_name, relative_abundance)
        data["Gene"].append(gene_name)
        data["Expression Change"].append(float(relative_abundance))


# # Creating DataFrame
# df = pd.DataFrame(data)

# # Setting up the matplotlib figure
# plt.figure(figsize=(10, 8))

# # Creating a heatmap
# sns.heatmap(df.pivot("Gene", "Gene", "Expression Change"), cmap="coolwarm", cbar_kws={'label': 'Expression Change'})

# # Adjustments and labels
# plt.title('Gene Expression Change Heatmap\n(Negative values: Lower in HIV, Positive values: Higher in HIV)')
# plt.xlabel('')
# plt.ylabel('')
# plt.xticks([])  # Remove x-axis labels for clarity since they are redundant
# plt.yticks(rotation=0)  # Ensure gene names are readable

# plt.tight_layout()
# plt.show()

# import matplotlib.pyplot as plt
# import seaborn as sns
# import pandas as pd

# Gene expression data
data = {
    "Gene": ["RAP1A", "RAPGEF4", "ROCK2", "NCF2", "NCF1", "CXCR4", "SIPA1", "ITGAL", "MYL7", "MSN", "PLCG2", "ITGA4", 
             "RASSF5", "MYL9", "CYBA", "VCAM1", "MYLPF", "PRKCG", "ICAM1", "PLCG1", "ITGB2", "ITGB1", "RAP1B", "RAC1", 
             "JAM3", "ROCK1", "ARHGAP35", "CDH5", "ESAM", "RAPGEF3", "PECAM1", "MYL5", "NCF4", "EZR", "ARHGAP5", 
             "BCAR1", "MYL2", "PRKCA", "PXN", "CTNNB1", "MMP9", "MYL10", "OCLN", "PTK2B", "PRKCB", "JAM2", "MMP2", 
             "F11R", "PTK2", "ITGAM", "CXCL12", "CTNND1", "THY1", "RHOA"],
    "Expression Change": [0.0, -0.06522410193887954, -0.28496530986929103, 0.0, -0.15957032663346504, 0.0, 
                          -0.38738203172445507, 0.0, -0.01469067419733941, -0.08885180140445947, -0.5968610797892853, 
                          0.0, -0.24896531109061323, 0.0010674356111297872, -0.5874797985673279, 0.0, 0.0, 0.0, 
                          -0.12595129550226722, -0.6671762633585933, 0.0, 0.0, -0.3415232147398935, 
                          -0.0005621135469364812, -0.672421536917975, 0.0, 0.0, 0.0, -0.0005621135469364812, 0.0, 0.0, 
                          -0.0017052711350571685, -0.0005621135469364812, -0.5618494360086701, -0.3172024149204049, 
                          -0.16708023689986898, 0.0, -0.03940473976829686, -0.5795870037882362, 0.0, 
                          -0.07316941357298118, 0.0, -0.03492676188705274, -0.0455501277961027, -0.2430034267247894, 
                          -0.38714692719912047, -0.5130698269600161, 0.0, -0.13901242054154037, -0.30552382662232036, 
                          0.0, -0.5779079910806839, -0.0675672085978641, -0.24945170266055885]
}

# Create a DataFrame from the data dictionary
df = pd.DataFrame(data)

# Since we only have one metric (expression change), we create a single column heatmap
# We first need to set the Gene column as the index
df.set_index('Gene', inplace=True)

# Create a single column dataframe for the heatmap
heatmap_data = pd.DataFrame(df['Expression Change'])

# Create the heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(heatmap_data, annot=True, cmap='coolwarm', cbar_kws={'label': 'Expression Change'})

# Adjust the plot
plt.title('Gene Expression Change Heatmap\n(Negative values: Lower in HIV, Positive values: Higher in HIV)')
plt.xlabel('Expression Change')
plt.ylabel('Gene')
plt.yticks(rotation=0)  # Make sure the gene names are horizontal for readability

plt.tight_layout()
plt.show()