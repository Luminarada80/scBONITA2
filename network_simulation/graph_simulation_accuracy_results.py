import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

dataframes = []
for file in os.listdir("/home/emoeller/github/scBONITA2/"):
    if file.endswith(".log"):
        print(file)
        file_path = os.path.join("/home/emoeller/github/scBONITA2/", file)

        # 1) Load the CSV
        df = pd.read_csv(file_path)
        dataframes.append(df)
        
df = pd.concat(dataframes)

# 2) Compute accuracy = 100 - error_pct
df["accuracy"] = 100 - df["error_pct"]

# 3) Ensure num_rules is integer (so seaborn orders it correctly)
df["num_rules"] = df["num_rules"].astype(int)

# 4) Create the box‐and‐whisker plot
plt.figure(figsize=(10, 6))
ax = sns.boxplot(
    data=df,
    x="num_rules",
    y="accuracy",
    hue="rule_set",
    palette={"scBONITA": "#1682b1", "simulated": "#cb5f17"}
)

plt.title("Rule Recovery Accuracy by Number of Genes")
plt.xlabel("Number of Genes")
plt.ylabel("Accuracy (%)")
plt.ylim((0, 100))

# Move legend outside top‐right
#   loc='upper left' and bbox_to_anchor=(1, 1) places it just outside the axes on the right.
ax.legend(
    title="Rule Set",
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    borderaxespad=0
)

plt.tight_layout()
plt.show()