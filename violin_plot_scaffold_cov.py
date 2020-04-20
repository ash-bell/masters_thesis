import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("combined_scaffold_coverage.tsv", sep = "\t")

sns.set(style="whitegrid", rc={'figure.figsize':(18,12)}, font_scale = 1.5)

ax = sns.violinplot(data=df, orient="h", palette="colorblind", 
                    scale= "area", inner="quartile", 
                    order = ["1", "3", "5", "7", "9", "11", "2", "4", "6", "8", "10", "12"])
ax.set(xlabel="Scaffold Coverage", ylabel="Metagenome", xlim = (0,50))

plt.show()
