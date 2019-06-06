import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

ani = pd.read_csv("/home/minion/Desktop/fastANI.output", 
                  sep = "\t", header = None, 
                  names = ["Query", "Reference", "ANI", "Genes_1", "Genes_2"])
data = ani[["Query", "Reference", "ANI"]]
pivot_table = data.pivot("Query", "Reference", "ANI")
#pivot_table.fillna(value=80)
pivot_table.sort_index(level=0, ascending=True, inplace=True)
plt.figure(figsize=(20,20))
plt.xlabel("Query")
plt.ylabel("Reference")
plt.title("ANI")
sns.heatmap(pivot_table, square=True, cmap="PuBu")
plt.savefig("/home/minion/Desktop/fastANI_heatmap.png", dpi=150)
