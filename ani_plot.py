""" 

Script from https://github.com/nanoporetech/marine-phage-paper-scripts written by jbeaulaurier
I DO NOT claim this script as my own

"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

ani = pd.read_csv("/home/minion/projects/SAGs/fastani/trimmed_fastANI.tsv", 
                  sep = "\t", header = None, 
                  names = ["Query", "Reference", "ANI", "Genes_1", "Genes_2"])
data = ani[["Query", "Reference", "ANI"]]
pivot_table = data.pivot("Query", "Reference", "ANI")
pivot_table.drop("all_sags", inplace=True)
pivot_table.drop("all_sags", axis = 1, inplace=True)
plt.figure(figsize=(100,100))
plt.xlabel("Query")
plt.ylabel("Reference")
plt.title("ANI")
sns.heatmap(pivot_table, cmap="YlOrRd")
plt.savefig("/home/minion/Desktop/fastANI_heatmap.png", dpi=150)
