import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

df = pd.read_csv("bins.csv", names = ["Completeness", "Redundancy", "Taxa", "Size"])

colours = {"k__Archaea":"#1f77b4", 
    "k__Bacteria":"#ff7f0e", 
    "p__Proteobacteria":"#2ca02c", 
    "root":"#d62728", 
    "p__Cyanobacteria":"#9467bd", 
    "o__Rickettsiales":"#8c564b", 
    "f__Rhodobacteraceae":"#e377c2", 
    "c__Gammaproteobacteria":"#7f7f7f", 
    "p__Euryarchaeota":"#bcbd22", 
    "s__algicola":"#17becf",
    "p__Actinobacteria":"#9467bd"
    }

fig, ax1 = plt.subplots(figsize=(15, 15), facecolor='w', edgecolor='k')

ax1.scatter(df.Completeness, df.Redundancy, c = df.Taxa.apply(lambda x: colours[x]), s = df.Size*250)
ax1.scatter(x=-1, y=-1, c="#1f77b4", s=1*250, label = "Archaea")
ax1.scatter(x=-1, y=-1, c="#ff7f0e", s=1*250, label = "Bacteria")
ax1.scatter(x=-1, y=-1, c="#2ca02c", s=1*250, label = "Proteobacteria")
ax1.scatter(x=-1, y=-1, c="#d62728", s=1*250, label = "Unknown")
ax1.scatter(x=-1, y=-1, c="#9467bd", s=1*250, label = "Cyanobacteria")
ax1.scatter(x=-1, y=-1, c="#8c564b", s=1*250, label = "Rickettsiales")
ax1.scatter(x=-1, y=-1, c="#e377c2", s=1*250, label = "Rhodobacteraceae")
ax1.scatter(x=-1, y=-1, c="#7f7f7f", s=1*250, label = "Gammaproteobacteria")
ax1.scatter(x=-1, y=-1, c="#bcbd22", s=1*250, label = "Euryarchaeota")
ax1.scatter(x=-1, y=-1, c="#17becf", s=1*250, label = "algicola")
ax1.scatter(x=-1, y=-1, c="#9467bd", s=1*250, label = "Actinobacteria")

ax1.add_patch(patches.Rectangle((90,0),10,5,linewidth=1,edgecolor='none',facecolor='darkgreen', alpha = 0.50, label = "High Quaility Genomes"))
ax1.add_patch(patches.Rectangle((50,5),50,5,linewidth=1,edgecolor='none',facecolor='gold', alpha = 0.50, label = "Medium Quaility Draft Genomes"))
ax1.set_xlim(0, 100)
ax1.set_ylim(0, 100)
ax1.set_xlabel("Percentage Completeness", fontsize= 15)
ax1.set_ylabel("Percentage Redundancy", fontsize=15)
ax1.tick_params(labelsize=15)
ax1.legend(fontsize=15, loc="upper left")

plt.show()
