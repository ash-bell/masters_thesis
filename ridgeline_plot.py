import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

df = pd.read_csv("combined_scaffold_coverage.tsv", sep = ",")
df.columns = ("200m 10", "80m 11", "200m 12", "80m 1", "200m 2", "80m 3", "200m 4", "80m 5", "200m 6", "80m 7", "200m 8", "80m 9")
df = df[["80m 1","80m 3","80m 5","80m 7","80m 9","80m 11","200m 2","200m 4","200m 6","200m 8","200m 10","200m 12"]]

data = pd.melt(df[df<50], var_name="g", value_name="x")
data.dropna(inplace=True)

# Initialize the FacetGrid object
col_pal=['#1f77b4','#1f77b4','#1f77b4','#1f77b4','#1f77b4','#1f77b4','#2ca02c','#2ca02c','#2ca02c','#2ca02c','#2ca02c','#2ca02c',]
pal = sns.color_palette(col_pal)
g = sns.FacetGrid(data, row="g", hue="g", aspect=15, height=1, palette=pal)

# Draw the densities in a few steps
g.map(sns.kdeplot, "x", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
g.map(sns.kdeplot, "x", clip_on=False, color="w", lw=2, bw=.2)
g.map(plt.axhline, y=0, lw=2, clip_on=False)


# Define and use a simple function to label the plot in axes coordinates
def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .5, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)


g.map(label, "x")

# Set the subplots to overlap
g.fig.subplots_adjust(hspace=-0.25)

# Remove axes details that don't play well with overlap
g.set_titles("")
g.set(yticks=[])
g.despine(bottom=True, left=True)
