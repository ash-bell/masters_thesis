import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

df = pd.read_csv(<kaiju_relative_abundance>, sep ="\t")
pivot_df = df.pivot(index="file", columns="taxon_name", values="percent")
pivot_df.rename(columns={"cannot be assigned to a (non-viral) order": "Above Order Classification"}, inplace=True)
df = pivot_df.T.nlargest(10, list(pivot_df.T))
pivot_df = df.T

df = pivot_df.reindex(['BATS_1.kaiju.out', 'BATS_3.kaiju.out', 'BATS_5.kaiju.out',
                  'BATS_7.kaiju.out', 'BATS_9.kaiju.out', 'BATS_11.kaiju.out',
                  'BATS_2.kaiju.out', 'BATS_4.kaiju.out', 'BATS_6.kaiju.out',
                  'BATS_8.kaiju.out', 'BATS_10.kaiju.out', 'BATS_12.kaiju.out'])

df.loc[:,list(df)].plot.bar(stacked=True, figsize=(10,10), colormap=ListedColormap(sns.color_palette("colorblind", len(list(df)))))
plt.legend(bbox_to_anchor=(1.05, -0.07), ncol=4)
plt.ylim(0, 100)
plt.xlabel("Metagenome", fontsize=15)
plt.ylabel("Relative Order Percentage Abundance", fontsize=15)
plt.xticks(np.arange(12),('1', '3', '5', '7', '9', '11', '2', '4', '6', '8', '10', '12'), rotation='horizontal')
