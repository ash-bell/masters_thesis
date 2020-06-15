import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import seaborn as sns

df = pd.read_csv("/content/BATS_order_summary.tsv", sep ="\t")
pivot_df = df.pivot(index="file", columns="taxon_name", values="percent")
pivot_df.rename(columns={"cannot be assigned to a (non-viral) order": "Above Order Classification"}, inplace=True)
df = pivot_df.T.nlargest(10, list(pivot_df.T))
pivot_df = df.T

df = pivot_df.reindex(['BATS_1.kaiju.out', 'BATS_3.kaiju.out', 'BATS_5.kaiju.out',
                  'BATS_7.kaiju.out', 'BATS_9.kaiju.out', 'BATS_11.kaiju.out',
                  'BATS_2.kaiju.out', 'BATS_4.kaiju.out', 'BATS_6.kaiju.out',
                  'BATS_8.kaiju.out', 'BATS_10.kaiju.out', 'BATS_12.kaiju.out'])

df.loc[:,list(df)].plot.bar(stacked=True, figsize=(10,10), colormap=ListedColormap(sns.color_palette("colorblind", len(list(df)))))
plt.legend(bbox_to_anchor=(1.05, -0.15), ncol=4)
plt.ylim(0, 100)
plt.xlabel("Metagenome", fontsize=15)
plt.ylabel("Relative Order Percentage Abundance", fontsize=15)
plt.xticks(np.arange(12),('1\n18:30\n8 JUL\n80m', '3\n06:00\n9 JUL\n80m', '5\n19:00\n9 JUL\n80m', '7\n06:00\n10 JUL\n80m', '9\n19:00\n10 JUL\n80m', '11\n06:00\n11 JUL\n80m',
                          '2\n18:30\n8 JUL\n200m', '4\n06:00\n9 JUL\n200m', '6\n19:00\n9 JUL\n200m', '8\n06:00\n10 JUL\n200m', '10\n19:00\n10 JUL\n200m', '12\n06:00\n11 JUL\n200m'), rotation='horizontal')
