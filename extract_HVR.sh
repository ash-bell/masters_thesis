##get 23S-5S ITS region from barrnap
for i in <all_SAR11_genomes>; do echo $i; barrnap $i --threads 16 --quiet >> rRNA_pos.tsv;  done
egrep "23S|5S" rRNA_pos.tsv > 23S_5S_only.tsv
cut -f1 23S_5S_only.tsv | awk 'array[$0]++' | grep -f - 23S_5S_only.tsv | grep "23S" | cut -f1,4,5 > 23S_unique.tsv
cut -f1 23S_5S_only.tsv | awk 'array[$0]++' | grep -f - 23S_5S_only.tsv | grep "5S" | cut -f1,4,5 > 5S_unique.tsv

```python
import pandas as pd
rRNA23S = pd.read_csv("23S_unique.tsv", sep = "\t", header = None)
rRNA5S = pd.read_csv("5S_unique.tsv", sep = "\t", header = None)
test = rRNA23S.merge(rRNA5S, on = [0])
test["sstart"]=test.min(axis=1)
test["send"]=test.max(axis=1)
df = test.drop(["1_x","2_x","1_y","2_y"], axis = 1)
df.to_csv("reg.bed", sep = "\t", index = None, header = None)
```

seqtk subseq <concatiated_fasta_all_SAR11> reg.bed > ITS_reg.fasta
