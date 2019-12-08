#!/bin/python
import pandas as pd
from scipy.stats import hypergeom
ITS = pd.read_csv("<HVR_gene_to_COG.tsv>", sep = "\t", header = None, names = ("genes", "COG"))
ITS = ITS.fillna("S")
ITS = ITS.replace("_000000.*","",regex=True)
ITS = ITS.groupby(["genes"])
ITS = ITS["COG"].value_counts().unstack().fillna(0)

all_genomes = pd.read_csv("<whole_genome_to_COG.tsv>", sep = "\t", header = None, names = ("genes", "COG"))
all_genomes = all_genomes.fillna("S")
all_genomes = all_genomes.replace("_00000.*","",regex=True)
all_genomes = all_genomes.groupby(["genes"])
all_genomes = all_genomes["COG"].value_counts().unstack().fillna(0)

ITS["sample_size"] = ITS.sum(axis = 1)
all_genomes["pop_size"] = all_genomes.sum(axis = 1)


hypegeom = ITS.merge(all_genomes, on = ["genes"])
hypegeom["hypergeom"] = hypergeom.sf(hypegeom["M_x"]-1, hypegeom["pop_size"], hypegeom["M_y"], hypegeom["sample_size"])
test = hypegeom["hypergeom"]
