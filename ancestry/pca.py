#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 14:15:20 2022

@author: S.S.
"""

import pandas as pd

# Cluster path root -- update to match your own environment.
REF_DATA_ROOT = "/path/to/reference_data/pgta/glimpse"

#add popID to ref PCA
coord = pd.read_csv(f"{REF_DATA_ROOT}/laser/laser.RefPC.coord", sep="\t")
pop = pd.read_csv(f"{REF_DATA_ROOT}/1kg/1kg_grch38_phenotype.tsv", sep="\t")

pop_dict = dict(zip(pop["Sample name"], pop["Superpopulation code"]))
coord["popID"] = coord["indivID"].apply(lambda x: pop_dict[x])
coord.to_csv(f"{REF_DATA_ROOT}/laser/laser.RefPC.pop.coord", sep="\t", index=False)
