#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose: find shortest euclidean distance between subset of ten mutations and data set with all analyzed mutations
@author: katerinapratsinis
Email: pratsink@ethz.ch
"""

import pandas as pd
import numpy as np


#adjust for files at hand
path = 'data/3_phyloTF_out/'
top = pd.read_csv(path + "phyloTF_SA_L4_ten_a.csv", index_col=0)
factor1 = pd.read_csv(path + "f1_phyloTF_SA_L4_comp_hf_a.csv", index_col=0)
factor3 = pd.read_csv(path + "f3_phyloTF_SA_L4_comp_hf_a.csv", index_col=0)
factor5 = pd.read_csv(path + "f5_phyloTF_SA_L4_comp_hf_a.csv", index_col=0)

column_names = top.columns

factor1 = factor1[column_names]
factor3 = factor3[column_names]
factor5 = factor5[column_names]

top_values = top.iloc[0, :]
factor1_values = factor1.iloc[0, :]
factor3_values = factor3.iloc[0, :]
factor5_values = factor5.iloc[0, :]

dist_factor1 = np.linalg.norm(top_values - factor1_values)
dist_factor3 = np.linalg.norm(top_values - factor3_values)
dist_factor5 = np.linalg.norm(top_values - factor5_values)

closest_factor = np.argmin([dist_factor1, dist_factor3, dist_factor5])

print(f"closest factor to ten: factor{closest_factor*2+1}")
