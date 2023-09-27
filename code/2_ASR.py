#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose: reconstruct ancestral state of mutations
@author: katerinapratsinis
Email: pratsink@ethz.ch
"""
import os
import time
import pandas as pd
import numpy as np
from pastml.ml import JOINT
from pastml.acr import pastml_pipeline

BASE_PATH = 'data/'  # path to the directory with the data
bs_nr = 10 # only necessary for comp_hf


def reconstruct(tree_file, align_file, work_dir, ds, lineage, size):

    df = pd.read_csv(align_file, index_col=0)
    traits = [column for column in df.columns]

    html_compressed = os.path.join(work_dir, f'pypastml_{ds}_L{lineage}_{size}_map.html')    # Path to the output compressed map visualisation
    html = os.path.join(work_dir, f'pypastml_{ds}_L{lineage}_{size}_tree.html')    # (Optional) path to the output tree visualisation

    pastml_pipeline(tree=tree_file, data=align_file, data_sep=',', columns=traits, name_column=traits[0],
                    # , prediction_method=pastml.ml.MAP, forced_joint=True
                    html_compressed=html_compressed, html=html, work_dir=work_dir, verbose=True) #, prediction_method=pastml.ml.MAP)#, forced_joint=True)


def unique_ancestral_states(pastml_path, addon):
    """ creates a df of unique ancestral states from marginal probs """
    cas_file = os.path.join(pastml_path, 'combined_ancestral_states.tab')
    cas = pd.read_table(cas_file)
    unique_vals = cas.iloc[:, 0].unique()
    oas = pd.DataFrame(columns=cas.columns)

    for val in unique_vals:
        df_filtered = cas[cas.iloc[:, 0] == val]
        row = df_filtered.iloc[0]
        oas = oas.append(row, ignore_index=True)

    for col in oas.columns[1:]:
        mp_file = os.path.join(pastml_path, f'marginal_probabilities.character_{col}.model_F81.tab')
        mp = pd.read_table(mp_file)
        oas[col] = np.where(mp['0'] < mp['1'], 1, 0)

    oas.to_csv(os.path.join(pastml_path, f'ancestral_states_{addon}.tab'), sep='\t', index=False)
    return oas


# run
def main():
    start = time.time()
    for ds in ['SA', 'GA']:
        for lineage in ['2', '4']:
            size_options = ['ten', 'comp_hf']
            for size in size_options:
                iteration_start = time.time()
                print(ds, lineage, size)
                path = 'data/'
                tree_file = os.path.join(BASE_PATH, f'1_preprocessing_out/dated_trees/dated_{ds}_L{lineage}.nwk')
                csv_file = os.path.join(BASE_PATH, f'1_preprocessing_out/binary/binary_{ds}_L{lineage}_{size}.csv')
                work_dir = os.path.join(BASE_PATH, f'2_ASR_out/pypastml_{ds}_L{lineage}_{size}/')

                for i in range(bs_nr):
                    bt_file = os.path.join(BASE_PATH, f'1_preprocessing_out/dated_trees/dated_{ds}_L{lineage}_bt{i+1}.tre')
                    reconstruct(bt_file, csv_file, work_dir, ds, lineage, size)
                    unique_ancestral_states(work_dir, addon=f'bt{i+1}')

                reconstruct(tree_file, csv_file, work_dir, ds, lineage, size)
                unique_ancestral_states(work_dir, addon='main')

                iteration_end = time.time()
                iteration_time = iteration_end - iteration_start
                print(f"Time taken for {ds}, {lineage}, {size}: {iteration_time} seconds")

    end = time.time()
    seconds = end - start
    seconds = seconds % (24 * 3600)
    hours, rem = divmod(seconds, 3600)
    minutes, seconds = divmod(rem, 60)
    print(f"Total runtime: {hours}:{minutes}:{seconds}")


if __name__ == "__main__":
    main()





