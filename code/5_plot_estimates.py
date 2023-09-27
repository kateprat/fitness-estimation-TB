#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose: plot sensitivity analyses and estimates
@author: katerinapratsinis
Email: pratsink@ethz.ch
"""

## import necessary packages
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import regex as re
from matplotlib.ticker import ScalarFormatter

#relabel
gene_dict = {
    'Rv0005': 'gyrB',
    'Rv0006': 'gyrA',
    'Rv0667': 'rpoB',
    'Rv0668': 'rpoC',
    'Rv0682': 'rpsL',
    'Rv1483': 'fabG1',
    'Rv1484': 'inhA',
    'Rv1908c': 'katG',
    'Rv2043c': 'pncA',
    'Rv2416c': 'eis',
    'Rv3794': 'embA',
    'Rv3795': 'embB',
    'Rv3919c': 'gidB'
}

def reformat_label(label):
    if lineage == '2':
        corr = "${{\\it{{pncA}}}}$ V139A"
    elif lineage == '4':
        corr = "${{\\it{{pncA}}}}$ K96T"

    appendage = ''
    if label.endswith('_R_L'):
        label = label[:-4]
        appendage = ', ' + corr + ' (L)'
        print(appendage)
    elif label.endswith('_L'):
        label = label[:-2]
        appendage = ' (L)'
    elif label.endswith('_R'):
        label = label[:-2]
        appendage = '; ' + corr
        print(appendage)

    if label.startswith('rrs_'):
        return f"${{\\it{{rrs}}}}$ {label.split('_')[1]}{appendage}"

    if label.startswith('rpo'):
        if label.endswith('104NA'):
            label = label[:-2]
            label = label + '0R'
        if label.endswith('125NA'):
            label = label[:-2]
            label = label + '2L'
        return f"${{\\it{{{label.split('_')[0]}}}}}$ {label.split('_')[2]} (C)"

    gene_code = label.split('_')[0]
    if gene_code in gene_dict:
        return f"${{\\it{{{gene_dict[gene_code]}}}}}$ {label.split('_')[1]}{appendage}"

    return label  # if no transformation found

for ds in ['SA', 'GA']:
    for lineage in ['2', '4']:
        if ds == 'SA':
            s = 0.36
        elif ds == 'GA':
            s = 0.4
        for size in ['comp_hf']:
            main_file = "data/4_chisq_out/fit_" + ds + "_L" + lineage + "_" + size + "_" + str(0.01) + ".csv"
            main_df = pd.read_csv(main_file, header=0, index_col=0, dtype={'site': str})
            main_df = main_df[['site', 'main', 'lower', 'upper']]

            # if values for P1040R will not print in SA L2, uncomment the following:
            #if ds == "SA" and lineage=="2":
            #    main_df.loc[main_df['site'] == "rpoC_Rv0668_P104NA", 'site'] = "rpoC_Rv0668_P1040R"

            path = 'data/3_phyloTF_out/'
            common = "phyloTF_" + ds + "_L" + lineage + "_" + size

            b_path = path + common + "_b.csv"
            b_df = pd.read_csv(b_path, header=0, index_col=0)
            b_df = b_df.rename(index={'main': 'lambda = 2'})
            new_columns = {col: col[:-2] if col.endswith("_c") else col for col in b_df.columns}
            b_df = b_df.rename(columns=new_columns)

            c_path = path + common + "_c.csv"
            c_df = pd.read_csv(c_path, header=0, index_col=0)
            c_df = c_df.rename(index={'main': 'lambda = 0.5'})
            new_columns = {col: col[:-2] if col.endswith("_c") else col for col in c_df.columns}
            c_df = c_df.rename(columns=new_columns)

            lambda_path = path + "est_lambda_" + common + "_a.csv"
            lambda_df = pd.read_csv(lambda_path, header=0, index_col=0)
            lambda_df = lambda_df.rename(index={'main': 'estimated lambda'})
            new_columns = {col: col[:-2] if col.endswith("_c") else col for col in lambda_df.columns}
            lambda_df = lambda_df.rename(columns=new_columns)

            d_path = path + common + "_d.csv"
            d_df = pd.read_csv(d_path, header=0, index_col=0)
            d_df = d_df.rename(index={'main': 'mu = 2'})
            new_columns = {col: col[:-2] if col.endswith("_c") else col for col in d_df.columns}
            d_df = d_df.rename(columns=new_columns)

            e_path = path + common + "_e.csv"
            e_df = pd.read_csv(e_path, header=0, index_col=0)
            e_df = e_df.rename(index={'main': 'mu = 0.5'})
            new_columns = {col: col[:-2] if col.endswith("_c") else col for col in e_df.columns}
            e_df = e_df.rename(columns=new_columns)

            mu_path = path + "est_mu_" + common + "_a.csv"
            mu_df = pd.read_csv(mu_path, header=0, index_col=0)
            mu_df = mu_df.rename(index={'main': 'estimated mu'})
            new_columns = {col: col[:-2] if col.endswith("_c") else col for col in mu_df.columns}
            mu_df = mu_df.rename(columns=new_columns)

            combined_df = pd.concat([b_df, c_df, lambda_df, d_df, e_df, mu_df])
            combined_df = combined_df.T
            combined_df = combined_df.reset_index().rename(columns={'index': 'site'})

            main_df = pd.merge(main_df, combined_df, on='site', how='left')

            for parameter in ["lambda", "mu"]:
                if len(main_df['site']) > 90:
                    fig, ax = plt.subplots(figsize=(18, len(main_df['site']) * 0.2))  
                elif len(main_df['site']) > 25:
                    fig, ax = plt.subplots(
                        figsize=(len(main_df['site']) * 0.25, len(main_df['site']) * 0.2)) 
                else:
                    fig, ax = plt.subplots(figsize=(10, 8))  # width: 10 inches, height: 20 inches
                ax.set_xlabel('Estimated Fitness Effect')
                ax.set_ylabel('Mutation')
                # dashed line at fitness value = 1
                ax.axvline(1, color='lightgrey', linestyle='--', zorder=0)

                colors = ['blue', 'orange', 'green']
                markers = ['D', 'p', '*']

                mu_path = path + "est_mu_" + common + "a.csv"

                est_path = path + "est_" + parameter + "_params_" + ds + "_L" + lineage + "_" + size + "_a.csv"
                est_df = pd.read_csv(est_path, header=0, index_col=0)
                est_value = str(round(est_df.loc[0, parameter],4))

                columns_to_plot = [parameter + " = 2", parameter+" = 0.5", "estimated " + parameter]

                # bars to plot
                for i, row in main_df.iterrows():
                    ax.barh(row["site"], width=row["upper"] - row["lower"], left=row["lower"], height=0.5,
                            color='darkgrey', zorder = 1)
                    # colored dots
                    for j, column in enumerate(columns_to_plot):
                        value = row[column]
                        color = colors[j]
                        marker = markers[j]
                        if not pd.isna(value):
                            ax.scatter(value, i, color=color, marker=marker, alpha=0.8, zorder=2)

                # MLE as black line
                ax.vlines(main_df["main"], ymin=np.array(range(len(main_df)))+0.25, ymax=np.array(range(len(main_df)))-0.25,
                          color='k', zorder=3)

                ax.set_yticks(range(len(main_df)))
                legend_labels = [reformat_label(label) for label in main_df["site"]]
                ax.set_yticklabels(legend_labels)

                if parameter == "lambda":
                    legend_labels = ["$\\" + parameter + "_0$ = 2, ", "$\\" + parameter + "_0$ = 0.5",
                                     "$\hat{\\" + parameter + "}_0$ = " + est_value]
                else:
                    legend_labels = ["$\\" + parameter + " = 2, \\sigma$ = " + str(s/2),
                                     "$\\" + parameter + " = 0.5, \\sigma$ = " + str(s*2),
                                     "$\hat{\\" + parameter + "} = " + est_value  + ", \\sigma$ = " + str(s)]

                legend_handles = [plt.Line2D([0], [0], marker=marker, color='w', markerfacecolor=color, markersize=10)
                                  for color, marker in zip(colors, markers)]

                ax.legend(legend_handles, legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=14)

                plt.subplots_adjust(right=2)  # adjust to squeeze in legend
                #plt.title(ds + ", L" + lineage + ", variations for $\\" + parameter + "_0$")
                ax.invert_yaxis()
                plt.tight_layout()
                plt.savefig('data/5_sensitivity-analysis_out/comparison_' + ds + "_L" + lineage + '_' + parameter +'.png')
                #plt.show()
