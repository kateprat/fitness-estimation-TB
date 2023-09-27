#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose: plot MLI with CI and bootstrap trees
@author: katerinapratsinis
Email: pratsink@ethz.ch
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import tensorflow as tf
assert tf.__version__>="2.0"

import os
#os.chdir('/Users/katerinapratsinis/Documents/CBB/2023FS/LS/code/')

import time
start = time.time()

printplot = True
sensitivity = 'a'
bs = 10 #number of bootstrap trees to iterate through (--> plotted as grey circles)


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
        if label.endswith('104NA_c'):
            label = label[:-4]
            label = label + '0R'
        if label.endswith('125NA_c'):
            label = label[:-4]
            label = label + '2L'
        return f"${{\\it{{{label.split('_')[0]}}}}}$ {label.split('_')[2]} (C)"

    gene_code = label.split('_')[0]
    if gene_code in gene_dict:
        return f"${{\\it{{{gene_dict[gene_code]}}}}}$ {label.split('_')[1]}{appendage}"

    return label  # if no transformation found

for ds in ['SA', 'GA']:
    for lineage in ['2', '4']:
        if ds == 'SA':
            size_options = ['comp_hf']
        elif lineage == '2':
            size_options = ['comp_hf']
        elif lineage == '4':
            size_options = ['comp_hf']
        for size in size_options:
            path = 'data/'
            csv_file = path + "4_chisq_out/fit_" + ds + "_L" + lineage + "_" + size + ".csv"
            site_effects_df = pd.read_csv(csv_file, header=0, index_col=0)

            # Extract the numeric part of the 'site' column and convert to integers
            # Extract 4-digit numbers
            site_effects_df['numeric_part'] = site_effects_df['site'].str.extract('(\d{4})')

            # Convert to integer, and for NaN (not found) replace with a large number
            site_effects_df['numeric_part'] = site_effects_df['numeric_part'].fillna(99999).astype(int)

            df = site_effects_df.sort_values('numeric_part')

            # Delete the auxiliary column
            df = df.drop('numeric_part', axis=1)
            df['site'] = df['site'].apply(reformat_label)

            if len(df['site']) > 90:
                fig, ax = plt.subplots(figsize=(18, len(df['site']) * 0.2))  # width: 10 inches, height: 20 inches
            elif len(df['site']) > 25:
                fig, ax = plt.subplots(
                    figsize=(len(df['site']) * 0.25, len(df['site']) * 0.2))  # width: 10 inches, height: 20 inches
            else:
                fig, ax = plt.subplots(figsize=(10,8))  # width: 10 inches, height: 20 inches

            # color palette
            colors = sns.color_palette("hls", len(df))
            # CI bars
            bars = ax.barh(df['site'], df['upper']-df['lower'], left=df['lower'], color=colors, edgecolor='black')
            # black lines
            for bar, value in zip(bars, df['main']):
                plt.plot([value, value], [bar.get_y(), bar.get_y() + bar.get_height()], color='black')
            if printplot:
                # Plot the bootstrap values
                if df.shape[1] > 4:
                    for i, bar in enumerate(bars):
                        for bs in df.columns[2:-3]:  # assuming the bootstrap columns start from index 4
                            plt.scatter(df.iloc[i][bs], bar.get_y() + 0.5 * bar.get_height(), s=75, edgecolor='grey',
                                        facecolor='none', zorder=2)

                bs_min = df[df.columns[2:-3]].min().min()
                bs_max = df[df.columns[2:-3]].max().max()

                x_min = min(df['lower'].min(), bs_min)
                x_max = min(4, max(df['upper'].max(), bs_max))
                x_padding = 0.05 * (df['upper'].max() - df['lower'].min())
                ax.set_xlim(x_min - x_padding, x_max + x_padding)

                ax.set_xlabel('Transmission Fitness Effect', fontsize=14)
                ax.set_ylabel('Mutation', fontsize=14)

                ax.tick_params(axis='x', labelsize=14)
                ax.tick_params(axis='y', labelsize=14)

                # vertical line at x = 1
                plt.axvline(x=1, color='grey', linestyle='--')

                # reverse the y-axis to get horizontal barplots
                ax.invert_yaxis()
                plt.tight_layout(pad=0.3)
                plt.savefig('data/4_chisq_out/sorted_' + ds + '_L' + lineage + '.png')
                if ds=='GA' and lineage == '4':
                    plt.savefig('data/4_chisq_out/sorted_' + ds + '_L' + lineage + '.png', dpi=220)
                #plt.show()

            compact = df.drop(df.iloc[:,2:-3], axis = 1)
            # remove insignificant sites
            to_drop = []
            compact['MLE'] = round(df['main'],3)
            compact['mut_freq'] = round(df['mut_freq'],3)
            compact['95CI'] = round(df['lower'],3).astype(str) + ' - ' + round(df['upper'],3).astype(str)
            for i in range(compact.shape[0]):
                if compact.loc[i,'lower']<1 and compact.loc[i,'upper']>1:
                    to_drop.append(i)
            compact = compact.drop(to_drop)
            compact = compact.drop(['lower','upper','main'], axis = 1)

            compact.to_csv('data/4_chisq_out/sorted_' + ds + '_L' + lineage + '_compact.csv')
            df.to_csv('data/4_chisq_out/sorted_' + ds + '_L' + lineage + '.csv')



end = time.time()
seconds = end - start
seconds = seconds % (24 * 3600)
hour = seconds // 3600
seconds %= 3600
minutes = seconds // 60
seconds %= 60
print("Runtime: %d:%02d:%02d" % (hour, minutes, seconds))
