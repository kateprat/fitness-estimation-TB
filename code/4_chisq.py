#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose: find upper and lower limit for CI
@author: katerinapratsinis
Email: pratsink@ethz.ch
"""

from ete3 import Tree          # to import a Newick tree (.tre file)
import numpy as np
import pandas as pd
from scipy.stats import chi2
import seaborn as sns
import matplotlib.pyplot as plt

import tensorflow as tf
assert tf.__version__>="2.0"

from phyloTF2.FitModelKeras import FitModelKeras
#import phyloTF2.RandomEffectsFitModel as RandomEffectsFitModel

from phyloTF2.TreeLikeLoss import TreeLikeLoss
from phyloTF2.TensorTree import TensorTree
import phyloTF2.TreeUtils as TreeUtils

import os
# adjust path as needed
#os.chdir('/Users/katerinapratsinis/Documents/CBB/2023FS/LS/code/')

import time
start = time.time()
printplot = True

sensitivity = 'a' # standard settings for birth & death rate =1
bs = 10 # number of bootstrap trees to iterate through
LR=0.001 # learning rate

for ds in ['SA', 'GA']:
    for lineage in ['2', '4']:
        for size in ['comp_hf']:
            path = 'data/'
            tree_file = path + '2_ASR_out/pypastml_' + ds + '_L' + lineage + '_' + size + '/named.tree_dated_' + ds + '_L' + lineage + '.nwk'        # tree file with the phylogeny we will fit our model to
            tab_file = path + '2_ASR_out/pypastml_' + ds + '_L' + lineage + '_' + size + '/ancestral_states' + '_main.tab'
            fit_file = path + "3_phyloTF_out/phyloTF_" + ds + "_L" + lineage + "_" + size + "_" +sensitivity + ".csv"
            freq_file = path + '1_preprocessing_out/' + ds + '_L' + lineage + '_comp-freq.csv'

            tree = Tree(tree_file, format=1)
            tree, tree_times = TreeUtils.add_tree_times(tree)
            tree = TreeUtils.index_branches(tree)

            df = pd.read_csv(tab_file, sep="\t", header=0, index_col=0)
            site_effects_df = pd.read_csv(fit_file, header=0, index_col=0)
            site_effects_df.loc['lower']=np.nan
            site_effects_df.loc['upper']=np.nan

            freqs = pd.read_csv(freq_file, header = 0, index_col=0)
            freqs.index = freqs.index.astype(str)
            intersect_cols = freqs.index.intersection(site_effects_df.columns)
            new_row = freqs.loc[intersect_cols, 'mut_freq']
            site_effects_df = site_effects_df.append(new_row)

            fitness_effects = site_effects_df.iloc[0].to_numpy()
            sites = df.shape[1]

            features_dic = {}

            for row_name, row_data in df.iterrows():
                row_array = np.array(row_data)
                row_array = row_array.astype(int)
                features_dic[row_name] = row_array

            beta = 1
            d = 1

            gamma = 0
            if ds == 'SA':
                s = 0.36
            elif ds == 'GA':
                s = 0.4
            rho = s
            dt = 1.0

            params = {'beta': beta, 'd': d, 'gamma': gamma, 's': s, 'rho': rho, 'dt': dt,
                      'time_intervals': 0, 'sites': sites, 'fitness_effects': fitness_effects}
            est_params = {'site_effects': False, 'beta': False, 'd': False, 'gamma': False,
                          's': False, 'rho': False}

            absolute_time = max(tree_times)  # final absolute time of sampling
            time_intervals = np.array([max(tree_times)])  # time intervals in absolute time

            final_time = max(tree_times)  # final time in tree time
            root_time = absolute_time - final_time  # time of root in absolute time
            time_intervals = time_intervals - root_time  # time intervals in tree time w/ root at t = 0
            params.update(time_intervals=time_intervals)

            def likelihood(params, est_params, tree, features_dic):
                model = FitModelKeras(params, est_params)
                tt = TensorTree(tree,features_dic,**params)
                tt.check_line_time_steps()
                fit_vals, bdm_params = model.call(tt)
                like_loss =TreeLikeLoss(tt, params)
                return float(like_loss.call(fit_vals,bdm_params,tt))

            max_ll = likelihood(params, est_params,tree,features_dic)

            def log_ratio(max_ll, alt_ll):
                return 2*(max_ll-alt_ll)

            test_stat = chi2.ppf(1 - 0.05, df=1)

            for i in range(sites):
                print(i)
                # get maximum likelihood of site
                print("fit: ", fitness_effects)

                lower_fitness = np.copy(fitness_effects)
                lower_fitness[i] = np.floor(lower_fitness[i]*1000)/1000
                params.update(fitness_effects=lower_fitness)
                lower_ll = likelihood(params, est_params, tree, features_dic)

                while log_ratio(max_ll,lower_ll) < test_stat and lower_fitness[i]>0:
                    print("lower log ratio for site ", i, " with fitness", lower_fitness[i],
                          ": ", log_ratio(max_ll, lower_ll))
                    lower_fitness[i]=lower_fitness[i] - LR
                    params.update(fitness_effects=lower_fitness)
                    lower_ll = likelihood(params, est_params,tree,features_dic)
                site_effects_df.xs('lower')[i] = lower_fitness[i] + LR

                upper_fitness = np.copy(fitness_effects)
                upper_fitness[i] = np.ceil(fitness_effects[i]*1000)/1000
                params.update(fitness_effects=upper_fitness)
                upper_ll = likelihood(params, est_params, tree, features_dic)

                #upper_axis_limit = 7

                while log_ratio(max_ll,upper_ll) < test_stat: #and upper_fitness[i] < upper_axis_limit:
                    print("upper log ratio for site ", i, " with fitness", upper_fitness[i],
                          ": ", log_ratio(max_ll, upper_ll))
                    upper_fitness[i]=upper_fitness[i] + LR
                    params.update(fitness_effects=upper_fitness)
                    upper_ll = likelihood(params, est_params,tree,features_dic)
                site_effects_df.xs('upper')[i] = upper_fitness[i] - LR

            site_effects_df.head()


            df=site_effects_df.T
            df.index.name = 'site'
            df.reset_index(inplace=True)

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
                # plot the bootstrap values
                if df.shape[1]>4:
                    for i, bar in enumerate(bars):
                        for bs in df.columns[2:-3]:  # assuming the bootstrap columns start from index 4
                            plt.scatter(df.loc[i, bs], bar.get_y() + 0.5*bar.get_height(), s=75, edgecolor='grey',
                            facecolor='none', zorder=2)

                bs_min = df[df.columns[2:-3]].min().min()
                bs_max = df[df.columns[2:-3]].max().max()

                x_min = min(df['lower'].min(), bs_min)
                x_max = max(df['upper'].max(), bs_max)
                x_padding = 0.05 * (df['upper'].max() - df['lower'].min())
                ax.set_xlim(x_min - x_padding, x_max + x_padding)

                ax.set_xlabel('Transmission Fitness Effect')
                ax.set_ylabel('Mutation Site')

                # vertical line at x = 1
                plt.axvline(x=1, color='grey', linestyle='--')

                # reverse the y-axis to get horizontal barplots
                ax.invert_yaxis()
                plt.tight_layout(pad=0.3)
                #plt.savefig('data/4_chisq_out/fit_' + ds + '_L' + lineage + '_' + size + '.png')
                #plt.show()

            compact = df.drop(df.iloc[:, 2:-3], axis=1)
            to_drop = []
            compact['MLE'] = round(df['main'], 3)
            compact['95CI'] = round(df['lower'], 3).astype(str) + ' - ' + round(df['upper'],3).astype(str)
            for i in range(compact.shape[0]):
                if compact.loc[i,'lower']<1 and compact.loc[i,'upper']>1:
                    to_drop.append(i)
            compact = compact.drop(to_drop)
            compact = compact.drop(['lower', 'upper', 'main'], axis=1)

            df.to_csv('data/4_chisq_out/fit_' + ds + '_L' + lineage + '_' + size + '.csv')
            #compact.to_csv('data/4_chisq_out/fit_' + ds + '_L' + lineage + '_' + size + "_" + str(LR) + '_compact.csv')



end = time.time()
seconds = end - start
seconds = seconds % (24 * 3600)
hour = seconds // 3600
seconds %= 3600
minutes = seconds // 60
seconds %= 60
print("Runtime: %d:%02d:%02d" % (hour, minutes, seconds))
