#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose: phyloTF
@author: katerinapratsinis
Email: pratsink@ethz.ch
"""

# created with template from https://davidrasm.github.io/phyloTF2/

import time
import numpy as np
import pandas as pd
import seaborn as sns
import tensorflow as tf
from ete3 import Tree
from tensorflow import keras
from phyloTF2.FitModelKeras import FitModelKeras
from phyloTF2.RandomEffectsFitModel import RandomEffectsFitModel
from phyloTF2.TreeLikeLoss import TreeLikeLoss
from phyloTF2.TensorTree import TensorTree
import phyloTF2.TreeUtils as TreeUtils
from phyloTF2.L1Regularizer import L1Regularizer
import matplotlib.pyplot as plt

assert tf.__version__ >= "2.0"

start = time.time()

bs_nr = 0  # =10 when analyzing bootstrapped trees
plot_convergence_plot = True

lr = 0.001
addon = "est_sigma_"  # in case of alternative preceding output label

# relabel locus --> gene
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

#################### SETTING UP BIRTH-DEATH-SAMPLING MODEL ####################
for ds in ['SA', 'GA']:
    for lineage in ['2', '4']:
        size_options = ['comp_hf'] # options: "ten" (subset), "comp_hf" (entire data set including compensatory mutations)
        for size in size_options:
            if size == 'ten':
                regularizer = "0"
            else:
                regularizer = "1"
            for sensitivity in ['a']:  # a:stat.quo; b: b=2; c: b=0.5; d: d=2; e: d=0.5
                beta = 1  # =1/year, transmission / birth rate
                d = 1  # =1/year death / removal rate

                if "b" in sensitivity:
                    beta = 2
                if "c" in sensitivity:
                    beta = 0.5
                if "d" in sensitivity:
                    d = 2
                if "e" in sensitivity:
                    d = 0.5

                factor = 0
                n_epochs = 10000
                if ds == 'SA':
                    s = 0.36 / d #sampling rate at removal
                    if lineage == '2':
                        factor = 5.9
                    elif lineage == '4':
                        n_epochs = 16000
                elif ds == 'GA':
                    s = 0.4/d
                    if lineage == '2':
                        factor = 3.1


                gamma = 0.0  # migration / transition rate (none here)

                rho = s  # same assumption as s, sampling fraction at present
                dt = 1.0  # status quo. time step interval for update pE0s along branch

                convergence_plot = False
                path = 'data/2_ASR_out/'
                tree_base = path + 'pypastml_' + ds + '_L' + lineage + '_' + size + '/named.tree_dated_' + ds + '_L' + lineage  # tree file with the phylogeny we will fit our model to
                tab_base = path + 'pypastml_' + ds + '_L' + lineage + '_' + size + '/ancestral_states'

                file_dict = {'main': {'tree_file': tree_base + '.nwk', 'tab_file': tab_base + '_main.tab'}}

                for i in range(bs_nr):
                    key = 'bs' + str(i + 1)
                    file_dict[key] = {'tree_file': tree_base + '_bt' + str(i + 1) + '.nwk',
                                      'tab_file': tab_base + '_bt' + str(i + 1) + '.tab'}

                features_dic = {}
                sites = {}
                feature_names = {}

                idx = 0
                for key in file_dict.keys():
                    tree_file = file_dict[key]['tree_file']
                    tab_file = file_dict[key]['tab_file']

                    tree = Tree(tree_file, format=1)
                    tree, tree_times = TreeUtils.add_tree_times(tree)
                    tree = TreeUtils.index_branches(tree)

                    df = pd.read_csv(tab_file, sep="\t", header=0, index_col=0)
                    features_dic[key] = {}
                    for row_name, row_data in df.iterrows():
                        row_array = np.array(row_data)
                        row_array = row_array.astype(int)
                        features_dic[key][row_name] = row_array

                    if idx == 0:
                        sites = df.shape[1]
                        site_effects_df = pd.DataFrame(columns=df.columns.values)
                        params = {'beta': beta, 'd': d, 'gamma': gamma, 's': s, 'rho': s, 'dt': dt,
                                  'time_intervals': 0, 'sites': sites}
                        # adjust in case of additional parameter estimation
                        est_params = {'site_effects': True, 'beta': False, 'd': False, 'gamma': False,
                                      's': False, 'rho': False}
                    idx = idx + 1

                    ########################## SETTING UP TIME INTERVALS ##########################
                    absolute_time = max(tree_times)  # final absolute time of sampling
                    time_intervals = np.array([absolute_time])

                    final_time = max(tree_times)  # final time in tree time
                    root_time = absolute_time - final_time  # time of root in absolute time
                    time_intervals = time_intervals - root_time  # time intervals in tree time w/ root at t = 0
                    params.update(time_intervals=time_intervals)

                    ########################## SETTING UP THE TensorTree ##########################
                    tt = TensorTree(tree, features_dic[key], **params)

                    tt.check_line_time_steps()

                    ######################### BUILDING THE FITNESS MODEL ##########################
                    branches = np.max(tt.birth_branch_indexes) + 1  # num of branches: plus one b/c indexed from 0
                    model = FitModelKeras(params, est_params)
                    fit_vals, bdm_params = model.call(tt)
                    like_loss = TreeLikeLoss(tt, params)

                    ############################ ADDING REGULARIZATION ############################
                    if regularizer == "1":
                        reg = L1Regularizer(factor, offset=1.0)  # Lasso regularization penalizes sum of absolute values of weights

                    ############################## FITTING THE MODEL ##############################
                    # epoch and learning rate of 0.001: adjusted to problem at hand to ensure convergence
                    optimizer = keras.optimizers.Adam(learning_rate=lr)
                    if key == 'main' and sensitivity == 'a':
                        convergence_plot = True
                        legend_labels = [reformat_label(label) for label in df.columns.values]
                        parameters_history = []
                    for epoch in range(1, n_epochs + 1):
                        with tf.GradientTape() as tape:
                            fit_vals, bdm_params = model(tt)
                            if regularizer != "0":
                                reg_penalty = reg.call(model.site_effects)
                            else:
                                reg_penalty = 0
                            penalty = 0
                            loss = -(like_loss.call(fit_vals, bdm_params, tt) + penalty) + reg_penalty
                        gradients = tape.gradient(loss, model.trainable_variables)
                        optimizer.apply_gradients(zip(gradients, model.trainable_variables))

                        if key == 'main' and sensitivity == 'a':
                            parameters = model.site_effects.numpy()
                            parameters_history.append(parameters)

                        if epoch % 10 == 0:
                            print("Epoch", epoch, "loss =", str(-loss.numpy()),
                                  "params =",
                                  str(model.trainable_variables[0].numpy()))

                    # dictionary with keys as the column names and values as the corresponding lists
                    data = {'lambda': model.beta_series.numpy(), 'mu': model.d_series.numpy(),
                            'sigma': model.s_series.numpy(), 'rho': model.rho.numpy()}
                    est_df = pd.DataFrame(data)

                    ############################ SAVING THE ESTIMATES #############################
                    site_effects_ests = model.trainable_variables[0].numpy().reshape((1, sites))
                    print(site_effects_ests)
                    site_effects_df.loc[key] = site_effects_ests[0]

                path_out = 'data/3_phyloTF_out/'
                file_name = path_out + addon + "phyloTF_" + ds + "_L" + lineage + "_" + size + "_" + sensitivity + ".csv"
                site_effects_df.to_csv(file_name)
                est_path = path_out + addon + "params_" + ds + "_L" + lineage + "_" + size + "_" + sensitivity + ".csv"
                # est_df.to_csv(est_path) #uncomment to save estimated parameters for transmission, recovery or sampling rate

                if convergence_plot and plot_convergence_plot:
                    parameters_history = np.array(parameters_history)
                    fig, ax = plt.subplots(figsize=(10, 8))
                    for i in range(parameters_history.shape[1]):
                        plt.plot(parameters_history[:, i], label=legend_labels[i])
                    plt.xlabel('Epoch')
                    plt.ylabel('Fitness Estimate')
                    #plt.title(ds + ', L' + lineage + ', LR = ' + str(lr) + ', L1 factor = ' + str(factor) + ', sensitivity: ' + sensitivity)
                    if len(legend_labels)>25:
                        spalte = 2
                    else:
                        spalte = 1
                    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', prop={'size': 8}, ncol=spalte)
                    plt.subplots_adjust(right=0.7)
                    # plt.show()
                    fig_name = path_out + addon + "phyloTF_" + ds + "_L" + lineage + "_" + size + "_convergence.png"
                    fig.savefig(fig_name, bbox_inches='tight')

end = time.time()
seconds = end - start
seconds = seconds % (24 * 3600)
hour = seconds // 3600
seconds %= 3600
minutes = seconds // 60
seconds %= 60
print("Runtime: %d:%02d:%02d" % (hour, minutes, seconds))

