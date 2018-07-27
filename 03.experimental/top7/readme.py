import os

import rstoolbox

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


def plot_step1(df):
    fig = plt.figure(figsize=(10, 5))
    grid = (1, 2)
    ax00 = plt.subplot2grid(grid, (0, 0), fig=fig)
    sns.boxplot(y = df['rmsd_drift'], data=df, ax=ax00, showfliers=False)
    ax00.axhline(df[df['description'] == '87_ffl_twostrand_1qys_Lys_0001'].iloc[0]['rmsd_drift'], color='red', linewidth=3)
    rstoolbox.utils.add_top_title(ax00, 'RMSD Drift from TOP7')
    ax00.set_ylabel('RMSD')
    ax01 = plt.subplot2grid(grid, (0, 1), fig=fig)
    sns.boxplot(y = df['hbond_lr_bb'], data=df, ax=ax01, showfliers=False)
    ax01.axhline(df[df['description'] == '87_ffl_twostrand_1qys_Lys_0001'].iloc[0]['hbond_lr_bb'], color='red', linewidth=3)
    rstoolbox.utils.add_top_title(ax01, 'Long-Range Hydrogen Bonds')
    ax01.set_ylabel('REU')
    plt.tight_layout()
    plt.show()


def plot_full(df):
    fig = plt.figure(figsize=(15, 5))
    grid = (1, 3)
    ax00 = plt.subplot2grid(grid, (0, 0), fig=fig)
    sns.boxplot(y = df['score'], data=df, ax=ax00, showfliers=False)
    ax00.axhline(df[df['description'] == 'sd_top7_12top7_graft_0001'].iloc[0]['score'], color='red', linewidth=3)
    rstoolbox.utils.add_top_title(ax00, 'RosettaScore')
    ax00.set_ylabel('REU')
    ax01 = plt.subplot2grid(grid, (0, 1), fig=fig)
    sns.boxplot(y = df['hbond_lr_bb'], data=df, ax=ax01, showfliers=False)
    ax01.axhline(df[df['description'] == 'sd_top7_12top7_graft_0001'].iloc[0]['hbond_lr_bb'], color='red', linewidth=3)
    rstoolbox.utils.add_top_title(ax01, 'Long-Range Hydrogen Bonds')
    ax01.set_ylabel('REU')
    ax02 = plt.subplot2grid(grid, (0, 2), fig=fig)
    sns.boxplot(y = df['pack'], data=df, ax=ax02, showfliers=False)
    ax02.axhline(df[df['description'] == 'sd_top7_12top7_graft_0001'].iloc[0]['pack'], color='red', linewidth=3)
    rstoolbox.utils.add_top_title(ax02, 'Packing')
    ax02.set_ylabel('Percentage')
    plt.tight_layout()
    plt.show()


def plot_part(df):
    fig = plt.figure(figsize=(15, 5))
    grid = (1, 3)
    ax00 = plt.subplot2grid(grid, (0, 0), fig=fig)
    sns.boxplot(y = df['score'], data=df, ax=ax00, showfliers=False)
    ax00.axhline(df[df['description'] == 'sd_top7_nterm2_45top7_graft_0002'].iloc[0]['score'], color='red', linewidth=3)
    rstoolbox.utils.add_top_title(ax00, 'RosettaScore')
    ax00.set_ylabel('REU')
    ax01 = plt.subplot2grid(grid, (0, 1), fig=fig)
    sns.boxplot(y = df['hbond_lr_bb'], data=df, ax=ax01, showfliers=False)
    ax01.axhline(df[df['description'] == 'sd_top7_nterm2_45top7_graft_0002'].iloc[0]['hbond_lr_bb'], color='red', linewidth=3)
    rstoolbox.utils.add_top_title(ax01, 'Long-Range Hydrogen Bonds')
    ax01.set_ylabel('REU')
    ax02 = plt.subplot2grid(grid, (0, 2), fig=fig)
    sns.boxplot(y = df['pack'], data=df, ax=ax02, showfliers=False)
    ax02.axhline(df[df['description'] == 'sd_top7_nterm2_45top7_graft_0002'].iloc[0]['pack'], color='red', linewidth=3)
    rstoolbox.utils.add_top_title(ax02, 'Packing')
    ax02.set_ylabel('Percentage')
    plt.tight_layout()
    plt.show()


def load_experimental():
    designs = {'top7_full': {}, 'top7_partial': {}}
    designs['top7_full'].setdefault('MALS', pd.read_csv(os.path.join('selection', 'top7_full', 'experimental', 'MALS', 'top7_full_peak.csv')))
    designs['top7_partial'].setdefault('MALS', rstoolbox.io.read_MALS(os.path.join('selection', 'top7_partial', 'experimental', 'MALS', 'top7_partial_peak.txt')))

    designs['top7_full'].setdefault('CD', pd.read_csv(os.path.join('selection', 'top7_full', 'experimental', 'CD', 'top7_cd.csv.gz') ))
    designs['top7_partial'].setdefault('CD', rstoolbox.io.read_CD(os.path.join('selection', 'top7_partial', 'experimental', 'CD')))
    designs['top7_partial']['CD'].drop(columns=['voltage', 'title', 'bin', 'minw'], inplace=True)
    designs['top7_partial']['CD'] = pd.DataFrame(designs['top7_partial']['CD'])

    designs['top7_full'].setdefault('SPR', rstoolbox.io.read_SPR(os.path.join('selection', 'top7_full', 'experimental', 'SPR', 'top7_full.csv')))
    designs['top7_partial'].setdefault('SPR', rstoolbox.io.read_SPR(os.path.join('selection', 'top7_partial', 'experimental', 'SPR', 'top7_partial.txt')))
    return designs


def data_plot(data, melt_wave=220):
    fig = plt.figure(figsize=(20, 25))
    grid = (7, 4)
    for i, decoy in enumerate(['top7_full', 'top7_partial']):
        dfd = data[decoy]
        # MALS
        ax = plt.subplot2grid(grid, (i, 0), fig=fig)
        rstoolbox.utils.add_left_title(ax, decoy)
        if i == 0:
            rstoolbox.utils.add_top_title(ax, 'MALS')
        df = dfd['MALS']
        rstoolbox.plot.plot_MALS(df, ax, lscolor=False, mwcolor=False)
        # CD
        ax = plt.subplot2grid(grid, (i, 1), fig=fig)
        if i == 0:
            rstoolbox.utils.add_top_title(ax, 'CD')
        df = dfd['CD']
        rstoolbox.plot.plot_CD(df, ax)
        # SPR
        ax = plt.subplot2grid(grid, (i, 2), fig=fig)
        if i == 0:
            rstoolbox.utils.add_top_title(ax, 'SPR')
        df = dfd['SPR']
        rstoolbox.plot.plot_SPR(df, ax, datacolor='black', fitcolor='red')

    plt.tight_layout()
    #plt.savefig('top7_compare2.svg')
    plt.show()

