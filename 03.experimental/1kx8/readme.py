import os

import rstoolbox

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def plot_step1(df):
    fig = plt.figure(figsize=(10, 5))
    grid = (1, 2)
    ax00 = plt.subplot2grid(grid, (0, 0), fig=fig)
    sns.boxplot(y = df['score'], data=df, ax=ax00, showfliers=False)
    ax00.axhline(df[df['description'] == '188_1kx8_0033'].iloc[0]['score'], color='red', linewidth=3)
    rstoolbox.utils.add_top_title(ax00, 'RosettaScore')
    ax00.set_ylabel('REU')
    ax01 = plt.subplot2grid(grid, (0, 1), fig=fig)
    sns.boxplot(y = df['hbond_sr_bb'], data=df, ax=ax01, showfliers=False)
    ax01.axhline(df[df['description'] == '188_1kx8_0033'].iloc[0]['hbond_sr_bb'], color='red', linewidth=3)
    rstoolbox.utils.add_top_title(ax01, 'Short-Range Hydrogen Bonds')
    ax01.set_ylabel('REU')
    plt.tight_layout()
    plt.show()


def plot_step2(df):
    fig = plt.figure(figsize=(10, 5))
    grid = (1, 2)
    ax00 = plt.subplot2grid(grid, (0, 0), fig=fig)
    sns.boxplot(y = df['score'], data=df, ax=ax00, showfliers=False)
    rstoolbox.utils.add_top_title(ax00, 'RosettaScore')
    ax00.set_ylabel('REU')
    ax01 = plt.subplot2grid(grid, (0, 1), fig=fig)
    sns.boxplot(y = df['hbond_sr_bb'], data=df, ax=ax01, showfliers=False)
    rstoolbox.utils.add_top_title(ax01, 'Short-Range Hydrogen Bonds')
    ax01.set_ylabel('REU')
    plt.tight_layout()
    plt.show()


def plot_step3(df):
    fig = plt.figure(figsize=(15, 5))
    grid = (1, 3)
    ax00 = plt.subplot2grid(grid, (0, 0), fig=fig)
    sns.boxplot(y = df['cav_vol'], data=df, ax=ax00, showfliers=False)
    ax00.axhline(df[df['description'] == '138_188_1kx8_0033_0018_0001'].iloc[0]['cav_vol'], color='red', linewidth=2)
    ax00.axhline(df[df['description'] == '158_188_1kx8_0033_0028_0001'].iloc[0]['cav_vol'], color='red', linewidth=2)
    ax00.axhline(df[df['description'] == '158_188_1kx8_0033_0044_0001'].iloc[0]['cav_vol'], color='red', linewidth=2)
    ax00.axhline(df[df['description'] == '215_188_1kx8_0033_0024_0001'].iloc[0]['cav_vol'], color='red', linewidth=2)
    ax00.axhline(df[df['description'] == '72_188_1kx8_0033_0018_0001'].iloc[0]['cav_vol'], color='red', linewidth=2)
    ax00.axhline(df[df['description'] == '85_188_1kx8_0033_0016_0001'].iloc[0]['cav_vol'], color='red', linewidth=2)
    rstoolbox.utils.add_top_title(ax00, 'Cavity Volume')
    ax00.set_ylabel('vol')
    ax01 = plt.subplot2grid(grid, (0, 1), fig=fig)
    sns.boxplot(y = df['BUNS'], data=df, ax=ax01, showfliers=False)
    ax01.axhline(df[df['description'] == '138_188_1kx8_0033_0018_0001'].iloc[0]['BUNS'], color='red', linewidth=2)
    ax01.axhline(df[df['description'] == '158_188_1kx8_0033_0028_0001'].iloc[0]['BUNS'], color='red', linewidth=2)
    ax01.axhline(df[df['description'] == '158_188_1kx8_0033_0044_0001'].iloc[0]['BUNS'], color='red', linewidth=2)
    ax01.axhline(df[df['description'] == '215_188_1kx8_0033_0024_0001'].iloc[0]['BUNS'], color='red', linewidth=2)
    ax01.axhline(df[df['description'] == '72_188_1kx8_0033_0018_0001'].iloc[0]['BUNS'], color='red', linewidth=2)
    ax01.axhline(df[df['description'] == '85_188_1kx8_0033_0016_0001'].iloc[0]['BUNS'], color='red', linewidth=2)
    rstoolbox.utils.add_top_title(ax01, 'Buried-UNsatisfied')
    ax01.set_ylabel('BUNS')
    ax02 = plt.subplot2grid(grid, (0, 2), fig=fig)
    sns.boxplot(y = df['rmsd'], data=df, ax=ax02, showfliers=False)
    ax02.axhline(df[df['description'] == '138_188_1kx8_0033_0018_0001'].iloc[0]['rmsd'], color='red', linewidth=2)
    ax02.axhline(df[df['description'] == '158_188_1kx8_0033_0028_0001'].iloc[0]['rmsd'], color='red', linewidth=2)
    ax02.axhline(df[df['description'] == '158_188_1kx8_0033_0044_0001'].iloc[0]['rmsd'], color='red', linewidth=2)
    ax02.axhline(df[df['description'] == '215_188_1kx8_0033_0024_0001'].iloc[0]['rmsd'], color='red', linewidth=2)
    ax02.axhline(df[df['description'] == '72_188_1kx8_0033_0018_0001'].iloc[0]['rmsd'], color='red', linewidth=2)
    ax02.axhline(df[df['description'] == '85_188_1kx8_0033_0016_0001'].iloc[0]['rmsd'], color='red', linewidth=2)
    rstoolbox.utils.add_top_title(ax02, 'RMSD drift')
    ax02.set_ylabel('RMSD')
    plt.tight_layout()
    plt.show()


def load_experimental():
    designs = {}
    for r in range(1, 3):
        r = 'selection/round{}'.format(r)
        for c in range(1, 8):
            for f in ['mota_1kx8_0{}', 'mota_1kx8_03_d{}', 'mota_1kx8_d{}']:
                f = f.format(c)
                cddir = os.path.join(r, f, 'experimental', 'CD')
                if os.path.isdir(cddir):
                    cd = rstoolbox.io.read_CD(cddir)
                    designs.setdefault(f, {}).setdefault('CD', cd)
                cddir = os.path.join(r, f, 'experimental', 'CDTCEP')
                if os.path.isdir(cddir):
                    cd = rstoolbox.io.read_CD(cddir)
                    designs.setdefault(f, {}).setdefault('CDTCEP', cd)
                sprfile = os.path.join(r, f, 'experimental', 'SPR', f + '.txt')
                if os.path.isfile(sprfile):
                    spr = rstoolbox.io.read_SPR(sprfile)
                    designs.setdefault(f, {}).setdefault('SPR', spr)
                malsfile = os.path.join(r, f, 'experimental', 'MALS', f + '.csv')
                malsfile2 = os.path.join(r, f, 'experimental', 'MALS', f + '_mm.csv')
                if os.path.isfile(malsfile):
                    malsfile2 = None if not os.path.isfile(malsfile2) else malsfile2
                    mals = rstoolbox.io.read_MALS(malsfile, malsfile2)
                    designs.setdefault(f, {}).setdefault('MALS', mals)
    return designs


def firstgen_plot(data, melt_wave=220):
    fig = plt.figure(figsize=(15, 22))
    grid = (6, 3)
    for i, ii in enumerate([1,2,3,4,6,7]):
        decoy = 'mota_1kx8_0{}'.format(ii)
        dfd = data[decoy]
        # CD
        ax = plt.subplot2grid(grid, (i, 0))
        rstoolbox.utils.add_left_title(ax, decoy)
        if i == 1:
            rstoolbox.utils.add_top_title(ax, 'CD')
        if 'CD' in dfd:
            df = dfd['CD']
            rstoolbox.plot.plot_CD(df, ax, sample=7)
        # Melting Point
        ax = plt.subplot2grid(grid, (i, 1))
        if i == 1:
            rstoolbox.utils.add_top_title(ax, 'Thermal Melt ({} nm)'.format(melt_wave))
        if 'CD' in dfd:
            df = dfd['CD']
            df = df[df['Wavelength'] == melt_wave]
            rstoolbox.plot.plot_thermal_melt(df, ax)
        # SPR
        ax = plt.subplot2grid(grid, (i, 2))
        if i == 1:
            rstoolbox.utils.add_top_title(ax, 'SPR')
        if 'SPR' in dfd:
            df = dfd['SPR']
            rstoolbox.plot.plot_SPR(df, ax, datacolor='black', fitcolor='red')
    plt.tight_layout()
    #plt.savefig('1kx8_gen1.svg')
    plt.show()


def secondgen_plot(data, melt_wave=220):
    fig = plt.figure(figsize=(15, 22))
    grid = (6, 3)
    for i, ii in enumerate([2,3,5,6,7]):
        decoy = 'mota_1kx8_d{}'.format(ii)
        dfd = data[decoy]
        # CD
        ax = plt.subplot2grid(grid, (i, 0))
        rstoolbox.utils.add_left_title(ax, decoy)
        if i == 1:
            rstoolbox.utils.add_top_title(ax, 'CD')
        if 'CD' in dfd:
            df = dfd['CD']
            rstoolbox.plot.plot_CD(df, ax, sample=7)
        # Melting Point
        ax = plt.subplot2grid(grid, (i, 1))
        if i == 1:
            rstoolbox.utils.add_top_title(ax, 'Thermal Melt ({} nm)'.format(melt_wave))
        if 'CD' in dfd:
            df = dfd['CD']
            df = df[df['Wavelength'] == melt_wave]
            rstoolbox.plot.plot_thermal_melt(df, ax)
        # SPR
        ax = plt.subplot2grid(grid, (i, 2))
        if i == 1:
            rstoolbox.utils.add_top_title(ax, 'SPR')
        if 'SPR' in dfd:
            df = dfd['SPR']
            rstoolbox.plot.plot_SPR(df, ax, datacolor='black', fitcolor='red')

    decoy = 'mota_1kx8_03_d1'
    i = 5
    dfd = data[decoy]
    # CD
    ax = plt.subplot2grid(grid, (i, 0))
    rstoolbox.utils.add_left_title(ax, decoy)
    if 'CD' in dfd:
        df = dfd['CD']
        rstoolbox.plot.plot_CD(df, ax, sample=7)
    # Melting Point
    ax = plt.subplot2grid(grid, (i, 1))
    if 'CD' in dfd:
        df = dfd['CD']
        df = df[df['Wavelength'] == melt_wave]
        rstoolbox.plot.plot_thermal_melt(df, ax)
    # SPR
    ax = plt.subplot2grid(grid, (i, 2))
    if 'SPR' in dfd:
        df = dfd['SPR']
        rstoolbox.plot.plot_SPR(df, ax, datacolor='black', fitcolor='red')

    plt.tight_layout()
    #plt.savefig('1kx8_gen2.svg')
    plt.show()

