import gzip

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

import rstoolbox

def parse_pfam( f ):
    """
    Read and process the PFAM file
    """
    data = {"pdb":[], "chain":[], "pfamrange":[], "pfam":[]}
    with gzip.open(f) as fd:
        for line in fd:
            l = [x.strip() for x in line.strip().split(";")]
            data["pdb"].append(l[0].lower())
            data["chain"].append(l[1])
            data["pfamrange"].append(l[-2])
            data["pfam"].append(l[3])
    return pd.DataFrame(data)

def pfam_length( pfam_range ):
    """
    Measure the length of the PFAM domain
    """
    r = [int(x) for x in pfam_range.split("-")]
    return r[1] - r[0] + 1

def fit_in_range( motif_range, pfam_range ):
    """
    As a single protein can have multiple PFAM domains, this function evaluates if the motif match falls
    inside the assigned PFAM domain.
    """
    if pfam_range == np.nan or isinstance(pfam_range, float): return False
    pfam_range = [int(x) for x in pfam_range.split("-")]
    motif_range = motif_range.replace(")(", "),(").replace("(", "").replace(")", "").strip("[]")
    cmatch = 0
    for r in motif_range.split(","):
        r = [int(x) for x in r.split("-")]
        if r[0] >= pfam_range[0] and r[1] <= pfam_range[1]:
           cmatch += 1
        else: return False
    return cmatch == len(motif_range.split(","))

def data_plot( df, ax, maxim, rmsd_lim, linestyle, color ):
    """
    Calculate true cumulative curves
    (seaborn does an aproximation that does not work well with very small recovery)
    """
    allvalues = df.sort("rmsd").groupby(["pdb","chain"]).head(1)[["rmsd"]].values
    raw, y, x = rstoolbox.analysis.cumulative(allvalues, max_count=maxim, upper_limit=rmsd_lim)
    ax.plot(x, y, color=sns.color_palette()[color], lw=4, linestyle=linestyle)
    return raw, x, y

def first_marker_plot( raw, x, ax, top, shape, color, linestyle ):
    """
    Top x marker
    """
    idx = (np.abs(np.array(raw)-top)).argmin()
    ax.axvline(x=x[idx], ymin=0, ymax=0.2, c=sns.color_palette()[color], linewidth=4, linestyle=linestyle, zorder=10)
    ax.plot([x[idx]], [0.2], shape, c=sns.color_palette()[color], markersize=12)


def plot( df, selfname, maxim, title, ax, known=None, knames=None, kcolor=None, annotate=False, top=10, rmsd_lim=5, max_length=180 ):
    """
    Full plot
    """
    # 1. Filter self-hits and error assignments
    df = df[(df["pdb"] != "eeee") & (df["pdb"] != selfname)]

    # 2. Make cumulative plot for best-hit/pdb-chain for all dataset
    raw, x, y = data_plot(df, ax, maxim, rmsd_lim, "solid", 0)
    if annotate:
        first_marker_plot( raw, x, ax, top, 'o', 0, "solid")

    # 3. Make cumulative plot for best-hit/pdb-chain for proteins of size "max_length" or smaller
    sdf = df[(df["length"] <= max_length) & (df["inrange"] == True)]
    raw, x, y = data_plot(sdf, ax, maxim, rmsd_lim, "dashed", 0)
    if annotate:
        first_marker_plot( raw, x, ax, top, 's', 0, "dashed")

    # 4. Add clashes data
    if "ddg" in df:
        # 4.1. Make cumulative plot for best-hit/pdb-chain for all non-clashing dataset
        raw, x, y = data_plot(df[(df["ddg"] <= 0) ], ax, maxim, rmsd_lim, "solid", 1)
        ax.fill_between(x, 0, y, color=sns.color_palette()[1], alpha=0.3)
        if annotate:
            first_marker_plot( raw, x, ax, top, 'o', 1, "solid")

        # 4.2. Make cumulative plot for best-hit/pdb-chain for all non-clashing dataset for proteins of size "max_length" or smaller
        sdf = df[(df["ddg"] <= 0) & (df["length"] <= max_length) & (df["inrange"] == True)]
        raw, x, y = data_plot(sdf, ax, maxim, rmsd_lim, "dashed", 1)
        if annotate:
            first_marker_plot( raw, x, ax, top, 's', 1, "dashed")

    # 5. Label where known designs match
    if known is not None:
        for i, k in enumerate(known):
            ax.axvline(x=k, ymin=0, ymax=1, c=kcolor[i], linewidth=2, linestyle="dashed", zorder=10)

    # 6. Label where own designs mathc
    if knames is not None:
        for i, k in enumerate(knames):
           ax.text(known[i] - 0.3, 0.8 - (0.06 * i), k, color=kcolor[i], horizontalalignment='center', size='x-small', weight='semibold')

    # 7. Show the max expected number to retrieve defined as "maxim"; value = 1 in y scale
    ax.text(4.02, 0.15, "n = {}".format(maxim), color="black", horizontalalignment='left', size='x-small', weight='semibold')

    # 8. Format axis
    ax.set_ylim(0, 1)
    ax.set_xlim(0, 5)
    ax.set_xlabel("RMSD")
    ax.set_ylabel("")
    ax.annotate(title, xy=(0.5, 1), xytext=(0, 5), xycoords='axes fraction',
                textcoords='offset points', ha='center', va='baseline')
