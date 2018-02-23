import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.patheffects
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch
import matplotlib.patches as mpatches
from matplotlib import transforms
from matplotlib.font_manager import FontProperties
from matplotlib.text import TextPath

from IPython.display import display_javascript, display_html, display

import rstoolbox

def show_info( info ):
    display_html(
    """<div style="width:100%;">
    <table style="width:60%;margin-left:20%">
        <tr style="border-top: 3px solid black;border-bottom: 1px solid black;">
            <th colspan="6" style="text-align:center;" >Benchmark</th>
        </tr>
        <tr>
            <th colspan="2" style="text-align:center;" >ID</th>
            <th colspan="2" style="text-align:center;" >CATH</th>
            <th colspan="2" style="text-align:center;" ># motifs</th>
        </tr>
        <tr>
            <td colspan="2" style="text-align:center;" >{benchmark[id]}</td>
            <td colspan="2" style="text-align:center;" >{benchmark[cath]}</td>
            <td colspan="2" style="text-align:center;" >{benchmark[motifs]}</td>
        </tr>
        <tr style="border-top: 3px solid black;border-bottom: 1px solid black;">
            <th colspan="6" style="text-align:center;" >Structures</th>
        </tr>
        <tr>
            <th colspan="3" style="text-align:center;" >Query</th>
            <th colspan="3" style="text-align:center;" >Reference</th>
        </tr>
        <tr>
            <td colspan="3" style="text-align:center;" >{structures[query]}</td>
            <td colspan="3" style="text-align:center;" >{structures[reference]}</td>
        </tr>
        <tr style="border-top: 3px solid black;border-bottom: 1px solid black;">
            <th colspan="6" style="text-align:center;" >Design</th>
        </tr>
        <tr>
            <th colspan="2" style="text-align:center;" >motif</th>
            <th colspan="2" style="text-align:center;" >chain</th>
            <th colspan="2" style="text-align:center;" >sequence shift</th>
        </tr>
        <tr>
            <td colspan="2" style="text-align:center;" >{design[motif]}</td>
            <td colspan="2" style="text-align:center;" >{design[chain]}</td>
            <td colspan="2" style="text-align:center;" >{design[shift]}</td>
        </tr>
        <tr style="border-top: 3px solid black;">
            <th colspan="2" style="text-align:center;border-right: 1px solid black;" >Experiments</th>
            <td colspan="2" style="text-align:center;" >{experiments[0]}</td>
            <td colspan="2" style="text-align:center;" >{experiments[1]}</td>
        </tr>
        <tr style="border-top: 3px solid black;border-bottom: 1px solid black;">
            <th colspan="6" style="text-align:center;" >Fragment types</th>
        </tr>
        <tr>
            <th colspan="1" style="text-align:center;" >{fragments[0][0]}</th>
            <td colspan="5" style="text-align:left;" >{fragments[0][1]}</td>
        </tr>
        <tr>
            <th colspan="1" style="text-align:center;" >{fragments[1][0]}</th>
            <td colspan="5" style="text-align:left;" >{fragments[1][1]}</td>
        </tr>
        <tr>
            <th colspan="1" style="text-align:center;" >{fragments[2][0]}</th>
            <td colspan="5" style="text-align:left;" >{fragments[2][1]}</td>
        </tr>
    </table>
    </div>""".format(**info), raw=True)

def plot_fragments( fragments, info, base ):

    for f in ["picker", "wauto"]:
        tst3 = fragments[(fragments["fragments"] == f) & (fragments["size"] == 3)]
        tst9 = fragments[(fragments["fragments"] == f) & (fragments["size"] == 9)]

        fig  = plt.figure(figsize=(60, 15))
        rstoolbox.plot.plot_fragment_profiles(fig, tst3, tst9, base.get_sequence(info["design"]["chain"])[0],
                                              base.get_structure(info["design"]["chain"])[0])
        plt.suptitle(f)
        plt.show()

def plot_main_summary( df, toprmsd=10, top=None ):
    fig = plt.figure(figsize=(40, 10))
    grid = (1, 4)
    ax00 = plt.subplot2grid(grid, (0, 0), fig=fig)
    ax01 = plt.subplot2grid(grid, (0, 1), fig=fig, sharey=ax00)
    ax02 = plt.subplot2grid(grid, (0, 2), fig=fig, sharey=ax00)
    ax03 = plt.subplot2grid(grid, (0, 3), fig=fig, sharey=ax00)

    ax = [ax00, ax01, ax02, ax03]
    tl = ["Global RMSD", "Alignment RMSD", "Motif Segment RMSD", "Query Segment RMSD"]
    for x, r in enumerate(["finalRMSD", "ALIGNRMSD", "MOTIFRMSD", "COMPRRMSD"]):
        sns.boxplot(x="experiment", y=r, hue="fragments", ax=ax[x],
                    data=df if top is None else df.sort_values(["score"]).groupby(["experiment", "fragments"]).head(top),
                    showfliers=False, hue_order=["wauto", "picker", "auto"], order=["nubinitio", "abinitio"])
        ax[x].legend_.remove()
        ax[x].set_xticklabels(["FFL", "abinitio"])
        rstoolbox.utils.add_top_title(ax[x], tl[x])
        if x == 0:
            rstoolbox.utils.add_right_title(ax00, df["benchmark"].values[0])
            ax[x].set_ylabel("RMSD")
        else:
            plt.setp(ax[x].get_yticklabels(), visible=False)
            ax[x].set_ylabel("")

    fig.subplots_adjust(wspace=0, hspace=0)

    fig.legend(handles=[
            mpatches.Patch(color=sns.color_palette()[0], label="wauto"),
            mpatches.Patch(color=sns.color_palette()[1], label="picker"),
            mpatches.Patch(color=sns.color_palette()[2], label="auto")
        ], ncol=3, loc='lower center', borderaxespad=-0.3)

    ax00.set_ylim(0, toprmsd)

    plt.show()

def line_plot( df, ax, bins, linestyle, color, area ):
    """
    Calculate true distribution curves
    (seaborn does an aproximation that does not work well with very small recovery)
    """
    raw, y, x = rstoolbox.analysis.cumulative(df, bins=bins, cumulative=0)
    ax.plot(x, y, color=sns.color_palette()[color], lw=4, linestyle=linestyle)
    if area:
        ax.fill_between(x, 0, y, color=sns.color_palette()[color], alpha=0.3)
    return raw, x, y

def plot_main_distributions( df, toprmsd=10 ):
    fig  = plt.figure(figsize=(30, 10))
    grid = (1, 3)

    ax00 = plt.subplot2grid(grid, (0, 0))
    ax01 = plt.subplot2grid(grid, (0, 1), sharey=ax00)
    ax02 = plt.subplot2grid(grid, (0, 2), sharey=ax00)
    ax   = [ax00, ax01, ax02]
    maxf = 0
    for _, f in enumerate(["wauto", "picker", "auto"]):
        for c, e in enumerate(["abinitio", "nubinitio"]):
            dfs = df[(df["fragments"] == f) & (df["experiment"] == e)]
            r, x, y = line_plot(dfs["COMPRRMSD"], ax[_], 100, "solid", c, False)
            maxf = np.max(y) if np.max(y) > maxf else maxf
            ax[_].set_xlim(0, toprmsd)
        rstoolbox.utils.add_top_title(ax[_], f)
        ax[_].set_xlabel("RMSD")
        if _ == 0:
            rstoolbox.utils.add_right_title(ax00, df["benchmark"].values[0])
            ax[_].set_ylabel("Frequency")
        else:
            plt.setp(ax[_].get_yticklabels(), visible=False)

    ax00.set_ylim(0, maxf + 0.01 if maxf < 0.1 else maxf + 0.05)

    fig.subplots_adjust(wspace=0.08, hspace=0)

    fig.legend(handles=[
            mpatches.Patch(color=sns.color_palette()[0], label="abinitio"),
            mpatches.Patch(color=sns.color_palette()[1], label="FFL")
        ], ncol=2, loc='lower center', borderaxespad=-0.3)

    plt.suptitle("Query Segment RMSD")
    plt.show()

def plot_aa_heatmaps( df, info, base, top=None ):

    core    = rstoolbox.components.Selection(base.get_label("core"))
    motif   = rstoolbox.components.Selection(base.get_label("motif"))
    core    = core - motif
    query   = rstoolbox.components.Selection(base.get_label("query"))
    contact = rstoolbox.components.Selection(base.get_label("contacts"))
    contact = contact - motif
    manual  = rstoolbox.components.Selection(base.get_label("picked"))
    aa = [core, query, contact, manual]

    fig = plt.figure(figsize=(50, 20))
    grid = (1, 4)
    ax00 = plt.subplot2grid(grid, (0, 0), fig=fig)
    ax01 = plt.subplot2grid(grid, (0, 1), fig=fig, sharey=ax00)
    ax02 = plt.subplot2grid(grid, (0, 2), fig=fig, sharey=ax00)
    ax03 = plt.subplot2grid(grid, (0, 3), fig=fig, sharey=ax00)
    ax = [ax00, ax01, ax02, ax03]

    df = df[(df["fragments"] == "wauto") & (df["experiment"] == "nubinitio")].sort_values("score")
    count = df.shape[0] if top is None else int(float(df.shape[0]) * top)
    for _, a in enumerate(aa):
        rstoolbox.plot.sequence_frequency_plot( df.head(count), info["design"]["chain"], ax[_], key_residues=a )

    plt.suptitle("FFL-wauto Sequence Recovery")
    plt.show()

def plot_aa_similarities( df, info, base, matrix="BLOSUM62", top=None ):
    df = df[(df["fragments"] == "wauto") & (df["experiment"] == "nubinitio")].sort_values("score")
    count = df.shape[0] if top is None else int(float(df.shape[0]) * top)
    df = df.head(count)
    stats = rstoolbox.analysis.positional_sequence_similarity(df, info["design"]["chain"])
    fig = plt.figure(figsize=(50, 10))
    grid = (1, 4)
    ax00 = plt.subplot2grid(grid, (0, 0), fig=fig, colspan=4)

    rstoolbox.plot.positional_sequence_similarity_plot(stats, ax00)

    fig.legend(handles=[
            mpatches.Patch(color="green", label="identities"),
            mpatches.Patch(color="orange", label="similarities")
        ], ncol=2, loc='lower center', borderaxespad=-0.3)

    plt.suptitle("FFL-wauto Sequence Recovery (BLOSUM62)")
    plt.show()

def check_success( df, info, base, top, max_rmsd, matrix="BLOSUM62" ):

    data = {"experiment": [], "fragments": [], "value": [], "evaluation": []}
    motif = rstoolbox.components.Selection(base.get_label("motif")) - 1
    for x in info["experiments"]:
        for f in info["fragments"]:
            dft = df[(df["fragments"] == f[0]) & (df["experiment"] == x)].sort_values("score")
            count = int(float(dft.shape[0]) * top)
            dft = dft.head(count)
            # rmsd
            data["experiment"].append(x)
            data["fragments"].append(f[0])
            data["evaluation"].append("rmsd")
            data["value"].append(float(dft[dft["COMPRRMSD"] <= max_rmsd].shape[0])/count)
            # identity
            data["experiment"].append(x)
            data["fragments"].append(f[0])
            data["evaluation"].append("identity")
            seqd = rstoolbox.analysis.positional_sequence_similarity(dft, info["design"]["chain"])
            seqd = seqd.drop(seqd.index[motif.to_list()]).mean()
            data["value"].append(seqd["identity_perc"])
            # similarity
            data["experiment"].append(x)
            data["fragments"].append(f[0])
            data["evaluation"].append("similarity")
            data["value"].append(seqd["positive_perc"])

    data = pd.DataFrame(data)

    fig = plt.figure(figsize=(50, 10))
    grid = (1, 3)
    ax00 = plt.subplot2grid(grid, (0, 0), fig=fig)
    ax01 = plt.subplot2grid(grid, (0, 1), fig=fig)
    ax02 = plt.subplot2grid(grid, (0, 2), fig=fig)

    sns.barplot(x="experiment", y="value", hue="fragments", data=data[data["evaluation"] == "rmsd"],
                ax=ax00, hue_order=["wauto", "picker", "auto"])
    rstoolbox.utils.add_top_title(ax00, "Query RMSD <= {}".format(max_rmsd))
    ax00.set_ylabel("percentage")
    ax00.legend_.remove()

    sns.barplot(x="experiment", y="value", hue="fragments", data=data[data["evaluation"] == "identity"],
                ax=ax01, hue_order=["wauto", "picker", "auto"])
    rstoolbox.utils.add_top_title(ax01, "Mean sequence identity")
    ax01.set_ylabel("percentage")
    ax01.legend_.remove()

    sns.barplot(x="experiment", y="value", hue="fragments", data=data[data["evaluation"] == "similarity"],
                ax=ax02, hue_order=["wauto", "picker", "auto"])
    rstoolbox.utils.add_top_title(ax02, "Mean sequence similarity")
    ax02.set_ylabel("percentage")
    ax02.legend_.remove()

    fig.legend(handles=[
            mpatches.Patch(color=sns.color_palette()[0], label="wauto"),
            mpatches.Patch(color=sns.color_palette()[1], label="picker"),
            mpatches.Patch(color=sns.color_palette()[2], label="auto")
        ], ncol=3, loc='lower center', borderaxespad=-0.3)

    plt.suptitle("Success at top {}%".format(top*100))
    plt.show()












