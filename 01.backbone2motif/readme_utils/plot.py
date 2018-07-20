import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from IPython.display import display_html

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
            <td colspan="2" style="text-align:center;" >FunFolDes</td>
        </tr>
        <tr style="border-top: 3px solid black;border-bottom: 1px solid black;">
            <th colspan="6" style="text-align:center;" >Fragment types</th>
        </tr>
        <tr>
            <th colspan="2" style="text-align:center;" >{fragments[0][2]}</th>
            <td colspan="4" style="text-align:left;" >{fragments[0][1]}</td>
        </tr>
        <tr>
            <th colspan="2" style="text-align:center;" >{fragments[1][2]}</th>
            <td colspan="4" style="text-align:left;" >{fragments[1][1]}</td>
        </tr>
    </table>
    </div>""".format(**info), raw=True)

def plot_fragments( fragments, info, base ):

    for f in info['fragments']:
        tst3 = fragments[(fragments["fragments"] == f[0]) & (fragments["size"] == 3)]
        tst9 = fragments[(fragments["fragments"] == f[0]) & (fragments["size"] == 9)]

        fig = plt.figure(figsize=(60, 15))
        rstoolbox.plot.plot_fragment_profiles(fig, tst3, tst9, base.get_sequence(info["design"]["chain"])[0],
                                              base.get_structure(info["design"]["chain"])[0])
        plt.suptitle(f[2])
        plt.show()


def plot_coverage( df ):
    gf = df.groupby(["experiment", "fragments"]).count()["description"].reset_index()
    gf = gf.replace('nubinitio', 'FunFolDes')
    matplotlib.rcParams.update({'font.size': 10})
    g = sns.factorplot(x="fragments", y="description", col="experiment",
                       order=["picker", "wauto"], data=gf, kind="bar",
                       size=4, aspect=.7)
    for axlist in g.axes:
        for ax in axlist:
            ax.set_ylim(0, 10296)
            ax.set_xticklabels(['seq-based', 'str-based'])
        axlist[0].set_ylabel('Number of Decoys')
    plt.show()


def plot_main_summary( df, toprmsd=10, top=None ):
    matplotlib.rcParams.update({'font.size': 30})
    df = df.replace('nubinitio', 'FunFolDes')
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
                    showfliers=False, hue_order=["wauto", "picker"], order=["FunFolDes", "abinitio"])
        ax[x].legend_.remove()
        ax[x].set_xticklabels(["FunFolDes", "abinitio"])
        rstoolbox.utils.add_top_title(ax[x], tl[x])
        if x == 0:
            rstoolbox.utils.add_right_title(ax00, df["benchmark"].values[0])
            ax[x].set_ylabel("RMSD")
        else:
            plt.setp(ax[x].get_yticklabels(), visible=False)
            ax[x].set_ylabel("")

    fig.subplots_adjust(wspace=0, hspace=0)

    fig.legend(handles=[
            mpatches.Patch(color=sns.color_palette()[0], label="str-based fragments"),
            mpatches.Patch(color=sns.color_palette()[1], label="seq-based fragments")
        ], ncol=2, loc='lower center', borderaxespad=-0.3)
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
    matplotlib.rcParams.update({'font.size': 25})
    df = df.replace('nubinitio', 'FunFolDes')
    fig  = plt.figure(figsize=(30, 10))
    grid = (1, 2)

    ax00 = plt.subplot2grid(grid, (0, 0), fig=fig)
    ax01 = plt.subplot2grid(grid, (0, 1), sharey=ax00, fig=fig)
    ax   = [ax00, ax01]
    maxf = 0
    for _, f in enumerate(["wauto", "picker"]):
        for c, e in enumerate(["abinitio", "FunFolDes"]):
            dfs = df[(df["fragments"] == f) & (df["experiment"] == e)]
            r, x, y = line_plot(dfs["COMPRRMSD"], ax[_], 100, "solid", c, False)
            maxf = np.max(y) if np.max(y) > maxf else maxf
            ax[_].set_xlim(0, toprmsd)
        ss = 'seq' if f == 'picker' else 'str'
        rstoolbox.utils.add_top_title(ax[_], '{}-based fragments'.format(ss))
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
            mpatches.Patch(color=sns.color_palette()[1], label="FunFolDes")
        ], ncol=2, loc='lower center', borderaxespad=-0.3)

    plt.suptitle("Query Segment RMSD")
    plt.show()


def plot_aa_heatmaps( df, info, base, top=None ):
    matplotlib.rcParams.update({'font.size': 15})
    core    = rstoolbox.components.Selection(base.iloc[0].get_label("core"))
    motif   = rstoolbox.components.Selection(base.iloc[0].get_label("motif"))
    core    = core - motif
    query   = rstoolbox.components.Selection(base.iloc[0].get_label("query"))
    contact = rstoolbox.components.Selection(base.iloc[0].get_label("contacts"))
    contact = contact - motif
    manual  = rstoolbox.components.Selection(base.iloc[0].get_label("picked"))
    aa = [core, query, contact, manual]
    tl = ["CORE", "QUERY REGION", "MOTIF CONTACT", "MANUAL"]

    fig = plt.figure(figsize=(15, 25))
    grid = (4, 1)
    ax00 = plt.subplot2grid(grid, (0, 0), fig=fig)
    ax01 = plt.subplot2grid(grid, (1, 0), fig=fig, sharey=ax00)
    ax02 = plt.subplot2grid(grid, (2, 0), fig=fig, sharey=ax00)
    ax03 = plt.subplot2grid(grid, (3, 0), fig=fig, sharey=ax00)
    ax = [ax00, ax01, ax02, ax03]

    df = df[(df["fragments"] == "wauto") & (df["experiment"] == "nubinitio")].sort_values("score")
    count = df.shape[0] if top is None else int(float(df.shape[0]) * top)
    for _, a in enumerate(aa):
        if not a.is_empty():
            rstoolbox.plot.sequence_frequency_plot( df.head(count), info["design"]["chain"],
                                                   ax[_], key_residues=a, xrotation=90,
                                                   cbar_kws={"orientation": "vertical"} )
        rstoolbox.utils.add_top_title(ax[_], tl[_])

    plt.suptitle("FunFolDes & Seq-based Fragments Sequence Recovery")
    plt.show()


def plot_aa_similarities( df, info, base, matrix="BLOSUM62", top=None ):
    matplotlib.rcParams.update({'font.size': 30})
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

    plt.suptitle("FunFolDes & Seq-based Fragments Sequence Recovery (BLOSUM62)")
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

def plot_global_main_distributions( df, info, base, toprmsd=10, topseq=0.5, top=0.1 ):
    fig  = plt.figure(figsize=(50, 140))
    grid = (14, 5)
    tl   = ["Global RMSD", "Alignment RMSD", "Motif Segment RMSD", "Query Segment RMSD"]

    for _i, data in enumerate(df):
        ax = []
        ax.append(plt.subplot2grid(grid, (_i, 0), fig=fig))
        ax.append(plt.subplot2grid(grid, (_i, 1), sharey=ax[0], fig=fig))
        ax.append(plt.subplot2grid(grid, (_i, 2), sharey=ax[0], fig=fig))
        ax.append(plt.subplot2grid(grid, (_i, 3), sharey=ax[0], fig=fig))
        ax.append(plt.subplot2grid(grid, (_i, 4), fig=fig))

        for x, r in enumerate(["finalRMSD", "ALIGNRMSD", "MOTIFRMSD", "COMPRRMSD"]):
            sns.boxplot(x="experiment", y=r, hue="fragments", ax=ax[x], data=data,
                        showfliers=False, hue_order=["wauto", "picker"], order=["nubinitio", "abinitio"])
            ax[x].legend_.remove()
            if _i == len(df) - 1:
                ax[x].set_xticklabels(["FFL", "abinitio"])
            else:
                ax[x].set_xticklabels([])
            if _i == 0:
                rstoolbox.utils.add_top_title(ax[x], tl[x])
            if x == 0:
                rstoolbox.utils.add_right_title(ax[x], data["benchmark"].values[0])
                ax[x].set_ylabel("RMSD")
            else:
                plt.setp(ax[x].get_yticklabels(), visible=False)
                ax[x].set_ylabel("")
        ax[0].set_ylim(0, toprmsd)

        ax[-1].yaxis.tick_right()
        seqlen = len(base[_i].get_sequence(info[_i]["design"]["chain"])[0])
        allres = rstoolbox.components.Selection("1-{}".format(seqlen))
        motif  = rstoolbox.components.Selection(base[_i].get_label("motif"))
        selection = allres - motif + info[_i]["design"]["shift"] - 1
        mdata = data[(data["fragments"] == "wauto") & (data["experiment"] == "nubinitio")]
        mdata = rstoolbox.analysis.sequence_similarity(data, info[_i]["design"]["chain"], selection.to_list() )
        identity = "blosum62_{}_identity".format(info[_i]["design"]["chain"])
        positive = "blosum62_{}_positive".format(info[_i]["design"]["chain"])
        mdata["identity_perc"] = mdata[identity]/float(seqlen)
        mdata["positive_perc"] = mdata[positive]/float(seqlen)
        mdata = rstoolbox.utils.split_values(mdata, {"split": [("identity_perc", "identity"), ("positive_perc", "positive")],
                                                     "names": ["recovery", "simtype"], "keep":["description"]})

        sns.boxplot(x="simtype", y="recovery", ax=ax[-1], data=mdata, showfliers=True, order=["identity", "positive"],
                    color=sns.color_palette()[2])
        if _i == 0:
            rstoolbox.utils.add_top_title(ax[-1], "Sequence Recovery FFL-autofrags")
        ax[-1].set_ylim(0, topseq)

    fig.subplots_adjust(wspace=0.01, hspace=0.07)

    fig.legend(handles=[
            mpatches.Patch(color=sns.color_palette()[0], label="automatic fragments"),
            mpatches.Patch(color=sns.color_palette()[1], label="picker fragments"),
            mpatches.Patch(color=sns.color_palette()[2], label="sequence recovery")
        ], ncol=3, loc='lower center')

    plt.savefig("SFig1.benchmark_overview.preview.png")
    plt.savefig("SFig1.benchmark_overview.preview.svg")
    plt.show()


def plot_global_query_distributions(df, info, base, toprmsd=10):
    fig  = plt.figure(figsize=(40, 40))
    grid = (4, 4)

    counter = 0
    for _i in range(4):
        for _j in range(4):
            if counter >= len(df):
                break
            ax = plt.subplot2grid(grid, (_i, _j), fig=fig)
            maxf = 0
            for c, f in enumerate(["wauto", "picker"]):
                dfs = df[counter][(df[counter]["experiment"] == "nubinitio") & (df[counter]["fragments"] == f)]
                r, x, y = line_plot(dfs["COMPRRMSD"], ax, 100, "solid", c, False)
                maxf = np.max(y) if np.max(y) > maxf else maxf
                ax.set_xlim(0, toprmsd)
            rstoolbox.utils.add_top_title(ax, info[counter]["benchmark"]["id"])
            ax.set_xlabel("RMSD")
            counter += 1

    fig.legend(handles=[
            mpatches.Patch(color=sns.color_palette()[0], label="automatic fragments"),
            mpatches.Patch(color=sns.color_palette()[1], label="picker fragments")
        ], ncol=2, loc='lower center', borderaxespad=-0.3)

    plt.suptitle("Query RMSD distribution")
    plt.savefig("SFig1B.benchmark_overview.queryRMSD.auto_vs_picker.png")
    plt.savefig("SFig1B.benchmark_overview.queryRMSD.auto_vs_picker.svg")
    plt.tight_layout()
    plt.show()


def plot_hmm(df):
    df = df.replace('nubinitio', 'FunFolDes')
    matplotlib.rcParams.update({'font.size': 10})
    g = sns.factorplot(x="fragments", y="count", col="experiment", order=["wauto", "picker"],
                       data=df, kind="bar", size=4, aspect=.7)
    for axlist in g.axes:
        for ax in axlist:
            ax.axhline(df.iloc[0]['family'])
            ax.set_xticklabels(['str-based', 'seq-based'])
    plt.show()
