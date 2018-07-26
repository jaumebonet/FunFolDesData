import os
import glob

import pandas as pd
import rstoolbox

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


def make_paths():
    alldirs = glob.glob("output/*binder*")
    format1 = []
    format2 = []
    for dirs in alldirs:
        if "_min" in dirs or "_pack" in dirs:
            format2.append(dirs)
        else:
            format1.append(dirs)
    return format1, format2


def get_data():
    dirs1, dirs2 = make_paths()
    dirs  = [dirs1, dirs2]
    merge = [
        ["condition", "motif", "helix", "cluster", "decoy"],
        ["condition", "motif", "helix", "binder_state", "cluster", "decoy"]
    ]
    gdf, ddf, cdf, mdf, rdf, df = [], [], [], [], [], []
    for i in range(2):
        gdf.append(rstoolbox.io.parse_rosetta_file(
            [os.path.join(x, "classicFFL_run.minisilent.gz") for x in dirs[i]],
            "output/_definitions/definition_run.{0}.json".format(i+1), True
        ))
        ddf.append(rstoolbox.io.parse_rosetta_file(
            [os.path.join(x, "classicFFL_binder_evaluate.minisilent.gz") for x in dirs[i]],
            "output/_definitions/definition_binder_evaluate.{0}.json".format(i+1), True
        ))
        cdf.append(rstoolbox.io.parse_rosetta_file(
            [os.path.join(x, "classicFFL_binder_clash.score.gz") for x in dirs[i]],
            "output/_definitions/definition_binder_clash.{0}.json".format(i+1), True
        ))
        mdf.append(rstoolbox.io.parse_rosetta_file(
            [os.path.join(x, "classicFFL_minimize_evaluate.score.gz") for x in dirs[i]],
            "output/_definitions/definition_minimize_evaluate.{0}.json".format(i+1), True
        ))
        rdf.append(rstoolbox.io.parse_rosetta_file(
            [os.path.join(x, "classicFFL_extra_rmsd.score.gz") for x in dirs[i]],
            "output/_definitions/definition_extra_rmsd.{0}.json".format(i+1), True
        ))
        df.append(gdf[i].merge(ddf[i], how="left", on=merge[i])
                        .merge(cdf[i], how="left", on=merge[i])
                        .merge(mdf[i], how="left", on=merge[i])
                        .merge(rdf[i], how="left", on=merge[i]))
        if i == 0:
            df[0] = rstoolbox.utils.add_column(df[0], "binder_state", "static")
            df[0].loc[df[0]["condition"] == "nobinder", 'binder_state'] = "nobinder"
    dfall = pd.concat(df, ignore_index=True, sort=True)
    return dfall.replace({'nobinder': 'no_target', 'binder': 'target', 'min': 'packmin'})


def plot_preview(df):
    fig  = plt.figure(figsize=(25, 10))
    rstoolbox.plot.multiple_distributions(df, fig, (2, 5), ['GlobalRMSD', 'LocalRMSDH', 'BUNS', 'bb_clash',
                                                            'ddg_binder_rmsd_drift', 'ddg_design_rmsd_drift',
                                                            'minimize_design_rmsd_drift', 'pre_ddg', 'post_ddg',
                                                            'score'],
                                         x='binder_state', order=['no_target', 'static', 'pack', 'packmin'],
                                         showfliers=False)
    plt.tight_layout()
    plt.show()
