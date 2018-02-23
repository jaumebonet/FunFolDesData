import os

import yaml
import pandas as pd

import rstoolbox

def load_benchmark_info( filename ):
    data = yaml.load(open(filename))
    return data

def load_baseline( data, folder="." ):
    score = {"scores": ["score",],
             "sequence": data["design"]["chain"],
             "structure": data["design"]["chain"],
             "labels": ["CORE", "PICKED", "MOTIF", "QUERY"]}
    df = rstoolbox.io.parse_rosetta_file("baseline/baseline.sc.gz", score)
    ctfile = "{0}/baseline/contacts_{1}_0001.csv.gz".format(folder, data["structures"]["query"].replace(".pdb", ""))
    ct, r1, r2 = rstoolbox.io.parse_rosetta_contacts(ctfile)
    selection = rstoolbox.components.Selection(data["design"]["motif"])
    motif = [x - data["design"]["shift"] + 1 for x in selection]
    nomotif = sorted(list(set(ct.index.values) - set(motif)))
    ct = ct[motif].loc[nomotif]
    contacts = rstoolbox.components.Selection(list(ct[(ct.T != 0).any()].index.values))
    return rstoolbox.utils.add_column(df, "lbl_CONTACTS", str(contacts))

def load_fragments( data, folder="." ):
    name = data["structures"]["query"].replace(".pdb", "")
    df   = []

    for f in data["fragments"]:
        for x in [3, 9]:
            df.append(rstoolbox.io.parse_rosetta_fragments("{0}/fragments/{1}.{2}.200.{3}mers".format(folder, name, f[0], x)))
            df[-1] = rstoolbox.utils.add_column(df[-1], "fragments", f[0])
            df[-1] = df[-1].add_quality_measure("{0}/fragments/{1}.{2}.200.{3}mers.qual".format(folder, name, f[0], x))

    return rstoolbox.components.FragmentFrame(pd.concat(df))

def load_main_data( data, base, folder="." ):
    scores = {"scores_ignore": ["loop*", "niccd_*", "fa_*", "hbond_*", "pro_close", "dslf_fa13", "atom_pair_constraint",
                                "rama", "omega", "p_aa_pp",  "yhh_planarity", "ref", "aa_composition", "time"],
              "sequence": data["design"]["chain"],
              "naming": ["experiment", "fragments"]
             }
    df = []
    for e in data["experiments"]:
        for f in data["fragments"]:
            infile = "{0}/output/{1}_{2}.minisilent.gz".format(folder, e, f[0])
            if os.path.isfile(infile):
                df.append(rstoolbox.io.parse_rosetta_file(infile, scores))
    df = pd.concat(df)
    df.add_reference( data["design"]["chain"],
                      base.get_sequence(data["design"]["chain"])[0],
                      shift=data["design"]["shift"])
    return rstoolbox.utils.add_column( df, "benchmark", data["benchmark"]["id"])
