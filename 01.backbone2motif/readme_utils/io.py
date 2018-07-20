import os

import yaml
import pandas as pd

import rstoolbox


def load_benchmark_info(filename):
    data = yaml.load(open(filename))
    return data


def load_baseline(data, folder="."):
    score = {"scores": ["score", ],
             "sequence": data["design"]["chain"],
             "structure": data["design"]["chain"],
             "labels": ["CORE", "PICKED", "MOTIF", "QUERY"]}
    sc_file = "{}/baseline/baseline.sc.gz"
    ct_file = "{0}/baseline/contacts_{1}_0001.csv.gz"
    df = rstoolbox.io.parse_rosetta_file(sc_file.format(folder), score)
    ctfile = ct_file.format(folder, data["structures"]["query"].replace(".pdb", ""))
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


def load_main_data( data, base, folder=".", ignore=None ):
    scores = {"scores_ignore": ["loop*", "niccd_*", "fa_*", "hbond_*", "pro_close", "dslf_fa13", "atom_pair_constraint",
                                "rama", "omega", "p_aa_pp",  "yhh_planarity", "ref", "aa_composition", "time"],
              "sequence": data["design"]["chain"],
              "naming": ["experiment", "fragments"]
             }
    df = []
    ignore = ignore if ignore is not None else []
    for e in data["experiments"]:
        for f in data["fragments"]:
            if f[0] not in ignore:
                infile = "{0}/output/{1}_{2}.minisilent.gz".format(folder, e, f[0])
                if os.path.isfile(infile):
                    df.append(rstoolbox.io.parse_rosetta_file(infile, scores))
    df = pd.concat(df)
    df.add_reference( data["design"]["chain"],
                      base.get_sequence(data["design"]["chain"])[0],
                      shift=data["design"]["shift"])
    return rstoolbox.utils.add_column( df, "benchmark", data["benchmark"]["id"])


def load_hmm_data( df, info, coverage, evalue ):
    import re
    import pandas as pd
    import rstoolbox

    # Family Matches
    cath_fam = 'hmm/cath-superfamily-seqs-{}.fa'.format(info['benchmark']['cath'].replace('CATH.', ''))
    cath_siz = rstoolbox.io.read_fasta(cath_fam).shape[0]
    cath_hmm = rstoolbox.io.read_hmmsearch('hmm/cath_family.hmm')
    cath_hmm = cath_hmm[(cath_hmm['full-e-value'] <= evalue) & (cath_hmm['acc'] > coverage)].drop_duplicates('description')
    cath_cov = float(cath_hmm.shape[0])/cath_siz

    # Design Matches
    data = {"experiment": [], "fragments": [], "count": [], "benchmark": [], "family": []}
    for g in df.groupby(["experiment", "fragments"]):
        fafile = os.path.join("hmm", "sequence_{}_{}.fa".format(g[0][0], g[0][1]))
        if not os.path.isfile(fafile):
            rstoolbox.io.write_fasta(g[1], seqID=info['design']['chain'], filename=fafile)
        facount = rstoolbox.io.read_fasta(fafile).shape[0]
        hmmfile = os.path.join("hmm", "search_{}_{}.hmm".format(g[0][0], g[0][1]))
        if not os.path.isfile(hmmfile):
            os.system("hmmsearch --max --noali -E 100 hmm/{2} {0} > {1}".format(fafile, hmmfile, info['benchmark']['hmm']))
        dhmm = rstoolbox.io.read_hmmsearch(hmmfile)
        dhmm = dhmm[(dhmm['full-e-value'] <= evalue) & (dhmm['acc'] > coverage)].drop_duplicates('description')
        data['experiment'].append(g[0][0])
        data['fragments'].append(g[0][1])
        data['benchmark'].append(info['benchmark']['id'])
        data['count'].append(float(dhmm.shape[0])/facount)
        data['family'].append(cath_cov)
    return pd.DataFrame(data)
