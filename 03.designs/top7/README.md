

```python
# Default Libraries
import os
import StringIO

# External Libraries
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from IPython.display import FileLink, FileLinks

# Own Libraries
import rstoolbox  # Learn about this library at: http://lpdi-epfl.github.io/RosettaSilentToolbox

import readme

# Global Parameters
plt.rcParams['svg.fonttype'] = 'none'
rstoolbox.utils.format_Ipython()
```


# Design of a RSVF site IV carrier from TOP7

> This summary does not track all attempts and variations, just the pathway towards obtaining the designs that were experimentally characterized.


```python
# Direct information obtain from scoring TOP7 with Rosetta
baseline = rstoolbox.io.get_sequence_and_structure(os.path.join('selection', 'baseline', '1qys.pdb'))
baseline
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>score</th>
      <th>fa_atr</th>
      <th>fa_rep</th>
      <th>fa_sol</th>
      <th>fa_intra_rep</th>
      <th>fa_intra_sol_xover4</th>
      <th>lk_ball_wtd</th>
      <th>fa_elec</th>
      <th>pro_close</th>
      <th>hbond_sr_bb</th>
      <th>hbond_lr_bb</th>
      <th>hbond_bb_sc</th>
      <th>hbond_sc</th>
      <th>dslf_fa13</th>
      <th>omega</th>
      <th>fa_dun</th>
      <th>p_aa_pp</th>
      <th>yhh_planarity</th>
      <th>ref</th>
      <th>rama_prepro</th>
      <th>time</th>
      <th>description</th>
      <th>sequence_A</th>
      <th>structure_A</th>
      <th>phi_A</th>
      <th>psi_A</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>191.625</td>
      <td>-498.177</td>
      <td>194.188</td>
      <td>320.719</td>
      <td>1.444</td>
      <td>19.121</td>
      <td>-7.686</td>
      <td>-112.63</td>
      <td>0.0</td>
      <td>-29.587</td>
      <td>-31.008</td>
      <td>-2.424</td>
      <td>-4.577</td>
      <td>0.0</td>
      <td>0.935</td>
      <td>296.765</td>
      <td>-8.478</td>
      <td>0.0</td>
      <td>22.146</td>
      <td>30.876</td>
      <td>2.0</td>
      <td>1qys_0001</td>
      <td>DIQVQVNIDDNGKNFDYTYTVTTESELQKVLNELKDYIKKQGAKRVRISITARTKKEAEKFAAILIKVFAELGYNDINVTFDGDTVTVEGQL</td>
      <td>LEEEEEEEELLLLEEEEEEEELLLLHHHHHHHHHHHHHHHHLLLEEEEEEELLLHHHHHHHHHHHHHHHHHLLLLEEEEEEELLEEEEEEEL</td>
      <td>[0.0, -91.4576, -100.939, ...]</td>
      <td>[-165.465, 111.061, 116.221, ...]</td>
    </tr>
  </tbody>
</table>
</div>



## Generation of fragments to guide the folding process

Structure-based fragments to guide the folding were generated applying the RosettaScript [mkfragments](design/mkfragments.xml) to `1QYS`, the PDB containing **TOP7**.

## Making the designs

The construction of the **TOP7**-derived designs implies the grafting to RSVF's site IV epitope with the presence of its neutralizing antibody 101F [motif.pdb](design/motif.pdb).

The process is split in several steps.

### Folding

Site IV was grafted into **TOP7** by means of the _NubInitioMover_ with the RosettaScript [ffl.xml](design/ffl.xml) into a [SLURM](https://slurm.schedmd.com/)-based cluster through the submiter script [submit_ffl.sh](design/submit_ffl.sh).

A total of 800 decoys were generated without allowing for any sequence mutation to the template residues. The aim of this firts run was to explore and find a folded conformation that would:

* carry and stabilize at backbone level site IV.
* fold back to a conformation close to **TOP7** keeping the proper beta pairing.
* shift enought the conformation to be productive in presenting site IV to 101F (avoid clashes).

From all the decoys, [87_ffl_twostrand_1qys_Lys_0001](design/87_ffl_twostrand_1qys_Lys_0001.pdb) was selected as seed for the next step of the process. The values of this decoy in the context of the design population are marked with a red line.


```python
step1 = rstoolbox.io.parse_rosetta_file('design/top7_ffl_folding.minisilent.gz')
readme.plot_step1(step1)
```


![png](README_files/README_4_0.png)


### Full Design

Full design was carried out applying to [87_ffl_twostrand_1qys_Lys_0001](design/87_ffl_twostrand_1qys_Lys_0001.pdb) the RosettaScript [design_full.xml](design/design_full.xml) submited via the SLURM script [submit_design_full.sh](design/submit_design_full.sh), generating a total of 92 designs from which [sd_top7_12top7_graft_0001](design/sd_top7_12top7_graft_0001.pdb) (highlighted in red in the distributions) was selected.


```python
fulldes = rstoolbox.io.parse_rosetta_file('design/sd_top7_', {'sequence': 'P'})
readme.plot_full(fulldes)
```


![png](README_files/README_6_0.png)


### Part Design

Partial design was carried out applying to [87_ffl_twostrand_1qys_Lys_0001](design/87_ffl_twostrand_1qys_Lys_0001.pdb) the RosettaScript [design_part.xml](design/design_part.xml) submited via the SLURM script [submit_design_part.sh](design/submit_design_part.sh), generating a total of 100 designs from which [sd_top7_nterm2_45top7_graft_0002](design/sd_top7_nterm2_45top7_graft_0002.pdb) (highlighted in red in the distributions) was selected.


```python
partdes = rstoolbox.io.parse_rosetta_file('design/sd_top7_nterm2_', {'sequence': 'P'})
readme.plot_part(partdes)
```


![png](README_files/README_8_0.png)


### Remodel and Resurface

Both selected designs were submited to REMODEL to shorten the first loop of the structure and fix surface residues in order to improve the quality of the decoys, yielding the final two designs as described in the paper.

# Experimental Data


```python
designs = readme.load_experimental()
readme.data_plot(designs)

```


![png](README_files/README_10_0.png)

