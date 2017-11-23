

```python
%load_ext autoreload
%reload_ext autoreload
%autoreload 2
```


```python
import warnings; warnings.simplefilter('ignore')
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import rstoolbox
import pandas as pd
import seaborn as sns
import numpy as np
import copy
import glob
import gzip
import os

sns.set(font_scale=3)
matplotlib.rcParams['svg.fonttype'] = 'none'
pd.set_option('display.max_columns', None)
sns.set_style("whitegrid")
```

# Scaffold Search

For each motif of interest, we perform a [_master_](http://www.grigoryanlab.org/master/) search over the __PDBnr30__ database provided by _master_ and then we use _Rosetta_ to place the binder of the motif according to the motif alignment to the matched scaffold, clean the sidechains, and evaluate clashes.

|__INFO__|__H__|__E__|__HLH__|__ELE__|__HLE__|__HxL__|
|:------:|:---:|:---:|:-----:|:-----:|:-----:|:-----:|
|Query   |2fx7 |3o41 |3ixt   |3vtt   |5tpn   |4jhw   |
|Match   |1z6n |     |3lhp   |3wei   |       |       |
|Design  |3lef |     |4l8i   |4x42   |       |       |
||![2fx7](images/2fx7.png)|![3o41](images/3o41.png)|![3ixt](images/3ixt.png)|![3vtt](images/3vtt.png)|![5tpn](images/5tpn.png)|![4jhw](images/4jhw.png)|


```python
dfs ={"5tpn": 1, "2fx7": 1, "3ixt": 1, "3vtt": 1, "3o41": 1, "4jhw": 1}
for k in dfs:
    if not isinstance(dfs[k], int):
        continue
    dfs[k] = pd.read_csv(os.path.join(k, "master_search.csv.gz"))
    if os.path.isfile(os.path.join(k, "ddg_match.csv.gz")):
        tmp = pd.read_csv(os.path.join(k, "ddg_match.csv.gz"))
        dfs[k] = dfs[k].merge(tmp, how="left", on=["cluster", "str"])
total_master_list = 17539
```

## Pfam assignations

Used to check domain length. The file ```pdbmap.gz```, is downloaded from: [```ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/```](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/)


```python
def parse_pfam( f ):
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
    r = [int(x) for x in pfam_range.split("-")]
    return r[1] - r[0] + 1

pfam = parse_pfam("pdbmap.gz")
pfam["length"] = pfam.apply(lambda row: pfam_length(row["pfamrange"]), axis=1)

```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>chain</th>
      <th>pdb</th>
      <th>pfam</th>
      <th>pfamrange</th>
      <th>length</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>B</td>
      <td>3h3m</td>
      <td>PF05400</td>
      <td>34-106</td>
      <td>73</td>
    </tr>
    <tr>
      <th>1</th>
      <td>A</td>
      <td>3h3m</td>
      <td>PF05400</td>
      <td>34-106</td>
      <td>73</td>
    </tr>
    <tr>
      <th>2</th>
      <td>C</td>
      <td>3nkz</td>
      <td>PF05400</td>
      <td>21-92</td>
      <td>72</td>
    </tr>
    <tr>
      <th>3</th>
      <td>B</td>
      <td>3nkz</td>
      <td>PF05400</td>
      <td>21-94</td>
      <td>74</td>
    </tr>
    <tr>
      <th>4</th>
      <td>D</td>
      <td>3nkz</td>
      <td>PF05400</td>
      <td>21-93</td>
      <td>73</td>
    </tr>
  </tbody>
</table>
</div>




```python
def fit_in_range( motif_range, pfam_range ):
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
    
for k in dfs:
    if isinstance(dfs[k], pd.DataFrame):
        dfs[k] = dfs[k].merge( pfam, how="left", on=["pdb", "chain"])
        dfs[k]["inrange"] = dfs[k].apply(lambda row: fit_in_range(row["range"], row["pfamrange"]), axis=1)
```


```python
tmp = dfs["4jhw"]
tmp["tt"] = tmp.apply(lambda row: fit_in_range(row["range"], row["pfamrange"]), axis=1)
tmp
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>chain</th>
      <th>cluster</th>
      <th>pdb</th>
      <th>range</th>
      <th>rmsd</th>
      <th>str</th>
      <th>ddg</th>
      <th>pfam</th>
      <th>pfamrange</th>
      <th>length</th>
      <th>inrange</th>
      <th>tt</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>e</td>
      <td>e</td>
      <td>eeee</td>
      <td>none</td>
      <td>0.0000</td>
      <td>error</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>1</th>
      <td>e</td>
      <td>e</td>
      <td>eeee</td>
      <td>none</td>
      <td>0.0000</td>
      <td>error</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>2</th>
      <td>e</td>
      <td>e</td>
      <td>eeee</td>
      <td>none</td>
      <td>0.0000</td>
      <td>error</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>3</th>
      <td>e</td>
      <td>e</td>
      <td>eeee</td>
      <td>none</td>
      <td>0.0000</td>
      <td>error</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>4</th>
      <td>e</td>
      <td>e</td>
      <td>eeee</td>
      <td>none</td>
      <td>0.0000</td>
      <td>error</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>5</th>
      <td>e</td>
      <td>e</td>
      <td>eeee</td>
      <td>none</td>
      <td>0.0000</td>
      <td>error</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>6</th>
      <td>e</td>
      <td>e</td>
      <td>eeee</td>
      <td>none</td>
      <td>0.0000</td>
      <td>error</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>7</th>
      <td>e</td>
      <td>e</td>
      <td>eeee</td>
      <td>none</td>
      <td>0.0000</td>
      <td>error</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>8</th>
      <td>e</td>
      <td>e</td>
      <td>eeee</td>
      <td>none</td>
      <td>0.0000</td>
      <td>error</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>9</th>
      <td>A</td>
      <td>x417</td>
      <td>3kgd</td>
      <td>[(164-171)(129-145)]</td>
      <td>1.5985</td>
      <td>full0001</td>
      <td>31472.695</td>
      <td>PF05189</td>
      <td>183-275</td>
      <td>93.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>10</th>
      <td>A</td>
      <td>x417</td>
      <td>3kgd</td>
      <td>[(164-171)(129-145)]</td>
      <td>1.5985</td>
      <td>full0001</td>
      <td>31472.695</td>
      <td>PF01137</td>
      <td>11-326</td>
      <td>316.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>11</th>
      <td>A</td>
      <td>x606</td>
      <td>3pqv</td>
      <td>[(160-167)(124-140)]</td>
      <td>1.8467</td>
      <td>full0001</td>
      <td>221.653</td>
      <td>PF05189</td>
      <td>182-282</td>
      <td>101.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>12</th>
      <td>A</td>
      <td>x606</td>
      <td>3pqv</td>
      <td>[(160-167)(124-140)]</td>
      <td>1.8467</td>
      <td>full0001</td>
      <td>221.653</td>
      <td>PF01137</td>
      <td>8-337</td>
      <td>330.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>13</th>
      <td>A</td>
      <td>x100</td>
      <td>3chm</td>
      <td>[(146-153)(120-136)]</td>
      <td>1.9302</td>
      <td>full0001</td>
      <td>8757.513</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>14</th>
      <td>A</td>
      <td>x868</td>
      <td>3zph</td>
      <td>[(191-198)(156-172)]</td>
      <td>1.9517</td>
      <td>full00001</td>
      <td>17415.512</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>15</th>
      <td>A</td>
      <td>x613</td>
      <td>1pyt</td>
      <td>[(6-13)(51-67)]</td>
      <td>1.9875</td>
      <td>full0001</td>
      <td>3539.164</td>
      <td>PF02244</td>
      <td>26-100</td>
      <td>75.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>16</th>
      <td>A</td>
      <td>x859</td>
      <td>2zfh</td>
      <td>[(58-65)(16-32)]</td>
      <td>1.9937</td>
      <td>full00001</td>
      <td>12759.057</td>
      <td>PF03091</td>
      <td>69-167</td>
      <td>99.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>17</th>
      <td>P</td>
      <td>x807</td>
      <td>2wzp</td>
      <td>[(123-130)(74-90)]</td>
      <td>1.9970</td>
      <td>full00001</td>
      <td>5952.324</td>
      <td>PF16774</td>
      <td>2-160</td>
      <td>159.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>18</th>
      <td>A</td>
      <td>x513</td>
      <td>3ni8</td>
      <td>[(4-11)(119-135)]</td>
      <td>2.0018</td>
      <td>full0001</td>
      <td>1881.634</td>
      <td>PF08327</td>
      <td>10-135</td>
      <td>126.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>19</th>
      <td>A</td>
      <td>x570</td>
      <td>1osc</td>
      <td>[(56-63)(14-30)]</td>
      <td>2.0147</td>
      <td>full0001</td>
      <td>12682.082</td>
      <td>PF03091</td>
      <td>67-165</td>
      <td>99.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>20</th>
      <td>A</td>
      <td>x136</td>
      <td>3dhx</td>
      <td>[(19-26)(72-88)]</td>
      <td>2.0173</td>
      <td>full00001</td>
      <td>5232.319</td>
      <td>PF09383</td>
      <td>267-339</td>
      <td>73.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>21</th>
      <td>A</td>
      <td>x406</td>
      <td>1k47</td>
      <td>[(192-199)(297-313)]</td>
      <td>2.0296</td>
      <td>full00001</td>
      <td>4157.014</td>
      <td>PF00288</td>
      <td>89-157</td>
      <td>69.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>22</th>
      <td>A</td>
      <td>x818</td>
      <td>1xk8</td>
      <td>[(56-63)(14-30)]</td>
      <td>2.0302</td>
      <td>full0001</td>
      <td>12702.917</td>
      <td>PF03091</td>
      <td>69-167</td>
      <td>99.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>23</th>
      <td>A</td>
      <td>x395</td>
      <td>1jqg</td>
      <td>[(6-13)(48-64)]</td>
      <td>2.0352</td>
      <td>full0001</td>
      <td>4115.408</td>
      <td>Propep_M14</td>
      <td>27-98</td>
      <td>72.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>24</th>
      <td>A</td>
      <td>x395</td>
      <td>1jqg</td>
      <td>[(6-13)(48-64)]</td>
      <td>2.0352</td>
      <td>full0001</td>
      <td>4115.408</td>
      <td>PF00246</td>
      <td>127-412</td>
      <td>286.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>25</th>
      <td>A</td>
      <td>x602</td>
      <td>3pm9</td>
      <td>[(220-227)(290-306)]</td>
      <td>2.0513</td>
      <td>full00001</td>
      <td>8157.496</td>
      <td>PF01565</td>
      <td>51-190</td>
      <td>140.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>26</th>
      <td>A</td>
      <td>x602</td>
      <td>3pm9</td>
      <td>[(220-227)(290-306)]</td>
      <td>2.0513</td>
      <td>full00001</td>
      <td>8157.496</td>
      <td>PF02913</td>
      <td>226-474</td>
      <td>249.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>27</th>
      <td>A</td>
      <td>x766</td>
      <td>1vk3</td>
      <td>[(265-272)(196-212)]</td>
      <td>2.0674</td>
      <td>full0001</td>
      <td>70773.969</td>
      <td>PF00586</td>
      <td>52-158</td>
      <td>107.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>28</th>
      <td>A</td>
      <td>x766</td>
      <td>1vk3</td>
      <td>[(265-272)(196-212)]</td>
      <td>2.0674</td>
      <td>full0001</td>
      <td>70773.969</td>
      <td>PF00586</td>
      <td>389-502</td>
      <td>114.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>29</th>
      <td>A</td>
      <td>x766</td>
      <td>1vk3</td>
      <td>[(265-272)(196-212)]</td>
      <td>2.0674</td>
      <td>full0001</td>
      <td>70773.969</td>
      <td>PF02769</td>
      <td>170-325</td>
      <td>156.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>11498027</th>
      <td>A</td>
      <td>x408</td>
      <td>4k5l</td>
      <td>[(798-805)(836-852)]</td>
      <td>5.1000</td>
      <td>full17638</td>
      <td>62665.434</td>
      <td>PF17432</td>
      <td>748-1083</td>
      <td>336.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>11498028</th>
      <td>A</td>
      <td>x483</td>
      <td>1mi1</td>
      <td>[(187-194)(333-349)]</td>
      <td>5.1000</td>
      <td>full10474</td>
      <td>132523.656</td>
      <td>PF02138</td>
      <td>2287-2563</td>
      <td>277.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498029</th>
      <td>A</td>
      <td>x483</td>
      <td>1mi1</td>
      <td>[(187-194)(333-349)]</td>
      <td>5.1000</td>
      <td>full10474</td>
      <td>132523.656</td>
      <td>PF14844</td>
      <td>2158-2255</td>
      <td>98.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498030</th>
      <td>A</td>
      <td>x718</td>
      <td>4tq3</td>
      <td>[(112-119)(157-173)]</td>
      <td>5.1000</td>
      <td>full16205</td>
      <td>31464.035</td>
      <td>PF01040</td>
      <td>29-275</td>
      <td>247.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>11498031</th>
      <td>A</td>
      <td>x558</td>
      <td>1okj</td>
      <td>[(126-133)(93-109)]</td>
      <td>5.1000</td>
      <td>full7396</td>
      <td>53445.531</td>
      <td>PF00814</td>
      <td>26-155</td>
      <td>130.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>11498032</th>
      <td>A</td>
      <td>x083</td>
      <td>2c2j</td>
      <td>[(17-24)(67-83)]</td>
      <td>5.1000</td>
      <td>full17388</td>
      <td>2064.492</td>
      <td>PF00210</td>
      <td>79-222</td>
      <td>144.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498033</th>
      <td>A</td>
      <td>x761</td>
      <td>1vfq</td>
      <td>[(56-63)(115-131)]</td>
      <td>5.1000</td>
      <td>full5825</td>
      <td>14516.954</td>
      <td>PF00241</td>
      <td>10-128</td>
      <td>119.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498034</th>
      <td>A</td>
      <td>x610</td>
      <td>4put</td>
      <td>[(665-672)(497-513)]</td>
      <td>5.1000</td>
      <td>full10023</td>
      <td>130886.227</td>
      <td>PF01432</td>
      <td>235-696</td>
      <td>462.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>11498035</th>
      <td>A</td>
      <td>x428</td>
      <td>4ks9</td>
      <td>[(74-81)(19-35)]</td>
      <td>5.1000</td>
      <td>full13525</td>
      <td>5609.623</td>
      <td>PF17408</td>
      <td>104-191</td>
      <td>88.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498036</th>
      <td>A</td>
      <td>x428</td>
      <td>4ks9</td>
      <td>[(74-81)(19-35)]</td>
      <td>5.1000</td>
      <td>full13525</td>
      <td>5609.623</td>
      <td>PF05292</td>
      <td>194-447</td>
      <td>254.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498037</th>
      <td>A</td>
      <td>x433</td>
      <td>3kxp</td>
      <td>[(178-185)(144-160)]</td>
      <td>5.1000</td>
      <td>full12000</td>
      <td>35696.746</td>
      <td>PF00561</td>
      <td>33-146</td>
      <td>114.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498038</th>
      <td>A</td>
      <td>x433</td>
      <td>3kxp</td>
      <td>[(178-185)(144-160)]</td>
      <td>5.1000</td>
      <td>full12000</td>
      <td>35696.746</td>
      <td>PF12146</td>
      <td>31-262</td>
      <td>232.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>11498039</th>
      <td>A</td>
      <td>x433</td>
      <td>3kxp</td>
      <td>[(178-185)(144-160)]</td>
      <td>5.1000</td>
      <td>full12000</td>
      <td>35696.746</td>
      <td>PF12697</td>
      <td>35-271</td>
      <td>237.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>11498040</th>
      <td>A</td>
      <td>x661</td>
      <td>2r76</td>
      <td>[(8-15)(104-120)]</td>
      <td>5.1000</td>
      <td>full15529</td>
      <td>14645.959</td>
      <td>PF04390</td>
      <td>36-158</td>
      <td>123.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498041</th>
      <td>A</td>
      <td>x106</td>
      <td>3cmw</td>
      <td>[(840-847)(801-817)]</td>
      <td>5.1000</td>
      <td>full11170</td>
      <td>132705.125</td>
      <td>PF00154</td>
      <td>9-270</td>
      <td>262.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498042</th>
      <td>A</td>
      <td>x802</td>
      <td>1wr2</td>
      <td>[(131-138)(158-174)]</td>
      <td>5.1000</td>
      <td>full6645</td>
      <td>55741.613</td>
      <td>PF13549</td>
      <td>11-234</td>
      <td>224.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>11498043</th>
      <td>A</td>
      <td>x099</td>
      <td>4cgv</td>
      <td>[(85-92)(54-70)]</td>
      <td>5.1000</td>
      <td>full10178</td>
      <td>10947.659</td>
      <td>PF00515</td>
      <td>133-166</td>
      <td>34.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498044</th>
      <td>A</td>
      <td>x099</td>
      <td>4cgv</td>
      <td>[(85-92)(54-70)]</td>
      <td>5.1000</td>
      <td>full10178</td>
      <td>10947.659</td>
      <td>PF00515</td>
      <td>201-234</td>
      <td>34.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498045</th>
      <td>A</td>
      <td>x037</td>
      <td>3azo</td>
      <td>[(594-601)(505-521)]</td>
      <td>5.1000</td>
      <td>full8465</td>
      <td>87533.398</td>
      <td>PF01738</td>
      <td>469-648</td>
      <td>180.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>11498046</th>
      <td>A</td>
      <td>x037</td>
      <td>3azo</td>
      <td>[(594-601)(505-521)]</td>
      <td>5.1000</td>
      <td>full8465</td>
      <td>87533.398</td>
      <td>PF00326</td>
      <td>440-651</td>
      <td>212.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>11498047</th>
      <td>O</td>
      <td>x686</td>
      <td>1s72</td>
      <td>[(35-42)(2-18)]</td>
      <td>5.1000</td>
      <td>full5404</td>
      <td>10627.818</td>
      <td>PF00828</td>
      <td>23-99</td>
      <td>77.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498048</th>
      <td>A</td>
      <td>x352</td>
      <td>2ii2</td>
      <td>[(156-163)(195-211)]</td>
      <td>5.1000</td>
      <td>full15170</td>
      <td>143836.109</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498049</th>
      <td>A</td>
      <td>x823</td>
      <td>1xsi</td>
      <td>[(623-630)(484-500)]</td>
      <td>5.1000</td>
      <td>full14083</td>
      <td>78820.781</td>
      <td>PF13802</td>
      <td>159-219</td>
      <td>61.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498050</th>
      <td>A</td>
      <td>x823</td>
      <td>1xsi</td>
      <td>[(623-630)(484-500)]</td>
      <td>5.1000</td>
      <td>full14083</td>
      <td>78820.781</td>
      <td>PF01055</td>
      <td>240-669</td>
      <td>430.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>11498051</th>
      <td>A</td>
      <td>x063</td>
      <td>3bk1</td>
      <td>[(480-487)(515-531)]</td>
      <td>5.1000</td>
      <td>full11317</td>
      <td>40752.629</td>
      <td>PF07521</td>
      <td>386-439</td>
      <td>54.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498052</th>
      <td>A</td>
      <td>x063</td>
      <td>3bk1</td>
      <td>[(480-487)(515-531)]</td>
      <td>5.1000</td>
      <td>full11317</td>
      <td>40752.629</td>
      <td>PF12706</td>
      <td>67-232</td>
      <td>166.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498053</th>
      <td>A</td>
      <td>x114</td>
      <td>3cvj</td>
      <td>[(156-163)(92-108)]</td>
      <td>5.1000</td>
      <td>full9605</td>
      <td>13344.507</td>
      <td>PF13580</td>
      <td>7-142</td>
      <td>136.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498054</th>
      <td>A</td>
      <td>x664</td>
      <td>2rad</td>
      <td>[(254-261)(72-88)]</td>
      <td>5.1000</td>
      <td>full13904</td>
      <td>54156.688</td>
      <td>PF05139</td>
      <td>102-436</td>
      <td>335.0</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11498055</th>
      <td>A</td>
      <td>x249</td>
      <td>1g6i</td>
      <td>[(397-404)(308-324)]</td>
      <td>5.1000</td>
      <td>full11587</td>
      <td>68211.195</td>
      <td>PF01532</td>
      <td>45-532</td>
      <td>488.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>11498056</th>
      <td>A</td>
      <td>x429</td>
      <td>3kt9</td>
      <td>[(54-61)(78-94)]</td>
      <td>5.1000</td>
      <td>full4520</td>
      <td>19647.674</td>
      <td>PF00498</td>
      <td>39-104</td>
      <td>66.0</td>
      <td>False</td>
      <td>True</td>
    </tr>
  </tbody>
</table>
<p>11498057 rows × 12 columns</p>
</div>




```python
def data_plot( df, ax, maxim, rmsd_lim, linestyle, color ):
    allvalues = df.sort("rmsd").groupby(["pdb","chain"]).head(1)[["rmsd"]].values
    raw, y, x = rstoolbox.analysis.cumulative(allvalues, max_count=maxim, upper_limit=rmsd_lim)
    ax.plot(x, y, color=sns.color_palette()[color], lw=4, linestyle=linestyle)
    return raw, x, y

def first_marker_plot( raw, x, ax, top, shape, color ):
    idx = (np.abs(np.array(raw)-top)).argmin()
    ax.axvline(x=x[idx], ymin=0, ymax=0.2, c=sns.color_palette()[color], linewidth=4, linestyle="dashed", zorder=10)
    ax.plot([x[idx]], [0.2], shape, c=sns.color_palette()[color], markersize=12)
    
    
def plot( df, selfname, maxim, title, ax, known=None, knames=None, kcolor=None, annotate=False, top=10, rmsd_lim=5, max_length=180 ):
    # 1. Filter self-hits and error assignments
    df = df[(df["pdb"] != "eeee") & (df["pdb"] != selfname)]

    # 2. Make cumulative plot for best-hit/pdb-chain for all dataset
    raw, x, y = data_plot(df, ax, maxim, rmsd_lim, "solid", 0)
    if annotate:
        first_marker_plot( raw, x, ax, top, 'o', 0)
        
    # 3. Make cumulative plot for best-hit/pdb-chain for proteins of size "max_length" or smaller
    sdf = df[(df["length"] <= max_length) & (df["inrange"] == True)]
    raw, x, y = data_plot(sdf, ax, maxim, rmsd_lim, "dashed", 0)
    if annotate:
        first_marker_plot( raw, x, ax, top, 's', 0)

    # 4. Add clashes data
    if "ddg" in df:
        # 4.1. Make cumulative plot for best-hit/pdb-chain for all non-clashing dataset
        raw, x, y = data_plot(df[(df["ddg"] <= 0) ], ax, maxim, rmsd_lim, "solid", 1)
        ax.fill_between(x, 0, y, color=sns.color_palette()[1], alpha=0.3)
        if annotate:
            first_marker_plot( raw, x, ax, top, 'o', 1)
        
        # 4.2. Make cumulative plot for best-hit/pdb-chain for all non-clashing dataset for proteins of size "max_length" or smaller
        sdf = df[(df["ddg"] <= 0) & (df["length"] <= max_length) & (df["inrange"] == True)]
        raw, x, y = data_plot(sdf, ax, maxim, rmsd_lim, "dashed", 1)
        if annotate:
            first_marker_plot( raw, x, ax, top, 's', 1)

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
```


```python
fig  = plt.figure(figsize=(30, 20))
grid = (2, 3)
ax11 = plt.subplot2grid(grid, (0, 0))
ax12 = plt.subplot2grid(grid, (0, 1))
ax13 = plt.subplot2grid(grid, (0, 2))
ax21 = plt.subplot2grid(grid, (1, 0))
ax22 = plt.subplot2grid(grid, (1, 1))
ax23 = plt.subplot2grid(grid, (1, 2))
if isinstance(dfs["2fx7"], pd.DataFrame):
    plot( dfs["2fx7"], "2fx7", total_master_list, "2FX7: HIV 4e10 (H)",    ax11, [0.89], ["1Z6N"], kcolor=["black"], annotate=True, max_length=180 )
if isinstance(dfs["3ixt"], pd.DataFrame):
    plot( dfs["3ixt"], "3ixt", total_master_list, "3IXT: RSV mota (HLH)",  ax12, [2.31], ["3LHP"], kcolor=["black"], annotate=True, max_length=180 )
if isinstance(dfs["5tpn"], pd.DataFrame):
    plot( dfs["5tpn"], "5tpn", total_master_list, "5TPN: RSV RSV90 (HLE)", ax13, annotate=True, max_length=180 )
if isinstance(dfs["3o41"], pd.DataFrame):
    plot( dfs["3o41"], "3o41", total_master_list, "3O41: RSV 101F (E)",    ax21, annotate=True, max_length=180 )
if isinstance(dfs["3vtt"], pd.DataFrame):
    plot( dfs["3vtt"], "3vtt", total_master_list, "3VTT: DEN3 ED3 (ELE)",  ax22, [0.68], ["3WEI"], kcolor=["black"], annotate=True, max_length=180 )
if isinstance(dfs["4jhw"], pd.DataFrame):
    plot( dfs["4jhw"], "4jhw", total_master_list, "4JHW: RSV D25 (HxL)",   ax23, annotate=True, max_length=180 )

plt.tight_layout()
plt.savefig("images/motif_recovery.svg")
plt.savefig("images/motif_recovery.png")
sns.plt.show()
```


![png](README_files/README_9_0.png)



```python
fig  = plt.figure(figsize=(30, 20))
grid = (2, 3)
ax11 = plt.subplot2grid(grid, (0, 0))
ax12 = plt.subplot2grid(grid, (0, 1))
ax13 = plt.subplot2grid(grid, (0, 2))
ax21 = plt.subplot2grid(grid, (1, 0))
ax22 = plt.subplot2grid(grid, (1, 1))
ax23 = plt.subplot2grid(grid, (1, 2))
plot( dfs["2fx7"], "2fx7", total_master_list, "2FX7: HIV 4e10 (H)",    ax11, [0.89], ["1Z6N"], kcolor=["black"], max_length=180 )
plot( dfs["3ixt"], "3ixt", total_master_list, "3IXT: RSV mota (HLH)",  ax12, [2.31, 3.25, 3.6, 3.7], ["3LHP", "1TIL", "1NVI", "1KX8"], kcolor=["black", "red", "orange", "blue"], max_length=180 )
plot( dfs["5tpn"], "5tpn", total_master_list, "5TPN: RSV RSV90 (HLE)", ax13, max_length=180 )
plot( dfs["3o41"], "3o41", total_master_list, "3O41: RSV 101F (E)",    ax21, max_length=180 )
plot( dfs["3vtt"], "3vtt", total_master_list, "3VTT: DEN3 ED3 (ELE)",  ax22, [0.68], ["3WEI"], kcolor=["black"], max_length=180 )
plot( dfs["4jhw"], "4jhw", total_master_list, "4JHW: RSV D25 (HxL)",   ax23, max_length=180 )
plt.tight_layout()
plt.savefig("images/motif_recovery_with_designs.svg")
plt.savefig("images/motif_recovery_with_designs.png")
sns.plt.show()
```


![png](README_files/README_10_0.png)



```python

```
