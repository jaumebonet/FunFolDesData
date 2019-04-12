# Experimental Designs and Characterization

The folders in this section contain the scripts and data to mimic the grafting design experiments presented in the paper.

Each folder contains all the data related to that particular case. Inside it, there are three main folders:

* **design:** Contain all the data and script to make the design again (as far as Rosetta would allow it).

* **benchmark:** Contains a small benchmark testing the outcomes for different conditions of the design in the presence of the target antibody.

* **selection:** Contains the final predicted structures and obtained experimental data for all the designs that were selected and tested in experimentally.

| WARNING: Since the publication of this data, several features have changed in Rosetta. First, the default score function is not talaris2014 anymore. To execute the RosettaScripts as they are, use the flag `-restore_talaris_behavior`. Second, Rosetta now ignores empty occupancies by default. This will affect the 1KX8 design, for which no density is present for the CA of residue 100A. To solve this issue, add the flag `-ignore_zero_occupancy 0`|
| --- |
