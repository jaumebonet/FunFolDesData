# DATABASE

## Vall - Fragment Database

The Vall database is needed to create the fragments.  

Due to its size, it cannot be directly uploaded. Thus, once you have all the pieces run:  

```bash
cat vall.jul19.2011.gz-part-a* > vall.jul19.2011.gz
```

## PFAM HMM database

HMM are used in several assessments. **REQUIRES HMMER!**

The main HMM file for Pfam: ``Pfam-A.hmm.gz`` is over 5Gb, which makes it difficult to add
to the repo. It can be retrieved with:

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

Version ``31.0`` is the one used for the publication analysis as the latest release on
publication date.
