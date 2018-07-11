minisilent.py -in:files abinitio_strfrags_ -out:file abinitio_strfrags.minsilent.gz
minisilent.py -in:files abinitio_seqfrags_ -out:file abinitio_seqfrags.minsilent.gz
tar zcvf abinitio_seqfrags.tgz abinitio_seqfrags_*
tar zcvf abinitio_strfrags.tgz abinitio_strfrags_*
rm abinitio_seqfrags_*
rm abinitio_strfrags_*

