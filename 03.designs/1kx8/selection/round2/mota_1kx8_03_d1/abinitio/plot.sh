regplot_rosetta.py -in:file abinitio_strfrags.minsilent.gz -in:y score -in:x rms -out:silent -plot:color 1 -plot:title ${1}_str -out:file ${1}_str.png
regplot_rosetta.py -in:file abinitio_seqfrags.minsilent.gz -in:y score -in:x rms -out:silent -plot:color 0 -plot:title ${1}_seq -out:file ${1}_seq.png
