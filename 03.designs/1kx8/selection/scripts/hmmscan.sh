HMMDIR=/Volumes/MiniTwo/bin/external/hmmer-3.1b2-macosx-intel/binaries/
GITDIR=$(git rev-parse --show-toplevel)

${HMMDIR}/hmmscan ${GITDIR}/resources/database/Pfam-A.hmm ../${1}.fasta > ${1}.HMMmatch
