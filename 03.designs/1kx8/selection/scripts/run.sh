RSDIR=/Volumes/MiniTwo/bin/Rosetta/main/source/bin/
GITDIR=$(git rev-parse --show-toplevel)


${RSDIR}/rosetta_scripts.macosclangrelease -parser:protocol frags.xml -s ${1} -nooutput true -parser:script_vars fragsDB=${GITDIR}/resources/database/vall.jul19.2011.gz
