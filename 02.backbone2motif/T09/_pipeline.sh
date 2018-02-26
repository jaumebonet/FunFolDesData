# Fixed Paths (for local executions)
export ROSETTA=/Volumes/MiniTwo/bin/Rosetta/main/source/bin/
export COMPILE=macosclangrelease
export TOOLS=/Volumes/MiniTwo/bin/Rosetta/tools/
export RSBIN=/Volumes/MiniTwo/bin/python/RosettaSilentToolbox/rstoolbox/bin/

# Setup query pdb.
QUERYPDB=3phxB

# Run baseline data on our target
function baseline() {
  cd baseline
  $ROSETTA/rosetta_scripts.$COMPILE -parser:protocol baseline.xml -s ../structures/$QUERYPDB.pdb -out:file:silent baseline.sc
  gzip baseline.sc
  gzip contacts_*.csv
  cd -
}

# Generate fragments
function fragments() {
  cd fragments/auto
  # $ROSETTA/rosetta_scripts.$COMPILE -parser:protocol frags.xml -s ../../structures/$QUERYPDB.pdb -score_only tmp.sc
  # rm tmp.sc
  # cd -
  # cd fragments/wauto
  # $ROSETTA/rosetta_scripts.$COMPILE -parser:protocol frags.xml -s ../../structures/$QUERYPDB.pdb -score_only tmp.sc
  # rm tmp.sc
  # cd -
  # cd fragments/picker
  # $TOOLS/fragment_tools/make_fragments.pl ../../structures/$QUERYPDB.fa
  cd ..
  ln -s auto/auto.200.3mers $QUERYPDB.auto.200.3mers
  if [ ! -f $QUERYPDB.auto.200.3mers.qual ]; then
    $ROSETTA/r_frag_quality.$COMPILE -in:file:native ../structures/$QUERYPDB.pdb -f $QUERYPDB.auto.200.3mers -out:qual $QUERYPDB.auto.200.3mers.qual
  fi
  ln -s auto/auto.200.9mers $QUERYPDB.auto.200.9mers
  if [ ! -f $QUERYPDB.auto.200.9mers.qual ]; then
    $ROSETTA/r_frag_quality.$COMPILE -in:file:native ../structures/$QUERYPDB.pdb -f $QUERYPDB.auto.200.9mers -out:qual $QUERYPDB.auto.200.9mers.qual
  fi
  ln -s picker/$QUERYPDB.200.3mers $QUERYPDB.picker.200.3mers
  if [ ! -f $QUERYPDB.picker.200.3mers.qual ]; then
    $ROSETTA/r_frag_quality.$COMPILE -in:file:native ../structures/$QUERYPDB.pdb -f $QUERYPDB.picker.200.3mers -out:qual $QUERYPDB.picker.200.3mers.qual
  fi
  ln -s picker/$QUERYPDB.200.9mers $QUERYPDB.picker.200.9mers
  if [ ! -f $QUERYPDB.picker.200.9mers.qual ]; then
    $ROSETTA/r_frag_quality.$COMPILE -in:file:native ../structures/$QUERYPDB.pdb -f $QUERYPDB.picker.200.9mers -out:qual $QUERYPDB.picker.200.9mers.qual
  fi
  ln -s wauto/wauto.200.3mers $QUERYPDB.wauto.200.3mers
  if [ ! -f $QUERYPDB.wauto.200.3mers.qual ]; then
    $ROSETTA/r_frag_quality.$COMPILE -in:file:native ../structures/$QUERYPDB.pdb -f $QUERYPDB.wauto.200.3mers -out:qual $QUERYPDB.wauto.200.3mers.qual
  fi
  ln -s wauto/wauto.200.9mers $QUERYPDB.wauto.200.9mers
  if [ ! -f $QUERYPDB.wauto.200.9mers.qual ]; then
    $ROSETTA/r_frag_quality.$COMPILE -in:file:native ../structures/$QUERYPDB.pdb -f $QUERYPDB.wauto.200.9mers -out:qual $QUERYPDB.wauto.200.9mers.qual
  fi
  cd ..
}


# Run experiments
# In cscs (to write)

# Minimize results
function minimize() {
  cd output
  python $RSBIN/minisilent.py -in:files "abinitio_auto/*[0-9].gz"    -out:file abinitio_auto.minisilent.gz
  python $RSBIN/minisilent.py -in:files "abinitio_picker/*[0-9].gz"  -out:file abinitio_picker.minisilent.gz
  python $RSBIN/minisilent.py -in:files "abinitio_wauto/*[0-9].gz"   -out:file abinitio_wauto.minisilent.gz
  python $RSBIN/minisilent.py -in:files "nubinitio_auto/*[0-9].gz"   -out:file nubinitio_auto.minisilent.gz
  python $RSBIN/minisilent.py -in:files "nubinitio_picker/*[0-9].gz" -out:file nubinitio_picker.minisilent.gz
  python $RSBIN/minisilent.py -in:files "nubinitio_wauto/*[0-9].gz"  -out:file nubinitio_wauto.minisilent.gz
  cd -
}

# Execute
baseline
fragments
minimize

