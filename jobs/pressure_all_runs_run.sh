#! /bin/bash
path_slow_control=/eos/user/p/pcotte/311data/slow_control
path_analysis=/eos/user/p/pcotte/311analysis
for run in $path_slow_control/* ; do
  run=$(basename $run);
  if [ "$run" == "all" ] ; then
    continue; 
  fi;
  root -l -q -b $path_analysis/loadHeader.cc $path_analysis/loadLib.cc ''$path_analysis'/slow_control/pressure.cc("'$run'")'
done
