#! /bin/bash
for dir in /afs/cern.ch/user/p/pcotte/eos_link/311data/slow_control/* ; do 
#  if [ -e "$dir/pressure_and_temp" ] ; then 
#    if [ -e "$dir/pressure_and_temp/temperature_cryostat.root" ] ; then 
#      continue
#    fi
#  fi
  dir=$(basename $dir)
  if [ "$dir" == "0" ] || [[ $dir == "all"* ]] ; then 
    continue
  fi
  root -l -q -b ../loadHeader.cc ../loadLib.cc 'Pressure.cc("'$dir'")'
  echo $dir
done

