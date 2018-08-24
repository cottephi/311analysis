#! /bin/bash
if [ $# -eq 0 ] ; then 
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc /eos/user/p/pcotte/311analysis/utils/track_cuts.cc
else
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/utils/track_cuts.cc('$1',"'$2'","'$3'","'$4'","'$5'",'$6','$7')' #> '/eos/user/p/pcotte/311analysis/jobs/TC_'$1'_'$2'_'$3'_'$4'_'$5'.txt'
fi
