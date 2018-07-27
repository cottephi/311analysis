#! /bin/bash
if [ $# -eq 4 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/utils/track_cuts_MC.cc('$1','$2',"'$3'","'$4'")' #> MC_$1_$2.txt
else
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc /eos/user/p/pcotte/311analysis/utils/track_cuts_MC.cc #> MC.txt
fi
