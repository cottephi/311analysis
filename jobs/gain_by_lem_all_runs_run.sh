#! /bin/bash
if [ $# -eq 1 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/gain_by_lem_all_runs.cc('$1')' #> /eos/user/p/pcotte/311analysis/jobs/GAIN_LEM_$1.txt
elif [ $# -eq 2 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/gain_by_lem_all_runs.cc('$1',"'$2'")' #> /eos/user/p/pcotte/311analysis/jobs/GAIN_LEM_$1_$2.txt
elif [ $# -eq 3 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/gain_by_lem_all_runs.cc('$1',"'$2'","'$3'")' #> /eos/user/p/pcotte/311analysis/jobs/GAIN_LEM_$1_$2_$3.txt
elif [ $# -eq 4 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/gain_by_lem_all_runs.cc('$1',"'$2'","'$3',"'$4'"")' #> /eos/user/p/pcotte/311analysis/jobs/GAIN_LEM_$1_$2_$3.txt
elif [ $# -eq 5 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/gain_by_lem_all_runs.cc('$1',"'$2'","'$3'","'$4'","'$5'")' #> /eos/user/p/pcotte/311analysis/jobs/GAIN_LEM_$1_$2_$3.txt
else
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc /eos/user/p/pcotte/311analysis/gain/gain_by_lem_all_runs.cc #> /eos/user/p/pcotte/311analysis/jobs/GAIN_LEM.txt
fi
