#! /bin/bash
if [ $# -eq 1 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/gain_ana_dQds.cc("'$1'")' #> /eos/user/p/pcotte/311analysis/jobs/GAIN_DQDS_$1.txt
elif [ $# -eq 2 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/gain_ana_dQds.cc("'$1'","'$2'")' #> /eos/user/p/pcotte/311analysis/jobs/GAIN_DQDS_$1_$2.txt
else
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc /eos/user/p/pcotte/311analysis/gain/gain_ana_dQds.cc #> /eos/user/p/pcotte/311analysis/jobs/GAIN_DQDS.txt
fi
