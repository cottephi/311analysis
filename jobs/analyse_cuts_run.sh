#! /bin/bash
if [ $# -eq 1 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/cuts_analysis/cuts_analysis.cc('$1')' #> CA_$1.txt
elif [ $# -eq 2 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/cuts_analysis/cuts_analysis.cc('$1',"'$2'")' #> CA_$1_$2.txt
elif [ $# -eq 3 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/cuts_analysis/cuts_analysis.cc('$1',"'$2'","'$3'")' #> CA_$1_$2_$3.txt
elif [ $# -eq 4 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/cuts_analysis/cuts_analysis.cc('$1',"'$2'","'$3'","'$4'")' #> CA_$1_$2_$3.txt
elif [ $# -eq 5 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/cuts_analysis/cuts_analysis.cc('$1',"'$2'","'$3'","'$4'","'$5'")' #> CA_$1_$2_$3.txt
elif [ $# -eq 6 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/cuts_analysis/cuts_analysis.cc('$1',"'$2'","'$3'","'$4'","'$5'",'$6')' #> CA_$1_$2_$3.txt
else
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc /eos/user/p/pcotte/311analysis/cuts_analysis/cuts_analysis.cc #> CA.txt
fi
