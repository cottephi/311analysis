#! /bin/bash
if [ $# -eq 1 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/purity/purity.cc('$1')' #> PUR_$1.txt
elif [ $# -eq 2 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/purity/purity.cc('$1',"'$2'")' #> PUR_$1_$2.txt
elif [ $# -eq 3 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/purity/purity.cc('$1',"'$2'","'$3'")' #> PUR_$1_$2_$3.txt
else
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc /eos/user/p/pcotte/311analysis/purity/purity.cc #> PUR.txt
fi
