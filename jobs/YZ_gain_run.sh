#! /bin/bash
if [ $# -eq 1 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/YZ_gain.cc('$1')' #> YZ_$1.txt
elif [ $# -eq 2 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/YZ_gain.cc('$1',"'$2'")' #> YZ_$1_$2.txt
elif [ $# -eq 3 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/YZ_gain.cc('$1',"'$2'",'$3')' #> YZ_$1_$2_$3.txt
elif [ $# -eq 4 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/YZ_gain.cc('$1',"'$2'",'$3',"'$4'")' #> YZ_$1_$2_$3_$4.txt
elif [ $# -eq 5 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/YZ_gain.cc('$1',"'$2'","'$3'","'$4'","'$5'")' #> YZ_$1_$2_$3_$4.txt
else
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc /eos/user/p/pcotte/311analysis/gain/YZ_gain.cc #>YZ.txt
fi
