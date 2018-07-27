#! /bin/bash
if [ $# -eq 1 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/charging_up_all.cc('$1')' #> CHAR_ALL_$1.txt
elif [ $# -eq 2 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/charging_up_all.cc('$1',"'$2'")' #> CHAR_ALL_$1_$2.txt
elif [ $# -eq 3 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/charging_up_all.cc('$1',"'$2'","'$3'")' #> CHAR_ALL_$1_$2_$3.txt
else
  root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc /eos/user/p/pcotte/311analysis/gain/charging_up_all.cc #>CHAR_ALL.txt
fi
