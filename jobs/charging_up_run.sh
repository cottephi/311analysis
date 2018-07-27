#! /bin/bash
if [ $# -eq 1 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/charging_up.cc('$1')' #> CHAR_$1.txt
elif [ $# -eq 2 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/charging_up.cc('$1',"'$2'")' #> CHAR_$1_$2.txt
elif [ $# -eq 3 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/charging_up.cc('$1',"'$2'","'$3'")' #> CHAR_$1_$2_$3.txt
elif [ $# -eq 4 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/charging_up.cc('$1',"'$2'","'$3'","'$4'")' #> CHAR_$1_$2.txt
elif [ $# -eq 5 ] ; then
  root -l -b -q /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/gain/charging_up.cc('$1',"'$2'","'$3'","'$4'","'$5'")' #> CHAR_$1_$2.txt
else
  root -l -b -q /eos/user/p/pcotte/311analysis/loadLib.cc /eos/user/p/pcotte/311analysis/gain/charging_up.cc #>CHAR.txt
fi
