#! /bin/bash
root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/slow_control/pressure.cc("'$1'")' #> press_all.txt
