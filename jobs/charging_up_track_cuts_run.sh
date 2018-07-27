#! /bin/bash
root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/utils/charging_up_track_cuts.cc('$1',"'$2'","'$3'","'$4'","'$5'",'$6')' #> '/eos/user/p/pcotte/311analysis/jobs/CHAR_TC_'$1'_'$2'_'$3'.txt'
