#! /bin/bash
root -l -b -q /eos/user/p/pcotte/311analysis/loadHeader.cc /eos/user/p/pcotte/311analysis/loadLib.cc '/eos/user/p/pcotte/311analysis/utils/YZ_track_cuts.cc('$1',"'$2'","'$3'","'$4'","'$5'",'$6')' #> '/eos/user/p/pcotte/311analysis/jobs/YZ_TC_'$1'_'$2'_'$3'_'$4'_'$5'.txt'
