#! /bin/bash
outdir=/eos/user/p/pcotte/311analysis/jobs
path_wa105_311data="/eos/experiment/wa105/offline/LArSoft/Data/Parser/2018_Feb_05"

if [ $# -eq 2 ] ; then
  bsub -q 1nh -R "pool>30000" -J MPV_vs_stuff_$1_$2 $outdir/job.sh "MPV_vs_stuff_run.sh $1 $2"
else
  bsub -q 1nh -R "pool>30000" -J MPV_vs_stuff $outdir/job.sh
fi
