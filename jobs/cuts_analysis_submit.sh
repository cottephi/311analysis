#! /bin/bash
outdir=/eos/user/p/pcotte/311analysis/jobs

if [ $# -eq 2 ] ; then
  bsub -q 1nh -R "pool>30000" -J purity_$1_$2 $outdir/job.sh "cuts_analysis_run.sh $1 $2"
else
  bsub -q 1nh -R "pool>30000" -J purity $outdir/job.sh "cuts_analysis_run.sh"
fi
