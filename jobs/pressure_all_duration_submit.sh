#! /bin/bash
outdir=/eos/user/p/pcotte/311analysis/jobs
for direc in /eos/user/p/pcotte/311data/slow_control/all_* ; do
  bsub -q 8nh -R "pool>30000" -J pressure_all_duration $outdir/job.sh "pressure_all_duration_run.sh $(basename $direc)"
done
