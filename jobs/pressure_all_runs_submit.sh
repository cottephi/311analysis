#! /bin/bash
outdir=/eos/user/p/pcotte/311analysis/jobs

bsub -q 8nh -R "pool>30000" -J getdata_all_runs $outdir/job.sh "pressure_all_runs_run.sh"
