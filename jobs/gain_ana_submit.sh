outdir=/eos/user/p/pcotte/311analysis/jobs
bsub -q 8nh -R "pool>30000" -J gain_ana_$1 $outdir/job.sh "gain_ana_run.sh $1"
