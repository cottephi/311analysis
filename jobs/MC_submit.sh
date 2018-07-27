#! /bin/bash
outdir=/eos/user/p/pcotte/311analysis/jobs
let nfiles=$(ls /eos/experiment/wa105/offline/LArSoft/MC/MC6/ROOT/g4detsim/ | wc -l)
let maxfile=501
let step=20
let i=0
finish=false
while [ "$finish" = false ] ; do
  let file0=$i*$maxfile
  if [ $file0 -ge $nfiles ] ; then
    break
  fi
  let file1=$i+1
  let file1=$file1*$maxfile
  let file1=$file1-1
  if [ $file1 -ge $nfiles ] ; then
    let file1=$nfiles-1
    finish=true
  fi
  bsub -q 8nh -R "pool>30000" -J MC_${file0}_$file1 $outdir/job.sh "MC_run.sh $file0 $file1 $1 $2"
  echo $file0 $file1
  let i=$i+1
done

