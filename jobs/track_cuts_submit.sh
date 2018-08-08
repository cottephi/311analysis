#! /bin/bash
outdir=/eos/user/p/pcotte/311analysis/jobs
#runs='{1189,1190,1191,1194,1195,1196,1197,1182,1183,1187,1188,840}'
runs='{}'
cut_type='Ds'
version='July'
method_dQ='sum'
method_ds='local'
save_tracks='true'
need_help=false
while  [ -n "$1" ] ; do
  case "$1" in
    -r) runs=$2
        shift ;;
    -c) cut_type=$2 
        shift ;;
    -v) version=$2
        shift ;;
    -q) method_dQ=$2
        shift ;;
    -s) method_ds=$2
        shift ;;
    -t) save_tracks=$2
        shift ;;
    -h) need_help=true
        shift;;
    *) echo "Option $1 not recognized" ;;
  esac
  shift
done

if [ $need_help = true ] ; then
  echo "-r run -c cut -v version -q method_dQ -s method_ds -t save_tracks -h help"
  exit 1
fi
path_wa105_311data="/eos/experiment/wa105/offline/LArSoft/Data/Reco/2018_Feb_05/ROOT/recofast/"

if [ "$runs" = "{}" ] ; then
  let i=0
  array=()
  for run in $path_wa105_311data/* ; do
    if [ ! -d $run ]  ; then
      continue
    fi
    if [ "$(ls -A $run)" ] ; then
      run=$(basename $run)
      array+=($run)
      if [ $i -eq 4 ] ; then
        runs='{'
        for item in ${array[*]}  ; do
          runs+=$item
          runs+=','
        done
        array=()
        runs=${runs%?}
        runs+='}'
        bsub -q 8nh -R "pool>30000" -J track_cuts_${runs}_${cut_type}_${version}_${method_dQ}_${method_ds}_${save_tracks} $outdir/job.sh "track_cuts_run.sh $runs $cut_type $version $method_dQ $method_ds $save_tracks 'true'"
        let i=0
      else
        let i=i+1
      fi #if i=lim
    fi #if ls
  done # for run
  if [ ${#array[*]} -ne 0 ] ; then
    runs='{'
    for item in ${array[*]}  ; do
      runs+=$item
      runs+=','
    done
    array=()
    runs=${runs%?}
    runs+='}'
    bsub -q 8nh -R "pool>30000" -J track_cuts_${runs}_${cut_type}_${version}_${method_dQ}_${method_ds}_${save_tracks} $outdir/job.sh "track_cuts_run.sh $runs $cut_type $version $method_dQ $method_ds $save_tracks 'true'"
  fi #if array not empty
else
  bsub -q 8nh -R "pool>30000" -J track_cuts_${runs}_${cut_type}_${version}_${method_dQ}_${method_ds}_${save_tracks} $outdir/job.sh "track_cuts_run.sh $runs $cut_type $version $method_dQ $method_ds $save_tracks 'true'"
fi
