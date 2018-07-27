#! /bin/bash
if [ $# -ne 5 ] && [ $# -ne 3 ] ; then 
 echo "ERROR: need 5 arguments if using human time format: start and end moments, both must have date and time (like 2017/06/28 17:00:00 2017/06/29 17:00:00), and the run number. If using unix time, need  arguments."
 exit 1
fi

run=""
path="/eos/user/p/pcotte/311analysis/slow_control"
path_data="/eos/user/p/pcotte/311data/slow_control"
tstart=""
tend=""
if [ $# -eq 3 ] ; then 
  tstart=$(date "+%Y/%m/%d %H:%M:%S" -d @$1)
  tend=$(date "+%Y/%m/%d %H:%M:%S" -d @$2)
  run=$3
else
  tstart="$1 $2"
  tend="$3 $4"
  run=$5
fi
if ! date -d "$tstart" > /dev/null 2>&1 ; then
 echo "ERROR: invalide start time format"
 exit 1
fi
if ! date -d "$tend" > /dev/null 2>&1 ; then
 echo "ERROR: invalide end time format"
 exit 1
fi

$path/exp Vavu-Debi ssh wa105db@lxplus.cern.ch "~/tools/pcotte_getdata/Pcotte_getdata.sh $tstart $tend ~/tools/pcotte_getdata/LEMs_pressure_and_temperature_parameters.dat chien 1> out.txt 2> err.txt"
echo " "
$path/exp Vavu-Debi scp wa105db@lxplus.cern.ch:~/out.txt $path/$run.txt
echo " "
if [ ! -f $path/$run.txt ] ; then
  echo "ERROR: $run.txt not found"
  exit 1
fi
$path/exp Vavu-Debi scp wa105db@lxplus.cern.ch:~/err.txt $path/${run}_err.txt
echo " "
if [ ! -f $path/${run}_err.txt ] ; then
  echo "ERROR: ${run}_err.txt not found"
  exit 1
fi
$path/exp Vavu-Debi ssh wa105db@lxplus.cern.ch 'rm ~/out.txt'
echo " "
$path/exp Vavu-Debi ssh wa105db@lxplus.cern.ch 'rm ~/err.txt'
echo " "

if [ -s $path/err.txt ] ; then
  echo "Error occured on wa105db, check ${run}_err.txt for details"
  if [ -f $path/$run.txt ] ; then
    rm $path/$run.txt
  fi
  exit 1
else
  if [ -d $path_data/$run ] ; then
    rm -r $path_data/$run
  fi
  mkdir $path_data/$run
  mv $path/$run.txt $path_data/$run/
  rm $path/${run}_err.txt
fi

exit 0
