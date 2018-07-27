#! /bin/bash
#first, go to ssh wa105db@lxplus.cern.ch (passwd Dazake51)

if [ $# -ne 0 ] ; then
  echo "WARNING: does not take arguments"
fi

for Dir in ./dir_* ; do
  dir_data=$(basename $Dir)
  if ! [ -d "$dir_data" ] ; then
    continue
  fi
  if ! ls $dir_data/Data_* 1> /dev/null 2>&1; then
    continue
  fi
  cp $dir_data/Data_* ./
done

rsync -avh -e ssh wa105db@lxplus.cern.ch:~/tools/pcotte_getdata/Data_* ../../311data/slow_control/

if ! ls ../../311data/slow_control/Data_* 1> /dev/null 2>&1; then
  echo "No Data file transfered"
  exit 1
fi

for data in ../../311data/slow_control/Data_* ; do
  if [ -d "$data" ] ; then
    echo "$data is a directory. Skipping it."
    continue
  fi
  direc=$(echo $data | sed -r 's/Data_//g')
  if [ -d "$direc" ] ; then
    rm -r $direc
 fi
  mkdir $direc
  mv $data $direc/
done
