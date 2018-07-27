#! /bin/bash 
direc="/eos/user/p/pcotte/311data/slow_control/all"
file="$direc/0.txt"
if [ ! -d $direc ] ; then
  echo "can not find folder"
  exit 1
fi
if [ ! -f $file ] ; then
  echo "can not find folder first file"
  exit 1
fi
first_line=$(head -n 1 $file)
nfiles=$(ls $direc | wc -l)
let dfiles=20
let I=0
let i=0
first=true
while [ $i -lt $nfiles ] ; do
  if [ ! -d $direc/../all_$I ] ; then
    echo "Creating directory $direc/../all_$I.."
    mkdir $direc/../all_$I
    echo "...done"
  fi
  empty="true"
  if [ -f $direc/../all_$I/all_$I.dat ] ; then
    echo "  Found previous all_$I.dat"
    if [ -s $direc/../all_$I/all_$I.dat ] ; then
      echo "    all_$I.dat is not empty. Ignoring..."
      empty="false"
    fi
    if [ "$empty" == "true" ] ; then
      echo "  Removing empty all_$I.dat files..."
      rm $direc/../all_$I/all_$I.dat
      echo "  ...done"
    fi
  fi
  if [ "$empty" == "false" ] ; then
      let i=$i+$dfiles
      let I=$I+1
    continue
  fi
  echo "Writing header line in all_$I.dat..."
  echo $first_line > $direc/../all_$I/all_$I.dat
  echo "...done"
  let j=0
  echo "Copying file to all_$I..."
  while [ $j -lt $dfiles ] ; do
    let ifile=$i+$j
    if [ -f $direc/$ifile.txt ] ; then
      cp $direc/$ifile.txt $direc/../all_$I/
    fi
    let j=$j+1
  done
  echo "...done"
  let i=$i+$dfiles
  echo "Filling all_$I.dat"
  tail -q -n +2 $direc/../all_$I/*.txt >> $direc/../all_$I/all_$I.dat
  echo "...done"
  echo "Removing i.txt files..."
  rm $direc/../all_$I/*.txt
  echo "...done"
  let I=$I+1
done
echo "All done"

