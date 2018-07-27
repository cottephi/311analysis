#! /bin/bash
#let start=1499792254
let start=1504630654
let end=1510331674
let oneday=86400
let start2=$start+$oneday
#let i=0
let i=56
while [ $start -lt $end ] ; do
  if [ $start2 -gt $end ] ; then
    let start2=$end
  fi
  if ! /eos/user/p/pcotte/311analysis/slow_control/Pcotte_getdata.sh $start $start2 $i ; then
    exit 1
  fi
  if [ ! -f /eos/user/p/pcotte/311data/slow_control/$i/$i.txt ] ; then
    echo "ERROR: $i.txt not found"
    exit 1
  fi
  mv /eos/user/p/pcotte/311data/slow_control/$i/$i.txt /eos/user/p/pcotte/311data/slow_control/all/$i.txt
  rm -r /eos/user/p/pcotte/311data/slow_control/$i/
  
  echo "$start $start2 $i" >> /eos/user/p/pcotte/311analysis/slow_control/log.txt
  
  let start=$start2
  let start2=$start2+$oneday
  let i=$i+1
done

exit 0
