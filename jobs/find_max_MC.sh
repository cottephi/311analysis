#! /bin/bash

rm ./max*
max=""
let i=1
while [ $i -le 5 ] ; do
  for mydir in ./LSFJOB* ; do
    tail -n 7 $mydir/STDOUT | head -n 1 | cut -d " " -f $i >> column_$i.txt
  done
  max=${max}$(sort -n -k1 column_$i.txt | tail -n 1)" "
  let i=i+1
done
rm ./column_*
echo $max > max.txt

