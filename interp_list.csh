#!/bin/csh

set list=`ls ${4}*0`
echo $list
echo 
set n=`expr ${2} - ${1} + 1`
echo $n
set m=`sseq $1 $2`
echo $m
echo $m
foreach name ($list)
  swampi -n ${n} -m "${m}" $TREEHOME/interp/interp ${3} ${name} 
  $TREEHOME/sdftofits/sph2fits.p4icc grid_${name} ${name}.fits
end

echo "Done."
