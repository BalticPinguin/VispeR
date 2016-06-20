#!/bin/bash

# this might be necessary to re-set.
alias VispeR.py='~/bin/smallscript/VispeR.py'

function test () {
system=$1
cd $system
for files in $@
do
   if [[ $files == $system ]]
   then
      continue
   fi
   (VispeR.py $files\.inp
   diff $files\.log results/$files\.log > $files\.diff
   isdifferent=$(wc $files\.diff | awk '{print $1}')
   if [[ $isdifferent -ge 1 ]]
   then
      echo "some error in system $files ."
   fi
   )
   num_jobs=$( jobs | wc | awk '{print $1}')
   while [ $num_jobs -ge 4 ]
   do
      sleep 30
      num_jobs=$( jobs | wc | awk '{print $1}')
   done
done
cd ../
}

system=Acrolein
echo "test Acrolein:"
(test $system changed gradient shift)

system=Anisole 
echo "test Anisole:"
(test $system changed shift smsc)

system=Artificial
echo "test artifical system:"
(test $system MultiMode oneMode)

system=Benzene 
echo "test Benzene:"
(test $system DR_exc FC_exc FC_fluor)

system=Phenol 
echo "test Phenole:"
(test $system DR_abs FC_abs FC_abs2)

wait
