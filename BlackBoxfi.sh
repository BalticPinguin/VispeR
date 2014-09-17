#!/bin/bash

Threshold=6
if ! [ -x "$1" ]  #if the first argument is not an executable file
then
   echo $1 'is not executable or does not exist.'
fi
float='^-?[0-9]+([.][0-9]+)?$'
args=("${@: 2}")
echo -ne "args:"
echo ${args[*]}
result=$(./$1 ${*: 2})
echo $result
let len=$#-2
#test which parameters are numbers (and can be changed numerically) and which not:
for j in `seq 0 1 "$len"` ; do
   if  [[ ${args[$j]} =~ $float ]] #if args is a number:
   then
      result=$(./$1 ${args[*]}) #this might be unneccesary here!
      delta=`echo "${args[$j]} / 100" |bc -l`
      args[$j]=`echo "${args[$j]} - $delta" | bc -l ` #change args
      for i in `seq 1 1 "$Threshold"`; do
	 result1=$(./$1 ${args[*]}) 
	 echo $result1 $result
	 if [[ `echo "$result1 < $result"| bc -l` -eq 1 ]] ; then
	    args[$j]=`echo "${args[$j]} - $delta" | bc -l` #change args
	 elif [[ `echo "$result1 > $result"| bc -l` -eq 1 ]] ; then
	    delta=-$delta
	    args[$j]=`echo "${args[$j]} - $delta" | bc -l` #change args
	 else 
	    echo "foo"
	    break
	 fi
	 result=$result1
      done
   fi
   echo ""
done

echo -ne "new parameters are (in same order as input):"
echo ${args[*]}
