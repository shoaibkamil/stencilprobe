#!/bin/bash
# finds best block for 512^3 problem

for tj in 16 32 64 128 256 512 
do
	for ti in 16 32 64 128 256 512
	do
		echo "#define SIZE 512" > run.h
		echo "#define NUM_TRIALS 5" >> run.h
		echo "#define TI $ti" >> run.h
		echo "#define TJ $tj" >> run.h
		make blocked_probe
		./probe > out-$ti\-$tj\.out
		echo "done with $ti $tj"
	done
done

for x in 16 32 64 128 256 512
do 
	for y in 16 32 64 128 256 512
	do 
		echo -n $x x $y, 
		cat out-$x\-$y\.out | grep elapsed | awk '{print $4}' | sed 's/time://' | perl min.pl
	done
done
