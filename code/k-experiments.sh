#!/bin/bash
 
# k - experiments.sh
# Sets up EMS-GT.
#  Generates (l,d) planted motif problem datasets.
#  Performs r experimental runs (testing k-values 3,4,...,7,8 per run) and formats results.

# @author Aia Sia
# @version 1.0 9/09/2015

r=20 # number of experimental runs
javac -d bin src/*.java # compile java source files

# cd src/PMS8; make; mv Debug/PMS8  ../../bin/PMS8;  # uncomment these lines to compile PMS8 and qPMS9 from scratch.
# cd ../qPMS9; make; mv NoMpi/qpms9 ../../bin/qPMS9; # see Readme.md files in src/PMS8 and src/qPMS9 for dependency info.

d=1
cd datasets
# (l,d): (9,2) (11,3) (13,4) (15,5) 
for l in 9 11 13 15 17 # only include (17,6) for systems with > 10GB RAM
do
	((d++))
	for k in 3 4 6 7 8
		# print headers in EMS_GT result files
		echo "l,d,run,time(s),time(min),memuse,memuse after GC,motif,motifs found" > ../results/EMS_GT-$l,$d
		echo "l,d,run,time(s),time(min),memuse,memuse after GC,motif,motifs found" > ../results/EMS_GT_32-$l,$d
		echo "l,d,run,time(s),time(min),memuse,memuse after GC,motif,motifs found" > ../results/EMS_GT_64-$l,$d
		# start r runs
		for((i=1; i <= r; i++)) do
			# generate a unique (l,d) dataset for this run
			java -cp ../bin DatasetGenerator  $l $d $i
			# test all programs on this dataset
			java -cp ../bin EMS_GT    $l,$d,$i $k	>> ../results/EMS_GT-$l,$d--k$k
			java -cp ../bin EMS_GT_32 $l,$d,$i $k	>> ../results/EMS_GT_32-$l,$d--k$k
			java -cp ../bin EMS_GT_64 $l,$d,$i $k	>> ../results/EMS_GT_64-$l,$d--k$k
		done
	done
done
cd ..
