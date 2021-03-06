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
for l in 9 11 13 15 # (l,d): (9,2) (11,3) (13,4) (15,5), and (17,6) for systems with > 10GB RAM
do
	((d++))
	# print headers in EMS_GT result files
	
	# start r runs
	for((i=1; i <= r; i++)) do
		# generate a unique (l,d) dataset for this run
		# java -cp ../bin DatasetGenerator  $l $d $i

		# test all programs on this dataset with all values of k	
		for k in 3 4 5 6 7 8
		do
			echo "l,d,run,time(s),time(min),memuse,memuse after GC,motif,motifs found" > ../results/E32-$l,$d,k=$k
			java -cp ../bin EMS_GT_32 $l,$d,$i $k	>> ../results/E32-$l,$d,k=$k
			echo "k = $k : EMS_GT_32 finished $l, $d, $i"
		done
	done
done
l=17
d=6
# start r runs
for((i=1; i <= r; i++)) do
	# generate a unique (l,d) dataset for this run
	# java -cp ../bin DatasetGenerator  $l $d $i

	# test all programs on this dataset with all values of k	
	for k in 3 4 5 6
	do
		echo "l,d,run,time(s),time(min),memuse,memuse after GC,motif,motifs found" > ../results/E32-$l,$d,k=$k
		java -cp ../bin -Xmx14g EMS_GT_32 $l,$d,$i $k	>> ../results/E32-$l,$d,k=$k
		echo "k = $k : EMS_GT_32 finished $l, $d, $i"
	done
done
cd ..

git add --all .
git commit -m 'Server tests finished'
while [ true ]
do    
	if git push --repo=https://aiiasiia:ghh3lln0@github.com/aiiasiia/thesis
	then
		break
	else
		echo "retrying..."
	fi
done

