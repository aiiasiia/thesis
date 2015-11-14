#!/bin/bash

r=1
javac -d bin src/*.java

d=1
cd datasets
for((i=1; i <= r; i++)) do
	java -cp ../bin DatasetGenerator  9 2 $i
	java -cp ../bin DatasetGenerator 11 3 $i

	for k in 3 4
	do
		java -cp ../bin EMS_GT_32  9,2,$i $k	>> ../results/E32-9,2,k=$k
		echo "k = $k : finished  9,2,$i"
		java -cp ../bin EMS_GT_32 11,3,$i $k	>> ../results/E32-11,3,k=$k
		echo "k = $k : finished 11,3,$i"
	done
done

cd ..
git add --all .
git commit -m 'Shell script worked!'
s='git status --porcelain'
while[ -n s ] do
	echo "git push to aiiasiia/thesis"
	git push --repo https://aiiasiia:ghh3lln0@github.com/aiiasiia/thesis > log
done
