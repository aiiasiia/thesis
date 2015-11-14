#!/bin/bash
while true
do
 	if git push -q --repo=https://aiiasiia:ghh3lln0@github.com/aiiasiia/thesis
 	then
 		break
 	else
 		echo "retrying"
 	fi
done