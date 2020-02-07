#!/bin/bash

for f in `ls *.csv`
do
	echo $f
	python kontaminacje.py $f
done
