#!/bin/bash

prefix=$(timedatectl | awk '$1 == "Local" {tmp=$4"_"$5;gsub(/[-,:]/, "_", tmp);print "reses_"tmp}')

echo "copying files to $prefix"

#copy all of the plots
cp -r /scratch/sstromsw/junc_simuls results/$prefix

#move the log files into the archive
for fname in *.out; do
    mv "$fname" "results/$prefix/$fname"
done

#copy whatever parameters were used for the simulation
cp run.sh "results/$prefix/run.sh"
cp params.conf "results/$prefix/params.conf"
cp junc_template.geom "results/$prefix/junc_template.geom"

#make the archive
tar -cvf results/"$prefix".tar.gz results/$prefix
