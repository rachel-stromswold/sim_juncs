#!/bin/bash

prefix=$(timedatectl | awk '$1 == "Local" {tmp=$4"_"$5;gsub(/[-,:]/, "_", tmp);print "figures_"tmp}')

echo "copying files to $prefix"

#copy all of the plots
cp -r /scratch/sstromsw/junc_simuls/figures results/$prefix

#move the log files into the archive
for fname in *.out; do
    mv "$fname" "results/$prefix/$fname"
done

#copy whatever parameters were used for the simulation
cp params.conf "$prefix/params.conf"
cp junc_template.eps "$prefix/junc_template.eps"

#make the archive
tar -cvf results/"$prefix".tar.gz results/$prefix
