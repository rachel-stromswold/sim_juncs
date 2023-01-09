#!/bin/bash

res_dir="/scratch/$(whoami)/results"
prefix=$(timedatectl | awk '$1 == "Local" {tmp=$4"_"$5;gsub(/[-,:]/, "_", tmp);print "reses_"tmp}')

echo "copying files to $prefix"

#copy all of the plots
cp -r /scratch/sstromsw/junc_simuls "$res_dir/$prefix"

#move the log files into the archive
for fname in *.out; do
    mv "$fname" "$res_dir/$prefix/$fname"
done

#copy whatever parameters were used for the simulation
cp run.sh "$res_dir/$prefix/run.sh"
cp params.conf "$res_dir/$prefix/params.conf"
cp junc_template.geom "$res_dir/$prefix/junc_template.geom"

#make the archive
tar -cvf $res_dir/"$prefix".tar.gz $res_dir/$prefix
