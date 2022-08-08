#!/bin/bash
#SBATCH -p femto --time=02:00:00 --mem=244G
#SBATCH -a 1-4

run_local="f"
run_simuls="t"
make_movie="f"
SCALE=8.0

#set the output directory name
pname="."

for var in "$@"; do
    if [ $var == "-l" ]; then
        run_local="t"
    elif [ $var == "-m" ]; then
        make_movie="t"
    elif [ $var == "--dry" ]; then
        run_simuls="f"
    fi
done

#only load modules if we're running this on the cluster
if [ $run_local == "f" ]; then
    #pname="/local_scratch/$SLURM_JOB_ID"
    pname=/scratch/sstromsw/junc_simuls
    echo "using $pname"

    module load python3/3.8.8
    module load meep/b3
    module load gcc/11.2.0/b1       
    module load hdf5/1.12.1/b1      
    module load openmpi/1.10.7/b1
    module load eigen/3.3.7/b1
fi

if [ "$#" -ge 2 ]; then
    pname=$1
fi

echo "output name: $pname        run local? $run_local"

#clean up the working directory of old pngs
rm -rf $pname/figures
rm -f $pname/errors.txt

#the widths used for each job in the array. These are similar to those used in Boolakee et al.
widths=("0.001" "0.002" "0.005" "0.01")
thickness="0.03"
resolution="12.0"
h5dir="$pname/test_$SLURM_ARRAY_TASK_ID"
mkdir $h5dir
if [ $run_simuls == "t" ]; then
    rm -f $h5dir/*
    #valgrind --leak-check=full --track-origins=yes ./sim_geom --out-dir $h5dir --grid-res ${resolutions[$((SLURM_ARRAY_TASK_ID-1))]}
    python set_dims.py ${widths[$((SLURM_ARRAY_TASK_ID-1))]} $thickness > $h5dir/junc.eps
    ./sim_geom --out-dir $h5dir --grid-res $resolution --geom-file $h5dir/junc.eps
    eps_file=$(ls | grep eps)
    echo "eps_fname = $eps_file"
    cp $eps_file "$h5dir/$eps_file"
fi
#make the test plots in check_enes.py
mkdir "$pname"/figures
if [ $make_movie == "t" ] && [ $SLURM_ARRAY_TASK_ID == "1" ]; then
	python check_enes.py --prefix $h5dir --grid-res $resolution --gap-width ${widths[$((SLURM_ARRAY_TASK_ID-1))]} --gap-thick $thickness --movie
    rm -f "$pname"/figures/out.mp4
	#combine the images to make an animation
	ffmpeg -r 30 -f image2 -s 1920x1080 -i $h5dir/im_%d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p "$pname"/figures/out.mp4
fi

python check_enes.py --prefix $h5dir --grid-res $resolution --gap-width ${widths[$((SLURM_ARRAY_TASK_ID-1))]} --gap-thick $thickness | tee $h5dir/tmp.log
awk '$1 == "square_errors:" { $1=""; print $0; }' $h5dir/tmp.log >> $pname/errors.txt
python plot_convergence.py --in-prefix $pname --out-prefix $pname/figures

#move the plots into a folder where we can view them
cp "$h5dir"/space_plot.pdf "$pname"/figures/space_plot_${widths[$((SLURM_ARRAY_TASK_ID-1))]}.pdf
cp "$h5dir"/cross_plot.pdf "$pname"/figures/cross_plot_${widths[$((SLURM_ARRAY_TASK_ID-1))]}.pdf
