#!/bin/bash
#SBATCH -p femto --time=02:00:00 --mem=244G
#SBATCH -a 1-4

run_local="f"
run_simuls="t"
make_movie="t"
SCALE=8.0

#set the output directory name
pname="."
oname="."

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
    pname="/local_scratch/$SLURM_JOB_ID"
    #pname="/scratch/sstromsw/junc_simuls"
    oname="/scratch/sstromsw/junc_simuls"

    module load python3/3.8.8
    module load meep/b3
    module load gcc/11.2.0/b1
    module load hdf5/1.12.1/b1
    module load openmpi/1.10.7/b1
    module load eigen/3.3.7/b1
fi

echo "working directory name: $pname, output name: $oname, run local? $run_local"

#clean up the working directory of old pngs
rm -f $oname/errors.txt
mkdir $oname/figures

#the widths used for each job in the array. These are similar to those used in Schiffrin et al. For a SiO2 junction.
widths=("0.02" "0.05" "0.1" "0.5")
thickness="0.2"
resolution="12.0"
h5dir="$pname/test_$SLURM_ARRAY_TASK_ID"
mkdir $h5dir
if [ $run_simuls == "t" ]; then
    rm -f $h5dir/*
    #valgrind --leak-check=full --track-origins=yes ./sim_geom --out-dir $h5dir --grid-res ${resolutions[$((SLURM_ARRAY_TASK_ID-1))]}
    python set_dims.py ${widths[$((SLURM_ARRAY_TASK_ID-1))]} $thickness > $h5dir/junc.geom
    ./sim_geom --out-dir $h5dir --grid-res $resolution --geom-file "$h5dir"/junc.geom
    eps_file=$(ls | grep eps)
    echo "eps_fname = $eps_file"
    cp $eps_file "$h5dir/$eps_file"
fi
#make the test plots in check_enes.py
mkdir "$pname"/figures
if [ $make_movie == "t" ]; then
	python check_enes.py --prefix $h5dir --grid-res $resolution --gap-width ${widths[$((SLURM_ARRAY_TASK_ID-1))]} --gap-thick $thickness --movie cross
    rm -f "$pname"/figures/out.mp4 | tee $h5dir/tmp.log
	#combine the images to make an animation
    outname="$oname/figures/out_test_$SLURM_ARRAY_TASK_ID"
    echo "saving to $outname"
	ffmpeg -r 30 -f image2 -s 1920x1080 -i $h5dir/im_%d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p "$h5dir"/out.mp4
    mv "$h5dir"/out.mp4 "$outname".mp4
else
    python check_enes.py --prefix $h5dir --grid-res $resolution --gap-width ${widths[$((SLURM_ARRAY_TASK_ID-1))]} --gap-thick $thickness | tee $h5dir/tmp.log
fi

#make the frequency space plot
python time_space.py --prefix $h5dir --gap-width ${widths[$((SLURM_ARRAY_TASK_ID-1))]} --gap-thick $thickness

#move the plots into a folder where we can view them
cp "$h5dir"/eps.pdf "$oname"/figures/eps_"$SLURM_ARRAY_TASK_ID".pdf
cp "$h5dir"/fields.txt "$oname"/figures/fields_"$SLURM_ARRAY_TASK_ID".txt
cp "$h5dir"/field_samples.h5 "$oname"/figures/field_samples_"$SLURM_ARRAY_TASK_ID".h5
cp "$h5dir"/space_plot.pdf "$oname"/figures/space_plot_test_"$SLURM_ARRAY_TASK_ID".pdf
cp "$h5dir"/cross_plot.pdf "$oname"/figures/cross_plot_test_"$SLURM_ARRAY_TASK_ID".pdf
cp "$h5dir"/tdom_plot.pdf "$oname"/figures/tdom_plot_test_"$SLURM_ARRAY_TASK_ID".pdf
cp "$h5dir"/fdom_plot.pdf "$oname"/figures/fdom_plot_test_"$SLURM_ARRAY_TASK_ID".pdf
cp "$h5dir"/fdom_plot_2d.pdf "$oname"/figures/fdom_plot_2d_test_"$SLURM_ARRAY_TASK_ID".pdf
