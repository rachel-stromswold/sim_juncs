#!/bin/bash
#SBATCH -p standard --time=20:00:00 --mem=244G
#SBATCH -a 1-14

#this is an SBATCH argument when running on femto: -C Gold6330

run_local="f"
run_simuls="t"
make_movie="f"
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
    oname="/scratch/$(whoami)/junc_simuls"
    
    #module load python3/3.8.8 meep/b3 hdf5/1.12.1/b1 openmpi/1.10.7/b1
fi

mkdir $oname
echo "working directory name: $pname, output name: $oname, run local? $run_local"

# comparisons of lengths
# lambda = 0.76 um, length = 0.05,0.2 um
# plasma frequency of gold=5.8 eV, corresponds to wavelength 0.213 um
# conparisons of frequencies
# lambda = 0.157, 0.76, 10.6 um, length=0.05 um <- Herman P.R., Schiffrin, Pattel
# bowtie vs rectangular comparison
# lambda = 0.76 um, length=0.05, 0.2 um

#r_760_50 r_760_200 r_760_1500 r_157_50 r_10600_50 b_760_50 b_760_200
#same for infinite
wavelengths=("0.76" "0.76" "0.157" "10.6" "0.76" "0.76" "0.76" "0.76" "0.157" "10.6" "0.76" "0.76")
widths=("0.05" "0.20" "0.05" "0.05" "0.05" "0.20" "0.05" "0.20" "0.05" "0.05" "0.05" "0.20")
juncs=("junc_box.geom" "junc_box.geom" "junc_box.geom" "junc_box.geom" "junc_bowtie.geom" "junc_bowtie.geom" "junc_box.geom" "junc_box.geom" "junc_box.geom" "junc_box.geom" "junc_bowtie.geom" "junc_bowtie.geom")
use_infs=("0" "0" "0" "0" "0" "0" "1" "1" "1" "1" "1" "1")
names=("760_50_fin" "760_200_fin" "157_50_fin" "10600_50_fin" "760_50_fin_bow" "760_200_fin_bow" "760_50_inf" "760_200_inf" "157_50_inf" "10600_50_inf" "760_50_inf_bow" "760_200_inf_bow")

#the widths used for each job in the array. These are similar to those used in Schiffrin et al. For a SiO2 junction.
thickness="0.2"
n_cycles="0.5"
resolution="12.0"
h5dir="$pname/test_$SLURM_ARRAY_TASK_ID"

#get the value stored in the array for this job id
slind=$((SLURM_ARRAY_TASK_ID-1))
junc=${juncs[$slind]}
width=${widths[$slind]}
wavelength=${wavelengths[$slind]}
use_inf=${use_infs[$slind]}
suf_name=${names[$slind]}
this_n=${n_goups[$slind]}

#make the h5 directory and subdirectories
mkdir $h5dir
mkdir $h5dir/fit_figs

if [ $run_simuls == "t" ]; then
    rm -f $h5dir/*
    #valgrind --leak-check=full --track-origins=yes ./sim_geom --out-dir $h5dir --grid-res ${resolutions[$((SLURM_ARRAY_TASK_ID-1))]}
    if [ $wavelength == "0.157" ]; then
        ./sim_geom --out-dir $h5dir --grid-res $resolution --geom-file $junc --save-span 10 --opts "width=$width;thick=$thickness;inf_thick=$use_inf;wavelen=$wavelength;n_cycles=$n_cycles"
    else
        ./sim_geom --out-dir $h5dir --grid-res $resolution --geom-file $junc  --opts "width=$width;thick=$thickness;inf_thick=$use_inf;wavelen=$wavelength;n_cycles=$n_cycles"
    fi
fi
#make the test plots in check_enes.py
mkdir "$pname"/figures
if [ $make_movie == "t" ]; then
	python check_enes.py --prefix $h5dir --grid-res $resolution --gap-width $width --gap-thick $thickness --movie cross
    rm -f "$pname"/figures/out.mp4 | tee $h5dir/tmp.log
	#combine the images to make an animation
    outname="$oname/figures/out_test_$SLURM_ARRAY_TASK_ID"
    echo "saving to $outname"
	ffmpeg -r 30 -f image2 -s 1920x1080 -i $h5dir/im_%d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p "$h5dir"/out.mp4
    mv "$h5dir"/out.mp4 "$outname".mp4
else
        python check_enes.py --prefix $h5dir --grid-res $resolution --gap-width $width --gap-thick $thickness | tee $h5dir/tmp.log
fi

#plot samples of the fields
python time_space.py --fname $h5dir/field_samples.h5 --prefix $h5dir
#set the parameters used by the paper for figures
plotting_opts=""
if [ [ $width == "0.05" ] && [ $junc == *"box"* ] ] || [ $wavelength == "0.157" ]; then
    #this is a special case that appears twice in the figures
    python phase_plot.py --fname $h5dir/field_samples.h5 --prefix $h5dir --gap-width $width --gap-thick $thickness
    cp "$h5dir"/heatmap_grp0.svg "$oname"/figures/heatmap_center_"$suf_name"_no_y.svg
    cp "$h5dir"/heatmap_grp1.svg "$oname"/figures/heatmap_offset_"$suf_name"_no_y.svg
    plotting_opts=$plotting_opts + "--plot-y-labels"
elif [ [ $width == "0.20" ] && [ $junc == *"bowtie"* ] ] || [ $wavelength ==  "10.6" ]; then
     plotting_opts=$plotting_opts + "--plot-cbar"
fi
if [ $use_inf == "1" ]; then
    plotting_opts=$plotting_opts + "--plot-x-labels"
fi
python phase_plot.py --fname $h5dir/field_samples.h5 --prefix $h5dir --gap-width $width --gap-thick $thickness $plotting_opts

#move the plots into a folder where we can view them
mkdir $oname/figures
cp "$h5dir"/field_samples.h5 "$oname"/field_samples_"$suf_name".h5
cp "$h5dir"/space_plot.svg "$oname"/figures/space_plot_"$suf_name".svg
cp "$h5dir"/cross_plot.svg "$oname"/figures/cross_plot_"$suf_name".svg
cp "$h5dir"/tdom_plot.svg "$oname"/figures/tdom_plot_"$suf_name".svg
cp "$h5dir"/amps.svg "$oname"/figures/amps_"$suf_name".svg
cp "$h5dir"/avgs.svg "$oname"/figures/avgs_"$suf_name".svg
cp "$h5dir"/sigs.svg "$oname"/figures/sigs_"$suf_name".svg
cp "$h5dir"/heatmap_grp0.svg "$oname"/figures/heatmap_center_"$suf_name".svg
cp "$h5dir"/heatmap_grp1.svg "$oname"/figures/heatmap_offset_"$suf_name".svg
cp "$h5dir"/junc.pgm "$oname"/figures/junc_"$suf_name".pgm
