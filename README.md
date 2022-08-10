# sim_juncs
A tool for simulating molecular junctions illuminated by laser pulses using FDTD methods

# dependencies
-meep v1.24 (https://github.com/NanoComp/meep/)
-hdf5
-openmpi
-eigen3 (https://eigen.tuxfamily.org/)

To build the project you will also need a C++ compiler and cmake.

# building
To build the project first navigate to the desired installation directory. Then enter the following commands:
```
git clone https://github.com/rachel-stromswold/sim_juncs
cd sim_juncs
mkdir build
cd build
cmake .. && make
```

## cmake options
By default, configuration files for a graphene junction are copied into the build directory as per Boolakee, T. et al. Nature 605, 251–255 (2022). To use a diferent junction geometry, you can pass the `-DJUNCTION_TYPE` argument when calling cmake. Valid names for this option are given as directories in the `junctions` folder. For example, to use a silica junction you can call the following:
```
cmake .. -DJUNCTION_TYPE=Au_SiO2_box
```

# running
Once you've completed the previous steps, call `run.sh -l` inside the build directory. If you want to run jobs on a cluster using slurm call `sbatch run.sh`.
