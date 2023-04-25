# sim_juncs
A tool for simulating molecular junctions illuminated by laser pulses using FDTD methods

# dependencies
* python3
* open MPI (https://github.com/open-mpi/ompi)
* hdf5 (https://www.hdfgroup.org/solutions/hdf5/)
* meep v1.24 (https://github.com/NanoComp/meep/)

To build the project you will also need a C++ compiler and either cmake or meson.

# building
To build the project first navigate to the desired installation directory. Then enter the following commands:
```
git clone https://github.com/rachel-stromswold/sim_juncs
cd sim_juncs
mkdir build
```
If you are using cmake, run the following:
```
cd build
cmake .. && make
```
If you are using meson, run this instead:
```
meson setup build && cd build
ninja compile
```

## cmake options
By default, configuration files for a graphene junction (with parameters based on Boolakee, T. et al. Nature 605, 251–255 (2022)) are copied into the build directory. To use a diferent junction geometry, you can pass the `-DJUNCTION_TYPE` argument when calling cmake. Valid names for this option are given as directories in the `junctions` folder. For example, to use a silica junction you can call the following:
```
cmake .. -DJUNCTION_TYPE=Au_SiO2_box
```
The same options are available in meson.

# running
Once you've completed the previous steps, call `run.sh -l` inside the build directory. If you want to run jobs on a cluster using slurm call `sbatch run.sh`.

# thank you to
* meep for supplying FDTD electromagnetic simulations
* eigen3 for providing utilities for linear algebra and geometry
