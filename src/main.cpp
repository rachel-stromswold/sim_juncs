#include <meep.hpp>
#include <math.h>
#include <hdf5.h>

#include "argparse.h"
#include "disp.hpp"
#define PREFIX_LEN 3

const double epsilon = 0.01;
const double n_monitor_spheres = 2;

int main(int argc, char **argv) {
    Settings args;

    //load defaults from the configuration file
    /*std::string name = "params.conf";
    char* name_dup = strdup(name.c_str());
    int ret = parse_conf_file(&args, name_dup);
    free(name_dup);
    if (ret) return ret;
    //iterate over user specified arguments for the test program
    ret = parse_args(&args, &argc, argv);
    if (ret) return ret;*/
    //pass the rest of the arguments to meep_mpi initialization
    meep::initialize mpi(argc, argv);

    double z_center = args.len/2 + args.pml_thickness;
    double y_loc = z_center;

    //create the geometry object
    bound_geom geom(args);

    //create the EM-wave source at the specified location
    meep::gaussian_src_time src(args.freq, args.gauss_width, args.gauss_start_time, args.gauss_cutoff*args.gauss_width);
    //pulse continuously for a time of 2*n_devs standard deviations with an instantaneous turn on and off
    //meep::continuous_src_time src(args.freq, 0.0, 0.0, 2*args.gauss_width*args.gauss_cutoff);
    meep::vec source_loc(y_loc, args.source_z);
    if (args.n_dims == 3) {
	//initialize a source right next to the upper pml
	meep::volume source_vol(meep::vec(0, 0, args.pml_thickness), meep::vec(2*z_center, 2*z_center, args.pml_thickness));
	geom.add_volume_source(meep::Ex, src, source_vol, args.amp);
    } else {
	if (args.n_dims == 1) source_loc = meep::vec(args.source_z);
	geom.add_point_source(meep::Ex, src, source_loc, args.amp);
    }

    //setup four points for monitoring the Poynting vector as a function of time. Two are near the source to the left and right, the others are near the simulation boundaries. All energy must pass through either the left or right (assuming 1 dimension)
    std::vector<meep::vec> monitor_locs;
    if (args.n_dims == 1) {
        monitor_locs.emplace_back(source_loc.z() - args.src_mon_dist);
        monitor_locs.emplace_back(source_loc.z() + args.src_mon_dist);
        monitor_locs.emplace_back(z_center - epsilon);
        monitor_locs.emplace_back(z_center + epsilon);
        monitor_locs.emplace_back(args.pml_thickness);
        monitor_locs.emplace_back(2*z_center - args.pml_thickness);
    } else if (args.n_dims == 2) {
        monitor_locs.emplace_back(y_loc, source_loc.z() - args.src_mon_dist);
        monitor_locs.emplace_back(y_loc, source_loc.z() + args.src_mon_dist);
        monitor_locs.emplace_back(y_loc, z_center - epsilon);
        monitor_locs.emplace_back(y_loc, z_center + epsilon);
        monitor_locs.emplace_back(y_loc, args.pml_thickness);
        monitor_locs.emplace_back(y_loc, 2*z_center - args.pml_thickness);
    } else {
        monitor_locs.emplace_back(z_center, y_loc, source_loc.z() - args.src_mon_dist);
        monitor_locs.emplace_back(z_center, y_loc, source_loc.z() + args.src_mon_dist);
        monitor_locs.emplace_back(z_center, y_loc, z_center - epsilon);
        monitor_locs.emplace_back(z_center, y_loc, z_center + epsilon);
        monitor_locs.emplace_back(z_center, y_loc, args.pml_thickness);
        monitor_locs.emplace_back(z_center, y_loc, 2*z_center - args.pml_thickness);
    }

    geom.run(args.out_dir, monitor_locs);

    cleanup_settings(&args);
    return 0;
}
