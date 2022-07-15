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

    //iterate over user specified arguments for the test program
    int ret = parse_args(&args, &argc, argv);
    if (ret) return ret;
    //pass the rest of the arguments to meep_mpi initialization
    meep::initialize mpi(argc, argv);
    //load defaults from the configuration file
    char* name;
    //check if a name was supplied via command line argument
    if (args.conf_fname) {
        name = args.conf_fname;
    } else {
        std::string cpp_name("params.conf");
        name = strdup(cpp_name.c_str());
    }
    ret = parse_conf_file(&args, name);
    //free memory if it was allocated
    if (!args.conf_fname) free(name);
    if (ret) return ret;

    double z_center = args.len/2 + args.pml_thickness;
    double y_loc = z_center;

    //create the geometry object
    parse_ercode ercode;
    bound_geom geom(args, &ercode);
    if (ercode) return (int)ercode;

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
