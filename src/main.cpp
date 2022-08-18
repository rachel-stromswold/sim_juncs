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
    correct_defaults(&args);

    double z_center = args.len/2 + args.pml_thickness;
    double y_loc = z_center;

    //create the geometry object
    parse_ercode ercode = E_SUCCESS;
    bound_geom geom(args, &ercode);
    if (ercode) return (int)ercode;

    geom.run(args.out_dir);

    cleanup_settings(&args);
    return 0;
}
