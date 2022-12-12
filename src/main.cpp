#include <meep.hpp>
#include <math.h>
#include <hdf5.h>

#include "argparse.h"
#include "disp.hpp"
#define PREFIX_LEN 3

const double epsilon = 0.01;
const double n_monitor_spheres = 2;

/**
 * insertMember isn't available in all versions of hdf5, and there were some issues that popped up calling it after running the simulations. It's better to fail early if something is wrong
 */
int test_h5_funcs() {
    printf("Now trying HDF5 library files\n");
    int ret = 1;
    /*try {
        H5::CompType fieldtype(sizeof(complex));
        fieldtype.insertMember("Re", HOFFSET(complex, re), H5_float_type);
        fieldtype.insertMember("Im", HOFFSET(complex, im), H5_float_type);
    } catch (...) {
        printf("insertMember doesn't work!\n");
        ret = 0;
    }*/
    return ret;
}

int main(int argc, char **argv) {
    parse_settings args;

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
    meep::verbosity = args.verbosity;
    //free memory if it was allocated
    if (!args.conf_fname) free(name);
    if (ret) return ret;
    correct_defaults(&args);

    //test hdf5
    int has_insert_member = test_h5_funcs();

    //create the geometry object
    parse_ercode ercode = E_SUCCESS;
    bound_geom geom(args, &ercode);
    if (ercode) return (int)ercode;

    geom.run(args.out_dir);
    //save the time series for fields
    geom.save_field_times(args.out_dir);

    cleanup_settings(&args);
    return 0;
}
