#include <meep.hpp>
#include <math.h>
#include <hdf5.h>
#include <chrono>

#include "argparse.h"
#include "disp.hpp"
#define PREFIX_LEN 3

const double epsilon = 0.01;
const double n_monitor_spheres = 2;

<<<<<<< HEAD
<<<<<<< Updated upstream
=======
=======
>>>>>>> a24f98f4a6e0cd48991594afc0e6ca961d44d232
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

<<<<<<< HEAD
/**
 * Load variables specified in the parse_settings struct s into the context con
 */
context context_from_settings(parse_settings& args) {
    context con;
    con.emplace("pml_thickness", make_val_num(args.pml_thickness));
    con.emplace("sim_length", make_val_num(args.len));
    con.emplace("length", make_val_num(2*args.pml_thickness + args.len));
    con.emplace("l_per_um", make_val_num(args.um_scale));
    //TODO: use context to handle this more elegantly
    if (args.user_opts) {
	char* save;
	char* line = strtok_r(args.user_opts, ";", &save);
	while (line) {
	    char* eq_loc = strchr(line, '=');
	    if (eq_loc) {
		*eq_loc = 0;
		parse_ercode er = E_SUCCESS;
		value tmp = con.parse_value(eq_loc+1, er);
		if (er) con.emplace(line, tmp);
		cleanup_val(&tmp);
	    }
	    line = strtok_r(args.user_opts, ";", &save);
	}
    }
    return con;
}

=======
>>>>>>> a24f98f4a6e0cd48991594afc0e6ca961d44d232
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

    //initialize a context for the scene based on parameters
    parse_ercode ercode = E_SUCCESS;
    context scene_con = context_from_settings(args);

    //test hdf5
    int has_insert_member = test_h5_funcs();

    //create the geometry object
    auto start = std::chrono::steady_clock::now();
    bound_geom geom(args, scene_con, &ercode);
    if (ercode) return (int)ercode;
    //time benchmarking
    auto end_init = std::chrono::steady_clock::now();
    auto this_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_init-start).count();
    std::cout << "initialization completed in " << this_time << " ms" << std::endl;

    geom.run(args.out_dir);
    //time benchmarking
    auto end_run = std::chrono::steady_clock::now();
    this_time = std::chrono::duration_cast<std::chrono::seconds> (end_run-end_init).count();
    std::cout << "simulation completed in " << this_time << " s" << std::endl;
    //save the time series for fields
    geom.save_field_times(args.out_dir);
    //time benchmarking
    auto end_write = std::chrono::steady_clock::now();
    this_time = std::chrono::duration_cast<std::chrono::milliseconds> (end_write-end_run).count();
    std::cout << "saving timeseries completed in " << this_time << " ms" << std::endl;
    this_time = std::chrono::duration_cast<std::chrono::seconds> (end_write-start).count();
    std::cout << "total time: " << this_time << " s" << std::endl;

    cleanup_settings(&args);
    return 0;
}
