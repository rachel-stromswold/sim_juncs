#include <meep.hpp>
#include <math.h>
#include <hdf5.h>
#include <chrono>

#include "argparse.h"
#include "disp.hpp"
#define PREFIX_LEN 3

const double epsilon = 0.01;
const double n_monitor_spheres = 2;

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

    //create the geometry object
    auto start = std::chrono::steady_clock::now();
    bound_geom geom(args, &ercode);
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
