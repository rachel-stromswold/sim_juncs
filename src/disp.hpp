#ifndef DISP_H
#define DISP_H

#include <meep.hpp>
#include <vector>
#include <random>
#include <math.h>
#include <H5Cpp.h>
#include <iostream>

#include "argparse.h"
#include "data_utils.hpp"
#include "cgs.hpp"
#define PREFIX_LEN 3

//meep seems to be confused about how large fields are. On the cluster I consistently saw writes past the official meep bounds. This constant is here to help prevent data clobbering
#define FIELDS_PAD_SIZE 10

#define DEFAULT_WIDTH_N 5
#define THICK_SCALE 1.0

//simulation parameters
const double sharpness = 5;
const double DEF_SEED = 0xd9a28bf3;
/*double z_center;
double eps_scale;*/

#ifdef STO_PREC_32
const H5::DataType H5_float_type = H5::PredType::NATIVE_FLOAT;
#else
const H5::DataType H5_float_type = H5::PredType::NATIVE_DOUBLE;
#endif

typedef struct {
    CompositeObject* c;
    double s;
} region_scale_pair;

/**
 * This is identical to the simple_material_function, but each term is multiplied by a constant scalar
 */
class cgs_material_function : public meep::material_function {
    _uint n_regions;
    region_scale_pair* regions;
    
    //scaling and offset
    double def_ret = 1.0;

    //we use a stochastic smoothing process
    _uint smooth_n = 1;
    double* smooth_pts = NULL;
    void generate_smooth_pts(double smooth_rad, uint64_t seed);

public:
    cgs_material_function(double p_def_ret=1.0, _uint p_smooth_n=1, double p_smooth_rad=DEFAULT_SMOOTH_RAD);
    cgs_material_function(CompositeObject* p_volume, std::string type="eps", double p_def_ret=1.0, _uint p_smooth_n=1, double p_smooth_rad=DEFAULT_SMOOTH_RAD);
    cgs_material_function(region_scale_pair p_volume, double p_def_ret=1.0, _uint p_smooth_n=1, double p_smooth_rad=DEFAULT_SMOOTH_RAD);
    cgs_material_function(const cgs_material_function& o);
    cgs_material_function(cgs_material_function&& o);
    cgs_material_function& operator=(cgs_material_function&& o);
    ~cgs_material_function();

    void add_region(CompositeObject* p_reg, std::string type="eps");
    virtual double chi1p1(meep::field_type ft, const meep::vec &r);
    virtual double eps(const meep::vec &r);
    virtual double mu(const meep::vec &r);
    virtual double conductivity(meep::component c, const meep::vec &r);
    virtual void sigma_row(meep::component c, double sigrow[3], const meep::vec &r);
    virtual double chi3(meep::component c, const meep::vec &r);
    virtual double chi2(meep::component c, const meep::vec &r);
    //returns whether we're to the left or the right of the dielectric boundary
    double in_bound(const meep::vec &r);
};

typedef struct {
    double omega_0;
    double gamma;
    double sigma;
    bool use_denom;
} drude_suscept;

typedef enum { SRC_GAUSSIAN, SRC_CONTINUOUS } src_type;
class source_info {
public:
    src_type type;
    meep::component component;
    double freq;
    double width;
    double start_time;
    double end_time;
    //double cutoff;
    double amplitude;

    source_info(std::string spec_str, const Scene& problem, parse_ercode* ercode);
};

struct sto_vec {
    _ftype x;
    _ftype y;
    _ftype z;
    sto_vec(_ftype xx, _ftype yy, _ftype zz) { x = xx;y = yy;z = zz; }
};

class bound_geom {
public:
    bound_geom(const Settings& s, parse_ercode* ercode=NULL);
    ~bound_geom();
    std::vector<meep::vec> get_monitor_locs();
    std::vector<data_arr> get_field_times() { return field_times; }

    void add_point_source(meep::component, const meep::src_time &src, const meep::vec &, std::complex<double> amp = 1.0);
    void add_volume_source(meep::component c, const meep::src_time &src, const Box& vol, std::complex<double> amp = 1.0);
    void run(const char* fname_prefix);
    void save_field_times(const char* fname_prefix, size_t* save_pts=NULL, size_t n_save_pts=0);

    double fs_to_meep_time(double t) { return t*0.299792458*um_scale; }
    double meep_time_to_fs(double t) { return t/(0.299792458*um_scale); }

    Scene problem;

#ifdef DEBUG_INFO
    std::vector<source_info> sources;
#endif
    std::vector<drude_suscept> parse_susceptibilities(char* const str, int* er);
    meep::structure* structure_from_settings(const Settings& s, Scene& problem, parse_ercode* ercode);

private:
    //TODO: this is rather sketchy, but for some reason using an std vector randomly changed sizes. I assumed this was a memory error but valgrind made no complaints.
    size_t n_locs = 0;
    meep::vec* monitor_locs = NULL;
    std::vector<data_arr> field_times;

    //meep objects
    meep::grid_volume vol;
    meep::structure* strct = NULL;
    meep::fields fields;

    //meep:fields seems to be confused about how long it is, so we add a whole bunch of dummy variables so that it doesn't clobber other variables
    char dummy_vals[FIELDS_PAD_SIZE];

    double um_scale;
    double ttot = 0;
    _uint n_t_pts = 0;
    double post_source_t = 0.0;

    //the arguments supplied will alter the location of the dielectric
    double pml_thickness;
    double len;
    double z_center;
    double eps_scale;
    _uint dump_span = 20;
};

#endif //DISP_H
