#ifndef DISP_H
#define DISP_H

#include <meep.hpp>
#include <vector>
#include <math.h>
#include <hdf5.h>

#include "argparse.h"
#include "cgs.hpp"
#define PREFIX_LEN 3

//simulation parameters
const double sharpness = 5;
/*double z_center;
double eps_scale;*/

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

public:
    cgs_material_function(double p_def_ret=1.0);
    cgs_material_function(CompositeObject* p_volume, std::string type="eps", double p_def_ret=1.0);
    cgs_material_function(region_scale_pair p_volume, double p_def_ret=1.0);
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

std::vector<drude_suscept> parse_susceptibilities(char* const str, int* er);
meep::structure* structure_from_settings(const Settings& s, parse_ercode* ercode);

class bound_geom {
    public:
	bound_geom(const Settings& s, parse_ercode* ercode=NULL);
	~bound_geom();

	void add_point_source(meep::component, const meep::src_time &src, const meep::vec &, std::complex<double> amp = 1.0);
	void add_volume_source(meep::component c, const meep::src_time &src, const meep::volume &source_vol, std::complex<double> amp = 1.0);
	void run(const char* fname_prefix, std::vector<meep::vec> locs);

	Scene sc;

    private:
	//meep objects
	meep::grid_volume vol;
    meep::structure* strct = NULL;
	meep::fields fields;

	std::vector<meep::monitor_point*> monitor_locs;

	double ttot = 0;
	_uint n_t_pts = 0;	
    double post_source_t = 0.0;
};

#endif //DISP_H
