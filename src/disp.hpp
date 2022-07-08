#ifndef GEOM_H
#define GEOM_H

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

/**
 * This is identical to the simple_material_function, but each term is multiplied by a constant scalar
 */
class cgs_material_function : public meep::material_function {
    CompositeObject* volume;
    
    //to reduce numeric artifacts (maybe?) you can apply a smooth sigmoid scaling function
    //TODO: find a way to make this work with arbitrary shapes
    double bound_sharpness;
    double bound_eps_scale;
    //returns whether we're to the left or the right of the dielectric boundary
    double in_bound(const meep::vec &r);

    //scaling and offset
    double sc;
    double off;

public:
    cgs_material_function(CompositeObject* p_volume, double scale, double offset, double p_bound_eps_scale=0.1, double p_bound_sharpness=5);

    virtual ~cgs_material_function() {}

    virtual double chi1p1(meep::field_type ft, const meep::vec &r);
    virtual double eps(const meep::vec &r);
    virtual double mu(const meep::vec &r);
    virtual double conductivity(meep::component c, const meep::vec &r);
    virtual void sigma_row(meep::component c, double sigrow[3], const meep::vec &r);
    virtual double chi3(meep::component c, const meep::vec &r);
    virtual double chi2(meep::component c, const meep::vec &r);
};

class bound_geom {
    public:
	bound_geom(const Settings& s);
	~bound_geom();

	void add_point_source(meep::component, const meep::src_time &src, const meep::vec &, std::complex<double> amp = 1.0);
	void add_volume_source(meep::component c, const meep::src_time &src, const meep::volume &source_vol, std::complex<double> amp = 1.0);
	void run(const char* fname_prefix, std::vector<meep::vec> locs);

	Scene sc;

    private:
	//meep objects
	meep::grid_volume vol;
	meep::structure* strct = NULL;
	meep::fields* fields = NULL;

	std::vector<meep::monitor_point*> monitor_locs;

	double ttot = 0;
	_uint n_t_pts = 0;
};

#endif //GEOM_H
