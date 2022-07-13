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
    double def_ret;

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
    _uint n_susceptibilities = 0;
    double* eps_2_omega;
    double* eps_2_gamma;
    double* eps_2_sigma;
    bool* eps_2_use_denom;
} susceptibility_list;

inline void cleanup_susceptibility_list(susceptibility_list* s) {
    if (s->eps_2_omega) free(s->eps_2_omega);
    if (s->eps_2_gamma) free(s->eps_2_gamma);
    if (s->eps_2_sigma) free(s->eps_2_sigma);
    if (s->eps_2_use_denom) free(s->eps_2_use_denom);
}

/**
 * Add the list of susceptibilities to the Settings file s. The string should have the format (omega_0,gamma_0,sigma_0),(omega_1,gamma_1,sigma_1),...
 * returns: 0 on success or an error code
 * 		-1 null string
 * 		-2 invalid or empty string
 * 		-3 insufficient memory
 */
inline int parse_susceptibilities(susceptibility_list* s, char* const str) {
    //check that the Settings struct is valid and allocate memory
    if (!s) return -1;
    _uint buf_size = BUF_SIZE;
    s->eps_2_omega = (double*)malloc(buf_size*sizeof(double));
    s->eps_2_gamma = (double*)malloc(buf_size*sizeof(double));
    s->eps_2_sigma = (double*)malloc(buf_size*sizeof(double));
    s->eps_2_use_denom = (bool*)malloc(buf_size*sizeof(int));
    //TODO: get gud
    //if (!(s->eps_2_omega && s->eps_2_gamma && s->eps_2_sigma && s->eps_2_use_denom)) 

    double cur_omega = 0.0;
    double cur_gamma = 0.0;
    double cur_sigma = 0.0;
    int cur_use_denom = 0;

    //used for strtok_r
    char* save_str;
    char* tok;
    //find the first entry
    char* cur_entry = strchr(str, '(');
    char* end = strchr(str, ')');

    //only proceed if we have pointers to the start and end of the current entry
    while (cur_entry && end) {
	//resize buffers if necessary
	if (s->n_susceptibilities == buf_size) {
	    buf_size *= 2;
	    s->eps_2_omega = (double*)realloc(s->eps_2_omega, buf_size*sizeof(double));
	    s->eps_2_gamma = (double*)realloc(s->eps_2_gamma, buf_size*sizeof(double));
	    s->eps_2_sigma = (double*)realloc(s->eps_2_sigma, buf_size*sizeof(double));
	    s->eps_2_use_denom = (bool*)realloc(s->eps_2_use_denom, buf_size*sizeof(double));
	}
	cur_omega = 0.0;
	cur_gamma = 0.0;
	cur_use_denom = 1;
	//the open paren must occur before the end paren
	if (cur_entry > end) return -2;

	//null terminate the parenthesis and tokenize by commas
	end[0] = 0;
	//read the omega value
	tok = trim_whitespace( strtok_r(cur_entry+1, ",", &save_str), NULL );
	cur_omega = strtod(tok, NULL);
	if (errno) { errno = 0;return -2; }
	//read the gamma value
	tok = trim_whitespace( strtok_r(NULL, ",", &save_str), NULL );
	if (!tok) return -2;
	cur_gamma = strtod(tok, NULL);
	if (errno) { errno = 0;return -2; }
	//read the sigma value
	tok = trim_whitespace( strtok_r(NULL, ",", &save_str), NULL );
	if (!tok) return -2;
	cur_sigma = strtod(tok, NULL);
	if (errno) { errno = 0;return -2; }
	//read the (optional) use denom flag
	tok = strtok_r(NULL, ",", &save_str);
	if (tok) {
	    tok = trim_whitespace(tok, NULL);
	    cur_use_denom = strtol(tok, NULL, 10);
	    if (errno) {
		errno = 0;
		cur_use_denom=1;
		if (strcmp(tok, "drude") == 0 || strcmp(tok, "true") == 0) cur_use_denom = 0;
	    }
	    if (strcmp(tok, "lorentz") == 0 || strcmp(tok, "false") == 0) cur_use_denom = 1;
	}

	//save the information
	s->eps_2_omega[s->n_susceptibilities] = cur_omega;
	s->eps_2_gamma[s->n_susceptibilities] = cur_gamma;
	s->eps_2_sigma[s->n_susceptibilities] = cur_sigma;
	s->eps_2_use_denom[s->n_susceptibilities] = cur_use_denom;
	++s->n_susceptibilities;

	//advance to the next entry
	if (end[1] == 0) break;
	cur_entry = strchr(end+1, '(');
	end = strchr(cur_entry, ')');
    }

    //realloc arrays to exactly fit memory
    buf_size = s->n_susceptibilities;
    s->eps_2_omega = (double*)realloc(s->eps_2_omega, buf_size*sizeof(double));
    s->eps_2_gamma = (double*)realloc(s->eps_2_gamma, buf_size*sizeof(double));
    s->eps_2_sigma = (double*)realloc(s->eps_2_sigma, buf_size*sizeof(double));
    s->eps_2_use_denom = (bool*)realloc(s->eps_2_use_denom, buf_size*sizeof(int));

    return 0;
}

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

#endif //DISP_H
