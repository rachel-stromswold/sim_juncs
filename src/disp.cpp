#include "disp.hpp"

const double eps_cutoff = 10;

Settings args;

/**
 * This is identical to the simple_material_function, but each term is multiplied by a constant scalar
 */
cgs_material_function::cgs_material_function(double p_def_ret) {
    def_ret = p_def_ret;
}

/**
 * This is identical to the simple_material_function, but each term is multiplied by a constant scalar
 */
cgs_material_function::cgs_material_function(CompositeObject* p_volume, std::string type, double p_def_ret) {
    //set the regions which are used
    n_regions = 1;
    regions = (region_scale_pair*)malloc(sizeof(region_scale_pair));
    if (regions) {
	double scale = def_ret;
	//lookup the scaling constant by reading metadata from the object
	if (p_volume->has_metadata(type)) {
	    scale = std::stod(p_volume->fetch_metadata(type));
	}
	regions[0].s = scale;//by default set the scale to whatever the default return value is
	regions[0].c = p_volume;
    } else {
	n_regions = 0;
    }
    //set default return value
    def_ret = p_def_ret;
}

/**
 * This is identical to the simple_material_function, but each term is multiplied by a constant scalar
 */
cgs_material_function::cgs_material_function(region_scale_pair p_reg, double p_def_ret) {
    //set the regions which are used
    n_regions = 1;
    regions = (region_scale_pair*)malloc(sizeof(region_scale_pair));
    if (regions) {
	regions[0] = p_reg;
    } else {
	n_regions = 0;
    }
    //set default return value
    def_ret = p_def_ret;
}

//copy constructor
cgs_material_function::cgs_material_function(const cgs_material_function& o) {
    def_ret = o.def_ret;
    n_regions = o.n_regions;
    regions = (region_scale_pair*)malloc(sizeof(region_scale_pair)*n_regions);
    //only duplicate entries if successfully allocated
    if (regions) {
	for (_uint i = 0; i < n_regions; ++i) regions[i] = o.regions[i];
    } else {
	n_regions = 0;
    }
}

//move constructor
cgs_material_function::cgs_material_function(cgs_material_function&& o) {
    def_ret = o.def_ret;
    n_regions = o.n_regions;
    regions = o.regions;
    o.regions = NULL;
    o.n_regions = 0;
}

//= constructor
cgs_material_function& cgs_material_function::operator=(cgs_material_function&& o) {
    //swap the default return values
    double tmp_def_ret = def_ret;
    def_ret = o.def_ret;
    o.def_ret = tmp_def_ret;
    //swap the regions and their number
    _uint tmp_n_regions = n_regions;
    region_scale_pair* tmp_regions = regions;
    n_regions = o.n_regions;
    regions = o.regions;
    o.n_regions = tmp_n_regions;
    o.regions = tmp_regions;

    return *this;
}

cgs_material_function::~cgs_material_function() {
    if (regions) free(regions);
    n_regions = 0;
}

/**
 * Add the region based on the CompositeObject p_reg.
 * p_reg: the composite object specifying the region to add
 * type: lookup the matching metadata key from the composite object and set the scale based on its value.
 */
void cgs_material_function::add_region(CompositeObject* p_reg, std::string type) {
    _uint new_n_regions = n_regions+1;
    region_scale_pair* tmp_regs = (region_scale_pair*)realloc(regions, sizeof(region_scale_pair)*new_n_regions);
    //only proceed if allocation of memory was successful
    if (tmp_regs) {
	n_regions = new_n_regions;
	regions[n_regions-1].c = p_reg;
	double scale = def_ret;
	//lookup the scaling constant by reading metadata from the object
	if (p_reg->has_metadata(type)) {
	    scale = std::stod(p_reg->fetch_metadata(type));
	}
	regions[n_regions-1].s = scale;//by default set the scale to whatever the default return value is
    }
}

//returns whether we're to the left or the right of the dielectric boundary 1 indicates we are, 0 that we are not
double cgs_material_function::in_bound(const meep::vec &r) {
    if (!regions) return def_ret;
    //look through each region and return the first that contains the vector
    for (_uint i = 0; i < n_regions; ++i) {
	int ret = regions[i].c->in(evec3(r.x(), r.y(), r.z()));
	if (ret) return regions[i].s;
    }
    return def_ret;
}

double cgs_material_function::chi1p1(meep::field_type ft, const meep::vec &r) {
    (void)ft;
    return in_bound(r);
}

double cgs_material_function::eps(const meep::vec &r) {
    return in_bound(r);
}

double cgs_material_function::mu(const meep::vec &r) {
    return in_bound(r);
}

double cgs_material_function::conductivity(meep::component c, const meep::vec &r) {
    (void)c;
    return in_bound(r);
}

void cgs_material_function::sigma_row(meep::component c, double sigrow[3], const meep::vec &r) {
    sigrow[0] = sigrow[1] = sigrow[2] = 0.0;
    sigrow[meep::component_index(c)] = in_bound(r);
}

double cgs_material_function::chi3(meep::component c, const meep::vec &r) {
    (void)c;
    return in_bound(r);
}

double cgs_material_function::chi2(meep::component c, const meep::vec &r) {
    (void)c;
    return in_bound(r);
}

bound_geom::bound_geom(const Settings& s) : sc(s.geom_fname) {
    //the arguments supplied will alter the location of the dielectric
    double z_center = s.len/2 + s.pml_thickness;
    double eps_scale = 1 / (sharpness*args.resolution);

    //initialize the volume
    if (s.n_dims == 1) {
	vol = meep::vol1d(2*z_center, s.resolution);
    } else if (s.n_dims == 2) {
	vol = meep::vol2d(2*z_center, 2*z_center, s.resolution);
    } else {
	vol = meep::vol3d(2*z_center, 2*z_center, 2*z_center, s.resolution);
    }

    //iterate over all objects specified in the scene
    std::vector<CompositeObject*> roots = sc.get_roots();
    //setup the structure with the infinite frequency dielectric component
    cgs_material_function inf_eps_func(s.ambient_eps);

    for (auto it = roots.begin(); it != roots.end(); ++it) {
	inf_eps_func.add_region(*it);
	/** ============================ DEBUG ============================ **/
	meep::vec test_loc_1(0.5,0.5,2);
	meep::vec test_loc_2(0.5,0.1,6);
	meep::vec test_loc_3(0.5,8,6);
	meep::vec test_loc_4(0.5,12,6);
	double ret_1 = inf_eps_func.eps(test_loc_1);
	double ret_2 = inf_eps_func.eps(test_loc_2);
	double ret_3 = inf_eps_func.eps(test_loc_1);
	double ret_4 = inf_eps_func.eps(test_loc_2);
	printf("%f %f %f %f\n", ret_1, ret_2, ret_3, ret_4);
	/** ============================ DEBUG ============================ **/
    }
    strct = new meep::structure(vol, inf_eps_func, meep::pml(args.pml_thickness));
    //read susceptibilities if they are available
    for (auto it = roots.begin(); it != roots.end(); ++it) {
	susceptibility_list cur_sups;
	int res = 0;
	if ((*it)->has_metadata("susceptibilities")) {
	    char* dat = strdup((*it)->fetch_metadata("susceptibilities").c_str());
	    res = parse_susceptibilities(&cur_sups, dat);
	    free(dat);
	}
	//add frequency dependent susceptibility
	for (_uint i = 0; i < cur_sups.n_susceptibilities; ++i) {
	    meep::lorentzian_susceptibility suscept( cur_sups.eps_2_omega[i]/s.um_scale, cur_sups.eps_2_gamma[i]/s.um_scale, !(cur_sups.eps_2_use_denom[i]) );
	    region_scale_pair tmp_pair = {*it, cur_sups.eps_2_sigma[i]};
	    cgs_material_function scale_func(tmp_pair, 0.0);
	    strct->add_susceptibility(scale_func, meep::E_stuff, suscept);
	}
    }

    //create the fields
    fields = new meep::fields(strct);
}

bound_geom::~bound_geom() {
    if (strct) delete strct;
    if (fields) delete fields;
    //delete all monitor locations
    for (_uint i = 0; i < monitor_locs.size(); ++i) delete monitor_locs[i];
}

void bound_geom::add_point_source(meep::component c, const meep::src_time &src, const meep::vec& source_loc, std::complex<double> amp) {
    fields->add_point_source(c, src, source_loc, args.amp);

    //set the total timespan based on the added source
    ttot = fields->last_source_time() + args.post_source_t;
    n_t_pts = (_uint)(ttot / fields->dt);
}

void bound_geom::add_volume_source(meep::component c, const meep::src_time &src, const meep::volume &source_vol, std::complex<double> amp) {
    fields->add_volume_source(c, src, source_vol, args.amp);

    //set the total timespan based on the added source
    ttot = fields->last_source_time() + args.post_source_t;
    n_t_pts = (_uint)(ttot / fields->dt);
}

void bound_geom::run(const char* fname_prefix, std::vector<meep::vec> locs) {
    fields->set_output_directory(fname_prefix);

    //open the file which will store poynting vector fluxes
    char flux_name[BUF_SIZE];
    snprintf(flux_name, BUF_SIZE, "%s/Poynting_fluxes.txt", fname_prefix);
    FILE* fp = fopen(flux_name, "w");

    //initialize the array of monitor locations and write to the fluxes header
    /*fprintf(fp, "#time, ");
    monitor_locs.resize(locs.size());
    for (_uint i = 0; i < locs.size(); ++i) {
	monitor_locs[i] = fields->get_new_point(locs[i]);
	fprintf(fp, "(%f,%f,%f) ", locs[i].x(), locs[i].y(), locs[i].z());
    }
    fprintf(fp, "\n");*/

    //figure out the number of digits before the decimal and after
    int n_digits_a = (int)(ceil(log(ttot)/log(10)));
    double rat = -log((double)(fields->dt))/log(10.0);
    int n_digits_b = (int)ceil(rat)+1;
    char h5_fname[BUF_SIZE];
    strcpy(h5_fname, "ex-");

    _uint i = 0;
    for (; fields->time() < ttot; ++i) {
        //magnetic and electric fields are stored at different times, we need to synchronize
        /*fields->synchronize_magnetic_fields();

	//fetch monitor points
	for (_uint j = 0; j < locs.size(); ++j) {
	    fields->get_point(monitor_locs[j], locs[j]);
	    fprintf(fp, "%f ", fields->get_field(meep::Ex, locs[j]));
	}
	fprintf(fp, "\n");

        //restore the fields to the original state to allow for further stepping
        fields->restore_magnetic_fields();*/

        //open an hdf5 file with a reasonable name
        if (i % 4 == 0) {
        size_t n_written = make_dec_str(h5_fname+PREFIX_LEN, BUF_SIZE-PREFIX_LEN, fields->time(), n_digits_a, n_digits_b);
        meep::h5file* file = fields->open_h5file(h5_fname);

        fields->output_hdf5(meep::Ex, vol.surroundings(), file);
        fields->step();

        //we're done with the file
        delete file;
}
    }
    fclose(fp);
}
