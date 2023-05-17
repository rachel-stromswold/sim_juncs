#include "disp.hpp"

const double eps_cutoff = 10;

value um_to_l(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    if (f.n_args < 1) { er = E_LACK_TOKENS;return ret; }
    if (f.args[0].type != VAL_NUM) { er = E_BAD_TOKEN;return ret; }
    //lookup the length scale
    value l_per_um = c.lookup("l_per_um");
    if (l_per_um.type != VAL_NUM) { er = E_NOT_DEFINED;return ret; }
    return make_val_num( (l_per_um.val.x)*(f.args[0].val.x) );
}
value fs_to_t(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    if (f.n_args < 1) { er = E_LACK_TOKENS;return ret; }
    if (f.args[0].type != VAL_NUM) { er = E_BAD_TOKEN;return ret; }
    //lookup the length scale
    value l_per_um = c.lookup("l_per_um");
    if (l_per_um.type != VAL_NUM) { er = E_NOT_DEFINED;return ret; }
    return make_val_num( LIGHT_SPEED*(l_per_um.val.x)*(f.args[0].val.x) );
}

/**
 * Load variables specified in the parse_settings struct s into the context con
 * TODO: place this in context to handle this more elegantly
 */
context context_from_settings(const parse_settings& args) {
    context con;
    con.emplace("pi", make_val_num(M_PI));
    con.emplace("pml_thickness", make_val_num(args.pml_thickness));
    con.emplace("sim_length", make_val_num(args.len));
    con.emplace("length", make_val_num(2*args.pml_thickness + args.len));
    con.emplace("l_per_um", make_val_num(args.um_scale));
    value tmp_out = make_val_str(args.out_dir);
    con.emplace("out_dir", tmp_out);
    cleanup_val(&tmp_out);
    //helpful functions
    value tmp_f = make_val_func("um_to_l", 1, &um_to_l);
    con.emplace("um_to_l", tmp_f);
    cleanup_val(&tmp_f);
    tmp_f = make_val_func("fs_to_t", 1, &fs_to_t);
    con.emplace("fs_to_t", tmp_f);
    cleanup_val(&tmp_f);
    //add user defined options
    if (args.user_opts) {
	line_buffer lb(args.user_opts, ';');
	con.read_from_lines(lb);
    }
    return con;
}



/**
 * Helper function which initializes an array of smooth points
 */
void cgs_material_function::generate_smooth_pts(double smooth_rad, uint64_t seed) {
    //generate a seed and set up Mersenne twister
    std::seed_seq seeder{(uint32_t)(seed & 0x0000ffff), (uint32_t)((seed & 0xffff0000) >> 32)};
    std::mt19937 gen(seeder);

    std::normal_distribution<double> gaussian(0, smooth_rad);//negative values are fine since we include all reflections anyway
    std::uniform_real_distribution<double> unif(0.0,1.0);

    //avoid leaking old tables
    if (smooth_pts) free(smooth_pts);
    smooth_pts = NULL;

    //iteratively decrease the number of smoothing points until we can successfully allocate. Note that we malloc this array even though it's constant because its size is unknown at compile time.
    if (smooth_n > 0)
	smooth_pts = (double*)malloc(24*sizeof(double)*smooth_n);
    while (smooth_pts == NULL && smooth_n > 1) {
	smooth_n /= 2;
	//NOTE: size is multiplied by the number of spatial dimensions times 8 so that each point is reflected through the origin
	smooth_pts = (double*)malloc(24*sizeof(double)*smooth_n);
    }

    //we want to only use the exact point if there is only one sample
    if (smooth_pts) {
	//fill up the array
	for (_uint i = 0; i < smooth_n; ++i) {
	    _uint j = 0;

	    //generate a random direction on the sphere
	    double theta_inv = unif(gen);
	    double cos_theta = 1-2*theta_inv;
	    double sin_theta = 2*sqrt(theta_inv*(1-theta_inv));
	    double phi = 2*M_PI*unif(gen);
	    double r = gaussian(gen);

	    //convert into cartesian coordinates
	    double x = r*sin_theta*cos(phi);
	    double y = r*sin_theta*cos(phi);
	    double z = r*cos_theta;

	    //iterate over each of the eight possible reflections
	    for (int xf = -1; xf < 2; xf += 2) {
		for (int yf = -1; yf < 2; yf += 2) {
		    for (int zf = -1; zf < 2; zf += 2) {
			smooth_pts[24*i+3*j]   = x*xf;
			smooth_pts[24*i+3*j+1] = y*yf;
			smooth_pts[24*i+3*j+2] = z*zf;
			++j;
		    }
		}
	    }
	}
	//now that we're done, include all the reflections
	smooth_n *= 8;
    } else {
	smooth_n = 0;
    }
}

/**
 * This is identical to the simple_material_function, but each term is multiplied by a constant scalar
 */
cgs_material_function::cgs_material_function(double p_def_ret, _uint p_smooth_n, double p_smooth_rad) {
    def_ret = p_def_ret;
    n_regions = 0;
    regions = NULL;

    //initialize smoothing
    smooth_n = p_smooth_n;
    //TODO: find a way to read arguments for the seed instead of using a constant at compile time. Note, we don't care about entropy or statistical quality. We mostly care about getting arbitrary as opposed to random values.
    generate_smooth_pts(p_smooth_rad, DEF_SEED);
}

/**
 * This is identical to the simple_material_function, but each term is multiplied by a constant scalar
 */
cgs_material_function::cgs_material_function(composite_object* p_volume, std::string type, double p_def_ret, _uint p_smooth_n, double p_smooth_rad) {
    //set the regions which are used
    n_regions = 1;
    regions = (region_scale_pair*)malloc(sizeof(region_scale_pair));
    if (regions) {
	double scale = def_ret;
	//lookup the scaling constant by reading metadata from the object
	if (p_volume->has_metadata(type)) {
	    value tmp_val = p_volume->fetch_metadata(type);
	    if (tmp_val.get_type() == VAL_NUM) scale = tmp_val.to_float();
	}
	regions[0].s = scale;//by default set the scale to whatever the default return value is
	regions[0].c = p_volume;
    } else {
	n_regions = 0;
    }
    //set default return value
    def_ret = p_def_ret;

    //initialize smoothing
    smooth_n = p_smooth_n;
    generate_smooth_pts(p_smooth_rad, DEF_SEED);
}

/**
 * This is identical to the simple_material_function, but each term is multiplied by a constant scalar
 */
cgs_material_function::cgs_material_function(region_scale_pair p_reg, double p_def_ret, _uint p_smooth_n, double p_smooth_rad) {
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

    //initialize smoothing
    smooth_n = p_smooth_n;
    generate_smooth_pts(p_smooth_rad, DEF_SEED);
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
    //copy smoothing data points
    smooth_n = o.smooth_n;
    smooth_pts = (double*)malloc(3*sizeof(double)*smooth_n);
    if (smooth_pts) {
	for (_uint i = 0; i < 3*smooth_n; ++i) smooth_pts[i] = o.smooth_pts[i];
    }
}

//move constructor
cgs_material_function::cgs_material_function(cgs_material_function&& o) {
    def_ret = o.def_ret;
    n_regions = o.n_regions;
    regions = o.regions;
    smooth_pts = o.smooth_pts;
    smooth_n = o.smooth_n;
    o.regions = NULL;
    o.n_regions = 0;
    o.smooth_pts = NULL;
    o.smooth_n = 0;
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
    //swap smoothing point data
    double* tmp_smooth_pts = smooth_pts;
    smooth_pts = o.smooth_pts;
    o.smooth_pts = tmp_smooth_pts;
    _uint tmp_smooth_n = smooth_n;
    smooth_n = o.smooth_n;
    o.smooth_n = tmp_smooth_n;

    return *this;
}

cgs_material_function::~cgs_material_function() {
    if (regions) free(regions);
    regions = NULL;
    n_regions = 0;
    if (smooth_pts) free(smooth_pts);
    smooth_pts = NULL;
    smooth_n = 0;
}

/**
 * Add the region based on the composite_object p_reg.
 * p_reg: the composite object specifying the region to add
 * type: lookup the matching metadata key from the composite object and set the scale based on its value.
 */
void cgs_material_function::add_region(composite_object* p_reg, std::string type) {
    _uint new_n_regions = n_regions+1;
    region_scale_pair* tmp_regs = (region_scale_pair*)realloc(regions, sizeof(region_scale_pair)*new_n_regions);
    //only proceed if allocation of memory was successful
    if (tmp_regs) {
	n_regions = new_n_regions;
	regions = tmp_regs;
	regions[n_regions-1].c = p_reg;
	double scale = def_ret;
	//lookup the scaling constant by reading metadata from the object
	if (p_reg->has_metadata(type)) {
	    value tmp_val = p_reg->fetch_metadata(type);
	    if (tmp_val.get_type() == VAL_NUM) scale = tmp_val.to_float();
	}
	regions[n_regions-1].s = scale;//by default set the scale to whatever the default return value is
    }
}

//returns whether we're to the left or the right of the dielectric boundary 1 indicates we are, 0 that we are not
double cgs_material_function::in_bound(const meep::vec &r) {
    if (!regions) return def_ret;

    //we stochastically smooth boundaries by taking smooth_n samples from a 3d Gaussian distribution with variance smooth_rad^2 centered around r
    double r_x = r.x();double r_y = r.y();double r_z = r.z();

    //look through each region and return the first that contains the vector
    double ret = 0;
    for (_uint i = 0; i < n_regions; ++i) {
	//initialize such that we always include the origin
	double this_ret = regions[i].c->in(vec3(r_x, r_y, r_z));
	if (smooth_pts) {
	    for (_uint j = 0; j < smooth_n; ++j) {
		this_ret += regions[i].c->in(vec3(r_x + smooth_pts[3*j], r_y + smooth_pts[3*j+1], r_z + smooth_pts[3*j+2]));
	    }
	}
	ret += def_ret + (regions[i].s - def_ret)*this_ret/(smooth_n+1);
    }
    return ret;
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

source_info::source_info(value info, parse_ercode& ercode) {
    ercode = E_SUCCESS;
    type = SRC_GAUSSIAN;
    component = meep::Ex;
    wavelen = 0.7;
    width = 1;
    phase = 0;
    start_time = 0;
    end_time = 0;
    amplitude = 1.0;
    parse_ercode tmp_er = E_SUCCESS;

    bool is_gauss = is_type(info, "Gaussian_source");
    bool is_contin = is_type(info, "CW_source");
    if (is_contin || is_gauss) {
	value tmp = info.val.c->lookup("component");
	if (tmp.type == VAL_NUM) {
	    switch ((basis_comp_vectors)(tmp.val.x)) {
		case C_EY: component = meep::Ey;break;
		case C_EZ: component = meep::Ez;break;
		case C_HX: component = meep::Hx;break;
		case C_HY: component = meep::Hy;break;
		case C_HZ: component = meep::Hz;break;
		default:break;
	    }
	}
	tmp = info.val.c->lookup("wavelength");
	if (tmp.type == VAL_NUM)
	    wavelen = tmp.val.x;
	tmp = info.val.c->lookup("amplitude");
	if (tmp.type == VAL_NUM)
	    amplitude = tmp.val.x;
	tmp = info.val.c->lookup("start_time");
	if (tmp.type == VAL_NUM)
	    start_time = tmp.val.x;
	tmp = info.val.c->lookup("end_time");
	if (tmp.type == VAL_NUM)
	    end_time = tmp.val.x;
	tmp = info.val.c->lookup("slowness");
	if (tmp.type == VAL_NUM)
	    width = tmp.val.x;
    } else {
	ercode = E_BAD_TYPE;
    }
    if (is_gauss) {
	value tmp = info.val.c->lookup("width");
	if (tmp.type == VAL_NUM)
	    width = tmp.val.x;
	tmp = info.val.c->lookup("phase");
	if (tmp.type == VAL_NUM)
	    phase = tmp.val.x;
	double cutoff = 5;
	tmp = info.val.c->lookup("cutoff");
	if (tmp.type == VAL_NUM)
	    cutoff = 6;
	end_time = start_time + 2*cutoff*width;
    }
    if (is_contin) type = SRC_CONTINUOUS;
}

/*source_info::source_info(value info, parse_ercode* ercode) {
    if (ercode) *ercode = E_SUCCESS;
    //initialize default values
    type = SRC_GAUSSIAN;
    component = meep::Ex;
    wavelen = 0.7;
    width = 1;
    phase = 0;
    start_time = 0;
    end_time = 0;
    amplitude = 1.0;
    parse_ercode tmp_er = E_SUCCESS;
    if (info.type == VAL_LIST) {
	if (info.n_els > 3) {
	    //read the component
	    if (info.val.l[1].type == VAL_STR && info.val.l[1].val.s) {
		if (strcmp(info.val.l[1].val.s, "Ex") == 0) {
		    component = meep::Ex;
		} else if (strcmp(info.val.l[1].val.s, "Ey") == 0) {
		    component = meep::Ey;
		} else if (strcmp(info.val.l[1].val.s, "Ez") == 0) {
		    component = meep::Ez;
		} else if (strcmp(info.val.l[1].val.s, "Hx") == 0) {
		    component = meep::Hx;
		} else if (strcmp(info.val.l[1].val.s, "Hy") == 0) {
		    component = meep::Hy;
		} else if (strcmp(info.val.l[1].val.s, "Hz") == 0) {
		    component = meep::Hz;
		}
	    } else {
		tmp_er = E_BAD_TOKEN;
	    }
	    //read the wavelength
	    if (info.val.l[2].get_type() == VAL_NUM) {
		wavelen = info.val.l[2].get_val().x;
	    } else {
		tmp_er = E_BAD_VALUE;
	    }
	    //read the type
	    if (tmp_er == E_SUCCESS) {
		//figure out what to do depending on what type of pulse envelope this is
		if (info.val.l[0].type == VAL_STR && (strcmp(info.val.l[0].val.s, "gaussian") == 0 || strcmp(info.val.l[0].val.s, "Gaussian") == 0)) {
		    type = SRC_GAUSSIAN;
		    //set default values for the envelope
		    start_time = 0;
		    double cutoff = 5;
		    if (info.n_els < 4) {
			tmp_er = E_LACK_TOKENS;
		    } else {
			if (info.val.l[3].get_type() == VAL_NUM)
			    width = info.val.l[3].get_val().x;
			else
			    tmp_er = E_BAD_TOKEN;
			//set default values
			amplitude = 1;
			start_time = 0;
			cutoff = DEFAULT_WIDTH_N;
			//read the phase if specified
			if (info.n_els > 4) {
			    if (info.val.l[4].get_type() == VAL_NUM)
				phase = info.val.l[4].get_val().x;
			    else
				tmp_er = E_BAD_TOKEN;
			    if (errno) tmp_er = E_BAD_TOKEN;
			}
			//read the start time if supplied
			if (info.n_els > 5) {
			    if (info.val.l[5].get_type() == VAL_NUM)
				start_time = info.val.l[5].get_val().x;
			    else
				tmp_er = E_BAD_TOKEN;
			    if (errno) tmp_er = E_BAD_TOKEN;
			}
			//read the cutoff if supplied
			if (info.n_els > 6) {
			    if (info.val.l[6].get_type() == VAL_NUM)
				cutoff = info.val.l[6].get_val().x;
			    else
				tmp_er = E_BAD_TOKEN;
			    if (errno) tmp_er = E_BAD_TOKEN;
			}
			//read the amplitude if supplied
			if (info.n_els > 7) {
			    if (info.val.l[7].get_type() == VAL_NUM)
				amplitude = info.val.l[7].get_val().x;
			    else
				tmp_er = E_BAD_TOKEN;
			    if (errno) tmp_er = E_BAD_TOKEN;
			}
			end_time = start_time + 2*cutoff*width;
		    }
		} else if (strcmp(info.val.l[0].val.s, "continuous") == 0) {
		    type = SRC_CONTINUOUS;
		    phase = 0;
		    if (info.n_els < 4) {
			tmp_er = E_LACK_TOKENS;
		    } else {
			if (info.val.l[3].get_type() == VAL_NUM)
			    start_time = info.val.l[3].get_val().x;
			else
			    tmp_er = E_BAD_TOKEN;
			if (errno) tmp_er = E_BAD_TOKEN;
			//read the end time if supplied
			if (info.n_els > 4) {
			    if (info.val.l[4].get_type() == VAL_NUM)
				end_time = info.val.l[4].get_val().x;
			    else
				tmp_er = E_BAD_TOKEN;
			    if (errno) tmp_er = E_BAD_TOKEN;
			}
			//read the width if supplied
			if (info.n_els > 5) {
			    if (info.val.l[5].get_type() == VAL_NUM)
				width = info.val.l[5].get_val().x;
			    else
				tmp_er = E_BAD_TOKEN;
			    if (errno) tmp_er = E_BAD_TOKEN;
			}
			//read the amplitude if supplied
			if (info.n_els > 6) {
			    if (info.val.l[6].get_type() == VAL_NUM)
				amplitude = info.val.l[6].get_val().x;
			    else
				tmp_er = E_BAD_TOKEN;
			    if (errno) tmp_er = E_BAD_TOKEN;
			}
		    }
		}
	    }
	} else {
	    tmp_er = E_LACK_TOKENS;
	}
    }

    //set the error code if there was one and the caller wants it
    if (ercode) *ercode = tmp_er;
}*/

gaussian_src_time_phase::gaussian_src_time_phase(double f, double w, double phase, double st, double et) {
    omega = 2*M_PI*f;
    width = w;
    phi = phase+M_PI;
    peak_time = 0.5*(st+et);
    cutoff = (et-st)*0.5;
    //fwidth = gaussian_bandwidth(width);
    // correction factor so that current amplitude (= d(dipole)/dt) is
    // ~ 1 near the peak of the Gaussian.
    amp = 1.0 / std::complex<double>(0, -omega);

    // this is to make last_source_time as small as possible
    while (exp(-cutoff * cutoff / (2 * width * width)) < 1e-100)
        cutoff *= 0.9;
    cutoff = float(cutoff); // don't make cutoff sensitive to roundoff error
}

std::complex<double> gaussian_src_time_phase::dipole(double time) const {
    double tt = time - peak_time;
    if (float(fabs(tt)) > cutoff) return 0.0;

    return exp(-tt*tt / (2*width*width)) * std::polar(1.0, -omega*tt - phi) * amp;
}

bool gaussian_src_time_phase::is_equal(const src_time &t) const {
    const gaussian_src_time_phase *tp = dynamic_cast<const gaussian_src_time_phase *>(&t);
    if (tp)
        return (tp->omega == omega && tp->width == width && tp->peak_time == peak_time && tp->cutoff == cutoff);
    else
        return 0;
}

/**
 * Add the list of susceptibilities to the parse_settings file s. The string should have the format (omega_0,gamma_0,sigma_0),(omega_1,gamma_1,sigma_1),...
 * returns: a vector list of susceptibilities
 * If an error is encountered, a code will be saved to er if it is not NULL
 *  0 on success
 *  -1 not a list or insufficient arguments
 *  -2 bad argument supplied
 */
std::vector<drude_suscept> bound_geom::parse_susceptibilities(value val, int* er) {
    std::vector<drude_suscept> ret;
    if (val.type == VAL_LIST) {
	drude_suscept cur_sus;
	for (size_t i = 0; i < val.n_els; ++i) {
	    value cur = val.val.l[i];
	    if (cur.type != VAL_LIST || cur.n_els < 2) {
		if (er) *er = -1;
		return ret;
	    }
	    if (cur.val.l[0].type != VAL_NUM || cur.val.l[1].type != VAL_NUM || cur.val.l[2].type != VAL_NUM) {
		if (er) *er = -2;
		return ret;
	    }
	    cur_sus.omega_0 = cur.val.l[0].val.x;
	    cur_sus.gamma = cur.val.l[1].val.x;
	    cur_sus.sigma = cur.val.l[2].val.x;
	    cur_sus.use_denom = 1;
	    //read whether this is a Drude or Lorentz susceptibility
	    if (cur.n_els > 3) {
		if (cur.val.l[3].type == VAL_NUM) {
		    cur_sus.use_denom = (cur.val.l[3].val.x != 0.0);
		} else if (cur.val.l[3].type == VAL_STR) {
		    char* tok = cur.val.l[3].val.s;
		    if (strcmp(tok, "drude") == 0 || strcmp(tok, "true") == 0) cur_sus.use_denom = 0;
		    if (strcmp(tok, "lorentz") == 0 || strcmp(tok, "false") == 0) cur_sus.use_denom = 1;
		} else {
		    if (er) *er = -2;
		    return ret;
		}
	    }
	    ret.push_back(cur_sus);
	}
    }
    if (er) *er = 0;
    return ret;
}

/**
 * Read a composite_object specifying monitor locations into a list of monitor locations
 * returns: an error code if one was encountered or E_SUCCESS
 */
parse_ercode bound_geom::parse_monitors(value vl) {
    if (vl.type == VAL_LIST) {
	parse_ercode er = E_SUCCESS;
	for (size_t i = 0; i < vl.n_els; ++i) {
	    //see if we can cast to a vector and check for errors
	    value vec_cast = vl.val.l[i].cast_to(VAL_3VEC, er);
	    if (er != E_SUCCESS) { cleanup_val(&vec_cast);return er; }
	    //convert to a meep vector and append
	    meep::vec tmp_vec(vec_cast.val.v->x(), vec_cast.val.v->y(), vec_cast.val.v->z());
	    monitor_locs.push_back(tmp_vec);
	    cleanup_val(&vec_cast);
	}
	return er;
    }
    return E_BAD_TYPE;
}

double dummy_eps(const meep::vec& r) { (void)r;return 1.0; }

/**
 * This is a helper function for the bound_geom constructor. Meep doesn't implement copy or move constructors so we have to initialize the structure immediately so that the fields can be initialized in turn.
 */
meep::structure* bound_geom::structure_from_settings(const parse_settings& s, scene& problem) {
    pml_thickness = s.pml_thickness;
    len = s.len;

    //read information about the problem size from the geometry file
    std::vector<composite_object*> data = problem.get_data();
    for (size_t i = 0; i < data.size(); ++i) {
	if (data[i]->has_metadata("pml_thickness")) {
	    value tmp_val = data[i]->fetch_metadata("pml_thickness");
	    if (tmp_val.get_type() == VAL_NUM) len = tmp_val.val.x;
	}
	if (data[i]->has_metadata("length")) {
	    value tmp_val = data[i]->fetch_metadata("length");
	    if (tmp_val.get_type() == VAL_NUM) len = tmp_val.val.x;
	}
    }

    //the arguments supplied will alter the location of the dielectric
    z_center = len/2 + pml_thickness;
    eps_scale = 1 / (sharpness*s.resolution);

    //initialize the volume
    if (s.n_dims == 1) {
	vol = meep::vol1d(2*z_center, s.resolution);
    } else if (s.n_dims == 2) {
	vol = meep::vol2d(2*z_center, 2*z_center, s.resolution);
    } else {
	vol = meep::vol3d(2*z_center, 2*z_center, 2*z_center, s.resolution);
    }

    //iterate over all objects specified in the scene
    std::vector<composite_object*> roots = problem.get_roots();
    //setup the structure with the infinite frequency dielectric component
    cgs_material_function inf_eps_func(s.ambient_eps, s.smooth_n, s.smooth_rad);

    std::vector<double> thicknesses(roots.size());
    for (size_t i = 0; i < roots.size(); ++i) {
	//we need to adjust the thickness and conductivity of the sample if it is 2D
	thicknesses[i] = 1.0;
	if (roots[i]->has_metadata("make_2d") && roots[i]->fetch_metadata("make_2d").val.x != 0) {
	    thicknesses[i] = THICK_SCALE / s.resolution;
	    roots[i]->rescale(vec3(1.0, 1.0, thicknesses[i]));
	}
	inf_eps_func.add_region(roots[i]);
    }
    meep::structure* strct = new meep::structure(vol, inf_eps_func, meep::pml(s.pml_thickness));
    //read susceptibilities if they are available
    for (size_t i = 0; i < roots.size(); ++i) {
	std::vector<drude_suscept> cur_sups;
	int res = 0;
	if (roots[i]->has_metadata("susceptibilities")) {
	    value sup_val = roots[i]->fetch_metadata("susceptibilities");
	    cur_sups = parse_susceptibilities(sup_val, &res);
	}
	//add frequency dependent susceptibility
	for (_uint j = 0; j < cur_sups.size(); ++j) {
	    //create the susceptibility, accounting for the scale factor correction
	    double omega_0 = cur_sups[j].omega_0 / s.um_scale;
	    double gamma = cur_sups[j].gamma / s.um_scale;
	    double sigma = cur_sups[j].sigma / thicknesses[i];
	    meep::lorentzian_susceptibility suscept( omega_0, gamma, !(cur_sups[j].use_denom) );
	    region_scale_pair tmp_pair = {roots[i], sigma};
	    //add the susceptibility to the appropriate region
	    cgs_material_function scale_func(tmp_pair, 0.0, s.smooth_n, s.smooth_rad);
	    strct->add_susceptibility(scale_func, meep::E_stuff, suscept);
	}
    }
    return strct;
}

/**
 * Constructor for the bound_geom file. parse_settings are read from s.
 * s: parse_settings object to read.
 * ercode: if an error occurs while parsing the .geom file and ercode is not NULL, a code for the error is saved there
 */
bound_geom::bound_geom(const parse_settings& s, parse_ercode* ercode) :
    problem(s.geom_fname, context_from_settings(s), ercode),
    strct(structure_from_settings(s, problem)),
    fields(strct)
{
    if (*ercode != E_SUCCESS) {
        printf("Scene parsing failed, exiting.\n");
        exit(1);
    }
    //TODO: actually fix the memory issues
    for (size_t i = 0; i < FIELDS_PAD_SIZE; ++i) dummy_vals[i] = 0;
    um_scale = s.um_scale;
    if (ercode) *ercode = E_SUCCESS;

    //we have to kludge it to get around the very f** annoying fact that meep doesn't have default constructors for fields, just read the structure_from_settings comment
    double z_center = s.len/2 + s.pml_thickness;
    printf("using simulation side length %f, resolution %f\n", s.len, s.resolution);
    if (s.n_dims == 1) {
	vol = meep::vol1d(2*z_center, s.resolution);
    } else if (s.n_dims == 2) {
	vol = meep::vol2d(2*z_center, 2*z_center, s.resolution);
    } else {
	vol = meep::vol3d(2*z_center, 2*z_center, 2*z_center, s.resolution);
    }

    post_source_t = s.post_source_t;

    //add fields specified in the problem
    context& c = problem.get_context();
    for (size_t i = c.size(); i > 0; --i) {
	parse_ercode tmp_er;
	value inst = c.peek_val(i);
	source_info cur_info(inst, tmp_er);
	//only inspect instances
	if (tmp_er == E_SUCCESS) {
	    value reg = inst.val.c->lookup("region");
	    parse_ercode tmp_er = E_BAD_TYPE;
	    if (is_type(reg, "Box")) {
		tmp_er = E_SUCCESS;
		parse_ercode tmp_er2 = E_SUCCESS;
		value p1 = reg.val.c->lookup("pt_1").cast_to(VAL_3VEC, tmp_er);
		value p2 = reg.val.c->lookup("pt_2").cast_to(VAL_3VEC, tmp_er2);
		if (tmp_er == E_SUCCESS && tmp_er2 == E_SUCCESS) {
		    double x_0 = p1.val.v->x();double x_1 = p2.val.v->x();
		    double y_0 = p1.val.v->y();double y_1 = p2.val.v->y();
		    double z_0 = p1.val.v->z();double z_1 = p2.val.v->z();
		    meep::volume source_vol(meep::vec(x_0, y_0, z_0), meep::vec(x_1, y_1, z_1));
		    double c_by_a = 0.299792458*s.um_scale;
		    //We specify the width of the pulse in units of the oscillation period
		    double frequency = 1 / (cur_info.wavelen*s.um_scale);
		    double width = cur_info.width*c_by_a;
		    double start_time = cur_info.start_time*c_by_a;
		    double end_time = cur_info.end_time*c_by_a;
		    if (cur_info.type == SRC_GAUSSIAN) {
			printf("Adding Gaussian envelope: f=%f, w=%f, t_0=%f, t_f=%f (meep units)\n",
				frequency, width, start_time, end_time);
			gaussian_src_time_phase src(frequency, width, cur_info.phase, start_time, end_time);
			fields.add_volume_source(cur_info.component, src, source_vol, cur_info.amplitude);
		    } else if (cur_info.type == SRC_CONTINUOUS) { 
			printf("Adding continuous wave: f=%f, w=%f, t_0=%f, t_f=%f (meep units)\n",
				frequency, width, start_time, end_time);
			meep::continuous_src_time src(frequency, width, start_time, end_time);
			fields.add_volume_source(cur_info.component, src, source_vol, cur_info.amplitude);
		    }
#ifdef DEBUG_INFO
		    sources.push_back(cur_info);
#endif
		    //set the total timespan based on the added source
		    ttot = fields.last_source_time() + post_source_t*0.299792458*s.um_scale;
		}
		cleanup_val(&p1);
		cleanup_val(&p2);
	    }
	} else {
	    if (is_type(inst, "monitor")) {
		value vl = inst.val.c->lookup("locations");
		parse_ercode tmp_er = parse_monitors(vl);
		if (tmp_er != E_SUCCESS)
		    printf("warning: Invalid monitor location encountered! (all monitors must be vectors or lists with three elements)\n");
		monitor_clusters.push_back(monitor_locs.size());
	    }
	}
    }

    std::vector<composite_object*> data = problem.get_data();
    /*for (size_t i = 0; i < data.size(); ++i) {
	if (data[i]->has_metadata("type")) {
        value type = data[i]->fetch_metadata("type");
	    if (type.type == VAL_STR && strcmp(type.val.s, "field_source") == 0) {
		//figure out the volume for the field source
		const object* l_child = data[i]->get_child_l();
		const object* r_child = data[i]->get_child_r();
		object_type l_type = data[i]->get_child_type_l();
		//WLOG fix the left child to be the one we care about
		if (r_child && !l_child) {
		    l_child = r_child;
		    l_type = data[i]->get_child_type_r();
		}
		//make sure that the object is a box so that we can figure out the corner and the offset
		if (l_type == CGS_BOX) {
		    vec3 center = ((box*)l_child)->get_center();
		    vec3 offset = ((box*)l_child)->get_offset();
		    double x_0 = center.x() - offset.x();double x_1 = center.x() + offset.x();
		    double y_0 = center.y() - offset.y();double y_1 = center.y() + offset.y();
		    double z_0 = center.z() - offset.z();double z_1 = center.z() + offset.z();
		    meep::volume source_vol(meep::vec(x_0, y_0, z_0), meep::vec(x_1, y_1, z_1));

		    //only continue if a shape for the pulse was specified
		    if (data[i]->has_metadata("envelope")) {
			*ercode = E_SUCCESS;
			source_info cur_info(data[i]->fetch_metadata("envelope"), ercode);
			//create the EM-wave source at the specified location only if everything was read successfully
			if (*ercode == E_SUCCESS) {
			    double c_by_a = 0.299792458*s.um_scale;
			    //We specify the width of the pulse in units of the oscillation period
			    double frequency = 1 / (cur_info.wavelen*s.um_scale);
			    double width = cur_info.width*c_by_a;
			    double start_time = cur_info.start_time*c_by_a;
			    double end_time = cur_info.end_time*c_by_a;
			    if (cur_info.type == SRC_GAUSSIAN) {
				printf("Adding Gaussian envelope: f=%f, w=%f, t_0=%f, t_f=%f (meep units)\n",
					frequency, width, start_time, end_time);
				gaussian_src_time_phase src(frequency, width, cur_info.phase, start_time, end_time);
				fields.add_volume_source(cur_info.component, src, source_vol, cur_info.amplitude);
			    } else if (cur_info.type == SRC_CONTINUOUS) { 
				printf("Adding continuous wave: f=%f, w=%f, t_0=%f, t_f=%f (meep units)\n",
					frequency, width, start_time, end_time);
				meep::continuous_src_time src(frequency, width, start_time, end_time);
				fields.add_volume_source(cur_info.component, src, source_vol, cur_info.amplitude);
			    }
#ifdef DEBUG_INFO
			    sources.push_back(cur_info);
#endif
			}
		    }

		    //set the total timespan based on the added source
		    ttot = fields.last_source_time() + post_source_t*0.299792458*s.um_scale;
		} else {
		    printf("Error: only boxes are currently supported for field volumes");
		    if (ercode) *ercode = E_BAD_VALUE;
		}
	    } else if (type.type == VAL_STR && strcmp(type.val.s, "monitor") == 0) {
		parse_ercode tmp_er = parse_monitors(data[i]);
		if (tmp_er != E_SUCCESS) printf("warning: Invalid monitor location encountered! (all monitors must be vectors or lists with three elements)\n");
		monitor_clusters.push_back(monitor_locs.size());
	    }
	}
    }*/
    write_settings(stdout);
    dump_span = s.field_dump_span;
    dump_raw = s.dump_raw;
}

bound_geom::~bound_geom() {
    if (strct) delete strct;
}

/**
 * Add a point source specified by src at the location source_loc.
 * c: the component of the source (e.g. Ex, Hz, etc.)
 * src: the source, see meep's documentation
 * source_loc: the location of the source
 * amp: the amplitude of the source
 */
void bound_geom::add_point_source(meep::component c, const meep::src_time &src, const meep::vec& source_loc, std::complex<double> amp) {
    fields.add_point_source(c, src, source_loc, amp);

    //set the total timespan based on the added source
    ttot = fields.last_source_time() + post_source_t;
}

/**
 * Add a volume source specified by src in source_vol.
 * c: the component of the source (e.g. Ex, Hz, etc.)
 * src: the source, see meep's documentation
 * vol: the region for the source
 * amp: the amplitude of the source
 */
void bound_geom::add_volume_source(meep::component c, const meep::src_time &src, const box& vol, std::complex<double> amp) {
    vec3 center = vol.get_center();
    vec3 offset = vol.get_offset();
    double x_0 = center.x() - offset.x();double x_1 = center.x() + offset.x();
    double y_0 = center.y() - offset.y();double y_1 = center.y() + offset.y();
    double z_0 = center.z() - offset.z();double z_1 = center.z() + offset.z();
    meep::volume source_vol(meep::vec(x_0, y_0, z_0), meep::vec(x_1, y_1, z_1));

    fields.add_volume_source(c, src, source_vol, amp);

    //set the total timespan based on the added source
    ttot = fields.last_source_time() + post_source_t;
}

/**
 * Run the simulation
 * fname_prefix: all files are saved to fname_prefix if specified
 */
void bound_geom::run(const char* fname_prefix) {
    fields.set_output_directory(fname_prefix);
    char h5_fname[BUF_SIZE];

    //save the dielectric used
    printf("Set output directory to %s\n", fname_prefix);
    fields.output_hdf5(meep::Dielectric, fields.total_volume());

    //avoid divisions by zero by defaulting to 1
    if (dump_span == 0) dump_span = 1;

    //make sure the time series corresponding to each monitor point is long enough to hold all of its information
    size_t n_locs = monitor_locs.size();
    field_times.resize(n_locs);
    n_t_pts = (_uint)( (ttot+fields.dt/2) / fields.dt );
    for (_uint j = 0; j < n_locs; ++j) {
        field_times[j].reserve( 1+ (n_t_pts/dump_span) );
        //make_data_arr(&(field_times[j]), n_t_pts/dump_span);
    }

    //figure out the number of digits before the decimal and after
    int n_digits_a = (int)(ceil(log(ttot)/log(10)));
    double rat = -log((double)(fields.dt))/log(10.0);
    int n_digits_b = (int)ceil(rat)+1;
    strcpy(h5_fname, "ex-");
    //run the simulation
    printf("starting simulations\n");
    _uint i = 0;
    try {
        for (; i < n_t_pts; ++i) {
            //write each of the monitor locations, but only if we've reached a savepoint
            if (i % dump_span == 0) {
                //fetch monitor points
                for (_uint j = 0; j < n_locs; ++j) {
                    std::complex<double> val = fields.get_field(meep::Ex, monitor_locs[j]);
                    complex tmp(val.real(), val.imag());
                    field_times[j].push_back(tmp);
                    if (tmp.re > 1000) {
                        printf("divergence in run at (i,j)=(%d,%d) (%f)\n", i, j, tmp.re);
                    }
                }
                //save the raw hdf5 files if requested
                if (dump_raw) {
                    make_dec_str(h5_fname+PREFIX_LEN, BUF_SIZE-PREFIX_LEN, fields.time(), n_digits_a, n_digits_b);
                    meep::h5file* file = fields.open_h5file(h5_fname);
                    fields.output_hdf5(meep::Ex, vol.surroundings(), file);
                    delete file;
                }
                printf("    %d%% complete\n", 100*i/n_t_pts);
            }
            fields.step();
        }
    } catch (std::runtime_error& err) {
        std::cout << "error on step " << i << ": " << err.what() << std::endl;
    }

    printf("values written to dummies:\n\t");
    for (size_t i = 0; i < FIELDS_PAD_SIZE; ++i) { printf("%x ", dummy_vals[i]); }
    printf("\nSimulations completed\n");
}

/**
 * Save time and wavelenuency domain data for each monitor_loc
 * fname: filename to use for saving information
 * save_pts: an array of indices of points that the user wants to save. If NULL (default), then all points are saved.
 * n_locs: the size of the array pointed to by save_pts
 * returns -1 on error, 0 otherwise
 */
void bound_geom::save_field_times(const char* fname_prefix) {
    char out_name[SMALL_BUF_SIZE];

    size_t n_locs = monitor_locs.size();
    //create the field type and specify members
    H5::CompType fieldtype(sizeof(complex));

    hid_t float_member_id = H5_float_type.getId();
    herr_t ret_val = H5Tinsert(fieldtype.getId(), "Re", HOFFSET(complex, re), float_member_id);
    if (ret_val < 0) return;
    ret_val = H5Tinsert(fieldtype.getId(), "Im", HOFFSET(complex, im), float_member_id);
    if (ret_val < 0) return;
    //use the space of rank 1 tensors with a dimension of n_t_pts
    hsize_t t_dim[1];
    t_dim[0] = {n_t_pts/dump_span};
    hsize_t f_dim[1];
    H5::DataSpace t_space(1, t_dim);

    //create the location type and specify members
    H5::CompType loctype(sizeof(sto_vec));
    ret_val = H5Tinsert(loctype.getId(), "x", HOFFSET(sto_vec, x), float_member_id);
    ret_val = H5Tinsert(loctype.getId(), "y", HOFFSET(sto_vec, y), float_member_id);
    ret_val = H5Tinsert(loctype.getId(), "z", HOFFSET(sto_vec, z), float_member_id);
    //use the space of rank 1 tensors with dimension of the number of monitor points
    hsize_t l_dim[1];
    l_dim[0] = {n_locs};
    H5::DataSpace l_space(1, l_dim);

    //create the source type and specify members
    H5::CompType srctype(sizeof(source_info));
    ret_val = H5Tinsert(srctype.getId(), "wavelen", HOFFSET(source_info, wavelen), float_member_id);
    ret_val = H5Tinsert(srctype.getId(), "width", HOFFSET(source_info, width), float_member_id);
    ret_val = H5Tinsert(srctype.getId(), "phase", HOFFSET(source_info, phase), float_member_id);
    ret_val = H5Tinsert(srctype.getId(), "start_time", HOFFSET(source_info, start_time), float_member_id);
    ret_val = H5Tinsert(srctype.getId(), "end_time", HOFFSET(source_info, end_time), float_member_id);
    ret_val = H5Tinsert(srctype.getId(), "amplitude", HOFFSET(source_info, amplitude), float_member_id);

    //take the fourier transform for each point
    std::vector<data_arr> tdom_fields(n_locs);
    std::vector<data_arr> fdom_fields(n_locs);
    for (_uint j = 0; j < n_locs; ++j) {
        //copy the time domain data
        make_data_arr(&(tdom_fields[j]), field_times[j].size());
        for (size_t i = 0; i < field_times[j].size(); ++i) tdom_fields[j].buf[i] = field_times[j][i];
	tdom_fields[j].size = field_times[j].size();
        //take the fourier transform
        fdom_fields[j] = fft(tdom_fields[j]);
        f_dim[0] = fdom_fields[j].size;
    }
    H5::DataSpace f_space(1, f_dim);

    //open the file which will store the fields as a function of time
    snprintf(out_name, BUF_SIZE, "%s/field_samples.h5", fname_prefix);
    printf("saving field output to %s\n", out_name);
    H5::H5File file(out_name, H5F_ACC_TRUNC);

    //save metadata
    H5::Group info_group = file.createGroup("info");
    //we need to convert meep vectors to sto_vecs
    sto_vec* tmp_vecs = (sto_vec*)malloc(sizeof(sto_vec)*n_locs);
    for (_uint i = 0; i < n_locs; ++i) {
        tmp_vecs[i].x = (_ftype)(monitor_locs[i].x());
        tmp_vecs[i].y = (_ftype)(monitor_locs[i].y());
        tmp_vecs[i].z = (_ftype)(monitor_locs[i].z());
    }
    //write the start and end times for each simulation
    hsize_t t_info_dim[1];
    t_info_dim[0] = {3};
    H5::DataSpace t_info_space(1, t_info_dim);
    H5::DataSet t_info_dataset(info_group.createDataSet("time_bounds", H5_float_type, t_info_space));
    //set the start and end times for the simulation TODO: is the first boundary necessary or can it be left implicit?
    _ftype time_boundaries[3];
    time_boundaries[0] = 0.0;time_boundaries[1] = meep_time_to_fs(ttot);            //start and end of simulation times
    time_boundaries[2] = (_ftype)(time_boundaries[1]-time_boundaries[0])*dump_span/n_t_pts;   //step size
    t_info_dataset.write(time_boundaries, H5_float_type);
    //write the number of clusters
    hsize_t n_info_dim[1];
    n_info_dim[0] = {1};
    H5::DataSpace n_info_space(1, n_info_dim);
    H5::DataSet n_c_info_dataset(info_group.createDataSet("n_clusters", H5::PredType::NATIVE_HSIZE, n_info_space));
    hsize_t hsize_sto = monitor_clusters.size();
    n_c_info_dataset.write(&hsize_sto, H5::PredType::NATIVE_HSIZE);
    //write the number of time points
    H5::DataSet n_t_info_dataset(info_group.createDataSet("n_time_points", H5::PredType::NATIVE_HSIZE, n_info_space));
    hsize_sto = n_t_pts/dump_span;
    n_t_info_dataset.write(&hsize_sto, H5::PredType::NATIVE_HSIZE);
    //write information about sources
    size_t n_srcs = sources.size();
    source_info* tmp_srcs = (source_info*)malloc(sizeof(source_info)*n_srcs);
    hsize_t src_dim[1];
    src_dim[0] = {n_srcs};
    H5::DataSpace src_space(1, src_dim);
    H5::DataSet src_dataset(info_group.createDataSet("sources", srctype, src_space));
    //copy the sources into the buffer we just allocated and convert to fs
    for (_uint i = 0; i < n_srcs; ++i) tmp_srcs[i] = sources[i];
    src_dataset.write(tmp_srcs, srctype);
    free(tmp_srcs);

    //iterate over each desired point and save its time series and fourier transform
    printf("found %lu monitor locations\n", n_locs);
    //iterate over monitor locations respecting groups, i will track monitor index, j will track group index
    _uint i = 0;
    size_t n_group_digits = (size_t)(log(monitor_clusters.size()) / log(10)) + 1;
    size_t n_pt_digits = (size_t)(log(n_locs) / log(10)) + 1;
    //we need to keep track of which monitor location each cluster starts with
    size_t mon_loc_offset = 0;
    for (_uint j = 0; j < monitor_clusters.size()+1; ++j) {
	_uint max_i;
	if (j >= monitor_clusters.size())
	    max_i = n_locs;
	else
	    max_i = monitor_clusters[j];
	//write the name for the current group
	strncpy(out_name, CLUSTER_NAME, SMALL_BUF_SIZE);
	write_number(out_name + strlen(CLUSTER_NAME), SMALL_BUF_SIZE-strlen(CLUSTER_NAME), j, n_group_digits);
	H5::Group clust_group = file.createGroup(out_name);
	//write the locations in this cluster
	//set the size of the space to only include however many points are in the current cluster
	size_t n_pts = max_i - mon_loc_offset;
	l_dim[0] = {n_pts};
	l_space = H5::DataSpace(1, l_dim);
	H5::DataSet l_dataset(clust_group.createDataSet("locations", loctype, l_space));
	l_dataset.write(tmp_vecs+mon_loc_offset, loctype);
	//update the monitor location offset
	mon_loc_offset = max_i;

	//write the (prefix for) the points in the group
	strncpy(out_name, POINT_NAME, SMALL_BUF_SIZE);
	for (; i < max_i; ++i) {
        //two is the minimum number of points needed to take np.diff which is used in data analysis
	    if (tdom_fields[i].size < 2) {
		printf("Error: monitor location %d has insufficient points\n", i);
		break;
	    }
	    //create a unique name for the point that respects alphabetic sorting
	    write_number(out_name + strlen(POINT_NAME), SMALL_BUF_SIZE-strlen(POINT_NAME), i, n_pt_digits);
	    printf("%d ", i);
	    H5::Group cur_group = clust_group.createGroup(out_name);

	    //write the time and wavelenuency domain data to the file
	    H5::DataSet t_dataset(cur_group.createDataSet("time", fieldtype, t_space));
	    H5::DataSet f_dataset(cur_group.createDataSet("frequency", fieldtype, f_space));
	    t_dataset.write(tdom_fields[i].buf, fieldtype);
	    f_dataset.write(fdom_fields[i].buf, fieldtype);
	}
    printf("\nfinished writing cluster %d\n", j);
    }
    file.close();
    printf("finished writing hdf5 file!\n");
    free(tmp_vecs);
}

void bound_geom::write_settings(FILE* fp) {
    fprintf(fp, "source end time = %f\ntotal simulation time = %f\n", fields.last_source_time(), ttot);
    //print information about sources
    fprintf(fp, "sources {\n\tn = %lu\n", sources.size());
    size_t jj = 0;
    for (auto it = sources.begin(); it != sources.end(); ++it) {
        fprintf(fp, "\tsource_%lu {\n", jj);
        if (it->type == SRC_GAUSSIAN)
            fprintf(fp, "\t\ttype = Gaussian\n");
        else
            fprintf(fp, "\t\ttype = continuous\n");
        fprintf(fp, "\t\tcomponent = %d\n", it->component);
        fprintf(fp, "\t\tvacuum wavelength = %f (meep length)\n", it->wavelen);
        fprintf(fp, "\t\tpulse length = %f (meep time)\n", it->width);
        fprintf(fp, "\t\tstart time = %f\n", it->start_time);
        fprintf(fp, "\t\tend time = %f\n", it->end_time);
        fprintf(fp, "\t\tphase = %f\n", it->phase);
        fprintf(fp, "\t\tamplitude = %f\n\t}", it->amplitude);
        ++jj;
    }
    fprintf(fp, "}\n");
    //print information about monitor locations
    fprintf(fp, "monitor_clusters {\n\tn = %lu\n", monitor_clusters.size());
    jj = 0;
    for (auto it = monitor_clusters.begin(); it != monitor_clusters.end(); ++it) {
        size_t end_ind = *it;
        fprintf(fp, "\t%s%lu {\n\t\tn = %lu\n", CLUSTER_NAME, std::distance(monitor_clusters.begin(), it), end_ind-jj);
        for (; jj < end_ind; ++jj) {
            fprintf(fp, "\t\t(%f, %f, %f)\n", monitor_locs[jj].x(), monitor_locs[jj].y(), monitor_locs[jj].z());
        }
        fprintf(fp, "\t}\n");
    }
    fprintf(fp, "}\n");
}
