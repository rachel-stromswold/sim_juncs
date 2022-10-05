#include "disp.hpp"

const double eps_cutoff = 10;

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
cgs_material_function::cgs_material_function(CompositeObject* p_volume, std::string type, double p_def_ret, _uint p_smooth_n, double p_smooth_rad) {
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
	regions = tmp_regs;
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

    //we stochastically smooth boundaries by taking smooth_n samples from a 3d Gaussian distribution with variance smooth_rad^2 centered around r
    double r_x = r.x();double r_y = r.y();double r_z = r.z();

    //look through each region and return the first that contains the vector
    double ret = 0;
    for (_uint i = 0; i < n_regions; ++i) {
	//initialize such that we always include the origin
	double this_ret = regions[i].c->in(evec3(r_x, r_y, r_z));
	if (smooth_pts) {
	    for (_uint j = 0; j < smooth_n; ++j) {
		this_ret += regions[i].c->in(evec3(r_x + smooth_pts[3*j], r_y + smooth_pts[3*j+1], r_z + smooth_pts[3*j+2]));
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

source_info::source_info(std::string spec_str, const Scene& problem, parse_ercode* ercode) {
    //read the specification for the pulse as a function
    char* spec = strdup( spec_str.c_str() );
    cgs_func env_func;
    parse_ercode tmp_er = problem.parse_func(spec, -1, env_func, NULL);

    //figure out the field component that the user wants to add (Default to electric field polarized in x direction)
    component = meep::Ex;
    //all sources require that a component and frequency be specified
    if (tmp_er == E_SUCCESS) {
	if (env_func.n_args > 2) {
	    if (strcmp(env_func.args[0], "Ex") == 0) {
		component = meep::Ex;
	    } else if (strcmp(env_func.args[0], "Ey") == 0) {
		component = meep::Ey;
	    } else if (strcmp(env_func.args[0], "Ez") == 0) {
		component = meep::Ez;
	    } else if (strcmp(env_func.args[0], "Hx") == 0) {
		component = meep::Hx;
	    } else if (strcmp(env_func.args[0], "Hy") == 0) {
		component = meep::Hy;
	    } else if (strcmp(env_func.args[0], "Hz") == 0) {
		component = meep::Hz;
	    } else {
		tmp_er = E_BAD_TOKEN;
	    }

	    if (tmp_er == E_SUCCESS) {
		wavelen = strtod(env_func.args[1], NULL);
		if (errno) tmp_er = E_BAD_TOKEN;
		//figure out what to do depending on what type of pulse envelope this is
		if (strcmp(env_func.name, "gaussian") == 0 || strcmp(env_func.name, "Gaussian") == 0) {
		    type = SRC_GAUSSIAN;
		    //set default values for the envelope
		    start_time = 0;
		    double cutoff = 5;
		    if (env_func.n_args < 3) {
			tmp_er = E_LACK_TOKENS;
		    } else {
			width = strtod(env_func.args[2], NULL);
			if (errno) tmp_er = E_BAD_TOKEN;
			//set default values
			amplitude = 1;
			start_time = 0;
			cutoff = DEFAULT_WIDTH_N;
			//read the start time if supplied
			if (env_func.n_args > 3) {
			    start_time = strtod(env_func.args[3], NULL);
			    if (errno) tmp_er = E_BAD_TOKEN;
			}
			//read the start time if supplied
			if (env_func.n_args > 3) {
			    phase = strtod(env_func.args[3], NULL);
			    if (errno) tmp_er = E_BAD_TOKEN;
			}
			//read the cutoff if supplied
			if (env_func.n_args > 5) {
			    cutoff = strtod(env_func.args[5], NULL);
			    if (errno) tmp_er = E_BAD_TOKEN;
			}
			//read the amplitude if supplied
			if (env_func.n_args > 6) {
			    amplitude = strtod(env_func.args[6], NULL);
			    if (errno) tmp_er = E_BAD_TOKEN;
			}
                end_time = start_time + 2*cutoff*width;
		    }
		} else if (strcmp(env_func.name, "continuous") == 0) {
		    type = SRC_CONTINUOUS;
		    if (env_func.n_args < 3) {
			tmp_er = E_LACK_TOKENS;
		    } else {
			start_time = strtod(env_func.args[2], NULL);
			if (errno) tmp_er = E_BAD_TOKEN;
			//read the end time if supplied
			if (env_func.n_args > 3) {
			    end_time = strtod(env_func.args[3], NULL);
			    if (errno) tmp_er = E_BAD_TOKEN;
			}
			//read the width if supplied
			if (env_func.n_args > 4) {
			    width = strtod(env_func.args[4], NULL);
			    if (errno) tmp_er = E_BAD_TOKEN;
			}
			//read the amplitude if supplied
			if (env_func.n_args > 5) {
			    amplitude = strtod(env_func.args[5], NULL);
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

    free(spec);
}
// bandwidth (in frequency units, not angular frequency) of the
// continuous Fourier transform of the Gaussian source function
// when it has decayed by a tolerance tol below its peak value
static double gaussian_bandwidth(double width) {
    double tol = 1e-7;
    return sqrt(-2.0 * log(tol)) / (width * M_PI);
}

//TODO: remove?
gaussian_src_time_phase::gaussian_src_time_phase(double f, double my_fwidth, double phase, double s) {
    omega = 2*M_PI*f;
    width = 1.0 / my_fwidth;
    phi = phase;
    peak_time = width * s;
    cutoff = width * s;
    fwidth = gaussian_bandwidth(width);
    // correction factor so that current amplitude (= d(dipole)/dt) is
    // ~ 1 near the peak of the Gaussian.
    amp = 1.0 / std::complex<double>(0, -omega);

    // this is to make last_source_time as small as possible
    while (exp(-cutoff * cutoff / (2 * width * width)) < 1e-100)
        cutoff *= 0.9;
    cutoff = float(cutoff); // don't make cutoff sensitive to roundoff error
}

gaussian_src_time_phase::gaussian_src_time_phase(double f, double w, double phase, double st, double et) {
    omega = 2*M_PI*f;
    width = w;
    phi = phase;
    peak_time = 0.5*(st+et);
    cutoff = (et-st)*0.5;
    fwidth = gaussian_bandwidth(width);
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
 * Add the list of susceptibilities to the Settings file s. The string should have the format (omega_0,gamma_0,sigma_0),(omega_1,gamma_1,sigma_1),...
 * returns: a vector list of susceptibilities
 * If an error is encountered, a code will be saved to er if it is not NULL
 *  0 on success
 *  -1 null string
 *  -2 invalid or empty string
 *  -3 insufficient memory
 */
std::vector<drude_suscept> bound_geom::parse_susceptibilities(char* const str, int* er) {
    std::vector<drude_suscept> ret;
    //check that the Settings struct is valid and allocate memory
    if (!str) {
	set_ercode(er, -1);
	return ret;
    }

    drude_suscept cur_sus;

    //used for strtok_r
    char* save_str;
    char* tok;
    //find the first entry
    char* cur_entry = strchr(str, '(');
    char* end = strchr(str, ')');

    //only proceed if we have pointers to the start and end of the current entry
    while (cur_entry && end) {
	cur_sus.omega_0 = 0.0;
	cur_sus.gamma = 0.0;
	cur_sus.sigma = 0.0;
	cur_sus.use_denom = 1;
	//the open paren must occur before the end paren
	if (cur_entry > end) {
	    set_ercode(er, -2);
	    return ret;
	}

	//null terminate the parenthesis and tokenize by commas
	end[0] = 0;
	//read the omega value
	tok = trim_whitespace( strtok_r(cur_entry+1, ",", &save_str), NULL );
	cur_sus.omega_0 = strtod(tok, NULL);
	if (errno) {
	    errno = 0;
	    set_ercode(er, -2);
	    return ret;
	}
	//read the gamma value
	tok = trim_whitespace( strtok_r(NULL, ",", &save_str), NULL );
	if (!tok) {
	    set_ercode(er, -2);
	    return ret;
	}
	cur_sus.gamma = strtod(tok, NULL);
	if (errno) {
	    errno = 0;
	    set_ercode(er, -2);
	    return ret;
	}
	//read the sigma value
	tok = trim_whitespace( strtok_r(NULL, ",", &save_str), NULL );
	if (!tok) {
	    set_ercode(er, -2);
	    return ret;
	}
	cur_sus.sigma = strtod(tok, NULL);
	if (errno) {
	    errno = 0;
	    set_ercode(er, -2);
	    return ret;
	}
	//read the (optional) use denom flag
	tok = strtok_r(NULL, ",", &save_str);
	if (tok) {
	    tok = trim_whitespace(tok, NULL);
	    cur_sus.use_denom = strtol(tok, NULL, 10);
	    if (errno) {
		errno = 0;
		cur_sus.use_denom = 1;
		if (strcmp(tok, "drude") == 0 || strcmp(tok, "true") == 0) cur_sus.use_denom = 0;
	    }
	    if (strcmp(tok, "lorentz") == 0 || strcmp(tok, "false") == 0) cur_sus.use_denom = 1;
	}

	//save the information
	ret.push_back(cur_sus);

	//advance to the next entry
	if (end[1] == 0) break;
	cur_entry = strchr(end+1, '(');
	if (!cur_entry) break;
	end = strchr(cur_entry, ')');
    }

    return ret;
}

/**
 * Read a CompositeObject specifying monitor locations into a list of monitor locations
 */
void bound_geom::parse_monitors(CompositeObject* comp) {
    //only continue if the user specified locations
    if (comp->has_metadata("locations")) {
	char* loc_str = strdup(comp->fetch_metadata("locations").c_str());
	//break up the string by end parens
	char* save_str;
	char* tok = NULL;
	//find the first entry
	char* cur_entry = strchr(loc_str, '(');
	char* end = strchr(loc_str, ')');

	double x = 0.0; double y = 0.0; double z = 0.0;

	//initialize the array of monitor locs and deallocate memory if necessary
	size_t cur_size = 2*SMALL_BUF_SIZE;
	//only proceed if we have pointers to the start and end of the current entry
	while (cur_entry && end) {
	    x = 0.0;
	    y = 0.0;
	    z = 0.0;
	    //the open paren must occur before the end paren
	    if (cur_entry > end) {
		std::cout << "Error: Invalid locations string" << comp->fetch_metadata("locations") << std::endl;
		break;
	    }

	    //null terminate the parenthesis and tokenize by commas
	    end[0] = 0;
	    //read the x value
	    tok = trim_whitespace( strtok_r(cur_entry+1, ",", &save_str), NULL );
	    x = strtod(tok, NULL);
	    if (errno) {
		printf("Error: Invalid token: %s", tok);
		errno = 0;
		break;
	    }
	    //read the y value
	    tok = trim_whitespace( strtok_r(NULL, ",", &save_str), NULL );
	    if (!tok) {
		std::cout << "Error: Invalid locations string" << comp->fetch_metadata("locations") << std::endl;
		break;
	    }
	    y = strtod(tok, NULL);
	    if (errno) {
		printf("Error: Invalid token: %s", tok);
		break;
	    }
	    //read the z value
	    tok = trim_whitespace( strtok_r(NULL, ",", &save_str), NULL );
	    if (!tok) {
		std::cout << "Error: Invalid locations string" << comp->fetch_metadata("locations") << std::endl;
		break;
	    }
	    z = strtod(tok, NULL);
	    if (errno) {
		printf("Error: Invalid token: %s", tok);
		break;
	    }

	    //construct the meep vector in place inside the buffer
	    monitor_locs.push_back(meep::vec(x, y, z));

	    //advance to the next entry
	    if (end[1] == 0) break;
	    cur_entry = strchr(end+1, '(');
	    if (!cur_entry) break;
	    end = strchr(cur_entry, ')');
	}
	free(loc_str);
    } else if (comp->has_metadata("spacing")) {
	double spacing = std::stod(comp->fetch_metadata("spacing"));
	//users may also specify monitor locations by providing a Box object and a spacing between each grid point
	const Object* l_child = comp->get_child_l();
	if (l_child && comp->get_child_type_l() == CGS_BOX) {
	    evec3 center = ((Box*)l_child)->get_center();
	    evec3 offset = ((Box*)l_child)->get_offset();
	    //fix the box such that each offset is positive
	    if (offset.x() < 0) offset.x() = -1*offset.x();
	    if (offset.y() < 0) offset.y() = -1*offset.y();
	    if (offset.z() < 0) offset.z() = -1*offset.z();
	    //find the two corners of the box
	    evec3 min = center - offset;
	    evec3 max = center + offset;
	    //iterate over the volume with step size=spacing
	    for (double x = min.x(); x <= max.x(); x += spacing) {
		for (double y = min.y(); y <= max.y(); y += spacing) {
		    for (double z = min.y(); z <= max.z(); z += spacing) {
			monitor_locs.push_back(meep::vec(x, y, z));
		    }
		}
	    }
	} 
    }
}

double dummy_eps(const meep::vec& r) { return 1.0; }

/**
 * This is a helper function for the bound_geom constructor. Meep doesn't implement copy or move constructors so we have to initialize the structure immediately so that the fields can be initialized in turn.
 */
meep::structure* bound_geom::structure_from_settings(const Settings& s, Scene& problem, parse_ercode* ercode) {
    pml_thickness = s.pml_thickness;
    len = s.len;

    //read information about the problem size from the geometry file
    std::vector<CompositeObject*> data = problem.get_data();
    for (size_t i = 0; i < data.size(); ++i) {
	if (data[i]->has_metadata("pml_thickness")) {
	    pml_thickness = std::stod(data[i]->fetch_metadata("pml_thickness"));
	}
	if (data[i]->has_metadata("length")) {
	    len = std::stod(data[i]->fetch_metadata("length"));
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
    std::vector<CompositeObject*> roots = problem.get_roots();
    //setup the structure with the infinite frequency dielectric component
    cgs_material_function inf_eps_func(s.ambient_eps, s.smooth_n, s.smooth_rad);

    std::vector<double> thicknesses(roots.size());
    for (size_t i = 0; i < roots.size(); ++i) {
	//we need to adjust the thickness and conductivity of the sample if it is 2D
	thicknesses[i] = 1.0;
	if (roots[i]->has_metadata("make_2d") && roots[i]->fetch_metadata("make_2d") != "false" && roots[i]->fetch_metadata("make_2d") != "0") {
	    thicknesses[i] = THICK_SCALE / s.resolution;
	    roots[i]->rescale(evec3(1.0, 1.0, thicknesses[i]));
	}
	inf_eps_func.add_region(roots[i]);
	/** ============================ DEBUG ============================ **/
	meep::vec test_loc_1(0.5,0.5,2);
	meep::vec test_loc_2(0.1,0.5,6);
	meep::vec test_loc_3(0.5,0.1,6);
	meep::vec test_loc_4(8,0.5,6);
	meep::vec test_loc_5(0.5,8,6);
	meep::vec test_loc_6(12,0.5,6);
	meep::vec test_loc_7(0.5,12,6);
	double ret_1 = inf_eps_func.eps(test_loc_1);
	double ret_2 = inf_eps_func.eps(test_loc_2);
	double ret_3 = inf_eps_func.eps(test_loc_3);
	double ret_4 = inf_eps_func.eps(test_loc_4);
	double ret_5 = inf_eps_func.eps(test_loc_5);
	double ret_6 = inf_eps_func.eps(test_loc_6);
	double ret_7 = inf_eps_func.eps(test_loc_7);
	printf("%f %f %f %f %f %f %f\n", ret_1, ret_2, ret_3, ret_4, ret_5, ret_6, ret_7);
	/** ============================ DEBUG ============================ **/
    }
    meep::structure* strct = new meep::structure(vol, inf_eps_func, meep::pml(s.pml_thickness));
    //read susceptibilities if they are available
    for (size_t i = 0; i < roots.size(); ++i) {
	std::vector<drude_suscept> cur_sups;
	int res = 0;
	if (roots[i]->has_metadata("susceptibilities")) {
	    char* dat = strdup(roots[i]->fetch_metadata("susceptibilities").c_str());
	    cur_sups = parse_susceptibilities(dat, &res);
	    free(dat);
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
 * Constructor for the bound_geom file. Settings are read from s.
 * s: Settings object to read.
 * ercode: if an error occurs while parsing the .geom file and ercode is not NULL, a code for the error is saved there
 */
bound_geom::bound_geom(const Settings& s, parse_ercode* ercode) :
    problem(s.geom_fname, ercode),
    strct(structure_from_settings(s, problem, ercode)),
    fields(strct)
{
    um_scale = s.um_scale;
    if (ercode) *ercode = E_SUCCESS;

    //we have to kludge it to get around the very f** annoying fact that meep doesn't have default constructors for fields, just read the structure_from_settings comment
    double z_center = s.len/2 + s.pml_thickness;
    double eps_scale = 1 / (sharpness*s.resolution);
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
    std::vector<CompositeObject*> data = problem.get_data();
    for (size_t i = 0; i < data.size(); ++i) {
	if (data[i]->has_metadata("type")) {
	    if (data[i]->fetch_metadata("type") == "field_source") {
		//figure out the volume for the field source
		const Object* l_child = data[i]->get_child_l();
		const Object* r_child = data[i]->get_child_r();
		object_type l_type = data[i]->get_child_type_l();
		//WLOG fix the left child to be the one we care about
		if (r_child && !l_child) {
		    l_child = r_child;
		    l_type = data[i]->get_child_type_r();
		}
		//make sure that the object is a box so that we can figure out the corner and the offset
		if (l_type == CGS_BOX) {
		    evec3 center = ((Box*)l_child)->get_center();
		    evec3 offset = ((Box*)l_child)->get_offset();
		    double x_0 = center.x() - offset.x();double x_1 = center.x() + offset.x();
		    double y_0 = center.y() - offset.y();double y_1 = center.y() + offset.y();
		    double z_0 = center.z() - offset.z();double z_1 = center.z() + offset.z();
		    meep::volume source_vol(meep::vec(x_0, y_0, z_0), meep::vec(x_1, y_1, z_1));

		    //only continue if a shape for the pulse was specified
		    if (data[i]->has_metadata("envelope")) {
			source_info cur_info(data[i]->fetch_metadata("envelope"), problem, ercode);
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
	    } else if (data[i]->fetch_metadata("type") == "monitor") {
		parse_monitors(data[i]);
		monitor_clusters.push_back(monitor_locs.size());
	    }
	}
    }
    printf("source end time: %f, total simulation time: %f\n", fields.last_source_time(), ttot);

    dump_span = s.field_dump_span;
}

bound_geom::~bound_geom() {
    if (strct) delete strct;
}

/*std::vector<meep::vec> bound_geom::get_monitor_locs() {
    std::vector<meep::vec> ret(n_locs);
    for(size_t i = 0; i < n_locs; ++i) ret[i] = monitor_locs[i];
    return ret;
}*/

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
void bound_geom::add_volume_source(meep::component c, const meep::src_time &src, const Box& vol, std::complex<double> amp) {
    evec3 center = vol.get_center();
    evec3 offset = vol.get_offset();
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

    //save the dielectric used
    printf("Set output directory to %s\n", fname_prefix);
    fields.output_hdf5(meep::Dielectric, fields.total_volume());

    //make sure the time series corresponding to each monitor point is long enough to hold all of its information
    size_t n_locs = monitor_locs.size();
    field_times.resize(n_locs);
    n_t_pts = (_uint)( (ttot+fields.dt/2) / fields.dt );
    for (_uint j = 0; j < n_locs; ++j) {
	make_data_arr(&(field_times[j]), n_t_pts);
    }

    //figure out the number of digits before the decimal and after
    int n_digits_a = (int)(ceil(log(ttot)/log(10)));
    double rat = -log((double)(fields.dt))/log(10.0);
    int n_digits_b = (int)ceil(rat)+1;
    char h5_fname[BUF_SIZE];
    strcpy(h5_fname, "ex-");

    //ensure that no field dumps are saved
    if (dump_span == 0) dump_span = n_t_pts+1;
    //run the simulation
    printf("starting simulations\n");
    for (_uint i = 0; i < n_t_pts; ++i) {
	//fetch monitor points
	for (_uint j = 0; j < n_locs; ++j) {
	    std::complex<double> val = fields.get_field(meep::Ex, monitor_locs[j]);
	    field_times[j].buf[i].re = val.real();
	    field_times[j].buf[i].im = val.imag();
	    if (field_times[j].buf[i].re > 1000) {
		printf("divergence in run at (i,j)=(%d,%d) (%f)\n", i, j, field_times[j].buf[i].re);
	    }
	}

        //open an hdf5 file with a reasonable name
        if (i % dump_span == 0) {
	    size_t n_written = make_dec_str(h5_fname+PREFIX_LEN, BUF_SIZE-PREFIX_LEN, fields.time(), n_digits_a, n_digits_b);
	    /*meep::h5file* file = fields.open_h5file(h5_fname);
	    fields.output_hdf5(meep::Ex, vol.surroundings(), file);
	    delete file;*/
	    fields.step();
	    printf("    %d%% complete\n", 100*i/n_t_pts);
	    //we're done with the file
	}
    }
    printf("Simulations completed\n");
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
    fieldtype.insertMember("Re", HOFFSET(complex, re), H5_float_type);
    fieldtype.insertMember("Im", HOFFSET(complex, im), H5_float_type);
    //use the space of rank 1 tensors with a dimension of n_t_pts
    hsize_t t_dim[1];
    t_dim[0] = {n_t_pts};
    hsize_t f_dim[1];
    H5::DataSpace t_space(1, t_dim);

    //create the location type and specify members
    H5::CompType loctype(sizeof(sto_vec));
    loctype.insertMember("x", HOFFSET(sto_vec, x), H5_float_type);
    loctype.insertMember("y", HOFFSET(sto_vec, y), H5_float_type);
    loctype.insertMember("z", HOFFSET(sto_vec, z), H5_float_type);
    //use the space of rank 1 tensors with dimension of the number of monitor points
    hsize_t l_dim[1];
    l_dim[0] = {n_locs};
    H5::DataSpace l_space(1, l_dim);

    //create the source type and specify members
    H5::CompType srctype(sizeof(source_info));
    /*H5::DataType enum_src_type(H5T_ENUM, sizeof(src_type));
    H5::DataType enum_component(H5T_ENUM, sizeof(meep::component));
    srctype.insertMember("type", HOFFSET(source_info, type), enum_src_type);
    srctype.insertMember("component", HOFFSET(source_info, component), enum_component);*/
    srctype.insertMember("wavelen", HOFFSET(source_info, wavelen), H5_float_type);
    srctype.insertMember("width", HOFFSET(source_info, width), H5_float_type);
    srctype.insertMember("start_time", HOFFSET(source_info, start_time), H5_float_type);
    srctype.insertMember("end_time", HOFFSET(source_info, end_time), H5_float_type);
    srctype.insertMember("amplitude", HOFFSET(source_info, amplitude), H5_float_type);

    //take the fourier transform for each point
    std::vector<data_arr> fours(n_locs);
    for (_uint j = 0; j < n_locs; ++j) {
        fours[j] = fft(field_times[j]);
        f_dim[0] = fours[j].size;
    }
    H5::DataSpace f_space(1, f_dim);

    //open the file which will store the fields as a function of time
    snprintf(out_name, BUF_SIZE, "%s/field_samples.h5", fname_prefix);
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
    t_info_dim[0] = {2};
    H5::DataSpace t_info_space(1, t_info_dim);
    H5::DataSet t_info_dataset(info_group.createDataSet("time_bounds", H5_float_type, t_info_space));
    //set the start and end times for the simulation TODO: is the first boundary necessary or can it be left implicit?
    _ftype time_boundaries[2];
    time_boundaries[0] = 0.0;time_boundaries[1] = meep_time_to_fs(ttot);
    t_info_dataset.write(time_boundaries, H5_float_type);
    //write the number of clusters
    hsize_t n_info_dim[1];
    n_info_dim[0] = {1};
    H5::DataSpace n_info_space(1, n_info_dim);
    H5::DataSet n_info_dataset(info_group.createDataSet("n_clusters", H5::PredType::NATIVE_HSIZE, n_info_space));
    hsize_t n_clusters = monitor_clusters.size();
    n_info_dataset.write(&n_clusters, H5::PredType::NATIVE_HSIZE);
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
    printf("found %d monitor locations\n", n_locs);
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
	    if (field_times[i].size < n_t_pts) {
		printf("Error: monitor location %d has insufficient points\n", i);
		break;
	    }
	    //create a unique name for the point that respects alphabetic sorting
	    write_number(out_name + strlen(POINT_NAME), SMALL_BUF_SIZE-strlen(POINT_NAME), i, n_pt_digits);
	    printf("saving point %d to group %s\n", i, out_name);
	    H5::Group cur_group = clust_group.createGroup(out_name);

	    //write the time and wavelenuency domain data to the file
	    H5::DataSet t_dataset(cur_group.createDataSet("time", fieldtype, t_space));
	    H5::DataSet f_dataset(cur_group.createDataSet("frequency", fieldtype, f_space));
	    t_dataset.write(field_times[i].buf, fieldtype);
	    f_dataset.write(fours[i].buf, fieldtype);
	}
    }
    free(tmp_vecs);
}
