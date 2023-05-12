#include "cgs_data.hpp"

/** ======================================================== scripting language functions ================================== **/

basis_comp_vectors read_component_string(const char* str) {
    basis_comp_vectors component = C_EX;
    if (strcmp(str, "Ex") == 0) {
	component = C_EX;
    } else if (strcmp(str, "Ey") == 0) {
	component = C_EY;
    } else if (strcmp(str, "Ez") == 0) {
	component = C_EZ;
    } else if (strcmp(str, "Hx") == 0) {
	component = C_HX;
    } else if (strcmp(str, "Hy") == 0) {
	component = C_HY;
    } else if (strcmp(str, "Hz") == 0) {
	component = C_HZ;
    }
    return component;
}
value cgs_gen_gaussian_source(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    if (f.n_args < 5) { er = E_LACK_TOKENS;return ret; }
    //check that entries have valid types
    if (f.args[0].type != VAL_STR) { er = E_BAD_TOKEN;return ret; }
    if (!f.args[0].val.s) { er = E_BAD_VALUE;return ret; }
    if (f.args[1].type != VAL_NUM) { er = E_BAD_TOKEN;return ret; }
    if (f.args[2].type != VAL_NUM) { er = E_BAD_TOKEN;return ret; }
    if (f.args[3].type != VAL_NUM) { er = E_BAD_TOKEN;return ret; }
    if (f.args[4].type != VAL_NUM) { er = E_BAD_TOKEN;return ret; }
    value vl = f.args[f.n_args-1].cast_to(VAL_LIST, er);
    if (er || vl.type != VAL_LIST) {
	er = E_BAD_VALUE;
	return ret;
    }
    ret.type = VAL_INST;ret.val.c = new context(&c);
    value tmp = make_val_str("Gaussian_source");
    ret.val.c->emplace( "type",  tmp);
    cleanup_val(&tmp);
    ret.val.c->emplace( "component", make_val_num((double)read_component_string(f.args[0].val.s)) );
    ret.val.c->emplace("wavelength", f.args[1]);
    ret.val.c->emplace("amplitude", f.args[2]);
    ret.val.c->emplace("t_0", f.args[3]);
    ret.val.c->emplace("width", f.args[4]);
    ret.val.c->emplace("phase", f.args[5]);
    //read additional parameters
    ret.val.c->emplace("cutoff", make_val_num(5));
    ret.val.c->emplace("start_time", make_val_num(5));
    ret.val.c->emplace("end_time", make_val_num(5));
    for (size_t i = 5; i < f.n_args; ++i) {
	if (f.arg_names[i]) {
	    if (strcmp(f.arg_names[i], "cutoff") == 0) ret.val.c->set_value("cutoff", f.args[i]);
	    else if (strcmp(f.arg_names[i], "start_time") == 0) ret.val.c->set_value("start_time", f.args[i]);
	    else if (strcmp(f.arg_names[i], "end_time") == 0) ret.val.c->set_value("end_time", f.args[i]);
	}
    }
    //the last argument is always the region
    ret.val.c->emplace("region", vl);
    cleanup_val(&vl);
    return ret;
}
value cgs_gen_continuous_source(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    if (f.n_args < 4) { er = E_LACK_TOKENS;return ret; }
    if (f.args[0].type != VAL_STR) { er = E_BAD_TOKEN;return ret; }
    if (!f.args[0].val.s) { er = E_BAD_VALUE;return ret; }
    if (f.args[1].type != VAL_NUM) { er = E_BAD_TOKEN;return ret; }
    if (f.args[2].type != VAL_NUM) { er = E_BAD_TOKEN;return ret; }
    if (f.args[3].type != VAL_NUM) { er = E_BAD_TOKEN;return ret; }
    value vl = f.args[f.n_args-1].cast_to(VAL_LIST, er);
    if (er || vl.type != VAL_LIST) {
	er = E_BAD_VALUE;
	return ret;
    }
    ret.type = VAL_INST;ret.val.c = new context(&c);
    value tmp = make_val_str("CW_source");
    ret.val.c->emplace("type", tmp);
    cleanup_val(&tmp);
    ret.val.c->emplace( "component", make_val_num((double)read_component_string(f.args[0].val.s)) );
    ret.val.c->emplace("wavelength", f.args[1]);
    ret.val.c->emplace("amplitude", f.args[2]);
    ret.val.c->emplace("start_time", f.args[3]);
    //read additional parameters
    ret.val.c->emplace("slowness", make_val_num(5));
    //use the fourth positional argument or set the end time to infinity
    if (f.n_args > 4 && f.args[4].type == VAL_NUM)
	ret.val.c->emplace("end_time", f.args[4]);
    else
	ret.val.c->emplace("end_time", make_val_num(std::numeric_limits<double>::max()));
    //read named parameters
    for (size_t i = 5; i < f.n_args; ++i) {
	if (f.arg_names[i]) {
	    if (strcmp(f.arg_names[i], "slowness") == 0) ret.val.c->set_value("slowness", f.args[i]);
	    else if (strcmp(f.arg_names[i], "end_time") == 0) ret.val.c->set_value("end_time", f.args[i]);
	}
    }
    //the last argument is always the region
    ret.val.c->emplace("region", vl);
    cleanup_val(&vl);
    return ret;
}
value cgs_gen_box(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    if (f.n_args < 2) { er = E_LACK_TOKENS;return ret; }
    value corn_1 = f.args[0].cast_to(VAL_3VEC, er);
    if (er != E_SUCCESS) return ret;
    value corn_2 = f.args[1].cast_to(VAL_3VEC, er);
    if (er != E_SUCCESS) return ret;
    ret.type = VAL_INST;ret.val.c = new context(&c);
    value tmp = make_val_str("Box");
    ret.val.c->emplace("type", tmp);
    cleanup_val(&tmp);
    ret.val.c->emplace("pt_1", corn_1);
    ret.val.c->emplace("pt_2", corn_2);
    cleanup_val(&corn_1);
    cleanup_val(&corn_2);
    return ret;
}
value cgs_gen_plane(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    if (f.n_args < 2) { er = E_LACK_TOKENS;return ret; }
    value normal;
    double offset = 0;
    if (f.n_args < 3) {
	//this means we have a normal vector and an offset
	normal = f.args[0].cast_to(VAL_3VEC, er);
	if (er != E_SUCCESS || f.args[1].type != VAL_NUM) { er = E_BAD_VALUE;return ret; }
	offset = f.args[1].val.x;
    } else {
	//this means we have three points defining the plane
	value point_1 = f.args[0].cast_to(VAL_3VEC, er);
	if (er != E_SUCCESS) return ret;
	value point_2 = f.args[1].cast_to(VAL_3VEC, er);
	if (er != E_SUCCESS) { cleanup_val(&point_1);return ret; }
	value point_3 = f.args[2].cast_to(VAL_3VEC, er);
	if (er != E_SUCCESS) { cleanup_val(&point_1);cleanup_val(&point_2);return ret; }
	//figure out the normal and offset from the information provided
	vec3 diff_2 = *(point_2.val.v) - *(point_1.val.v);
	vec3 diff_3 = *(point_3.val.v) - *(point_1.val.v);
	normal = make_val_vec3( diff_2.cross(diff_3).normalize() );
	offset = normal.val.v->dot(*(point_1.val.v));
	cleanup_val(&point_1);
	cleanup_val(&point_2);
	cleanup_val(&point_3);
    }
    ret.type = VAL_INST;ret.val.c = new context(&c);
    value tmp = make_val_str("Plane");
    ret.val.c->emplace("type", tmp);
    cleanup_val(&tmp);
    ret.val.c->emplace("normal", normal);
    ret.val.c->emplace("offset", make_val_num(offset));
    cleanup_val(&normal);
    return ret;
}
value cgs_gen_sphere(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    if (f.n_args < 2) { er = E_LACK_TOKENS;return ret; }
    //if we have enough tokens make sure we have both elements as vectors
    value cent = f.args[0].cast_to(VAL_3VEC, er);
    if (er != E_SUCCESS) return ret;
    if (f.args[1].type != VAL_NUM) { er = E_BAD_VALUE;return ret; }
    ret.type = VAL_INST;ret.val.c = new context(&c);
    value tmp = make_val_str("Sphere");
    ret.val.c->emplace("type", tmp);
    cleanup_val(&tmp);
    ret.val.c->emplace("center", cent);
    ret.val.c->emplace("radius", f.args[1]);
    //cleanup
    cleanup_val(&cent);
    return ret;
}
value cgs_gen_cylinder(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    if (f.n_args < 3) { er = E_LACK_TOKENS;return ret; }
    value cent = f.args[0].cast_to(VAL_3VEC, er);
    if (er != E_SUCCESS) return ret;
    if (f.args[1].type != VAL_NUM || f.args[2].type != VAL_NUM) { er = E_BAD_VALUE;return ret; }
    //by default set the top and bottom radius to be the same
    value r2 = f.args[2];
    if (f.n_args > 3 && f.args[3].type == VAL_NUM) r2 = f.args[3];
    ret.type = VAL_INST;ret.val.c = new context(&c);
    value tmp = make_val_str("Cylinder");
    ret.val.c->emplace("type", tmp);
    cleanup_val(&tmp);
    ret.val.c->emplace("center", cent);
    ret.val.c->emplace("h", f.args[1]);
    ret.val.c->emplace("r1", f.args[2]);
    ret.val.c->emplace("r2", r2);
    cleanup_val(&cent);
    return ret;
}
value cgs_gen_composite(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    if (f.n_args < 1) { er = E_LACK_TOKENS;return ret; }
    value vl = f.args[f.n_args-1].cast_to(VAL_LIST, er);
    if (er || vl.type != VAL_LIST) {
	er = E_BAD_VALUE;
	return ret;
    }
    ret.type = VAL_INST;ret.val.c = new context(&c);
    value tmp = make_val_str("Composite");
    ret.val.c->emplace("type", tmp);
    cleanup_val(&tmp);
    for (size_t i = 0; i < f.n_args-1; ++i) {
	if (f.arg_names[i]) ret.val.c->emplace(f.arg_names[i], f.args[i]);
    }
    ret.val.c->emplace("geometry", vl);
    cleanup_val(&vl);
    return ret;
}
value cgs_gen_union(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    if (f.n_args < 1) {
	er = E_LACK_TOKENS;
	return ret;
    }
    value vl = f.args[0].cast_to(VAL_LIST, er);
    if (er || vl.type != VAL_LIST) {
	er = E_BAD_VALUE;
	return ret;
    }
    ret.type = VAL_INST;ret.val.c = new context(&c);
    value tmp = make_val_str("Union");
    ret.val.c->emplace("type", tmp);
    cleanup_val(&tmp);
    ret.val.c->emplace("geometry", vl);
    cleanup_val(&vl);
    return ret;
}
value cgs_gen_intersect(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    if (f.n_args < 1) { er = E_LACK_TOKENS;return ret; }
    value vl = f.args[0].cast_to(VAL_LIST, er);
    if (er || vl.type != VAL_LIST) {
	cleanup_val(&vl);
	er = E_BAD_VALUE;
	return ret;
    }
    ret.type = VAL_INST;ret.val.c = new context(&c);
    value tmp = make_val_str("Intersect");
    ret.val.c->emplace("type", tmp);
    cleanup_val(&tmp);
    ret.val.c->emplace("geometry", vl);
    cleanup_val(&vl);
    return ret;
}
value cgs_gen_complement(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    if (f.n_args < 1) { er = E_LACK_TOKENS;return ret; }
    value vl = f.args[0].cast_to(VAL_LIST, er);
    if (er || vl.type != VAL_LIST) {
	cleanup_val(&vl);
	er = E_BAD_VALUE;
	return ret;
    }
    ret.type = VAL_INST;ret.val.c = new context(&c);
    value tmp = make_val_str("Complement");
    ret.val.c->emplace("type", tmp);
    cleanup_val(&tmp);
    ret.val.c->emplace("geometry", vl);
    cleanup_val(&vl);
    return ret;
}
value cgs_gen_difference(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    if (f.n_args < 1) {
	er = E_LACK_TOKENS;
	return ret;
    }
    value vl = f.args[0].cast_to(VAL_LIST, er);
    if (er || vl.type != VAL_LIST) {
	cleanup_val(&vl);
	er = E_BAD_VALUE;
	return ret;
    }
    if (vl.n_els < 2) {
	cleanup_val(&vl);
	er = E_BAD_VALUE;
	return ret;
    }
    ret.type = VAL_INST;ret.val.c = new context(&c);
    //setup the two sub elements, we subtract all elements after the first from the first
    value contents[2];
    contents[0] = vl.val.l[0];
    //get all of the elements which are subtracted
    value* complement_contents = (value*)malloc(sizeof(value)*(vl.n_els-1));
    for (size_t i = 1; i < vl.n_els; ++i) {
	complement_contents[i-1] = vl.val.l[i];
    }
    //turn that into a Complement instance
    contents[1].type = VAL_INST;
    contents[1].val.c = new context(&c);
    value tmp = make_val_str("Complement");
    contents[1].val.c->emplace("type", tmp);
    cleanup_val(&tmp);
    tmp = make_val_list(complement_contents, vl.n_els-1);
    contents[1].val.c->emplace("geometry", tmp);
    cleanup_val(&tmp);
    //actually put the difference in the instance
    tmp = make_val_str("Intersect");
    ret.val.c->emplace("type", tmp);
    cleanup_val(&tmp);
    tmp.type = VAL_LIST;
    tmp.n_els = 2;
    tmp.val.l = contents;
    ret.val.c->emplace("geometry", tmp);
    cleanup_val(&contents[1]);
    cleanup_val(&vl);
    return ret;
}
value cgs_gen_snapshot(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    //snapshot requires two arguments. The first is the filename that the picture should be saved to and the second is the camera location that should be used
    if (f.n_args < 2) { er = E_LACK_TOKENS;return ret; }
    //read the filename and the viewpoint
    if (f.args[0].type != VAL_STR) { er = E_BAD_VALUE;return ret; }
    value viewpoint = f.args[1].cast_to(VAL_3VEC, er);
    if (er != E_SUCCESS) return ret;
    //default values
    value res(make_val_num(DEF_IM_RES));
    value n_samples(make_val_num(DEF_TEST_N));
    value step(make_val_num(WALK_STEP));
    //read named arguments
    value cam_look = lookup_named(f, "look").cast_to(VAL_3VEC, er);
    if (cam_look.type != VAL_3VEC)
	cam_look = make_val_vec3( *(viewpoint.val.v)*-1 );
    value cam_up = lookup_named(f, "up").cast_to(VAL_3VEC, er);
    if (cam_up.type != VAL_3VEC)
	cam_up = make_val_vec3(vec3(0,0,1));
    //if lookup fails, the error code might be set, we use defaults so we don't care
    er = E_SUCCESS;
    value tmp = lookup_named(f, "resolution");
    if (tmp.type == VAL_NUM)
	res = tmp;
    tmp = lookup_named(f, "n_samples");
    if (tmp.type == VAL_NUM)
	n_samples = tmp;
    tmp = lookup_named(f, "step");
    if (tmp.type == VAL_NUM)
	step = tmp;
    value scale = lookup_named(f, "scale");
    if (scale.type != VAL_NUM)
	scale = make_val_num(0.5 / cam_look.val.v->norm());
    //setup the snapshot object struct
    ret.type = VAL_INST;ret.val.c = new context(&c);
    tmp = make_val_str("snapshot");
    ret.val.c->emplace("type", tmp);
    cleanup_val(&tmp);
    ret.val.c->emplace("fname", f.args[0]);
    ret.val.c->emplace("cam_v", viewpoint);
    ret.val.c->emplace("look_v", cam_look);
    ret.val.c->emplace("up_v", cam_up);
    ret.val.c->emplace("res", res);
    ret.val.c->emplace("n_samples", n_samples);
    ret.val.c->emplace("step", step);
    ret.val.c->emplace("scale", scale);
    cleanup_val(&viewpoint);
    cleanup_val(&cam_look);
    cleanup_val(&cam_up);
    return ret;
}

value cgs_gen_monitor(context& c, cgs_func f, parse_ercode& er) {
    value ret;
    //snapshot requires two arguments. The first is the filename that the picture should be saved to and the second is the camera location that should be used
    if (f.n_args < 1) { er = E_LACK_TOKENS;return ret; }
    if (f.args[0].type != VAL_LIST) { er = E_BAD_TYPE;return ret; }
    ret.type = VAL_INST;ret.val.c = new context(&c);
    value tmp = make_val_str("monitor");
    ret.val.c->emplace( "type",  tmp);
    cleanup_val(&tmp);
    ret.val.c->emplace( "locations", f.args[0] );
    //read the filename and the viewpoint
   return ret;
}

void setup_geometry_context(context& con) {
    //we have to set up the context with all of our functions
    value tmp_f = make_val_func("Gaussian_source", 5, &cgs_gen_gaussian_source);
    con.emplace("Gaussian_source", tmp_f);
    cleanup_val(&tmp_f);
    tmp_f = make_val_func("CW_source", 4, &cgs_gen_continuous_source);
    con.emplace("CW_source", tmp_f);
    cleanup_val(&tmp_f);
    tmp_f = make_val_func("Box", 2, &cgs_gen_box);
    con.emplace("Box", tmp_f);
    cleanup_val(&tmp_f);
    tmp_f = make_val_func("Plane", 2, &cgs_gen_plane);
    con.emplace("Plane", tmp_f);
    cleanup_val(&tmp_f);
    tmp_f = make_val_func("Sphere", 2, &cgs_gen_sphere);
    con.emplace("Sphere", tmp_f);
    cleanup_val(&tmp_f);
    tmp_f = make_val_func("Cylinder", 2, &cgs_gen_cylinder);
    con.emplace("Cylinder", tmp_f);
    cleanup_val(&tmp_f);
    tmp_f = make_val_func("Composite", 1, &cgs_gen_composite);
    con.emplace("Composite", tmp_f);
    cleanup_val(&tmp_f);
    tmp_f = make_val_func("Union", 1, &cgs_gen_union);
    con.emplace("Union", tmp_f);
    cleanup_val(&tmp_f);
    tmp_f = make_val_func("Intersect", 1, &cgs_gen_intersect);
    con.emplace("Intersect", tmp_f);
    cleanup_val(&tmp_f);
    tmp_f = make_val_func("Complement", 1, &cgs_gen_complement);
    con.emplace("Complement", tmp_f);
    cleanup_val(&tmp_f);
    tmp_f = make_val_func("Difference", 1, &cgs_gen_difference);
    con.emplace("Difference", tmp_f);
    cleanup_val(&tmp_f);
    tmp_f = make_val_func("snapshot", 2, &cgs_gen_snapshot);
    con.emplace("snapshot", tmp_f);
    cleanup_val(&tmp_f);
    tmp_f = make_val_func("monitors", 2, &cgs_gen_monitor);
    con.emplace("monitors", tmp_f);
    cleanup_val(&tmp_f);
}
