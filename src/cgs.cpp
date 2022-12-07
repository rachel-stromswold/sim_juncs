#include "cgs.hpp"

/** ======================================================== object ======================================================== **/

object::object(int p_invert) {
    invert = p_invert;
    trans_mat = emat3::Identity();
}

void object::set_inversion(int p_invert) {
    if (p_invert == 0)
	invert = 0;
    else
	invert = 1;
}

void object::rescale(const evec3& components) {
    trans_mat = trans_mat*(components.asDiagonal());
}

uint32_t lcg(uint32_t state) { return state*1664525 + 1013904223; }//LCG PRNG with parameters from Numerical Recipes
void object::draw(const char* out_fname, double scale, evec3 cam_pos, evec3 cam_look, size_t res, size_t n_samples) {
    double x,y,z;
    //initialize random state
    uint32_t state = lcg(lcg(1465432554));
    //cross the look direction with the z direction
    evec3 x_comp(cam_look.y(), -cam_look.x(), 0);
    evec3 y_comp = x_comp.cross(cam_look);
    emat3 view_mat;
    view_mat << x_comp/x_comp.norm(), y_comp/y_comp.norm(), cam_look/cam_look.norm();
    //store the image buffer in an array. We don't do anything fancy for shading, this is just a z-buffer. Thus, we initialize to white (the maximum distance).
    _uint8 im_arr[res*res];
    for (size_t i = 0; i < res*res; ++i) im_arr[i]=0xff;
    //sample random points. If a point is in the object, and it's closer than any other sampled point, then store it in the z-buffer.
    for (size_t i = 0; i < n_samples; ++i) {
	//get a random point
	state = lcg(state);
	x = (double)state / UINT32_MAX;
	state = lcg(state);
	y = (double)state / UINT32_MAX;
	state = lcg(state);
	z = (double)state / UINT32_MAX;
	state = lcg(state);
	//initialize a point from the random values and compute it's position projected onto the view matrix. The z-buffer will act as color intensity
	evec3 r(x, y, z);
	evec3 view_r = view_mat*(r*scale - cam_pos);
	//figure out the index and depth and store the result into the array
	size_t ind = floor( (1+view_r.x())*res/2 ) + res*floor( (1+view_r.y())*res/2 );
	if (ind >= res*res)
	    ind = 0;
	_uint8 depth = floor(IM_DEPTH*view_r.z());
	if (in(r))
	    if (depth < im_arr[ind]) im_arr[ind] = depth;
    }
    //open the file. We use a .pgm format since it's super simple. Then we write the header information.
    FILE* fp;
    if (out_fname) {
	fp = fopen(out_fname, "w");
	if (!fp) fp = stdout;
    } else {
	fp = stdout;
    }
    fprintf(fp, "P2\n%d %d\n%d\n", res, res, IM_DEPTH);
    //iterate through the image buffer and write
    for (size_t yy = 0; yy < res; ++yy) {
	for (size_t xx = 0; xx < res; ++xx) {
	    fprintf(fp, "%d ", im_arr[yy*res+xx]);
	}
	fprintf(fp, "\n");
    }
    fclose(fp);
}

/** ======================================================== sphere ======================================================== **/

sphere::sphere(evec3& p_center, double p_rad, int p_invert) : object(p_invert) {
    center = p_center;
    rad = p_rad;
}

sphere::sphere(evec3& p_center, evec3& p_rad, int p_invert) : object(p_invert) {
    center = p_center;
    rad = 1;
    rescale( evec3(1/p_rad.x(), 1/p_rad.y(), 1/p_rad.z()) );
}

int sphere::in(const evec3& r) {
    //find the relative offset between the input vector and the center of the box
    evec3 r_rel = trans_mat*(r - center);
    if (r_rel.norm() > rad) return invert;
    return 1 - invert;
}

box::box(evec3& p_corner_1, evec3& p_corner_2, int p_invert) : object(p_invert) {
    offset = (p_corner_2 - p_corner_1)/2;
    center = p_corner_1 + offset;
    //wlog, fix the offset to always be positive from the center
    if (offset.x() < 0) offset.x() *= -1;
    if (offset.y() < 0) offset.y() *= -1;
    if (offset.z() < 0) offset.z() *= -1;
}

/** ======================================================== box ======================================================== **/

int box::in(const evec3& r) {
    //find the relative offset between the input vector and the center of the box
    evec3 r_rel = trans_mat*(r - center);
    //check if we're outside the bounding box
    if (fabs(r_rel.x()) > offset.x()) {
        return invert;
    }
    if (fabs(r_rel.y()) > offset.y()) {
        return invert;
    }
    if (fabs(r_rel.z()) > offset.z()) {
        return invert;
    }
    //if we reach this point in execution we're inside the box
    return 1 - invert;
}

/**
 * The term plane is a bit misleading since we are describing the three dimensional volume. The plane has a fixed orientation defined by the normal vector and an offset. Any point with a component along the normal vector which is less than the offset will be included in the volume.
 */
plane::plane(evec3& p_normal, double p_offset, int p_invert) : object(0) {
    if (p_invert)
	offset = -p_offset;
    else
	offset = p_offset;
    normal = p_normal / p_normal.norm();
}

/**
 * Construct a plane defined by three points. The normal vector is fixed to be (point_2 - point_1)X(point_3 - point_1). Thus, flipping the order of the arguments point_1 and point_2 flips the direction of inclusion.
 */
plane::plane(evec3& point_1, evec3& point_2, evec3& point_3, int p_invert) : object(0) {
    evec3 diff_2 = point_2 - point_1;
    evec3 diff_3 = point_3 - point_1;
    normal = diff_2.cross(diff_3);
    normal /= normal.norm();
    offset = normal.dot(point_1);
    if (p_invert)
	offset *= -1;
}

int plane::in(const evec3& r) {
    double norm_comp = r.dot(normal);
    if (norm_comp > offset)
	return 0;
    else
	return 1;
}

/** ======================================================== cylinder ======================================================== **/

cylinder::cylinder(evec3& p_center, double p_height, double p_r1, double p_r2, int p_invert) : object(p_invert) {
    center = p_center;
    height = p_height;
    r1_sq = p_r1*p_r1;
    r1_sq_x_h = r1_sq*height;
    r2_sq = p_r2*p_r2;
}

int cylinder::in(const evec3& r) {
    //find the relative offset between the input vector and the center of the box
    evec3 r_rel = trans_mat*(r - center);
    if (r_rel.z() < 0 || r_rel.z() > height) return invert;
    double cent_rad_sq = r_rel.x()*r_rel.x() + r_rel.y()*r_rel.y();
    if (cent_rad_sq*height > (r2_sq - r1_sq)*r_rel.z() + r1_sq_x_h) return invert;

    return 1 - invert;
}

/** ======================================================== composite_object ======================================================== **/

composite_object::composite_object(combine_type p_cmb) {
    children[0] = NULL;
    children[1] = NULL;
    cmb = p_cmb;
}

composite_object::composite_object(combine_type p_cmb, const cgs_func& spec, int p_invert) : object(p_invert) {
    children[0] = NULL;
    children[1] = NULL;
    //iterate through arguments and add them to the metadata
    for (_uint i = 0; i < spec.n_args; ++i) {
	if (spec.arg_names[i]) {
	    //each piece of metadata should be separated by a named list
	    std::string tok_cpp(spec.arg_names[i]);
	    metadata[tok_cpp] = copy_val(spec.args[i]);
	}
    }
    cmb = p_cmb;
}

/**
 * Helper function which copies the child on the specified side. The caller is responsible for calling delete on the object and managing its lifetime.
 */
object* composite_object::copy_child(_uint side) const {
    //declare these outside of the switch statement so we don't bork the addresses
    composite_object* my_child = NULL;
    composite_object* comp = NULL;
    switch (child_types[side]) {
	//if the type is not a composite then copying is easy
	case CGS_BOX: return new box( *((box*)children[side]) );break;
	case CGS_SPHERE: return new sphere( *((sphere*)children[side]) );break;
	case CGS_CYLINDER: return new cylinder( *((cylinder*)children[side]) );break;
       //otherwise...
	case CGS_COMPOSITE: 
	    //allocate a new child and store a pointer of this child for convenience
	    my_child = (composite_object*)(children[side]);
	    comp = new composite_object(my_child->cmb);
	    //copy the metadata and child types
	    comp->metadata = my_child->metadata;
	    comp->child_types[0] = my_child->child_types[0];
	    comp->child_types[1] = my_child->child_types[1];
	    //Copying of the children must be done recursively. The base case is any child which is not a composite.
	    comp->children[0] = my_child->copy_child(0);
	    comp->children[1] = my_child->copy_child(1);
	    break;
	default: break;
    }
    return NULL;
}

//copy constructor
composite_object::composite_object(const composite_object& o) {
    //copy the "easy" stuff
    metadata = o.metadata;
    child_types[0] = o.child_types[0];
    child_types[1] = o.child_types[1];
    //copying children is potentially thorny. We have a helper function to do it for us
    children[0] = o.copy_child(0);
    children[1] = o.copy_child(1);
}

//move constructor
composite_object::composite_object(composite_object&& o) {
    metadata = o.metadata;
    //do this as a loop to avoid typing :p
    for (_uint side = 0; side < 2; ++side) {
	child_types[side] = o.child_types[side];
	children[side] = o.children[side];
	o.children[side] = NULL;
    }
}

//assignemnt operator
composite_object& composite_object::operator=(composite_object&& o) {
    //swap the easy stuff
    std::unordered_map<std::string, value> tmp_metadata = metadata;
    metadata = o.metadata;
    o.metadata = tmp_metadata;
    //do this as a loop to avoid typing :p
    for (_uint side = 0; side < 2; ++side) {
	//save the current
	object_type tmp_ctype = child_types[side];
	object* tmp_obj = children[side];
	//alter the current
	child_types[side] = o.child_types[side];
	children[side] = o.children[side];
	//store the current
	o.child_types[side] = tmp_ctype;
	o.children[side] = tmp_obj;
    }
    return *this;
}

//destructor
composite_object::~composite_object() {
    for (_uint i = 0; i < 2; ++i) {
	if (children[i]) {
	    if (child_types[i] == CGS_COMPOSITE)
		delete (composite_object*)(children[i]);
	    else
		delete children[i];
	}
    }
    for (auto it = metadata.begin(); it != metadata.end(); ++it) {
	cleanup_val(&(it->second));
    }
}

int composite_object::call_child_in(_uint side, const evec3& r) {
    if (child_types[side] == CGS_COMPOSITE) return ((composite_object*)children[side])->in(r);
    if (child_types[side] == CGS_SPHERE) return ((sphere*)children[side])->in(r);
    if (child_types[side] == CGS_BOX) return ((box*)children[side])->in(r);
    if (child_types[side] == CGS_PLANE) return ((plane*)children[side])->in(r);
    if (child_types[side] == CGS_CYLINDER) return ((cylinder*)children[side])->in(r);

    return invert;
}

void composite_object::add_child(_uint side, object* o, object_type p_type) {
    children[side] = o;
    child_types[side] = p_type;
    //differences can be expressed as intersections of negations
    if (cmb == CGS_DIFFERENCE && side == 1) {
	children[side]->set_inversion(1);	
    }
}

int composite_object::in(const evec3& r) {
    evec3 r_trans = trans_mat*r;
    //NOOP items should be ignored
    if (cmb == CGS_CMB_NOOP) return 0;
    //if the right child is null then we don't apply any operations
    if (children[1] == NULL) {
	if (children[0] == NULL)
	    return invert;
	else
	    return invert ^ call_child_in(0, r_trans);
    } else {
	if (children[0] == NULL) {
	    return invert ^ call_child_in(1, r_trans);
	} else {
	    int l_res = call_child_in(0, r_trans);
	    int r_res = call_child_in(1, r_trans);
	    if (cmb == CGS_UNION) {
		return invert ^ (l_res | r_res);
	    } else if (cmb == CGS_INTERSECT || cmb == CGS_DIFFERENCE) {
		return invert ^ (l_res & r_res);
	    }
	}
    }

    return invert;
}

/**
 * Based on the declaration syntax produce the appropriate geometric shape and store the result in ptr. Ptr must not be initialized before a call to this function to avoid a memory leak.
 * returns: 0 on success or an error code
 */
parse_ercode scene::make_object(const cgs_func& f, object** ptr, object_type* type, int p_invert) const {
    *ptr = NULL;
    parse_ercode er = E_SUCCESS;

    if (!f.name) return E_BAD_TOKEN;
    //switch between all potential types
    if (strcmp(f.name, "Box") == 0) {
	if (f.n_args < 2) return E_LACK_TOKENS;
	//if we have enough tokens make sure we have both elements as vectors
	value corn_1 = f.args[0].cast_to(VAL_3VEC, er);
	if (er != E_SUCCESS) return er;
	value corn_2 = f.args[1].cast_to(VAL_3VEC, er);
	if (er != E_SUCCESS) return er;
	//make the box
	if (ptr) *ptr = new box(*(corn_1.val.v), *(corn_2.val.v), p_invert);
	if (type) *type = CGS_BOX;
	//cleanup
	cleanup_val(&corn_1);
	cleanup_val(&corn_2);
    } else if (strcmp(f.name, "Plane") == 0) {
	if (f.n_args < 2) return E_LACK_TOKENS;
	if (f.n_args < 3) {
	    //this means we have a normal vector and an offset
	    value corn_1 = f.args[0].cast_to(VAL_3VEC, er);
	    if (er != E_SUCCESS || f.args[1].type != VAL_NUM) return E_BAD_VALUE;
	    //make the plane
	    if (ptr) *ptr = new plane(*(corn_1.val.v), f.args[1].val.x);
	    if (type) *type = CGS_PLANE;
	    cleanup_val(&corn_1);
	} else {
	    //this means we have three points defining the plane
	    value corn_1 = f.args[0].cast_to(VAL_3VEC, er);
	    if (er != E_SUCCESS) return er;
	    value corn_2 = f.args[1].cast_to(VAL_3VEC, er);
	    if (er != E_SUCCESS) return er;
	    value corn_3 = f.args[2].cast_to(VAL_3VEC, er);
	    if (er != E_SUCCESS) return er;
	    //make the plane
	    if (ptr) *ptr = new plane(*(corn_1.val.v), *(corn_2.val.v), *(corn_3.val.v));
	    if (type) *type = CGS_PLANE;
	    cleanup_val(&corn_1);
	    cleanup_val(&corn_2);
	    cleanup_val(&corn_3);
	}
    } else if (strcmp(f.name, "Sphere") == 0) {
	if (f.n_args < 2) return E_LACK_TOKENS;
	//if we have enough tokens make sure we have both elements as vectors
	value cent = f.args[0].cast_to(VAL_3VEC, er);
	if (er != E_SUCCESS || f.args[1].type != VAL_NUM) return E_BAD_VALUE;
	double rad = f.args[1].val.x;
	if (ptr) *ptr = new sphere(*(cent.val.v), rad, p_invert);
	if (type) *type = CGS_SPHERE;
	//cleanup
	cleanup_val(&cent);
    } else if (strcmp(f.name, "Cylinder") == 0) {
	if (f.n_args < 3) return E_LACK_TOKENS;
	value cent = f.args[0].cast_to(VAL_3VEC, er);
	if (er != E_SUCCESS || f.args[1].type != VAL_NUM || f.args[2].type != VAL_NUM) return E_BAD_VALUE;
	double h = f.args[1].val.x;
	double r1 = f.args[2].val.x;
	//by default assume that the radii are the same
	double r2 = r1;
	if (f.n_args > 3) {
	    if (f.args[3].type != VAL_NUM) return E_BAD_VALUE;
	    r2 = f.args[3].val.x;
	}
	if (ptr) *ptr = new cylinder(*(cent.val.v), h, r1, r2, p_invert);
	if (type) *type = CGS_CYLINDER;
	//cleanup
	cleanup_val(&cent);
    } else if (strcmp(f.name, "Composite") == 0) {
	if (ptr) *ptr = new composite_object(CGS_UNION, f, p_invert);
	if (type) *type = CGS_ROOT;
    } else if (strcmp(f.name, "union") == 0) {
	if (ptr) *ptr = new composite_object(CGS_UNION, f, p_invert);
	if (type) *type = CGS_COMPOSITE;
    } else if (strcmp(f.name, "intersect") == 0) {
	if (ptr) *ptr = new composite_object(CGS_INTERSECT, f, p_invert);
	if (type) *type = CGS_COMPOSITE;
    } else if (strcmp(f.name, "difference") == 0) {
	if (ptr) *ptr = new composite_object(CGS_DIFFERENCE, f, p_invert);
	if (type) *type = CGS_COMPOSITE;
    } else if (strcmp(f.name, "data") == 0) {
	if (ptr) *ptr = new composite_object(CGS_CMB_NOOP, f, p_invert);
	if (type) *type = CGS_DATA;
    } else {
	if (ptr) *ptr = NULL;
	if (type) *type = CGS_UNDEF;
    }

    return E_SUCCESS;
}

/**
 * Produce an appropriate transformation matrix
 */
parse_ercode scene::make_transformation(const cgs_func& f, emat3& res) const {
    parse_ercode er = E_SUCCESS;

    if (!f.name) return E_BAD_TOKEN;
    //switch between all potential types
    if (strcmp(f.name, "rotate") == 0) {
	if (f.n_args < 2) return E_LACK_TOKENS;
	if (f.args[0].type != VAL_NUM || f.args[1].type != VAL_3VEC) return E_BAD_VALUE;
	double angle = f.args[0].val.x;
	res = Eigen::AngleAxisd(angle, *(f.args[1].val.v));
    } else if (strcmp(f.name, "scale") == 0) {
	if (f.n_args < 1) return E_LACK_TOKENS;
	evec3 scale;
	if (f.args[0].type == VAL_3VEC) {
	    scale = *(f.args[0].val.v);
	} else if (f.args[0].type == VAL_NUM) {
	    double val = f.args[0].val.x;
	    scale.x() = val;
	    scale.y() = val;
	    scale.z() = val;
	}
	res = scale.asDiagonal();
    }

    return E_SUCCESS;
}

/** ============================ object_stack ============================ **/

/**
 * Insert an object onto the stack obeying the desired tree structure.
 */
parse_ercode object_stack::emplace_obj(object* obj, object_type p_type) {
    if (stack_ptr == 0) {
	if (p_type == CGS_COMPOSITE || p_type == CGS_ROOT || p_type == CGS_DATA) {
	    side_obj_pair cur(0, (composite_object*)obj);
	    return push(cur);
	} else {
	    return E_EMPTY_STACK;
	}
    } else {
	size_t ind = stack_ptr - 1;
	if (buf[ind].side < 2) {
	    buf[ind].obj->add_child(buf[ind].side, obj, p_type);
	    buf[ind].side += 1;
	} else {
	    //if we're at a leaf then we need to walk up the tree until we find a location where we can place the object
	    while (buf[ind].side > 1) {
		//if the stack pointer is at the root then this isn't a valid binary operation
		if (ind == 0) return E_NOT_BINARY;
		--ind;
	    }
	    buf[ind].obj->add_child(buf[ind].side, obj, p_type);
	    buf[ind].side += 1;
	    //we have to update the stack pointer appropriatly
	    stack_ptr = ind + 1;
	}

	//if it's a composite object then we'll need to add its children, so push it onto the stack
	if (p_type == CGS_COMPOSITE || p_type == CGS_ROOT) {
	    side_obj_pair cur(0, (composite_object*)obj);
	    return push(cur);
	}
    }

    return E_SUCCESS;
}

_uint object_stack::look_side() {
    if (stack_ptr == 0) return SIDE_UNDEF;
    return buf[stack_ptr - 1].side;
}

composite_object* object_stack::look_obj() {
    if (stack_ptr == 0) return NULL;
    return buf[stack_ptr].obj;
}

composite_object* object_stack::get_root() {
    if (stack_ptr == 0) return NULL;
    return buf[0].obj;
}

parse_ercode scene::fail_exit(parse_ercode er, FILE* fp) {
    fclose(fp);
    return er;
}

instance::instance(cgs_func decl) : type(decl.name) {
    fields.reserve(decl.n_args);
    value def_val;def_val.n_els = 0;def_val.type = VAL_UNDEF;def_val.val.x = 0;
    for (size_t i = 0; i < decl.n_args; ++i) {
	//if a default value was supplied, use that. Otherwise use an undefined initial value
	if (decl.args[i].type != VAL_UNDEF) {
	    fields.emplace_back(decl.arg_names[i], decl.args[i]);
	} else {
	    fields.emplace_back(decl.args[i].val.s, def_val);
	}
    }
}

parse_ercode scene::read_file(const char* p_fname) {
    parse_ercode er = E_SUCCESS;

    //generate the composite object on the fly
    object_stack tree_pos;
    stack<block_type> blk_stack;
    block_type last_type = BLK_UNDEF;
    composite_object* last_comp = NULL;
    int invert = 0;
    //keep track of transformations
    stack<emat3> transform_stack;
    emat3 cur_trans_mat;
    emat3 next_trans_mat;

    //open the file for reading and read individual lines. We need to remove whitespace. Handily, we know the line length once we're done.
    FILE* fp = NULL;
    if (p_fname) {
        fp = fopen(p_fname, "r");
    }
    if (fp) {
	char* buf = (char*)malloc(sizeof(char)*BUF_SIZE);
	size_t buf_size = BUF_SIZE;
	size_t lineno = 1;size_t next_lineno = 1;
	int line_len = read_cgs_line(&buf, &buf_size, fp, &next_lineno);
	char last_char = 0;
	//iterate over each line in the file
	//while (fgets(buf, BUF_SIZE, fp)) {//}
	while (line_len >= 0) {
	    //char* red_str = CGS_trim_whitespace(buf, &line_len);
	    for (int i = 0; i < line_len && buf[i]; ++i) {
		//check the most recent block pushed onto the stack
		block_type cur_type = BLK_UNDEF;
		if (blk_stack.size() > 0) cur_type = blk_stack.peek();
		//only interpret as normal code if we aren't in a comment or literal block
		if (cur_type != BLK_COMMENT && cur_type != BLK_LITERAL) {
		    //check for function calls
		    if (buf[i] == '(' && blk_stack.peek() != BLK_LITERAL) {
			//initialize a new cgs_func with the appropriate arguments
			cgs_func cur_func;
			char* endptr;
			cur_func = named_items.parse_func(buf, i, er, &endptr);
			switch (er) {
			    case E_BAD_TOKEN:
				printf("Error on line %d: Invalid function name \"%s\"\n", lineno, buf);
				free(buf);
				cleanup_func(&cur_func);
				return fail_exit(er, fp);
			    case E_BAD_SYNTAX:
				printf("Error on line %d: Invalid syntax\n", lineno);
				free(buf);
				cleanup_func(&cur_func);
				return fail_exit(er, fp);
			    default: break;
			}
			//check for class and function declarations
			char* dectype_start = token_block(cur_func.name, "def");
			if (dectype_start) {
			    cur_func.name = CGS_trim_whitespace(cur_func.name + KEY_DEF_LEN, NULL);

			    //TODO: finish
			} else if (dectype_start = token_block(cur_func.name, "class")) {
			    cur_func.name = CGS_trim_whitespace(cur_func.name + KEY_CLASS_LEN, NULL);
			    value tmp_val;
			    tmp_val.type = VAL_INST;
			    tmp_val.val.i = new instance(cur_func);
			    tmp_val.n_els = cur_func.n_args;
			}
			//try interpreting the function as a geometric object
			object* obj = NULL;
			object_type type;
			make_object(cur_func, &obj, &type, invert);
			if (!obj) {
			    emat3 tmp;
			    //if that failed try interpreting it as an operation (TODO: are there any useful things to put here?)
			    if (strcmp(cur_func.name, "invert") == 0) {
				last_type = BLK_INVERT;
			    } else if ((er = make_transformation(cur_func, tmp)) == E_SUCCESS) {
				last_type = BLK_TRANSFORM;
				next_trans_mat = tmp*cur_trans_mat;
			    }
			} else {
			    if (type == CGS_ROOT) {
				if (!tree_pos.is_empty()) {
				    printf("Error on line %d: Root composites may not be nested\n", lineno);
				} else {
				    last_comp = (composite_object*)obj;
				    last_type = BLK_ROOT;
				    //this is included so that we don't have to check whether something is a root or a composite every time
				    type = CGS_COMPOSITE;
				}
			    } else if (type == CGS_DATA) {
				data_objs.push_back((composite_object*)obj);
				last_type = BLK_DATA;
				type = CGS_COMPOSITE;
			    }
			    tree_pos.emplace_obj(obj, type);
			}
			//jump ahead until after the end of the function
			if (er == E_SUCCESS) i = endptr - buf;
			cleanup_func(&cur_func);
		    //check for comments
		    } else if (buf[i] == '/') {
			if (i < line_len-1) {
			    if (buf[i+1] == '/')
				break;
			    else if (buf[i+1] == '*')
				blk_stack.push(BLK_COMMENT);
			}
		    //check for blocks
		    } else if (buf[i] == '{') {
			switch (last_type) {
			    case BLK_INVERT: invert = 1 - invert;break;
			    case BLK_ROOT: roots.push_back(last_comp);break;
			    case BLK_TRANSFORM: cur_trans_mat=next_trans_mat;transform_stack.push(next_trans_mat);break;
			    default: break;
			}
			blk_stack.push(last_type);
			last_type = BLK_UNDEF;
		    } else if (buf[i] == '}') {
			block_type bt;
			if (blk_stack.pop(&bt) == E_EMPTY_STACK) printf("Error on line %d: unexpected '}'\n", lineno);
			switch (bt) {
			    case BLK_INVERT: invert = 1 - invert;break;
			    case BLK_ROOT: tree_pos.reset();
			    case BLK_TRANSFORM: transform_stack.pop(&cur_trans_mat);break;
			    default: break;
			}
		    //check for literal experessions enclosed in quotes
		    } else if (buf[i] == '\"') {
			if (blk_stack.peek() == BLK_LITERAL) {
			    blk_stack.pop(NULL);
			} else {
			    blk_stack.push(BLK_LITERAL);
			}
		    } else if (buf[i] == '=') {
			size_t tok_len;
			buf[i] = 0;
			char* tok = CGS_trim_whitespace(buf, &tok_len);
			size_t val_len;
			char* val = CGS_trim_whitespace(buf+i+1, &val_len);
			value v = named_items.parse_value(val, er);
			if (er == E_SUCCESS) named_items.set_value(tok, v);
			cleanup_val(&v);
		    }
		} else {
		    //check if we reached the end of a comment or string literal block
		    if (cur_type == BLK_COMMENT && (buf[i] == '*' && i < line_len-1 && buf[i+1] == '/')) {
			er = blk_stack.pop(NULL);
		    } else if (cur_type == BLK_LITERAL && (buf[i] == '\"' && last_char != '\\')) {
			er = blk_stack.pop(NULL);
		    }
		}
		last_char = buf[i];
	    }
	    //don't clutter up the tree if we have a global data object
	    if (blk_stack.is_empty()) {
		tree_pos.reset();
	    }
	    line_len = read_cgs_line(&buf, &buf_size, fp, &next_lineno);
	    lineno = next_lineno;
	}
	free(buf);
	fclose(fp);
    } else {
        printf("Error: couldn't open file %s for reading!\n", p_fname);
        return E_NOFILE;
    }

    return er;
}

scene::scene(const char* p_fname, parse_ercode* ercode) {
    if (ercode) *ercode = E_SUCCESS;
    *ercode = read_file(p_fname);
}

scene::scene(const scene& o) {
    roots.resize(o.roots.size());
    data_objs.resize(o.data_objs.size());
    for (_uint i = 0; i < roots.size(); ++i) {
	roots[i] = new composite_object( *(o.roots[i]) );
    }
    for (_uint i = 0; i < data_objs.size(); ++i) {
	data_objs[i] = new composite_object( *(o.data_objs[i]) );
    }
    named_items = o.named_items;
}

scene::scene(scene&& o) {
    roots = o.roots;
    data_objs = o.data_objs;
    named_items = o.named_items;
    o.roots.clear();
    o.data_objs.clear();
    o.named_items.reset();
}

scene& scene::operator=(scene& o) {
    roots = o.roots;
    data_objs = o.data_objs;
    named_items = o.named_items;

    return *this;
}

scene::~scene() {
    //TODO: double frees are bad lol
    for (_uint i = 0; i < roots.size(); ++i) {
	if (roots[i]) delete roots[i];
    }
    for (_uint i = 0; i < data_objs.size(); ++i) {
	if (data_objs[i]) delete data_objs[i];
    }
}
