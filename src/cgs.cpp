#include "cgs.hpp"

/** ======================================================== object ======================================================== **/

object::object(int p_invert) {
    invert = p_invert;
    trans_mat = mat3x3::id();
}

void object::set_inversion(int p_invert) {
    if (p_invert == 0)
	invert = 0;
    else
	invert = 1;
}

void object::rescale(const vec3& components) {
    trans_mat = trans_mat*(components.make_diagonal());
}

/** ======================================================== sphere ======================================================== **/

sphere::sphere(vec3& p_center, double p_rad, int p_invert) : object(p_invert) {
    center = p_center;
    rad = p_rad;
}

sphere::sphere(vec3& p_center, vec3& p_rad, int p_invert) : object(p_invert) {
    center = p_center;
    rad = 1;
    rescale( vec3(1/p_rad.x(), 1/p_rad.y(), 1/p_rad.z()) );
}

int sphere::in(const vec3& r) {
    //find the relative offset between the input vector and the center of the box
    vec3 r_rel = trans_mat*(r - center);
    if (r_rel.norm() > rad) return invert;
    return 1 - invert;
}

box::box(vec3& p_corner_1, vec3& p_corner_2, int p_invert) : object(p_invert) {
    offset = (p_corner_2 - p_corner_1)/2;
    center = p_corner_1 + offset;
    //wlog, fix the offset to always be positive from the center
    if (offset.x() < 0) offset.x() *= -1;
    if (offset.y() < 0) offset.y() *= -1;
    if (offset.z() < 0) offset.z() *= -1;
}

/** ======================================================== box ======================================================== **/

int box::in(const vec3& r) {
    //find the relative offset between the input vector and the center of the box
    vec3 r_rel = trans_mat*(r - center);
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
plane::plane(vec3& p_normal, double p_offset, int p_invert) : object(0) {
    if (p_invert)
	offset = -p_offset;
    else
	offset = p_offset;
    normal = p_normal / p_normal.norm();
}

/**
 * Construct a plane defined by three points. The normal vector is fixed to be (point_2 - point_1)X(point_3 - point_1). Thus, flipping the order of the arguments point_1 and point_2 flips the direction of inclusion.
 */
plane::plane(vec3& point_1, vec3& point_2, vec3& point_3, int p_invert) : object(0) {
    vec3 diff_2 = point_2 - point_1;
    vec3 diff_3 = point_3 - point_1;
    normal = diff_2.cross(diff_3).normalize();
    offset = normal.dot(point_1);
    if (p_invert)
	offset *= -1;
}

int plane::in(const vec3& r) {
    double norm_comp = r.dot(normal);
    if (norm_comp > offset)
	return 0;
    else
	return 1;
}

/** ======================================================== cylinder ======================================================== **/

cylinder::cylinder(vec3& p_center, double p_height, double p_r1, double p_r2, int p_invert) : object(p_invert) {
    center = p_center;
    height = p_height;
    r1_sq = p_r1*p_r1;
    r1_sq_x_h = r1_sq*height;
    r2_sq = p_r2*p_r2;
}

int cylinder::in(const vec3& r) {
    //find the relative offset between the input vector and the center of the box
    vec3 r_rel = trans_mat*(r - center);
    if (r_rel.z() < 0 || r_rel.z() > height) return invert;
    double cent_rad_sq = r_rel.x()*r_rel.x() + r_rel.y()*r_rel.y();
    if (cent_rad_sq*height > (r2_sq - r1_sq)*r_rel.z() + r1_sq_x_h) return invert;

    return 1 - invert;
}

/** ======================================================== composite_object ======================================================== **/

parse_ercode make_object(context* inst, object** ptr, object_type* type, int p_invert, stack<mat3x3>& transform_stack) {
    parse_ercode er = E_SUCCESS;
    value type_val = inst->lookup("__type__");
    if (type_val.type == VAL_STR) {
	if (strcmp(type_val.val.s, "Box") == 0) {
	    parse_ercode er1;
	    value pt_1 = inst->lookup("pt_1").cast_to(VAL_3VEC, er1);
	    value pt_2 = inst->lookup("pt_2").cast_to(VAL_3VEC, er);
	    if (er != E_SUCCESS || er1 != E_SUCCESS) {
		printf("Error: Box() supplied value that was not a valid vector");
		if (er == E_SUCCESS) er = er1;
		goto clean_box;
	    }
	    if (ptr) *ptr = new box(*(pt_1.val.v), *(pt_2.val.v), p_invert);
	    if (type) *type = CGS_BOX;
clean_box:
	    cleanup_val(&pt_1);
	    cleanup_val(&pt_2);
	    return er;
	} else if (strcmp(type_val.val.s, "Plane") == 0) {
	    value offset = inst->lookup("offset");
	    if (offset.type != VAL_NUM) { er = E_BAD_TYPE; return er; }
	    value normal = inst->lookup("normal").cast_to(VAL_3VEC, er);
	    if (er != E_SUCCESS)
		goto clean_plane;
	    if (ptr) *ptr = new plane(*(normal.val.v), offset.val.x, p_invert);
	    if (type) *type = CGS_PLANE;
clean_plane:
	    cleanup_val(&normal);
	    return er;
	} else if (strcmp(type_val.val.s, "Sphere") == 0) {
	    value radius = inst->lookup("radius");
	    if (radius.type != VAL_NUM) { er = E_BAD_TYPE; return er; }
	    value center = inst->lookup("center").cast_to(VAL_3VEC, er);
	    if (er != E_SUCCESS)
		goto clean_sphere;
	    if (ptr) *ptr = new sphere(*(center.val.v), radius.val.x, p_invert);
	    if (type) *type = CGS_SPHERE;
clean_sphere:
	    cleanup_val(&center);
	    return E_SUCCESS;
	} else if (strcmp(type_val.val.s, "Cylinder") == 0) {
	    value h = inst->lookup("h");
	    value r1 = inst->lookup("r1");
	    value r2 = inst->lookup("r2");
	    if (h.type != VAL_NUM || r1.type != VAL_NUM || r2.type != VAL_NUM) { er = E_BAD_TYPE; return er; }
	    value center = inst->lookup("center").cast_to(VAL_3VEC, er);
	    if (er != E_SUCCESS)
		goto clean_cylinder;
	    if (ptr) *ptr = new cylinder(*(center.val.v), h.val.x, r1.val.x, r2.val.x, p_invert);
	    if (type) *type = CGS_CYLINDER;
clean_cylinder:
	    cleanup_val(&center);
	    return E_SUCCESS;
	} else if (strcmp(type_val.val.s, "Union") == 0) {
	    value gv = inst->lookup("geometry");
	    if (gv.type != VAL_LIST) { er = E_BAD_TYPE; return er; }
	    if (ptr) *ptr = new composite_object(CGS_UNION, gv.val.l, gv.n_els, p_invert, transform_stack);
	    if (type) *type = CGS_COMPOSITE;
	    return E_SUCCESS;
	} else if (strcmp(type_val.val.s, "Intersect") == 0) {
	    value gv = inst->lookup("geometry");
	    if (gv.type != VAL_LIST) { er = E_BAD_TYPE; return er; }
	    if (ptr) *ptr = new composite_object(CGS_INTERSECT, gv.val.l, gv.n_els, p_invert, transform_stack);
	    if (type) *type = CGS_COMPOSITE;
	    return E_SUCCESS;
	} else if (strcmp(type_val.val.s, "Complement") == 0) {
	    if (ptr) *ptr = NULL;
	    if (type) *type = CGS_COMP_INVERT;
	    return E_SUCCESS;
	} else if (strcmp(type_val.val.s, "Rotate") == 0) {
	    if (ptr) *ptr = NULL;
	    if (type) *type = CGS_COMP_ROTATE;
	    return E_SUCCESS;
	}
    }
    if (ptr) *ptr = NULL;
    if (type) *type = CGS_UNDEF;
    return E_BAD_VALUE;
}

composite_object::composite_object(combine_type p_cmb) : object(false) {
    children[0] = NULL;
    children[1] = NULL;
    cmb = p_cmb;
}

void composite_object::append(composite_object** lc_ptr, object* obj, object_type type, bool is_last) {
    if (!lc_ptr) return;
    if ((*lc_ptr)->get_child_l() == NULL) {
	(*lc_ptr)->add_child(0, obj, type);
    } else if (is_last) {
	(*lc_ptr)->add_child(1, obj, type);
    } else {
	composite_object* ncobj = new composite_object(cmb);
	(*lc_ptr)->add_child(1, ncobj, CGS_COMPOSITE);
	*lc_ptr = ncobj;
	(*lc_ptr)->add_child(0, obj, type);
    }
}

void composite_object::init_from_list(value* l, size_t n, stack<mat3x3>& transform_stack) {
    //setup reading of the geometry
    object_stack tree_pos;
    stack<block_type> blk_stack;
    block_type last_type = BLK_UNDEF;
    composite_object* last_comp = this;
    int invert = 0;
    //keep track of transformations
    mat3x3 trans_mat = mat3x3::id();
    for (size_t i = 1; i <= transform_stack.size(); ++i) {
	trans_mat = trans_mat*transform_stack.peek(i);
    }
    for (size_t i = 0; i < n; ++i) {
	value g_obj = l[i];
	if (g_obj.type == VAL_INST) {
	    object* obj = NULL;
	    object_type type;
	    parse_ercode er = make_object(g_obj.val.c, &obj, &type, invert, transform_stack);
	    if (er) {
		printf("Encountered error generating composite object code %d.\n", er);
		break;
	    }
	    if (!obj) {
		//transformation objects are a special case that we handle here
		if (type) {
		    if (type == CGS_COMP_INVERT) {
			invert = 1 - invert;
			blk_stack.push(BLK_INVERT);
			obj = new composite_object(cmb, g_obj.val.c, 1-invert, transform_stack);
			append(&last_comp, obj, CGS_COMPOSITE, (i==n-1));
		    } else if (type  == CGS_COMP_ROTATE) {
			value va = g_obj.val.c->lookup("axis");
			value vt = g_obj.val.c->lookup("theta");
			if (va.type == VAL_3VEC && vt.type == VAL_NUM) {
			    mat3x3 tmp = make_rotation(vt.val.x, *(va.val.v));
			    transform_stack.push(tmp);
			    obj = new composite_object(cmb, g_obj.val.c, invert, transform_stack);
			    append(&last_comp, obj, CGS_COMPOSITE, (i==n-1));		    
			}
		    }
		}
	    } else {
		obj->set_trans_mat(trans_mat);
		append(&last_comp, obj, type, (i==n-1));
	    }
	}
    }
}

composite_object::composite_object(combine_type p_cmb, const cgs_func& spec, int p_invert) : object(p_invert) {
    children[0] = NULL;
    children[1] = NULL;
    cmb = p_cmb;
    //iterate through arguments and add them to the metadata
    for (_uint i = 0; i < spec.n_args; ++i) {
	if (spec.arg_names[i]) {
	    //each piece of metadata should be separated by a named list
	    std::string tok_cpp(spec.arg_names[i]);
	    metadata[tok_cpp] = copy_val(spec.args[i]);
	}
    }
}

composite_object::composite_object(combine_type p_cmb, value* list, size_t n_els, int p_invert, stack<mat3x3>& transform_stack) : object(p_invert) {
    children[0] = NULL;
    children[1] = NULL;
    child_types[0] = CGS_UNDEF;
    child_types[1] = CGS_UNDEF;
    cmb = p_cmb;
    init_from_list(list, n_els, transform_stack);
}

composite_object::composite_object(combine_type p_cmb, context* inst, int p_invert, stack<mat3x3>& transform_stack) : object(p_invert) {
    children[0] = NULL;
    children[1] = NULL;
    child_types[0] = CGS_UNDEF;
    child_types[1] = CGS_UNDEF;
    cmb = p_cmb;
    //read metadata
    size_t n_args = inst->size()-1;
    for (size_t i = 2; i <= n_args; ++i) {
	name_val_pair p = inst->peek(i);
	std::string tok_cpp(p.get_name());
	metadata[tok_cpp] = copy_val(p.get_val());
    }
    //get the geometry data
    value tmp_val = inst->peek_val();
    if (tmp_val.type == VAL_LIST) {
	init_from_list(tmp_val.val.l, tmp_val.n_els, transform_stack);
    } 
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

int composite_object::call_child_in(_uint side, const vec3& r) {
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

int composite_object::in(const vec3& r) {
    vec3 r_trans = trans_mat*r;
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

parse_ercode scene::read_file(const char* p_fname) {
    parse_ercode er = E_SUCCESS;

    size_t init_size = named_items.size();
    line_buffer lb(p_fname);
    er = named_items.read_from_lines(lb);
    if (er != E_SUCCESS) return er;

    for (size_t i = 1; i <= named_items.size() - init_size; ++i) {
	value inst = named_items.peek_val(i);
	//only inspect instances
	if (inst.type == VAL_INST) {
	    value type_val = inst.val.c->lookup("__type__");
	    //branch depending on instance type
	    if (type_val.type == VAL_STR) {
		if (strcmp(type_val.val.s, "Composite") == 0) {
		    //read the combination type and default to a union
		    combine_type cmb_type = CGS_UNION;
		    type_val = inst.val.c->lookup("combine_type");
		    if (type_val.type == VAL_STR) {
			if (strcmp(type_val.val.s, "intersect") == 0 || strcmp(type_val.val.s, "Intersect") == 0)
			    cmb_type = CGS_INTERSECT;
		    }
		    stack<mat3x3> transform_stack;
		    composite_object* ncobj = new composite_object(cmb_type, inst.val.c, false, transform_stack);
		    roots.push_back(ncobj);
		} else if (strcmp(type_val.val.s, "snapshot") == 0) {
		    value v_f = inst.val.c->lookup("fname");
		    //read vectors and cast them
		    value v_c = inst.val.c->lookup("cam_v").cast_to(VAL_3VEC, er);
		    if (er) { return er; }
		    value v_l = inst.val.c->lookup("look_v").cast_to(VAL_3VEC, er);
		    if (er) { cleanup_val(&v_c);return er; }
		    value v_u = inst.val.c->lookup("up_v").cast_to(VAL_3VEC, er);
		    if (er) { cleanup_val(&v_c);cleanup_val(&v_l);return er; }
		    //read numeric values
		    value v_res = inst.val.c->lookup("res");
		    value v_n = inst.val.c->lookup("n_samples");
		    value v_st = inst.val.c->lookup("step");
		    value v_sc = inst.val.c->lookup("scale");
		    if (v_f.type != VAL_STR) { er = E_BAD_TYPE;return er; }
		    if (v_res.type != VAL_NUM) { er = E_BAD_TYPE;return er; }
		    if (v_n.type != VAL_NUM) { er = E_BAD_TYPE;return er; }
		    if (v_st.type != VAL_NUM) { er = E_BAD_TYPE;return er; }
		    if (v_sc.type != VAL_NUM) { er = E_BAD_TYPE;return er; }
		    if (v_c.val.v->norm() == 0 || v_u.val.v->norm() == 0 || v_res.val.x == 0 || v_n.val.x == 0) {
			printf("invalid parameters passed to snapshot on line\n");
		    } else {
			//make the drawing
			rvector<2> scale_vec;
			scale_vec.el[0] = v_sc.val.x;
			scale_vec.el[1] = v_sc.val.x;
			size_t res = (size_t)(v_res.val.x);
			size_t n_samples = (size_t)(v_n.val.x);
			draw(v_f.val.s, *(v_c.val.v), *(v_l.val.v), *(v_u.val.v), scale_vec, res, res, n_samples, v_st.val.x);
		    }
		    cleanup_val(&v_c);
		    cleanup_val(&v_l);
		    cleanup_val(&v_u);
		}
	    }
	}
    }

    return er;
}

scene::scene(const char* p_fname, parse_ercode* ercode) {
    setup_geometry_context(named_items);
    if (ercode) *ercode = E_SUCCESS;
    *ercode = read_file(p_fname);
}

scene::scene(const char* p_fname, context con, parse_ercode* ercode) : named_items(con) {
    setup_geometry_context(named_items);
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

/*
 * convert from hsv to rgb
 * h,s,v hue satureation and value
 * rgb: a pointer to an array of three unsigned chars where the first second and third elements represent r,g, and b respectively
 */
uint32_t hsvtorgb(uint32_t hsv_col) {
    int g_start = 85;//floor(256/3)
    int b_start = 170;//floor(256*2/3)
    //int r_start = 255;
    uint32_t r,g,b;
    int h = get_r(hsv_col);
    int s = get_g(hsv_col);
    int v = get_b(hsv_col);
    int s_comp = 255-s;
    if (h < g_start) {
	r = (g_start-h)*v/g_start;
	g = h*v/g_start;
	r = (r*s_comp)/256 + 1;
	g = (g*s_comp)/256 + 1;
	b = s;
    } else if (h < b_start) {
	int h_rel = h-g_start;
	g = (g_start-h_rel)*v/g_start;
	b = h_rel*v/g_start;
	g = (g*s_comp)/256 + 1;
	b = (b*s_comp)/256 + 1;
	r = s;
    } else {
	int h_rel = h-b_start;
	b = (g_start-h_rel)*v/g_start;
	r = h_rel*v/g_start;
	b = (b*s_comp)/256 + 1;
	r = (r*s_comp)/256 + 1;
	g = s;
    }
    return (hsv_col & 0xff000000) | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
}
/**
 * Generate a list of hues (in hsv color space) associated with each object in the scene. This is a helper function for draw()
 */
std::vector<uint32_t> scene::get_cols() {
    std::vector<uint32_t> ret(roots.size());
    for (size_t k = 0; k < roots.size(); ++k) {
	uint32_t hue = (_uint8)((k*256)/roots.size());
	uint32_t alpha = 255;
	//check if the user specified a valid color/transparency for this object
	if (roots[k]->has_metadata("color")) {
	    value col_val = roots[k]->fetch_metadata("color");
	    if (col_val.type == VAL_NUM && col_val.val.x > 0 && col_val.val.x < 256) hue = (_uint8)col_val.val.x;
	}
	if (roots[k]->has_metadata("alpha")) {
	    value a_val = roots[k]->fetch_metadata("alpha");
	    if (a_val.type == VAL_NUM && a_val.val.x > 0 && a_val.val.x < 256) alpha = (_uint8)a_val.val.x;
	}
	ret[k] = ((alpha & 0xff) << 24) | ((hue & 0xff) << 16);
    }
    return ret;
}
/**
 * Given two 32 bit colors, c1 and c2, blend both of them with c1 on top of c2 and return the result.
 */
uint32_t blend(uint32_t c1, uint32_t c2) {
    //we use fixed point arithmetic with divisor 2^8, so a1=256->1.0 or a1=128->0.5.
    uint32_t a1 = get_a(c1);
    uint32_t a2 = get_a(c2);
    //special cases. If the top color has alpha=0, don't modify the bottom at all. If the top color has alpha=255 or the bottom color has alpha=0 then don't modify the top.
    if (a1 == 0)
	return c2;
    if (a1 == 255 || a2 == 0)
	return c1;
    //this is the expression for a1 + a2*(1-a1) = a1 + a2 - a1*a2 with divisor 2^16
    uint32_t p1 = (a1 << 8);
    uint32_t p2 = (a2 << 8) - a1*a2;
    uint32_t a0 = p1 + p2;
    //set r, g and b channels
    uint32_t ret = (get_b(c1)*p1 + get_b(c2)*p2) / a0;
    ret |= ( (get_g(c1)*p1 + get_g(c2)*p2) / a0 ) << 8;
    ret |= ( (get_r(c1)*p1 + get_r(c2)*p2) / a0 ) << 16;
    //set alpha channel
    ret |= ( a0 >> 8 ) << 24;
    return ret;
}
/*
 * Save the image specified by z_buf and c_buf to out_fname with the specified resolution
 * z_buf: the z buffer. This must have a size res*res
 * c_buf: the color index buffer. This must have a size res*res
 * res: the resolution of the image. Only square images are allowed
 */
void save_imbuf(const char* out_fname, uint32_t* c_buf, size_t res_x, size_t res_y) {
    //open the file. We use a .pgm format since it's super simple. Then we write the header information.
    FILE* fp;
    if (out_fname) {
	fp = fopen(out_fname, "w");
	if (!fp) fp = stdout;
    } else {
	fp = stdout;
    }
    uint32_t cur_col;
    fprintf(fp, "P3\n%lu %lu\n%d\n", res_x, res_y, IM_DEPTH);
    //iterate through the image buffer and write
    for (size_t yy = 0; yy < res_y; ++yy) {
	for (size_t xx = 0; xx < res_x; ++xx) {
	    uint32_t cur_col = c_buf[yy*res_y+xx];
	    //if this pixel is completely transparent, draw white. Otherwise draw the color
	    if (get_a(cur_col) == 0)
		fprintf(fp, "255 255 255 ");
	    else
		fprintf(fp, "%d %d %d ", get_r(cur_col), get_g(cur_col), get_b(cur_col));
	}
	fprintf(fp, "\n");
    }
    fclose(fp);
}
/**
 * Helper function for draw() which steps along the ray disp using the z buffer at indices n1 and n2 as a hint for the depth
 */
void scene::ray_step(_uint8* z_buf, _uint8* c_buf, size_t ind, size_t n1, size_t n2, vec3 disp) {

}
/**
 * Render the object and save the image to out_fname
 * out_fname: the filename for the picture
 * cam_pos: the position of the camera
 * cam_look: the direction the camera is looking.
 * cam_up: this determines the "up" direction for the camera. Note that the up vector is not directly used. Rather the vector in the cam_look,cam_up plane which is orthogonal to cam_look and has a positive dot product with cam_up is used.
 * transform: a transformation matrix applied to the points in the image.
 * res: the resolution of the saved image, only 1x1 aspect ratios are allowed
 * n_samples: the number of random samples to take, a higher value results in a cleaner image but takes longer to generate.
 */
void scene::draw(const char* out_fname, vec3 cam_pos, vec3 cam_look, vec3 cam_up, rvector<2> scale, size_t res_x, size_t res_y, size_t n_samples, double walk_step) {
    //depth scale factors and sampling region are determined by the magnitude of cam_look. Calculated these and then normalize cam_look
    double aa = 1/cam_look.norm();
    double max_draw_z = 2*cam_look.norm();
    cam_look = cam_look.normalize();
    //cross the look direction with the z directioncol_val.val.x > 0 
    vec3 x_comp = cam_look.cross(cam_up);
    //avoid divisions by zero
    if (x_comp.norm() == 0) {
	//if we're looking along the z axis, take the x component to be the x basis vector. Otherwise shift cam_up along the z axis and try again.
	if (cam_look.x() == 0 && cam_look.y() == 0) {
	    x_comp.x() = 1;
	} else {
	    cam_up.z() = cam_up.z() + 1;
	    x_comp = cam_look.cross(cam_up);
	}
    }
    x_comp = x_comp.normalize();
    vec3 y_comp = x_comp.cross(cam_look);
    scale.el[0] /= res_x;
    scale.el[1] /= res_y;
    //store the image buffer in an array. We don't do anything fancy for shading, this is just a z-buffer. Thus, we initialize to white (the maximum distance).
    _uint8 z_buf[res_x*res_y];
    uint32_t c_buf[res_x*res_y];
    for (size_t i = 0; i < res_x*res_y; ++i) {
	c_buf[i] = 0x00000000;
	z_buf[i] = 0xff;
    }
    //get object colors
    std::vector<uint32_t> cols = get_cols();

    double max_z = 0;
    double min_z = max_draw_z;

    //iterate across the pixels in the buffer
    for (long ii = 0; ii < res_x; ++ii) {
	for (long jj = res_y-1; jj >= 0; --jj) {
	    //orthographic projection
	    vec3 disp = scale.el[0]*(2*(double)ii-res_x)*x_comp + scale.el[1]*(2*(double)jj-res_y)*y_comp;
	    vec3 r0 = cam_pos+disp;
	    //take small steps until we hit the object
	    bool itz = true;
	    for (double z = 0; itz && z < max_draw_z; z += walk_step*max_draw_z) {
		vec3 r = r0 + z*cam_look;
		_uint8 depth = 255;
		uint32_t cur_col = 0;
		for (size_t k = 0; k < roots.size(); ++k) {
		    if (roots[k]->in(r)) {
			//if we hit the object then figure out the index in the image buffer and map the depth
			depth = floor( IM_DEPTH*(1 - exp(-z*aa)) );
			uint32_t this_col = hsvtorgb( cols[k] | (255 - depth) );
			//z_buf[ii + jj*res_y] = depth;	
			//in the case of a completely opaque object we may stop, otherwise we have to continue until we find the next object
			if (get_a(cols[k]) == 255) {
			    itz = false;
			    if (z > max_z)
				max_z = z;
			    if (z < min_z)
				min_z = z;
			    c_buf[ii + jj*res_y] = blend(cur_col, this_col);
			    break;
			} else {
			    cur_col = blend(cur_col, this_col);
			}
		    }
		}
		//if we reached the end and we only hit transparent objects, then write the color buffer
		if (itz && cur_col != 0) {
		    c_buf[ii + jj*res_y] = cur_col;
		}
	    }
	}
    }
    printf("max z: %f, min z: %f\n", max_z, min_z);
    save_imbuf(out_fname, c_buf, res_x, res_y);
}
void scene::draw(const char* out_fname, vec3 cam_pos) {
    vec3 cam_look = cam_pos*-1;
    vec3 cam_up(0,0,1);
    double xy_scale = 0.5 / cam_look.norm();
    rvector<2> scale;
    scale.el[0] = xy_scale;
    scale.el[1] = xy_scale;
    draw(out_fname, cam_pos, cam_look, cam_up, scale);
}
