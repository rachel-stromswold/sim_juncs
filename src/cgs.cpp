#include "cgs.hpp"

/**
 * Helper function for is_token which tests whether the character c is a token terminator
 */
bool is_char_sep(char c) {
    if (c == 0 || c == ' ' || c == '\t'  || c == '\n')
	return true;
    else
	return false;
}

/**
 * Helper function which looks at the string str at index i and tests whether it is a token with the matching name
 */
bool is_token(const char* str, size_t i, size_t len) {
    if (i > 0 && !is_char_sep(str[i-1])) return false;
    if (!is_char_sep(str[i+len])) return false;
    return true;
}

/**
 * Find the first index of the character c that isn't nested inside a block
 */
char* strchr_block(char* str, char c) {
    CGS_Stack<size_t> blk_stk;
    for (size_t i = 0; str[i] != 0; ++i) {
	if (str[i] == c && blk_stk.is_empty()) {
	    return str+i;
	}
	if (str[i] == '('/*)*/) {
	    blk_stk.push(i);
	} else if (str[i] == /*(*/')') {
	    blk_stk.pop(NULL);
	} else if (str[i] == '['/*]*/) {
	    blk_stk.push(i);
	} else if (str[i] == /*[*/']') {
	    blk_stk.pop(NULL);
	} else if (str[i] == '{'/*}*/) {
	    blk_stk.push(i);
	} else if (str[i] == /*{*/'}') {
	    blk_stk.pop(NULL);
	} else if (str[i] == '\"'/*"*/) {
	    //quotes are more complicated
	    if (!blk_stk.is_empty() && str[blk_stk.peek()] == '\"')
		blk_stk.pop(NULL);
	    else
		blk_stk.push(i);
	} else if (str[i] == '\''/*"*/) {
	    //quotes are more complicated
	    if (!blk_stk.is_empty() && str[blk_stk.peek()] == '\'')
		blk_stk.pop(NULL);
	    else
		blk_stk.push(i);
	}
    }
    return NULL;
}

/**
 * Find the first instance of a token (i.e. surrounded by whitespace) in the string str which matches comp
 */
char* token_block(char* str, const char* comp) {
    if (!str || !comp) return NULL;
    CGS_Stack<size_t> blk_stk;
    for (size_t i = 0; str[i] != 0; ++i) {
	if (blk_stk.is_empty()) {
	    size_t j = 0;
	    while (comp[j] && str[i+j] && str[i+j] == comp[j]) ++j;
	    if (comp[j] == 0 && is_token(str, i, j)) return str+i;
	}
	if (str[i] == '('/*)*/) {
	    blk_stk.push(i);
	} else if (str[i] == /*(*/')') {
	    blk_stk.pop(NULL);
	} else if (str[i] == '['/*]*/) {
	    blk_stk.push(i);
	} else if (str[i] == /*[*/']') {
	    blk_stk.pop(NULL);
	} else if (str[i] == '{'/*}*/) {
	    blk_stk.push(i);
	} else if (str[i] == /*{*/'}') {
	    blk_stk.pop(NULL);
	} else if (str[i] == '\"'/*"*/) {
	    //quotes are more complicated
	    if (!blk_stk.is_empty() && str[blk_stk.peek()] == '\"')
		blk_stk.pop(NULL);
	    else
		blk_stk.push(i);
	} else if (str[i] == '\''/*"*/) {
	    //quotes are more complicated
	    if (!blk_stk.is_empty() && str[blk_stk.peek()] == '\'')
		blk_stk.pop(NULL);
	    else
		blk_stk.push(i);
	}
    }
    return NULL;
}

/** ============================ Object ============================ **/

Object::Object(int p_invert) {
    invert = p_invert;
    trans_mat = emat3::Identity();
}

/*Object::Object(Eigen::Quaterniond& p_orientation, int p_invert) {
    invert = p_invert;
    trans_mat = p_orientation.toRotationMatrix();
}*/

void Object::set_inversion(int p_invert) {
    if (p_invert == 0)
	invert = 0;
    else
	invert = 1;
}

void Object::rescale(const evec3& components) {
    trans_mat = trans_mat*(components.asDiagonal());
}

Sphere::Sphere(evec3& p_center, double p_rad, int p_invert) : Object(p_invert) {
    center = p_center;
    rad = p_rad;
}

Sphere::Sphere(evec3& p_center, evec3& p_rad, int p_invert) : Object(p_invert) {
    center = p_center;
    rad = 1;
    rescale( evec3(1/p_rad.x(), 1/p_rad.y(), 1/p_rad.z()) );
}

int Sphere::in(const evec3& r) {
    //find the relative offset between the input vector and the center of the box
    evec3 r_rel = trans_mat*(r - center);
    if (r_rel.norm() > rad) return invert;
    return 1 - invert;
}

Box::Box(evec3& p_corner_1, evec3& p_corner_2, int p_invert) : Object(p_invert) {
    offset = (p_corner_2 - p_corner_1)/2;
    center = p_corner_1 + offset;
    //wlog, fix the offset to always be positive from the center
    if (offset.x() < 0) offset.x() *= -1;
    if (offset.y() < 0) offset.y() *= -1;
    if (offset.z() < 0) offset.z() *= -1;
}

/*Box::Box(evec3& p_corner_1, evec3& p_corner_2, Eigen::Quaterniond p_orientation, int p_invert) : Object(p_orientation, p_invert) {
    offset = (p_corner_2 - p_corner_1)/2;
    center = p_corner_1 + offset;
    //wlog, fix the offset to always be positive from the center
    if (offset.x() < 0) offset.x() *= -1;
    if (offset.y() < 0) offset.y() *= -1;
    if (offset.z() < 0) offset.z() *= -1;
}*/

int Box::in(const evec3& r) {
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
Plane::Plane(evec3& p_normal, double p_offset, int p_invert) : Object(0) {
    if (p_invert)
	offset = -p_offset;
    else
	offset = p_offset;
    normal = p_normal / p_normal.norm();
}

/**
 * Construct a plane defined by three points. The normal vector is fixed to be (point_2 - point_1)X(point_3 - point_1). Thus, flipping the order of the arguments point_1 and point_2 flips the direction of inclusion.
 */
Plane::Plane(evec3& point_1, evec3& point_2, evec3& point_3, int p_invert) : Object(0) {
    evec3 diff_2 = point_2 - point_1;
    evec3 diff_3 = point_3 - point_1;
    normal = diff_2.cross(diff_3);
    normal /= normal.norm();
    offset = normal.dot(point_1);
    if (p_invert)
	offset *= -1;
}

int Plane::in(const evec3& r) {
    double norm_comp = r.dot(normal);
    if (norm_comp > offset)
	return 0;
    else
	return 1;
}

Cylinder::Cylinder(evec3& p_center, double p_height, double p_r1, double p_r2, int p_invert) : Object(p_invert) {
    center = p_center;
    height = p_height;
    r1_sq = p_r1*p_r1;
    r1_sq_x_h = r1_sq*height;
    r2_sq = p_r2*p_r2;
}

int Cylinder::in(const evec3& r) {
    //find the relative offset between the input vector and the center of the box
    evec3 r_rel = trans_mat*(r - center);
    if (r_rel.z() < 0 || r_rel.z() > height) return invert;
    double cent_rad_sq = r_rel.x()*r_rel.x() + r_rel.y()*r_rel.y();
    if (cent_rad_sq*height > (r2_sq - r1_sq)*r_rel.z() + r1_sq_x_h) return invert;

    return 1 - invert;
}

Value make_val_str(const char* s) {
    Value v;
    v.type = VAL_STR;
    v.n_els = strlen(s) + 1;
    v.val.s = (char*)malloc(sizeof(char)*v.n_els);
    for (size_t i = 0; i < v.n_els; ++i) v.val.s[i] = s[i];
    return v;
}
Value make_val_std_str(std::string s) {
    Value v;
    v.type = VAL_STR;
    v.n_els = s.size();
    v.val.s = (char*)malloc(sizeof(char)*v.n_els);
    for (size_t i = 0; i < v.n_els; ++i) v.val.s[i] = s[i];
    return v;
}
Value make_val_list(const Value* vs, size_t n_vs) {
    Value v;
    v.type = VAL_LIST;
    v.n_els = n_vs;
    v.val.l = (Value*)malloc(sizeof(Value)*v.n_els);
    for (size_t i = 0; i < v.n_els; ++i) v.val.l[i] = vs[i];
    return v;
}
Value make_val_mat(Eigen::MatrixXd m) {
    Value v;
    v.type = VAL_MAT;
    v.n_els = 1;
    v.val.m = new Eigen::MatrixXd(m);
    return v;
}
Value make_val_vec3(evec3 vec) {
    Value v;
    v.type = VAL_3VEC;
    v.n_els = 1;
    v.val.v = new evec3(vec);
    return v;
}
void cleanup_val(Value* v) {
    if (v->type == VAL_STR && v->val.s) {
	free(v->val.s);
    } else if (v->type == VAL_LIST && v->val.l) {
	for (size_t i = 0; i < v->n_els; ++i) cleanup_val(v->val.l + i);
	free(v->val.l);
    } else if (v->type == VAL_MAT && v->val.m) {
	delete v->val.m;
    } else if (v->type == VAL_3VEC && v->val.v) {
	delete v->val.v;
    }
    v->val.x = 0;
}
Value copy_val(const Value o) {
    Value ret;
    ret.type = o.type;
    ret.n_els = o.n_els;
    //strings or lists must be copied
    if (o.type == VAL_STR) {
	ret.val.s = (char*)malloc(sizeof(char)*o.n_els);
	for (size_t i = 0; i < o.n_els; ++i) ret.val.s[i] = o.val.s[i];
    } else if (o.type == VAL_LIST) {
	ret.val.l = (Value*)calloc(o.n_els, sizeof(Value));
	for (size_t i = 0; i < o.n_els; ++i) ret.val.l[i] = copy_val(o.val.l[i]);
    } else if (o.type == VAL_MAT) {
	ret.val.m = new Eigen::MatrixXd(*(o.val.m));
    } else if (o.type == VAL_3VEC) {
	ret.val.v = new evec3(*(o.val.v));
    } else {
	ret.val.x = o.val.x;
    }
    return ret;
}
void swap_val(Value* a, Value* b) {
    //swap type and number of elements
    valtype tmp = a->type;
    a->type = b->type;
    b->type = tmp;
    size_t tmp_n = a->n_els;
    a->n_els = b->n_els;
    b->n_els = tmp_n;
    V tmp_v = a->val;
    a->val = b->val;
    b->val = tmp_v;
}
bool Value::operator==(std::string str) {
    if (type != VAL_STR) return false;
    size_t strlen = str.size();
    if (strlen != n_els-1) return false;
    for (size_t i = 0; i < strlen; ++i) {
	if (str[i] != val.s[i]) return false;
    }
    return true;
}
bool Value::operator!=(std::string str) {
    if (type != VAL_STR) return true;
    size_t strlen = str.size();
    if (strlen != n_els-1) return true;
    for (size_t i = 0; i < strlen; ++i) {
	if (str[i] != val.s[i]) return true;
    }
    return false;
}
char* Value::to_c_str() {
    if (type == VAL_STR)
	return val.s;
    return NULL;
}
double Value::to_float() {
    if (type == VAL_NUM)
	return val.x;
    return 0;
}
/**
 * Perform an in-place cast of the instance to the type t. An error is returned if a cast is impossible.
 */
Value Value::cast_to(valtype t, parse_ercode& er) const {
    er = E_SUCCESS;
    if (type == VAL_UNDEF) { er = E_BAD_VALUE;return copy_val(*this); }
    if (type == t) return copy_val(*this);
    Value ret;
    ret.type = t;
    ret.n_els = n_els;
    if (t == VAL_3VEC) {
	if (type == VAL_LIST) {
	    Value* tmp_lst = val.l;
	    //check that we have at least three numeric values
	    if (n_els < 3) { er = E_LACK_TOKENS;return copy_val(*this); }
	    if (tmp_lst[0].type != VAL_NUM || tmp_lst[1].type != VAL_NUM || tmp_lst[2].type != VAL_NUM) { er = E_BAD_TOKEN;return copy_val(*this); }
	    //actually change data and free the old
	    ret.val.v = new evec3(tmp_lst[0].val.x, tmp_lst[1].val.x, tmp_lst[2].val.x);
	    ret.n_els = 3;
	}
    } else if (t == VAL_LIST) {
	if (type == VAL_3VEC) {
	    ret.n_els = 3;
	    ret.val.l = (Value*)malloc(sizeof(Value)*ret.n_els);
	    for (size_t i = 0; i < ret.n_els; ++i) {
		ret.val.l[i].type = VAL_NUM;
		ret.val.l[i].n_els = 1;
		ret.val.l[i].val.x = (*val.v)(i);
	    }
	} else if (type == VAL_MAT) {
	    //TODO
	}
    } else if (type == VAL_STR) {
	//TODO
    }
    return ret;
}
void cleanup_func(cgs_func* f) {
    if (f) {
	for (size_t i = 0; i < f->n_args; ++i) cleanup_val(f->args + i);
    }
}
cgs_func copy_func(const cgs_func o) {
    cgs_func f;
    f.name = strdup(o.name);
    f.n_args = o.n_args;
    for (size_t i = 0; i < f.n_args; ++i) {
	f.args[i] = copy_val(o.args[i]);
	if (o.arg_names[i]) f.arg_names[i] = strdup(o.arg_names[i]);
    }
    return f;
}
void swap_func(cgs_func* a, cgs_func* b) {
    char* tmp_name = a->name;
    a->name = b->name;
    b->name = tmp_name;
    size_t tmp_n_args = a->n_args;
    a->n_args = b->n_args;
    b->n_args = tmp_n_args;
    //WLOG fix sf and lf to be the pointers to the functions smaller and larger number of arguments respectively
    cgs_func* sf = a;
    cgs_func* lf = b;
    size_t min_n_args = a->n_args;
    size_t max_n_args = b->n_args;
    if (b->n_args < a->n_args) {
	sf = b;
	lf = a;
	min_n_args = b->n_args;
	max_n_args = a->n_args;
    }
    //finally we can actually swap each element, a full swap only needs to be done for the minimal number of arguments
    size_t i = 0;
    for (; i < min_n_args; ++i) {
	Value tmp_val = sf->args[i];
	sf->args[i] = lf->args[i];
	lf->args[i] = tmp_val;
	tmp_name = sf->arg_names[i];
	sf->arg_names[i] = lf->arg_names[i];
	lf->arg_names[i] = tmp_name;
    }
    for (; i < max_n_args; ++i) {
	sf->args[i] = lf->args[i];
	sf->arg_names[i] = lf->arg_names[i];
    }
}

CompositeObject::CompositeObject(combine_type p_cmb) {
    children[0] = NULL;
    children[1] = NULL;
    cmb = p_cmb;
}

CompositeObject::CompositeObject(combine_type p_cmb, const cgs_func& spec, int p_invert) : Object(p_invert) {
    children[0] = NULL;
    children[1] = NULL;
    //iterate through arguments and add them to the metadata
    for (_uint i = 0; i < spec.n_args; ++i) {
	if (spec.arg_names[i]) {
	    //each piece of metadata should be separated by a named list
	    /*char* save_str;
	    char* tok = strtok_r(spec.args[i], "=", &save_str);
	    if (tok) {
		char* val = strtok_r(NULL, "=", &save_str);
		//only insert a metadata entry if all information was valid
		if (val) {
		    tok = CGS_trim_whitespace(tok, NULL);val = CGS_trim_whitespace(val, NULL);
		    std::string tok_cpp(tok);std::string val_cpp(val);
		    metadata[tok_cpp] = val_cpp;
		    free(tok);free(val);
		}
	    }*/
	    std::string tok_cpp(spec.arg_names[i]);
	    metadata[tok_cpp] = copy_val(spec.args[i]);
	}
    }
    cmb = p_cmb;
}

/**
 * Helper function which copies the child on the specified side. The caller is responsible for calling delete on the object and managing its lifetime.
 */
Object* CompositeObject::copy_child(_uint side) const {
    //declare these outside of the switch statement so we don't bork the addresses
    CompositeObject* my_child = NULL;
    CompositeObject* comp = NULL;
    switch (child_types[side]) {
	//if the type is not a composite then copying is easy
	case CGS_BOX: return new Box( *((Box*)children[side]) );break;
	case CGS_SPHERE: return new Sphere( *((Sphere*)children[side]) );break;
	case CGS_CYLINDER: return new Cylinder( *((Cylinder*)children[side]) );break;
       //otherwise...
	case CGS_COMPOSITE: 
	    //allocate a new child and store a pointer of this child for convenience
	    my_child = (CompositeObject*)(children[side]);
	    comp = new CompositeObject(my_child->cmb);
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
CompositeObject::CompositeObject(const CompositeObject& o) {
    //copy the "easy" stuff
    metadata = o.metadata;
    child_types[0] = o.child_types[0];
    child_types[1] = o.child_types[1];
    //copying children is potentially thorny. We have a helper function to do it for us
    children[0] = o.copy_child(0);
    children[1] = o.copy_child(1);
}

//move constructor
CompositeObject::CompositeObject(CompositeObject&& o) {
    metadata = o.metadata;
    //do this as a loop to avoid typing :p
    for (_uint side = 0; side < 2; ++side) {
	child_types[side] = o.child_types[side];
	children[side] = o.children[side];
	o.children[side] = NULL;
    }
}

//assignemnt operator
CompositeObject& CompositeObject::operator=(CompositeObject&& o) {
    //swap the easy stuff
    std::unordered_map<std::string, Value> tmp_metadata = metadata;
    metadata = o.metadata;
    o.metadata = tmp_metadata;
    //do this as a loop to avoid typing :p
    for (_uint side = 0; side < 2; ++side) {
	//save the current
	object_type tmp_ctype = child_types[side];
	Object* tmp_obj = children[side];
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
CompositeObject::~CompositeObject() {
    for (_uint i = 0; i < 2; ++i) {
	if (children[i]) {
	    if (child_types[i] == CGS_COMPOSITE)
		delete (CompositeObject*)(children[i]);
	    else
		delete children[i];
	}
    }
    for (auto it = metadata.begin(); it != metadata.end(); ++it) {
	cleanup_val(&(it->second));
    }
}

int CompositeObject::call_child_in(_uint side, const evec3& r) {
    if (child_types[side] == CGS_COMPOSITE) return ((CompositeObject*)children[side])->in(r);
    if (child_types[side] == CGS_SPHERE) return ((Sphere*)children[side])->in(r);
    if (child_types[side] == CGS_BOX) return ((Box*)children[side])->in(r);
    if (child_types[side] == CGS_PLANE) return ((Plane*)children[side])->in(r);
    if (child_types[side] == CGS_CYLINDER) return ((Cylinder*)children[side])->in(r);

    return invert;
}

void CompositeObject::add_child(_uint side, Object* o, object_type p_type) {
    children[side] = o;
    child_types[side] = p_type;
    //differences can be expressed as intersections of negations
    if (cmb == CGS_DIFFERENCE && side == 1) {
	children[side]->set_inversion(1);	
    }
}

int CompositeObject::in(const evec3& r) {
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
 * Convert a string separated by the character sep into a list of strings. For instance, if str == "arg1, arg2" and sep == ',' then the output will be a list of the form ["arg1", "arg2"]. If no sep character is found then the list will have one element which contains the whole string.
 * param str: string to parse into a list
 * param sep: separating character to use.
 * param listlen: location to save the length of the returned list to. Note that this pointer MUST be valid. It is not acceptable to call this function with listlen = NULL. (The returned list is not null terminated so this behavior ensures that callers are appropriately able to identify the length of the returned string.)
 * returns: list of strings separated by commas. This should be freed with a call to DTG_free(). In the event of an error, NULL is returned instead.
 * NOTE: The input string str is modified and the returned value uses memory allocated for str. Thus, the caller must ensure that str has a lifetime at least as long as the used string. Users should not try to call free() on any of the strings in the list, only the list itself.
 */
char** csv_to_list(char* str, char sep, size_t* listlen, parse_ercode& er) {
    er = E_SUCCESS;

    char** ret = NULL;/*(char**)malloc(sizeof(char*), &tmp_err);*/
    //we don't want to include separators that are in nested environments i.e. if the input is [[a,b],c] we should return "[a,b]","c" not "[a","b]","c"
    CGS_Stack<type_ind_pair> blk_stk;
    type_ind_pair tmp;

    //by default we ignore whitespace, only use it if we are in a block enclosed by quotes
    char* saveptr = str;
    size_t off = 0;
    size_t j = 0;
    size_t i = 0;
    bool verbatim = false;

    for (; str[i] != 0; ++i) {
	//if this is a separator then add the entry to the list
	if (str[i] == sep && blk_stk.is_empty()) {
	    ret = (char**)realloc(ret, sizeof(char*)*(off+1));

	    //append the element to the list
	    ret[off++] = saveptr;
	    //null terminate this string and increment j
	    str[j++] = 0;
	    saveptr = str + j;
	} else {
	    if (str[i] == '\\') {
		//check for escape sequences
		++i;
		bool placed_escape = false;
		switch (str[i]) {
		case 'n': str[j++] = '\n';placed_escape = true;break;
		case 't': str[j++] = '\t';placed_escape = true;break;
		case '\\': str[j++] = '\\';placed_escape = true;break;
		case '\"': str[j++] = '\"';placed_escape = true;break;
		default: er = E_BAD_SYNTAX;
		}
		if (placed_escape) continue;
	    } else if (str[i] == '\"') {
		tmp = blk_stk.peek();
		if (blk_stk.is_empty() || tmp.t != BLK_QUOTE) {
		    blk_stk.push(type_ind_pair(BLK_QUOTE, i));
		    verbatim = true;
		} else {
		    blk_stk.pop(&tmp);
		    verbatim = false;
		}
	    } else if (str[i] == '['/*]*/) {
		blk_stk.push(type_ind_pair(BLK_SQUARE, i));
	    } else if (str[i] == /*[*/']') {
		//don't fail if we reach the end of a block. This just means we've reached the end of the list
		if (blk_stk.pop(&tmp) != E_SUCCESS) break;
		if (tmp.t != BLK_SQUARE) { free(ret);er=E_BAD_SYNTAX;return NULL; }
	    } else if (str[i] == '('/*)*/) {
		blk_stk.push(type_ind_pair(BLK_PAREN, i));
	    } else if (str[i] == /*(*/')') {
		if (blk_stk.pop(&tmp) != E_SUCCESS) break;
		if (tmp.t != BLK_PAREN) { free(ret);er=E_BAD_SYNTAX;return NULL; }
	    } else if (str[i] == '{'/*}*/) {
		blk_stk.push(type_ind_pair(BLK_CURLY, i));
	    } else if (str[i] == /*{*/'}') {
		if (blk_stk.pop(&tmp) != E_SUCCESS) break;
		if (tmp.t != BLK_CURLY) { free(ret);er=E_BAD_SYNTAX;return NULL; }
	    }
	    if (!blk_stk.is_empty() || verbatim || (str[i] != ' ' && str[i] != '\t' && str[i] != '\n')) {
		//if this isn't whitespace then just copy the character
		str[j++] = str[i];
	    }
	}
    }
    //make sure that there weren't any unterminated blocks
    if (!blk_stk.is_empty()) {
	free(ret);
	er = E_BAD_SYNTAX;
	return NULL;
    }

    //make sure the string is null terminated
    //if (str[i] != 0) str[j] = 0;
    str[j] = 0;
    ret = (char**)realloc(ret, sizeof(char*)*(off+2));
    //add the last element to the list, but only if something was actually written, then set the length if requested
    if (j != 0) ret[off++] = saveptr;
    if (listlen) *listlen = off;
    ret[off] = NULL;
    return ret;
}

/**
 * Read a string of the format [x, y, z] into an Eigen::Vector3.
 * returns: 0 on success or an error code
 * 	-1: insufficient tokens
 * 	-2: one of the tokens supplied was invalid
 */
Value context::parse_list(char* str, parse_ercode& er) {
    Value sto;
    //find the start and the end of the vector
    char* start = strchr(str, '[');
    //make sure that the string is null terminated
    char* end = strchr_block(start+1, ']');
    if (!end) { er = E_BAD_SYNTAX;return sto; }
    *end = 0;

    //read the coordinates separated by spaces
    size_t n_els = 0;
    Value* lbuf;
    //check if this is a list interpretation
    char* for_start = token_block(start+1, KEY_FOR);
    size_t expr_len = 0;
    if (for_start) {
	//now look for a block labeled "in"
	char* in_start = token_block(for_start+KEY_FOR_LEN, "in");
	if (!in_start) { er = E_BAD_SYNTAX;return sto; }
	//the expression is the content before the for statement
	expr_len = for_start-(start+1);
	//the variable name is whatever is in between the "for" and the "in"
	in_start[0] = 0;
	char* var_name = CGS_trim_whitespace(for_start+KEY_FOR_LEN, &n_els);
	var_name = strndup(var_name, n_els+1);
	for_start[KEY_FOR_LEN+1+n_els] = ' ';
	//now parse the list we iterate over
	char* list_expr = strdup(in_start+KEY_IN_LEN);
	Value it_list = parse_value(list_expr, er);
	free(list_expr);
	if (er != E_SUCCESS || it_list.type != VAL_LIST) { er = E_BAD_VALUE;return sto; }
	//we now iterate through the list specified, substituting VAL in the expression with the current value
	for_start[0] = 0;
	emplace(var_name, sto);
	n_els = it_list.n_els;
	lbuf = (Value*)calloc(n_els, sizeof(Value));
	for (size_t i = 0; i < it_list.n_els; ++i) {
	    set_value(var_name, it_list.val.l[i]);
	    char* expr_name = strndup(start+1, expr_len);
	    lbuf[i] = parse_value(expr_name, er);
	    free(expr_name);
	    if (er != E_SUCCESS) { n_els = i;return sto; }
	}
	pop(NULL);
	//reset the string so that it can be parsed again
	for_start[0] = 'f';
	in_start[0] = 'i';
	//free the memory from the iteration list
	cleanup_val(&it_list);
	//free(expr_name);
	free(var_name);
    } else {
	char** list_els = csv_to_list(start+1, ',', &n_els, er);
	if (er != E_SUCCESS) { free(list_els);return sto; }
	lbuf = (Value*)calloc(n_els, sizeof(Value));
	for (size_t i = 0; list_els[i] && i < n_els; ++i) {
	    lbuf[i] = parse_value(list_els[i], er);
	    if (er != E_SUCCESS) { free(list_els);free(lbuf);return sto; }
	}
	free(list_els);
    }
    //cleanup and reset the string
    *end = ']';
    //set number of elements and type
    sto.type = VAL_LIST;
    sto.n_els = n_els;
    sto.val.l = lbuf;
    er = E_SUCCESS;
    return sto;
}

/**
 * Based on the declaration syntax produce the appropriate geometric shape and store the result in ptr. Ptr must not be initialized before a call to this function to avoid a memory leak.
 * returns: 0 on success or an error code
 */
parse_ercode Scene::make_object(const cgs_func& f, Object** ptr, object_type* type, int p_invert) const {
    *ptr = NULL;
    parse_ercode er = E_SUCCESS;

    if (!f.name) return E_BAD_TOKEN;
    //switch between all potential types
    if (strcmp(f.name, "Box") == 0) {
	if (f.n_args < 2) return E_LACK_TOKENS;
	//if we have enough tokens make sure we have both elements as vectors
	Value corn_1 = f.args[0].cast_to(VAL_3VEC, er);
	if (er != E_SUCCESS) return er;
	Value corn_2 = f.args[1].cast_to(VAL_3VEC, er);
	if (er != E_SUCCESS) return er;
	//make the box
	if (ptr) *ptr = new Box(*(corn_1.val.v), *(corn_2.val.v), p_invert);
	if (type) *type = CGS_BOX;
	//cleanup
	cleanup_val(&corn_1);
	cleanup_val(&corn_2);
    } else if (strcmp(f.name, "Plane") == 0) {
	if (f.n_args < 2) return E_LACK_TOKENS;
	if (f.n_args < 3) {
	    //this means we have a normal vector and an offset
	    Value corn_1 = f.args[0].cast_to(VAL_3VEC, er);
	    if (er != E_SUCCESS || f.args[1].type != VAL_NUM) return E_BAD_VALUE;
	    //make the plane
	    if (ptr) *ptr = new Plane(*(corn_1.val.v), f.args[1].val.x);
	    if (type) *type = CGS_PLANE;
	    cleanup_val(&corn_1);
	} else {
	    //this means we have three points defining the plane
	    Value corn_1 = f.args[0].cast_to(VAL_3VEC, er);
	    if (er != E_SUCCESS) return er;
	    Value corn_2 = f.args[1].cast_to(VAL_3VEC, er);
	    if (er != E_SUCCESS) return er;
	    Value corn_3 = f.args[2].cast_to(VAL_3VEC, er);
	    if (er != E_SUCCESS) return er;
	    //make the plane
	    if (ptr) *ptr = new Plane(*(corn_1.val.v), *(corn_2.val.v), *(corn_3.val.v));
	    if (type) *type = CGS_PLANE;
	    cleanup_val(&corn_1);
	    cleanup_val(&corn_2);
	    cleanup_val(&corn_3);
	}
    } else if (strcmp(f.name, "Sphere") == 0) {
	if (f.n_args < 2) return E_LACK_TOKENS;
	//if we have enough tokens make sure we have both elements as vectors
	Value cent = f.args[0].cast_to(VAL_3VEC, er);
	if (er != E_SUCCESS || f.args[1].type != VAL_NUM) return E_BAD_VALUE;
	double rad = f.args[1].val.x;
	if (ptr) *ptr = new Sphere(*(cent.val.v), rad, p_invert);
	if (type) *type = CGS_SPHERE;
	//cleanup
	cleanup_val(&cent);
    } else if (strcmp(f.name, "Cylinder") == 0) {
	if (f.n_args < 3) return E_LACK_TOKENS;
	Value cent = f.args[0].cast_to(VAL_3VEC, er);
	if (er != E_SUCCESS || f.args[1].type != VAL_NUM || f.args[2].type != VAL_NUM) return E_BAD_VALUE;
	double h = f.args[1].val.x;
	double r1 = f.args[2].val.x;
	//by default assume that the radii are the same
	double r2 = r1;
	if (f.n_args > 3) {
	    if (f.args[3].type != VAL_NUM) return E_BAD_VALUE;
	    r2 = f.args[3].val.x;
	}
	if (ptr) *ptr = new Cylinder(*(cent.val.v), h, r1, r2, p_invert);
	if (type) *type = CGS_CYLINDER;
	//cleanup
	cleanup_val(&cent);
    } else if (strcmp(f.name, "Composite") == 0) {
	if (ptr) *ptr = new CompositeObject(CGS_UNION, f, p_invert);
	if (type) *type = CGS_ROOT;
    } else if (strcmp(f.name, "union") == 0) {
	if (ptr) *ptr = new CompositeObject(CGS_UNION, f, p_invert);
	if (type) *type = CGS_COMPOSITE;
    } else if (strcmp(f.name, "intersect") == 0) {
	if (ptr) *ptr = new CompositeObject(CGS_INTERSECT, f, p_invert);
	if (type) *type = CGS_COMPOSITE;
    } else if (strcmp(f.name, "difference") == 0) {
	if (ptr) *ptr = new CompositeObject(CGS_DIFFERENCE, f, p_invert);
	if (type) *type = CGS_COMPOSITE;
    } else if (strcmp(f.name, "data") == 0) {
	if (ptr) *ptr = new CompositeObject(CGS_CMB_NOOP, f, p_invert);
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
parse_ercode Scene::make_transformation(const cgs_func& f, emat3& res) const {
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

/**
 * Given the string starting at token, and the index of an open paren parse the result into a cgs_func struct.
 * token: a c-string which is modified in place that contains the function
 * open_par_ind: the location of the open parenthesis
 * f: the cgs_func that information should be saved to
 * end: If not NULL, a pointer to the first character after the end of the string is stored here. If an error occurred during parsing end will be set to NULL.
 * returns: an errorcode if an invalid string was supplied.
 */
cgs_func context::parse_func(char* token, long open_par_ind, parse_ercode& er, char** end) {
    cgs_func f;

    //by default we want to indicate that we didn't get to the end
    if (end) *end = NULL;
    f.n_args = 0;
    //infer the location of the open paren index if the user didn't specify it
    if (open_par_ind < 0 || token[open_par_ind] != '('/*)*/) {
	char* par_char = strchr(token, '('/*)*/);
	//make sure there actually is an open paren
	if (par_char == NULL) { er = E_BAD_TOKEN;return f; }
	open_par_ind = par_char - token;
    }

    //break the string up at the parenthesis and remove surrounding whitespace
    token[open_par_ind] = 0;
    f.name = CGS_trim_whitespace(token, NULL);

    //now remove whitespace from the ends of the string
    char* arg_str = token+open_par_ind+1;
    //make sure that the string is null terminated
    char* term_ptr = strchr_block(arg_str, /*(*/')');
    if (!term_ptr) { er = E_BAD_SYNTAX;return f; }
    *term_ptr = 0;
    if (end) *end = term_ptr+1;
 
    //read the coordinates separated by spaces
    char** list_els = csv_to_list(arg_str, ',', &(f.n_args), er);
    if (er != E_SUCCESS) { free(list_els);return f; }
    //make sure that we don't go out of bounds, TODO: let functions accept arbitrarily many arguments?
    if (f.n_args >= ARGS_BUF_SIZE) {
	er = E_NOMEM;
	free(list_els);
	return f;
    }
    for (size_t i = 0; list_els[i] && i < f.n_args; ++i) {
	//handle named arguments
	char* eq_loc = strchr(list_els[i], '=');
	if (eq_loc) {
	    f.arg_names[i] = list_els[i];
	    *eq_loc = 0;
	    list_els[i] = eq_loc+1;
	} else {
	    f.arg_names[i] = NULL;
	}
	f.args[i] = parse_value(list_els[i], er);
	if (er != E_SUCCESS) { free(list_els);return f; }
    }

    //cleanup and reset the string
    //*term_ptr = /*(*/')';
    free(list_els);
    return f;
}

template <typename T>
parse_ercode CGS_Stack<T>::grow(size_t new_size) {
    size_t old_size = buf_len;
    buf_len = new_size;

    T* tmp_buf = (T*)realloc(buf, sizeof(T)*buf_len);
    if (tmp_buf) {
	buf = tmp_buf;
    } else {
	buf_len = old_size;
	return E_NOMEM;
    }

    return E_SUCCESS;
}

template <typename T>
CGS_Stack<T>::CGS_Stack() {
    stack_ptr = 0;
    buf_len = ARGS_BUF_SIZE;
    buf = (T*)malloc(sizeof(T)*buf_len);
    //check that allocation was successful
    if (!buf) {
	buf = NULL;
	buf_len = 0;
	stack_ptr = 0;
    }
}

template <typename T>
CGS_Stack<T>::~CGS_Stack<T>() {
    stack_ptr = 0;
    buf_len = 0;

    if (buf) {
	free(buf);
	buf = NULL;
    }
}

//swap
template <typename T>
void CGS_Stack<T>::swap(CGS_Stack<T>& o) {
    size_t tmp = buf_len;
    buf_len = o.buf_len;
    o.buf_len = tmp;
    tmp = stack_ptr;
    stack_ptr = o.stack_ptr;
    o.stack_ptr = tmp;
    T* tmp_buf = buf;
    buf = o.buf;
    o.buf = tmp_buf;
}

//copy constructor, assumes that the '=' operator does something sane
template <typename T>
CGS_Stack<T>::CGS_Stack(const CGS_Stack<T>& o) {
    stack_ptr = o.stack_ptr;
    //to save memory we'll only allocate however many entries the old object had
    buf_len = stack_ptr;
    buf = (T*)malloc(sizeof(T)*buf_len);
    if (!buf) {
	buf = NULL;
	buf_len = 0;
	stack_ptr = 0;
    }

    //copy the entries from the old object stack
    for (_uint i = 0; i < stack_ptr; ++i) buf[i] = o.buf[i];
}

//move constructor
template <typename T>
CGS_Stack<T>::CGS_Stack(CGS_Stack<T>&& o) {
    stack_ptr = o.stack_ptr;
    buf_len = o.buf_len;
    buf = o.buf;
    //invalidate the object that just got moved
    o.stack_ptr = 0;
    o.buf_len = 0;
    o.buf = NULL;
}

//assignment operator
template <typename T>
CGS_Stack<T>& CGS_Stack<T>::operator=(CGS_Stack<T> o) {
    swap(o);
    return *this;
}

/**
 * Push a side-object pair onto the stack.
 * returns: An error code if one occurred or E_SUCCESS if the operation was successful. It is possible for this function to fail if there wasn't sufficient space to push the object onto the stack.
 */
template <typename T>
parse_ercode CGS_Stack<T>::push(T val) {
    //we might need to allocate more memory
    parse_ercode res = E_SUCCESS;
    if (stack_ptr == buf_len) res = grow(2*buf_len);
    if (res == E_SUCCESS) buf[stack_ptr++] = val;

    return res;
}

/**
 * Pop a side-object pair into the values pointed to by side and obj respectively. This function is guaranteed to succeed, but if the stack is empty then SIDE_END will be assigned to side and NULL will be assigned to obj. The caller should check that the returned values are valid.
 */
template <typename T>
parse_ercode CGS_Stack<T>::pop(T* ptr) {
    //check if there is anything that can be popped
    if (stack_ptr == 0) return E_EMPTY_STACK;
    //assign the return values
    if (ptr) *ptr = buf[stack_ptr-1];
    //decrement the stack pointer
    --stack_ptr;

    return E_SUCCESS;
}

/**
 * Reset the stack to be empty
 */
template <typename T>
void CGS_Stack<T>::reset() { stack_ptr = 0; }

/**
 * Check wheter an entry matching key is already somewhere in the stack
 */
template <typename T>
bool CGS_Stack<T>::has(const T& key) {
    for (_uint i = 0; i < stack_ptr; ++i) {
	if (buf[i] == key) return true;
    }
    return false;
}

/**
 * Returns a copy of the object at the top of the stack. If the stack is empty then the object casted from 0 is returned.
 */
template <typename T>
T CGS_Stack<T>::peek(size_t ind) {
    if (stack_ptr >= ind && ind > 0)
	return buf[stack_ptr - ind];
    return (T)0;
}

/**
 * Insert an object onto the stack obeying the desired tree structure.
 */
parse_ercode ObjectStack::emplace_obj(Object* obj, object_type p_type) {
    if (stack_ptr == 0) {
	if (p_type == CGS_COMPOSITE || p_type == CGS_ROOT || p_type == CGS_DATA) {
	    side_obj_pair cur(0, (CompositeObject*)obj);
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
	    side_obj_pair cur(0, (CompositeObject*)obj);
	    return push(cur);
	}
    }

    return E_SUCCESS;
}

_uint ObjectStack::look_side() {
    if (stack_ptr == 0) return SIDE_UNDEF;
    return buf[stack_ptr - 1].side;
}

CompositeObject* ObjectStack::look_obj() {
    if (stack_ptr == 0) return NULL;
    return buf[stack_ptr].obj;
}

CompositeObject* ObjectStack::get_root() {
    if (stack_ptr == 0) return NULL;
    return buf[0].obj;
}

parse_ercode Scene::fail_exit(parse_ercode er, FILE* fp) {
    fclose(fp);
    return er;
}

/**
 * This acts similar to getline, but stops at a semicolon, newline (unless preceeded by a \), {, or }.
 * bufptr: a pointer to which the buffer is saved. If bufptr is NULL than a new buffer is allocated through malloc()
 * n: a pointer to a size_t with the number of characters in the buffer pointed to by bufptr. The call will return do nothing and return -1 if n is null but *bufptr is not.
 * fp: file pointer to read from
 * linecount: a pointer to an integer specifying the number of new line characters read.
 * Returns: -2 if an error occured, -1 if the end of the file was reached, or the number of characters read (including null termination) if successful
 */
int read_cgs_line(char** bufptr, size_t* n, FILE* fp, size_t* lineno) {
    //dereference pointers and interpret them correctly
    size_t n_lines = 0;
    if (lineno) n_lines = *lineno;
    size_t size = *n;
    char* buf = *bufptr;
    if (buf) {
	if (size == 0) return -2;
    } else {
	size = BUF_SIZE;
	buf = (char*)malloc(sizeof(char)*size);
    }

    int res = fgetc(fp);
    size_t i = 0;
    for (;; ++i) {
	if (i >= size) {
	    size *= 2;
	    buf = (char*)realloc(buf, sizeof(char)*size);
	}
	if (res == EOF) {
	    buf[i] = 0;
	    *n = size;
	    *bufptr = buf;
	    if (lineno) *lineno = n_lines;
	    return -1;
	} else if (res == ';') {
	    buf[i] = 0;
	    //only count one line end
	    res = fgetc(fp);
	    if (res == '\n')
		++n_lines;
	    else
		fseek(fp, -1, SEEK_CUR);
	    break;
	} else if (res == '{' || res == '}') {
	    buf[i] = (char)res;
	    res = (int)';';
	} else if (res == '\n') {
	    if (i == 0 || buf[i-1] != '\\') {
		long cur_pos = ftell(fp);
		int inc = 2;
		//we give this case it's own special logic because we want to roll the next open brace or semicolon onto this line
		while (res != EOF) {
		    //read until we hit an open curly brace or a semicolon or a non-whitespace character. If the character is non-whitespace return to the inital state. Otherwise write the semicolon/curly brace at the end of the current line
		    res = fgetc(fp);
		    if (res == '{'/*}*/) {
			buf[i] = (char)res;
			res = (int)';';
			break;
		    } else if (res == '\n') {
			n_lines += inc;//increment by 2 for the first newline we find
			inc = 1;
		    } else if (res != ' ' && res != '\t') {
			fseek(fp, cur_pos-1, SEEK_SET);//minus 1 to ensure the newline gets read when parsing the semicolon
			res = ';';
			--i;
			break;
		    }
		}
	    } else {
		i -= 2;
		res = fgetc(fp);
	    }
	} else {
	    buf[i] = (char)res;
	    res = fgetc(fp);
	}
    }
    *n = size;
    *bufptr = buf;
    if (lineno) *lineno = n_lines;
    return (int)(i+1);
}

/**
 * Iterate through the context and return a value to the variable with the matching name.
 * name: the name of the variable to set
 * returns: the matching value, no deep copies are performed
 */
Value context::lookup(const char* str) const {
    //iterate to find the item highest on the stack with a matching name
    for (long i = stack_ptr-1; i >= 0; --i) {
	if (strcmp(buf[i].name, str) == 0) return buf[i].val;
    }
    //return a default undefined value if nothing was found
    Value ret;
    ret.type = VAL_UNDEF;
    ret.val.x = 0;
    return ret;
}

/**
 * Iterate through the context and assign a new value to the variable with the matching name.
 * name: the name of the variable to set
 * new_val: the value to set the variable to
 * returns: E_SUCCESS if the variable with a matching name was found or E_NOT_DEFINED otherwise
 */
parse_ercode context::set_value(const char* name, Value new_val) {
    for (long i = stack_ptr-1; i >= 0; --i) {
	if (strcmp(buf[i].name, name) == 0) {
	    buf[i].val = new_val;
	    return E_SUCCESS;
	}
    }
    return E_NOT_DEFINED;
}

/**
 * remove the top n items from the stack
 */
parse_ercode context::pop_n(size_t n) {
    if (n > stack_ptr)
	return E_EMPTY_STACK;
    stack_ptr -= n;
    return E_SUCCESS;
}
/**
 * Make a vector argument with the x,y, and z coordinates supplied
 */
Value make_vec(cgs_func tmp_f, parse_ercode& er) {
    Value sto;sto.type = VAL_UNDEF;
    if (tmp_f.n_args < 3) { er = E_LACK_TOKENS;return sto; }
    if (tmp_f.args[0].type != VAL_NUM || tmp_f.args[1].type != VAL_NUM || tmp_f.args[2].type != VAL_NUM) { er = E_BAD_TOKEN;return sto; }
    sto.type = VAL_3VEC;
    sto.val.v = new evec3(tmp_f.args[0].val.x, tmp_f.args[1].val.x, tmp_f.args[2].val.x);
    sto.n_els = 3;
    er = E_SUCCESS;
    return sto;
}
/**
 * Make a range following python syntax. If one argument is supplied then a list with tmp_f.args[0] elements is created starting at index 0 and going up to (but not including) tmp_f.args[0]. If two arguments are supplied then the range is from (tmp_f.args[0], tmp_f.args[1]). If three arguments are supplied then the range (tmp_f.args[0], tmp_f.args[1]) is still returned, but now the spacing between successive elements is tmp_f.args[2].
 */
Value make_range(cgs_func tmp_f, parse_ercode& er) {
    Value sto;sto.type = VAL_UNDEF;
    if (tmp_f.n_args == 0) {
	er = E_LACK_TOKENS;
	return sto;
    } else if (tmp_f.n_args == 1) {
	if (tmp_f.args[0].type != VAL_NUM) { er = E_BAD_VALUE;return sto; }
	if (tmp_f.args[0].val.x < 0) { er = E_BAD_VALUE;return sto; }
	//interpret the argument as an upper bound starting from 0
	sto.type = VAL_LIST;
	sto.n_els = (size_t)(tmp_f.args[0].val.x);
	sto.val.l = (Value*)malloc(sizeof(Value)*sto.n_els);
	for (size_t i = 0; i < sto.n_els; ++i) {
	    sto.val.l[i].type = VAL_NUM;
	    sto.val.l[i].val.x = i;
	}
    } else if (tmp_f.n_args == 2) {
	if (tmp_f.args[0].type != VAL_NUM) { er = E_BAD_VALUE;return sto; }
	if (tmp_f.args[1].type != VAL_NUM) { er = E_BAD_VALUE;return sto; }
	//interpret the argument as an upper bound starting from 0
	sto.type = VAL_LIST;
	int list_start = (int)(tmp_f.args[0].val.x);
	int list_end = (int)(tmp_f.args[1].val.x);
	if (list_end < list_start) { er = E_BAD_VALUE;return sto; }
	sto.n_els = list_end - list_start;
	sto.val.l = (Value*)malloc(sizeof(Value)*sto.n_els);
	size_t j = 0;
	for (int i = list_start; i < list_end; ++i) {
	    sto.val.l[j].type = VAL_NUM;
	    sto.val.l[j].val.x = i;
	    ++j;
	}
    } else {
	if (tmp_f.args[0].type != VAL_NUM) { er = E_BAD_VALUE;return sto; }
	if (tmp_f.args[1].type != VAL_NUM) { er = E_BAD_VALUE;return sto; }
	if (tmp_f.args[2].type != VAL_NUM) { er = E_BAD_VALUE;return sto; }
	//interpret the argument as an upper bound starting from 0
	sto.type = VAL_LIST;
	double list_start = tmp_f.args[0].val.x;
	double list_end = tmp_f.args[1].val.x;
	double inc = tmp_f.args[2].val.x;
	if (list_end < list_start || inc == 0) { er = E_BAD_VALUE;return sto; }
	sto.n_els = (list_end - list_start) / inc;
	sto.val.l = (Value*)malloc(sizeof(Value)*sto.n_els);
	for (size_t i = 0; i < sto.n_els; ++i) {
	    sto.val.l[i].type = VAL_NUM;
	    sto.val.l[i].val.x = i*inc + list_start;
	}
    }
    er = E_SUCCESS;
    return sto;
}
/**
 * Take a list value and flatten it so that it has numpy dimensions (n) where n is the sum of the length of each list in the base list. Values are copied in order e.g flatten([0,1],[2,3]) -> [0,1,2,3]
 * cgs_func: the function with arguments passed
 */
Value flatten_list(cgs_func tmp_f, parse_ercode& er) {
    Value sto;sto.type = VAL_UNDEF;
    if (tmp_f.n_args < 1) { er = E_LACK_TOKENS;return sto; }
    if (tmp_f.args[0].type != VAL_LIST) { er = E_BAD_VALUE;return sto; }
    Value cur_list = tmp_f.args[0];
    size_t cur_st = 0;
    //this is used for estimating the size of the buffer we need. Take however many elements were needed for this list and assume each sub-list has the same number of elements
    size_t base_n_els = cur_list.n_els;
    //start with the number of elements in the lowest order of the list
    size_t buf_size = cur_list.n_els;
    sto.val.l = (Value*)malloc(sizeof(Value)*buf_size);
    //there may potentially be nested lists, we need to be able to find our way back to the parent and the index once we're done
    CGS_Stack<Value> lists;
    CGS_Stack<size_t> inds;
    //just make sure that there's a root level on the stack to be popped out
    lists.push(cur_list);
    inds.push(0);
    size_t start_depth = inds.size();
    size_t j = 0;
    do {
	size_t i = cur_st;
	start_depth = inds.size();
	for (; i < cur_list.n_els; ++i) {
	    if (cur_list.val.l[i].type == VAL_LIST) {
		lists.push(cur_list);
		inds.push(i+1);//push + 1 so that we start at the next index instead of reading the list again
		cur_list = cur_list.val.l[i];
		cur_st = 0;
		break;
	    }
	    if (j >= buf_size) {
		//-1 since we already have at least one element. no base_n_els=0 check is needed since that case will ensure the for loop is never evaluated
		buf_size += (base_n_els-1)*(i+1);
		Value* tmp_val = (Value*)realloc(sto.val.l, sizeof(Value)*buf_size);
		if (!tmp_val) {
		    free(sto.val.l);
		    cleanup_func(&tmp_f);
		    sto.type = VAL_UNDEF;
		    er = E_NOMEM;
		    return sto;
		}
		sto.val.l = tmp_val;
	    }
	    sto.val.l[j++] = copy_val(cur_list.val.l[i]);
	}
	//if we reached the end of a list without any sublists then we should return back to the parent list
	if (inds.size() <= start_depth) {
	    inds.pop(&cur_st);
	    lists.pop(&cur_list);
	}
    } while (!lists.is_empty());
    sto.type = VAL_LIST;
    sto.n_els = j;
    er = E_SUCCESS;
    return sto;
}
/**
 * Execute the mathematical operation in the string str at the location op_ind
 */
Value context::do_op(char* str, size_t i, parse_ercode& er) {
    Value sto;
    sto.type = VAL_UNDEF;
    sto.val.x = 0;
    sto.n_els = 0;
    //Store the operation before setting it to zero
    char term_char = str[i];
    str[i] = 0;
    //parse right and left values
    Value tmp_l = parse_value(str, er);
    if (er != E_SUCCESS) return sto;
    Value tmp_r = parse_value(str+i+1, er);
    if (er != E_SUCCESS) return sto;
    sto.type = VAL_NUM;

    //handle equality comparisons
    if (term_char == '=') {
	if (tmp_l.type != tmp_r.type) {
	    sto.val.x = 0;
	} else if (tmp_l.type == VAL_STR) {
	    if (tmp_l.size() == tmp_r.size() && strcmp(tmp_l.get_val().s, tmp_l.get_val().s) == 0)
		sto.val.x = 1;
	    else
		sto.val.x = 0;
	} else if (tmp_l.type == VAL_NUM) {
	    sto.val.x = 1 - (tmp_l.val.x - tmp_r.val.x);
	}
	++i;
    } else if (term_char == '>') {
	if (tmp_l.type != VAL_NUM || tmp_r.type != VAL_NUM) { er = E_BAD_VALUE;return sto; }
	if (str[i+1] == '=') {
	    sto.val.x = (tmp_l.val.x >= tmp_r.val.x)? 1: 0;
	    ++i;
	} else {
	    sto.val.x = (tmp_l.val.x > tmp_r.val.x)? 1: 0;
	}
    } else if (term_char == '<') {
	if (tmp_l.type != VAL_NUM || tmp_r.type != VAL_NUM) { er = E_BAD_VALUE;return sto; }
	if (str[i+1] == '=') {
	    sto.val.x = (tmp_l.val.x <= tmp_r.val.x)? 1: 0;
	    ++i;
	} else {
	    sto.val.x = (tmp_l.val.x < tmp_r.val.x)? 1: 0;
	}
    } else if (term_char == '+' && (i == 0 || str[i-1] != 'e')) {
	if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_NUM) {	
	    sto.val.x = tmp_l.val.x + tmp_r.val.x;
	} else if (tmp_l.type == VAL_MAT && tmp_r.type == VAL_MAT) {	
	    sto.type = VAL_MAT;
	    sto.val.m = new Eigen::MatrixXd(*tmp_l.val.m + *tmp_r.val.m);
	} else if (tmp_l.type == VAL_STR) {
	    //TODO: implement string concatenation
	}
    } else if (term_char == '-') {
	if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_NUM) {
	    sto.val.x = tmp_l.val.x - tmp_r.val.x;
	} else if (tmp_l.type == VAL_MAT && tmp_r.type == VAL_MAT) {	
	    sto.type = VAL_MAT;
	    sto.val.m = new Eigen::MatrixXd(*tmp_l.val.m - *tmp_r.val.m);
	}
    } else if (term_char == '*') {
	if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_NUM) {	
	    sto.val.x = tmp_l.val.x * tmp_r.val.x;
	} else if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_MAT) {	
	    sto.type = VAL_MAT;
	    sto.val.m = new Eigen::MatrixXd(tmp_l.val.x * *tmp_r.val.m);
	} else if (tmp_l.type == VAL_MAT && tmp_r.type == VAL_MAT) {	
	    sto.type = VAL_MAT;
	    sto.val.m = new Eigen::MatrixXd(*tmp_l.val.m * *tmp_r.val.m);
	}
    } else if (term_char == '/') {
	if (tmp_r.val.x == 0) { er = E_NAN;return sto; }//TODO: return a nan?
	if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_NUM) {
	    sto.val.x = tmp_l.val.x / tmp_r.val.x;
	} else if (tmp_r.type == VAL_NUM && tmp_l.type == VAL_MAT) {	
	    sto.type = VAL_MAT;
	    sto.val.m = new Eigen::MatrixXd(*tmp_r.val.m / tmp_r.val.x);
	}
    } else if (term_char == '^') {
	if (tmp_r.val.x == 0 || tmp_l.type != VAL_NUM || tmp_r.type != VAL_NUM) { er = E_NAN;return sto; }//TODO: return a nan?
	sto.type = VAL_NUM;
	sto.val.x = pow(tmp_l.val.x, tmp_r.val.x);
    }
    str[i] = term_char;
    return sto;
}

Value context::parse_value(char* str, parse_ercode& er) {
    Value sto;
    er = E_SUCCESS;

    //store locations of the first instance of different operators. We do this so we can quickly look up new operators if we didn't find any other operators of a lower precedence (such operators are placed in the tree first).
    int first_open_ind = -1;
    int last_close_ind = -1;

    //keeps track of open and close [], (), {}, and ""
    CGS_Stack<type_ind_pair> blk_stk;
    type_ind_pair start_ind;

    //keep track of the nesting level within parenthetical statements
    int nest_level = 0;

    //keep track of the precedence of the orders of operation (lower means executed later) ">,=,>=,==,<=,<"=4 "+,-"=3, "*,/"=2, "**"=1
    char op_prec = 0;
    size_t op_loc = 0;
    size_t reset_ind = 0;
    for (_uint i = 0; str[i] != 0; ++i) {
	if (str[i] == '['/*]*/) {
	    blk_stk.push({BLK_SQUARE, i});
	    ++nest_level;
	    if (first_open_ind == -1) { first_open_ind = i; }
	} else if (str[i] == /*[*/']') {
	    if (blk_stk.pop(&start_ind) != E_SUCCESS || start_ind.t != BLK_SQUARE) { er = E_BAD_SYNTAX;return sto; }
	    --nest_level;
	    //if (blk_stk.is_empty() && last_close_ind < 0) last_close_ind = i;
	    last_close_ind = i;
	} else if (str[i] == '('/*)*/) {
	    //keep track of open and close parenthesis, these will come in handy later
	    blk_stk.push({BLK_PAREN, i});
	    //TODO: decide if this is necessary
	    ++nest_level;
	    //only set the open index if this is the first match
	    if (first_open_ind == -1) { first_open_ind = i; }
	} else if (str[i] == /*(*/')') {
	    if (blk_stk.pop(&start_ind) != E_SUCCESS || start_ind.t != BLK_PAREN) { er = E_BAD_SYNTAX;return sto; }
	    --nest_level;
	    //only set the end paren location if it hasn't been set yet and the stack has no more parenthesis to remove, TODO: make this work with other block types inside a set of parenthesis
	    //if (blk_stk.is_empty() && last_close_ind < 0) last_close_ind = i;
	    last_close_ind = i;
	} else if (str[i] == '\"' && (i == 0 || str[i-1] != '\\')) {
	    start_ind = blk_stk.peek();
	    if (blk_stk.is_empty() || start_ind.t != BLK_QUOTE) {
		//quotes need to be handled in a special way
		blk_stk.push({BLK_QUOTE, i});
		//TODO: decide if this is necessary
		++nest_level;
		//only set the open index if this is the first match
		if (first_open_ind == -1) first_open_ind = i;
	    } else {
		blk_stk.pop(&start_ind);
		last_close_ind = i;
		--nest_level;
	    }
	}

	if (nest_level == 0) {
	    //keep track of the number of characters used by the operator
	    int code_n_chars = 1;
	    //check if we found a numeric operation symbol
	    if (((str[i] == '=' && str[i+1] == '=') || str[i] == '>' || str[i] == '<') && op_prec < 4) {
		op_prec = 4;
		op_loc = i;
	    } else if (i != 0 && (str[i] == '+' || str[i] == '-') && str[i-1] != 'e' && op_prec < 3) {
		//remember to recurse after we finish looping
		op_prec = 3;
		op_loc = i;
	    } else if (str[i] == '^' && op_prec < 2) {
		op_prec = 2;
		op_loc = i;
	    } else if ((str[i] == '*' || str[i] == '/') && op_prec < 1) {
		op_prec = 1;
		op_loc = i;
	    }
	}
    }
    if (nest_level > 0) {
	printf("expected close paren\n");
	er = E_BAD_SYNTAX;
	return sto;
    }

    //last try removing parenthesis
    if (op_prec <= 0) {
	//if there isn't a valid parenthetical expression, then we should interpret this as a variable
	if (first_open_ind < 0 || last_close_ind < 0) {
	    str = CGS_trim_whitespace(str, &reset_ind);
	    sto = lookup(str);
	    if (sto.type == VAL_UNDEF) {
		//try interpreting as a number
		errno = 0;
		sto.val.x = strtod(str, NULL);
		if (errno) {
		    er = E_NOT_DEFINED;
		    return sto;
		}
		er = E_SUCCESS;
		sto.type = VAL_NUM;
	    }
	    str[reset_ind] = ' ';//all whitespace is treated identically so it doesn't matter
	} else if (str[first_open_ind] == '\"' && str[last_close_ind] == '\"') {
	    //this is a string
	    sto.type = VAL_STR;
	    sto.n_els = last_close_ind-first_open_ind;
	    //allocate memory and copy
	    sto.val.s = (char*)malloc(sizeof(char)*sto.n_els);
	    for (size_t i = 0; i < sto.n_els-1; ++i) sto.val.s[i] = str[first_open_ind+i+1];
	    sto.val.s[sto.n_els-1] = 0;
	} else if (str[first_open_ind] == '[' && str[last_close_ind] == ']') {
	    //now parse the argument as a list
	    sto = parse_list(str+first_open_ind, er);
	    sto.type = VAL_LIST;
	    if (er != E_SUCCESS) return sto;
	} else if (str[first_open_ind] == '(' && str[last_close_ind] == ')') {
	    //check to see if this is a function call (it is if there are any non-whitespace characters before the open paren
	    int is_func = 0;
	    for (size_t i = 0; i < first_open_ind; ++i) {
		if (str[i] != ' ' && str[i] != '\t' && str[i] != '\n') {
		    //we can't leave this as zero in case the user needs to do some more operations
		    char term_char = str[last_close_ind+1];
		    str[last_close_ind+1] = 0;
		    cgs_func tmp_f;
		    char* f_end;
		    tmp_f = parse_func(str, first_open_ind, er, &f_end);
		    if (er != E_SUCCESS) { return sto; }
		    //TODO: allow user defined functions
		    er = E_BAD_TOKEN;
		    if (strcmp(tmp_f.name, "vec") == 0) {
			sto = make_vec(tmp_f, er);
		    } else if (strcmp(tmp_f.name, "range") == 0) {
			sto = make_range(tmp_f, er);
		    } else if (strcmp(tmp_f.name, "flatten") == 0) {
			sto = flatten_list(tmp_f, er);
		    }
		    str[last_close_ind+1] = term_char;
		    cleanup_func(&tmp_f);
		    return sto;
		}
	    }
	    //otherwise interpret this as a parenthetical expression
	    str[last_close_ind] = 0;
	    char* tmp_str = CGS_trim_whitespace(str+first_open_ind+1, &reset_ind);
	    sto = parse_value(str, er);
	    tmp_str[reset_ind] = ' ';
	    str[first_open_ind] = '(';
	    str[last_close_ind] = ')';
	}
    } else {
	sto = do_op(str, op_loc, er);
	if (er != E_SUCCESS) { sto.type=VAL_UNDEF;return sto; }
    }
    return sto;
}

user_func::user_func(cgs_func sig, char** bufptr, size_t* n, FILE* fp) {
    call_sig = copy_func(sig);
    //setup a buffer for code lines
    size_t cur_buf_size = FUNC_BUF_SIZE;
    code_lines = (char**)malloc(sizeof(char*)*cur_buf_size);
    //initialize tracking variables
    n_lines = 0;
    int nest_level = 0;
    int n_chars = read_cgs_line(bufptr, n, fp, &n_lines);
    //iterate until we hit an error, or the end of the function (n_chars=-1 or nest_level=-1)
    while (n_chars >= 0 && nest_level >= 0) {
	//expand the buffer if necessary
	if (n_lines == cur_buf_size) {
	    cur_buf_size *= 2;
	    code_lines = (char**)realloc(code_lines, sizeof(char*)*cur_buf_size);
	}
	char* cur_line = (char*)malloc(sizeof(char)*n_chars);
	char* buf = *bufptr;
	//copy the string and check to see if there was an open or close brace read
	for (size_t j = 0; j < n_chars; ++j) {
	    cur_line[j] = buf[j];
	    if (buf[j] == '{')
		++nest_level;
	    else if (buf[j] == '}')
		--nest_level;
	}
	cur_line[n_chars-1] = 0;
	code_lines[n_lines++] = CGS_trim_whitespace(cur_line, NULL);
	n_chars = read_cgs_line(bufptr, n, fp, &n_lines);
    }
}
user_func::~user_func() {
    cleanup_func(&call_sig);
    for (size_t j = 0; j < n_lines; ++j) free(code_lines[j]);
    free(code_lines);
}
user_func::user_func(const user_func& o) {
    call_sig = copy_func(o.call_sig);
    n_lines = o.n_lines;
    code_lines = (char**)malloc(sizeof(char*)*n_lines);
    for (size_t j = 0; j < n_lines; ++j) code_lines[j] = strdup(o.code_lines[j]);
}
user_func::user_func(user_func&& o) {
    call_sig = o.call_sig;
    n_lines = o.n_lines;
    code_lines = o.code_lines;
    o.call_sig.name = NULL;
    o.call_sig.n_args = 0;
    o.n_lines = 0;
    o.code_lines = NULL;
}
Value user_func::eval(context& c, cgs_func call, parse_ercode& er) {
    er = E_SUCCESS;
    std::unordered_map<std::string, Value> local_vars;
    size_t line = 0;
    Value bad_val;
    while (line < n_lines) {
	if (strncmp("return", code_lines[line], strlen("return")) == 0) {
	    char* tmp_buf = CGS_trim_whitespace(code_lines[line]+strlen("return")+1, NULL);
	    Value v = c.parse_value(tmp_buf, er);
	    if (er != E_SUCCESS) {
		std::string tmp(tmp_buf);
		if (local_vars.count(tmp) > 0)
		    return local_vars[tmp];
		else
		    return bad_val;
	    } else {
		return v;
	    }
	} else if (strncmp("if", code_lines[line], strlen("if")) == 0){
	    size_t len = 0;
	    char* tmp_buf = CGS_trim_whitespace(code_lines[line]+strlen("if")+1, &len);
	}
    }
    //return an undefined value by default
    Value ret;ret.type = VAL_UNDEF;ret.val.x = 0;return ret;
}

parse_ercode Scene::read_file(const char* p_fname) {
    parse_ercode er = E_SUCCESS;
    //char buf[BUF_SIZE];
    char* buf = (char*)malloc(sizeof(char)*BUF_SIZE);
    size_t buf_size = BUF_SIZE;
    char func_name[BUF_SIZE];
    size_t n_args = 0;

    //indicate that no names are stored in all of the specified fields
    func_name[0] = 0;
    //generate the composite object on the fly
    ObjectStack tree_pos;
    CGS_Stack<block_type> blk_stack;
    block_type last_type = BLK_UNDEF;
    CompositeObject* last_comp = NULL;
    int invert = 0;
    //keep track of transformations
    CGS_Stack<emat3> transform_stack;
    emat3 cur_trans_mat;
    emat3 next_trans_mat;

    //open the file for reading and read individual lines. We need to remove whitespace. Handily, we know the line length once we're done.
    FILE* fp = NULL;
    if (p_fname) {
        fp = fopen(p_fname, "r");
    }
    if (fp) {
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
		    if (buf[i] == '(' && blk_stack.peek() != BLK_LITERAL) {
			//initialize a new cgs_func with the appropriate arguments
			cgs_func cur_func;
			char* endptr;
			cur_func = named_items.parse_func(buf, i, er, &endptr);
			switch (er) {
			    case E_BAD_TOKEN:
				printf("Error on line %d: Invalid function name \"%s\"\n", lineno, buf);
				cleanup_func(&cur_func);
				return fail_exit(er, fp);
			    case E_BAD_SYNTAX:
				printf("Error on line %d: Invalid syntax\n", lineno);
				cleanup_func(&cur_func);
				return fail_exit(er, fp);
			    default: break;
			}
			//check to see if this is a user function declaration
			if (cur_func.name[0] == 'd' && cur_func.name[1] == 'e' && cur_func.name[2] == 'f' && cur_func.name[3] != 0) {
			    cur_func.name = CGS_trim_whitespace(cur_func.name + 3, NULL);

			    //TODO: finish
			}
			//try interpreting the function as a geometric object
			Object* obj = NULL;
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
				    last_comp = (CompositeObject*)obj;
				    last_type = BLK_ROOT;
				    //this is included so that we don't have to check whether something is a root or a composite every time
				    type = CGS_COMPOSITE;
				}
			    } else if (type == CGS_DATA) {
				data_objs.push_back((CompositeObject*)obj);
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
			char* tok = CGS_trim_whitespace(buf, &tok_len);
			size_t val_len;
			char* val = CGS_trim_whitespace(buf+i+1, &val_len);
			Value v = named_items.parse_value(val, er);
			if (er == E_SUCCESS) named_items.emplace(tok, v);
			tok[tok_len] = ' ';
			val[val_len] = ' ';
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
	fclose(fp);
    } else {
        printf("Error: couldn't open file %s for reading!\n", p_fname);
        return E_NOFILE;
    }
    free(buf);

    return er;
}

Scene::Scene(const char* p_fname, parse_ercode* ercode) {
    if (ercode) *ercode = E_SUCCESS;
    *ercode = read_file(p_fname);
}

Scene::Scene(const Scene& o) {
    roots.resize(o.roots.size());
    data_objs.resize(o.data_objs.size());
    for (_uint i = 0; i < roots.size(); ++i) {
	roots[i] = new CompositeObject( *(o.roots[i]) );
    }
    for (_uint i = 0; i < data_objs.size(); ++i) {
	data_objs[i] = new CompositeObject( *(o.data_objs[i]) );
    }
    named_items = o.named_items;
}

Scene::Scene(Scene&& o) {
    roots = o.roots;
    data_objs = o.data_objs;
    named_items = o.named_items;
    o.roots.clear();
    o.data_objs.clear();
    o.named_items.reset();
}

Scene& Scene::operator=(Scene& o) {
    roots = o.roots;
    data_objs = o.data_objs;
    named_items = o.named_items;

    return *this;
}

Scene::~Scene() {
    //TODO: double frees are bad lol
    for (_uint i = 0; i < roots.size(); ++i) {
	if (roots[i]) delete roots[i];
    }
    for (_uint i = 0; i < data_objs.size(); ++i) {
	if (data_objs[i]) delete data_objs[i];
    }
}
