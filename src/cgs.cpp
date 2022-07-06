#include "cgs.hpp"

Object::Object(int p_invert) {
    invert = p_invert;
    trans_mat = Eigen::Matrix3d::Identity();
}

Object::Object(Eigen::Quaterniond& p_orientation, int p_invert) {
    invert = p_invert;
    trans_mat = p_orientation.toRotationMatrix();
}

void Object::set_inversion(int p_invert) {
    if (p_invert == 0)
	invert = 0;
    else
	invert = 1;
}

Sphere::Sphere(evec3& p_center, double p_rad, int p_invert) : Object(p_invert) {
    center = p_center;
    rad = p_rad;
}

int Sphere::in(const evec3& r) {
    //find the relative offset between the input vector and the center of the box
    evec3 r_rel = r - center;
    if (r_rel.norm() > rad) return invert;
    return 1 - invert;
}

Box::Box(evec3& p_corner, evec3& p_offset, int p_invert) : Object(p_invert) {
    offset = p_offset/2;
    center = p_corner + offset;
    //wlog, fix the offset to always be positive from the center
    if (offset.x() < 0) offset.x() *= -1;
    if (offset.y() < 0) offset.y() *= -1;
    if (offset.z() < 0) offset.z() *= -1;
}

Box::Box(evec3& p_corner, evec3& p_offset, Eigen::Quaterniond p_orientation, int p_invert) : Object(p_orientation, p_invert) {
    offset = p_offset/2;
    center = p_corner + offset;
    //wlog, fix the offset to always be positive from the center
    if (offset.x() < 0) offset.x() *= -1;
    if (offset.y() < 0) offset.y() *= -1;
    if (offset.z() < 0) offset.z() *= -1;
}

int Box::in(const evec3& r) {
    //find the relative offset between the input vector and the center of the box
    evec3 r_rel = trans_mat*(r - center);
    //check if we're outside the bounding box
    if (abs(r_rel.x()) > offset.x()) return invert;
    if (abs(r_rel.y()) > offset.y()) return invert;
    if (abs(r_rel.z()) > offset.z()) return invert;
    //if we reach this point in execution we're inside the box
    return 1 - invert;
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
	//each piece of metadata should be separated by a named list
	char* save_str;
	char* tok = strtok_r(spec.args[i], "=", &save_str);
	if (tok) {
	    char* val = strtok_r(NULL, "=", &save_str);
	    //only insert a metadata entry if all information was valid
	    if (val) {
		std::string tok_cpp(tok);std::string val_cpp(val);
		metadata[tok_cpp] = val_cpp;
	    }
	}
    }
    cmb = p_cmb;
}

CompositeObject::~CompositeObject() {
    for (_uint i = 0; i < 2; ++i) {
	if (children[i]) {
	    if (child_types[i] == CGS_COMPOSITE)
		delete (CompositeObject*)(children[i]);
	    else
		delete children[i];
	}
    }
}

int CompositeObject::call_child_in(_uint side, const evec3& r) {
    if (child_types[side] == CGS_COMPOSITE) return ((CompositeObject*)children[side])->in(r);
    if (child_types[side] == CGS_SPHERE) return ((Sphere*)children[side])->in(r);
    if (child_types[side] == CGS_BOX) return ((Box*)children[side])->in(r);
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
    //if the right child is null then we don't apply any operations
    if (children[1] == NULL) {
	if (children[0] == NULL)
	    return invert;
	else
	    return invert ^ call_child_in(0, r);
    } else {
	if (children[0] == NULL) {
	    return invert ^ call_child_in(1, r);
	} else {
	    int l_res = call_child_in(0, r);
	    int r_res = call_child_in(1, r);
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
 * Read a string of the format [x, y, z] into an Eigen::Vector3.
 * returns: 0 on success or an error code
 * 	-1: insufficient tokens
 * 	-2: one of the tokens supplied was invalid
 */
parse_ercode parse_vector(char* str, evec3& sto) {
    //find the start and the end of the vector
    char* start = strchr(str, '[');
    char* end = strchr(str, ']');

    char* save_str;
    //read the coordinates separated by spaces
    //x
    char* tok = strtok_r(start+1, ",", &save_str);
    if (!tok) return E_LACK_TOKENS;
    tok = trim_whitespace(tok, NULL);
    sto.x() = strtod(tok, NULL);
    if (errno) { errno = 0;return E_BAD_TOKEN; }
    //y
    tok = strtok_r(NULL, ",", &save_str);
    if (!tok) return E_LACK_TOKENS;
    tok = trim_whitespace(tok, NULL);
    sto.y() = strtod(tok, NULL);
    if (errno) { errno = 0;return E_BAD_TOKEN; }
    //z
    tok = strtok_r(NULL, ",", &save_str);
    if (!tok) return E_LACK_TOKENS;
    tok = trim_whitespace(tok, NULL);
    sto.z() = strtod(tok, NULL);
    if (errno) { errno = 0;return E_BAD_TOKEN; }

    return E_SUCCESS;
}

/**
 * Based on the declaration syntax produce the appropriate geometric shape and store the result in ptr. Ptr must not be initialized before a call to this function to avoid a memory leak.
 * returns: 0 on success or an error code
 */
parse_ercode make_object(const cgs_func& f, Object** ptr, object_type* type, int p_invert) {
    *ptr = NULL;
    parse_ercode error = E_SUCCESS;

    if (!f.name) return E_BAD_TOKEN;
    //switch between all potential types
    if (strcmp(f.name, "Box") == 0) {
	if (f.n_args < 2) return E_LACK_TOKENS;
	evec3 corner;
	evec3 size;
	//read both vector arguments and throw errors if necessary
	if ((error = parse_vector(f.args[0], corner)) != E_SUCCESS) return error;
	if ((error = parse_vector(f.args[1], size)) != E_SUCCESS) return error;
	if (ptr) *ptr = new Box(corner, size, p_invert);
	if (type) *type = CGS_BOX;
    } else if (strcmp(f.name, "Sphere") == 0) {
	if (f.n_args < 2) return E_LACK_TOKENS;
	evec3 center;
	double rad = strtod(f.args[1], NULL);
	if (errno) { errno = 0;return E_BAD_TOKEN; }
	//read both vector arguments and throw errors if necessary
	if ((error = parse_vector(f.args[0], center)) != E_SUCCESS) return error;
	if (ptr) *ptr = new Sphere(center, rad, p_invert);
	if (type) *type = CGS_SPHERE;
    } else if (strcmp(f.name, "Cylinder") == 0) {
	if (f.n_args < 3) return E_LACK_TOKENS;
	evec3 center;
	double h = strtod(f.args[1], NULL);
	if (errno) { errno = 0;return E_BAD_TOKEN; }
	double r1 = strtod(f.args[2], NULL);
	if (errno) { errno = 0;return E_BAD_TOKEN; }
	//by default assume that the radii are the same
	double r2 = r1;
	if (f.n_args > 3) {
	    r2 = strtod(f.args[2], NULL);
	    if (errno) { errno = 0;return E_BAD_TOKEN; }
	}
	//read both vector arguments and throw errors if necessary
	if ((error = parse_vector(f.args[0], center)) != E_SUCCESS) return error;
	if (ptr) *ptr = new Cylinder(center, h, r1, r2, p_invert);
	if (type) *type = CGS_CYLINDER;
    } else if (strcmp(f.name, "Composite") == 0) {
	if (ptr) *ptr = new CompositeObject(CGS_UNION, f, p_invert);
	if (type) *type = CGS_COMPOSITE;
    } else if (strcmp(f.name, "union") == 0) {
	if (ptr) *ptr = new CompositeObject(CGS_UNION, f, p_invert);
	if (type) *type = CGS_COMPOSITE;
    } else if (strcmp(f.name, "intersect") == 0) {
	if (ptr) *ptr = new CompositeObject(CGS_INTERSECT, f, p_invert);
	if (type) *type = CGS_COMPOSITE;
    } else if (strcmp(f.name, "difference") == 0) {
	if (ptr) *ptr = new CompositeObject(CGS_DIFFERENCE, f, p_invert);
	if (type) *type = CGS_COMPOSITE;
    } else {
	if (ptr) *ptr = NULL;
	if (type) *type = CGS_UNDEF;
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
parse_ercode parse_func(char* token, size_t open_par_ind, cgs_func& f, char** end) {
    //by default we want to indicate that we didn't get to the end
    if (end) *end = NULL;
    f.n_args = 0;
    if (token[open_par_ind] != '(') return E_BAD_TOKEN;
    token[open_par_ind] = 0;
    f.name = trim_whitespace(token, NULL);

    //now remove whitespace from the ends of the string
    char* arg_str = token+open_par_ind+1;

    //keep track of information for individual arguments
    int s_block_in = 0;
    int p_block_in = 0;
    int q_block_in = 0;
    size_t j = 0;
    bool tok_start = false;
    char* cur_tok = arg_str;
    size_t last_non_space = 0;
    size_t i = 0;
    for (; arg_str[i] != 0; ++i) {
	//it is possible that there is an array or function inside one of the arguments, we should ignore field separators in those blocks
	if (arg_str[i] == '[') {
	    ++s_block_in;
	} else if (arg_str[i] == ']') {
	    --s_block_in;
	    if (s_block_in < 0) return E_BAD_SYNTAX;
	} else if (arg_str[i] == '(') {
	    ++p_block_in;
	} else if (arg_str[i] == ')') {
	    --p_block_in;
	    //when we reach an end paren without a matching open paren we should stop reading the function
	    if (p_block_in < 0) break;
	} else if (arg_str[i] == '\"') {
	    //TODO: handle single and double quotes separately
	    q_block_in = 1 - q_block_in;
	}
	
	//now if this character is a field separator we add it to the argument list
	if (arg_str[i] == ',' && s_block_in == 0 && p_block_in == 0 && q_block_in == 0) {
	    if (!tok_start) return E_LACK_TOKENS;
	    arg_str[last_non_space+1] = 0;
	    f.args[j++] = cur_tok;
	    //now reset the current token
	    cur_tok = arg_str + i + 1;
	    tok_start = false;
	    continue;
	}

	//we want to remove whitespace surrounding the arguments
	if (arg_str[i] == ' ' || arg_str[i] == '\t') {
	    if (!tok_start) {
		cur_tok += 1;
	    }
	} else {
	    last_non_space = i;
	    tok_start = true;
	}
    }
    //this means that there wasn't a closing parenthesis found
    if (arg_str[i] == 0) return E_BAD_SYNTAX;
    //add the last token
    if (j > 0) {
	if (!tok_start) return E_LACK_TOKENS;
	arg_str[last_non_space+1] = 0;
	f.args[j++] = cur_tok;
	arg_str[i] = 0;
    }
    f.n_args = j;
    //assign the end pointer if the caller wants it
    if (end) *end = arg_str + i + 1;

    return E_SUCCESS;
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
    swap(*this, o);
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
 * Insert an object onto the stack obeying the desired tree structure.
 */
parse_ercode ObjectStack::emplace_obj(Object* obj, object_type p_type) {
    if (stack_ptr == 0) {
	if (p_type == CGS_COMPOSITE) {
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
	if (p_type == CGS_COMPOSITE) {
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

Scene::Scene(const char* p_fname) {
    char buf[BUF_SIZE];
    char func_name[BUF_SIZE];
    size_t n_args = 0;
    bool in_comment = false;

    //we need to store whatever the current token is
    char* cur_token;

    //indicate that no names are stored in all of the specified fields
    func_name[0] = 0;
    //generate the composite object on the fly
    ObjectStack tree_pos;
    CGS_Stack<block_type> blk_stack;
    block_type last_type = BLK_UNDEF;
    int invert = 0;

    //open the file for reading and read individual lines. We need to remove whitespace. Handily, we know the line length once we're done.
    FILE* fp = fopen(p_fname, "r");
    size_t lineno = 1;
    size_t line_len = 0;
    while (fgets(buf, BUF_SIZE, fp)) {
	char* red_str = trim_whitespace(buf, &line_len);
	cur_token = buf;
	for (size_t i = 0; i < line_len; ++i) {
	    if (!in_comment) {
		if (buf[i] == '(') {
		    if (last_type != BLK_UNDEF) printf("Error on line %d: Expected '{' before function name", lineno);

		    //initialize a new cgs_func with the appropriate arguments
		    cgs_func cur_func;
		    char* endptr;
		    parse_ercode er = parse_func(cur_token, i, cur_func, &endptr);
		    switch (er) {
			case E_BAD_TOKEN: printf("Error on line %d: Invalid function name \"%s\"", lineno, cur_token);break;
			case E_BAD_SYNTAX: printf("Error on line %d: Invalid syntax", lineno);break;
			default: break;
		    }
		    //try interpreting the function as a geometric object
		    Object* obj = NULL;
		    object_type type;
		    make_object(cur_func, &obj, &type, invert);
		    if (!obj) {
			//if that failed try interpreting it as an operation (TODO: are there any useful things to put here?)
			if (strcmp(cur_func.name, "invert") == 0) {
			    last_type = BLK_INVERT;
			}
		    } else {
			tree_pos.emplace_obj(obj, type);
			objects.push_back(obj);
		    }
		//check for comments
		} else if (buf[i] == '/') {
		    if (i < line_len-1) {
			if (buf[i+1] == '/')
			    break;
			else if (buf[i+1] == '*')
			    in_comment = true;
		    }
		//check for blocks
		} else if (buf[i] == '{') {
		    switch (last_type) {
			case BLK_INVERT: invert = 1 - invert;break;
			default: break;
		    }
		    blk_stack.push(last_type);
		    last_type = BLK_UNDEF;
		} else if (buf[i] == '}') {
		    block_type bt;
		    if (blk_stack.pop(&bt) == E_EMPTY_STACK) printf("Error on line %d: unexpected '}'", lineno);
		    switch (bt) {
			case BLK_INVERT: invert = 1 - invert;break;
			default: break;
		    }
		}
	    } else if (buf[i] == '*' && i < line_len-1 && buf[i+1] == '/') {
		in_comment = false;
	    }
	}
	++lineno;
    }
}

Scene::~Scene() {
    //TODO: double frees are bad lol
    for (_uint i = 0; i < objects.size(); ++i) {
	if (objects[i]) delete objects[i];
    }
}

CompositeObject* Scene::get_root() {
    if (objects.size() == 0) return NULL;
    return (CompositeObject*)objects[0];
}
