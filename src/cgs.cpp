#include "cgs.hpp"

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

Value::Value(const char* s) {
    type = VAL_STR;
    n_els = strlen(s) + 1;
    val.s = (char*)malloc(sizeof(char)*n_els);
    for (size_t i = 0; i < n_els; ++i) val.s[i] = s[i];
}
Value::Value(std::string s) {
    type = VAL_STR;
    n_els = s.size();
    val.s = (char*)malloc(sizeof(char)*n_els);
    for (size_t i = 0; i < n_els; ++i) val.s[i] = s[i];
}
Value::Value(const Value* vs, size_t n_vs) {
    type = VAL_LIST;
    n_els = n_vs;
    val.l = (Value*)malloc(sizeof(Value)*n_els);
    for (size_t i = 0; i < n_els; ++i) val.l[i] = vs[i];
}
Value::Value(Eigen::MatrixXd m) {
    type = VAL_MAT;
    n_els = 1;
    val.m = new Eigen::MatrixXd(m);
}
Value::Value(evec3 vec) {
    type = VAL_3VEC;
    n_els = 1;
    val.v = new evec3(vec);
}
Value::~Value() {
    if (type == VAL_STR && val.s) {
	free(val.s);
    } else if (type == VAL_LIST && val.l) {
	free(val.l);
    } else if (type == VAL_MAT && val.m) {
	delete val.m;
    } else if (type == VAL_3VEC && val.v) {
	delete val.v;
    }
}
//copy 
Value::Value(const Value& o) {
    type = o.type;
    n_els = o.n_els;
    //strings or lists must be copied
    if (type == VAL_STR) {
	val.s = (char*)malloc(sizeof(char)*n_els);
	for (size_t i = 0; i < n_els; ++i) val.s[i] = o.val.s[i];
    } else if (type == VAL_LIST) {
	val.l = (Value*)malloc(sizeof(Value)*n_els);
	for (size_t i = 0; i < n_els; ++i) val.l[i] = o.val.l[i];
    } else if (type == VAL_MAT) {
	val.m = new Eigen::MatrixXd(*(o.val.m));
    } else if (type == VAL_3VEC) {
	val.v = new evec3(*(o.val.v));
    } else {
	val.x = o.val.x;
    }
}
//move
Value::Value(Value&& o) {
    type = o.type;
    n_els = o.n_els;
    if (type == VAL_STR) {
	val.s = o.val.s;
	o.val.s = NULL;
    } else if (type == VAL_LIST) {
	val.l = o.val.l;
	o.val.l = NULL;
    } else if (type == VAL_MAT) {
	val.m = o.val.m;
	o.val.m = NULL;
    } else if (type == VAL_3VEC) {
	val.v = o.val.v;
	o.val.v = NULL;
    } else {
	val.x = o.val.x;
    }
    o.type = VAL_UNDEF;
    o.n_els = 0;
}
//assign
/*Value& Value::operator=(const Value& o) {
    type = o.type;
    n_els = o.n_els;
    val = o.val;
    return *this;
}*/
Value& Value::operator=(Value o) {
    //swap type and number of elements
    valtype tmp = type;
    type = o.type;
    o.type = tmp;
    size_t tmp_n = n_els;
    n_els = o.n_els;
    o.n_els = tmp_n;
    V tmp_v = val;
    val = o.val;
    o.val = tmp_v;
    return *this;
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
    if (type == VAL_UNDEF) { er = E_BAD_VALUE;return Value(*this); }
    if (type == t) return Value(*this);
    Value ret;
    ret.type = t;
    if (type == VAL_LIST) {
	if (t == VAL_3VEC) {
	    Value* tmp_lst = val.l;
	    //check that we have at least three numeric values
	    if (n_els < 3) { er = E_LACK_TOKENS;return Value(*this); }
	    if (tmp_lst[0].type != VAL_NUM || tmp_lst[1].type != VAL_NUM || tmp_lst[2].type != VAL_NUM) return E_BAD_TOKEN;
	    //actually change data and free the old
	    ret.val.v = new evec3(tmp_lst[0].val.x, tmp_lst[1].val.x, tmp_lst[2].val.x);
	    ret.n_els = 3;
	}
    } else if (type == VAL_STR) {
	//*this = parse_value(val.s);
    }
    return ret;
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
		}
	    }*/
	    std::string tok_cpp(spec.arg_names[i]);
	    metadata[tok_cpp] = spec.args[i];
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

parse_ercode Scene::lookup_val(char* tok, Value& sto) const {
    errno = 0;
     std::string cpp_tok(tok);
    if (named_items.count(cpp_tok) == 0) return E_BAD_TOKEN;
    sto = named_items.at(cpp_tok);
    return E_SUCCESS;
}

parse_ercode Scene::parse_value(char* str, Value& sto) const {
    parse_ercode ret = E_SUCCESS;

    //store locations of the first instance of different operators. We do this so we can quickly look up new operators if we didn't find any other operators of a lower precedence (such operators are placed in the tree first).
    int first_open_ind = -1;
    int last_close_ind = -1;

    //keeps track of open and close [], (), {}, and ""
    CGS_Stack<type_ind_pair> blk_stk;
    type_ind_pair start_ind;

    //keep track of the nesting level within parenthetical statements
    int nest_level = 0;
    int found_valid_op = 0;//boolean

    //boolean (if true, then '+' and '-' operators are interpreted as signs instead of operations. We set to true to ensure this behaviour is performed for the first character
    int ignore_op = 1;

    //first try to find base scope addition and subtraction operations
    Value tmp_l, tmp_r;
    for (_uint i = 0; str[i] != 0; ++i) {
	if (str[i] == '['/*]*/) {
	    blk_stk.push({BLK_SQUARE, i});
	    if (first_open_ind == -1) { first_open_ind = i; }
	} else if (str[i] == /*[*/']') {
	    if (blk_stk.pop(&start_ind) != E_SUCCESS || start_ind.t != BLK_SQUARE) return E_BAD_SYNTAX;
	    //now parse the argument as a vector
	    parse_ercode tmp_er = parse_list(str+start_ind.i, sto);
	    if (tmp_er != E_SUCCESS) return tmp_er;
	    if (blk_stk.is_empty() && last_close_ind < 0) last_close_ind = i;
	} else if (str[i] == '('/*)*/) {
	    //keep track of open and close parenthesis, these will come in handy later
	    blk_stk.push({BLK_PAREN, i});
	    //TODO: decide if this is necessary
	    ++nest_level;
	    //only set the open index if this is the first match
	    if (first_open_ind == -1) { first_open_ind = i; }
	} else if (str[i] == /*(*/')') {
	    if (blk_stk.pop(&start_ind) != E_SUCCESS || start_ind.t != BLK_PAREN) return E_BAD_SYNTAX;
	    --nest_level;
	    //only set the end paren location if it hasn't been set yet and the stack has no more parenthesis to remove, TODO: make this work with other block types inside a set of parenthesis
	    if (blk_stk.is_empty() && last_close_ind < 0) last_close_ind = i;
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
	    //we need to store whether this was a valid operator separately in order to properly handle parenthetical groups. This valid_op can take on three values FLAG_NORMAL FLAG_OPER and FLAG_AND. The last flag is used to indicate that a logical operation of && or || was obtained. These operations should always be evaluated last (which means they are placed closer to the root of the tree).
	    int this_valid_op = FLAG_NORMAL;
	    //keep track of the number of characters used by the operator
	    int code_n_chars = 1;
	    //check if we found a numeric operation symbol
	    if ((str[i] == '=' && str[i+1] == '=')
	    || str[i] == '>' || str[i] == '<'
	    || str[i] == '+' || str[i] == '-' || str[i] == '*' || str[i] == '/') {
		//Store the operation before setting it to zero
		char term_char = str[i];
		str[i] = 0;
		//parse right and left values
		parse_ercode er = parse_value(str, tmp_l);
		if (er != E_SUCCESS) return er;
		er = parse_value(str+i+1, tmp_r);
		if (er != E_SUCCESS) return er;
		sto.type = VAL_NUM;
		//remember to recurse after we finish looping
		found_valid_op = 1;

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
		    if (tmp_l.type != VAL_NUM || tmp_r.type != VAL_NUM) return E_BAD_VALUE;
		    if (str[i+1] == '=') {
			sto.val.x = (tmp_l.val.x >= tmp_r.val.x)? 1: 0;
			++i;
		    } else {
			sto.val.x = (tmp_l.val.x > tmp_r.val.x)? 1: 0;
		    }
		} else if (term_char == '<') {
		    if (tmp_l.type != VAL_NUM || tmp_r.type != VAL_NUM) return E_BAD_VALUE;
		    if (str[i+1] == '=') {
			sto.val.x = (tmp_l.val.x <= tmp_r.val.x)? 1: 0;
			++i;
		    } else {
			sto.val.x = (tmp_l.val.x < tmp_r.val.x)? 1: 0;
		    }
		} else if (term_char == '+' && (i == 0 || str[i-1] != 'e')) {
		    if (i == 0 || str[i-1] == 'e') {
			found_valid_op = 0;
			str[i] = term_char;
		    } else if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_NUM) {	
			sto.val.x = tmp_l.val.x + tmp_r.val.x;
		    } else if (tmp_l.type == VAL_MAT && tmp_r.type == VAL_MAT) {	
			sto.type = VAL_MAT;
			sto.val.m = new Eigen::MatrixXd(*tmp_l.val.m + *tmp_r.val.m);
		    } else if (tmp_l.type == VAL_STR) {
			//TODO: implement string concatenation
		    }
		} else if (term_char == '-') {
		    //if this is a sign or exponent, then reset the string
		    if (i == 0 || str[i-1] == 'e') {
			found_valid_op = 0;
			str[i] = term_char;
		    } else if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_NUM) {
			sto.val.x = tmp_l.val.x - tmp_r.val.x;
		    } else if (tmp_l.type == VAL_MAT && tmp_r.type == VAL_MAT) {	
			sto.type = VAL_MAT;
			sto.val.m = new Eigen::MatrixXd(*tmp_l.val.m - *tmp_r.val.m);
		    }
		} else if (term_char == '*') {
		    if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_NUM) {	
			sto.val.x = tmp_l.val.x + tmp_r.val.x;
		    } else if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_MAT) {	
			sto.type = VAL_MAT;
			sto.val.m = new Eigen::MatrixXd(tmp_l.val.x * *tmp_r.val.m);
		    } else if (tmp_l.type == VAL_MAT && tmp_r.type == VAL_MAT) {	
			sto.type = VAL_MAT;
			sto.val.m = new Eigen::MatrixXd(*tmp_l.val.m * *tmp_r.val.m);
		    }
		} else if (term_char == '/') {
		    if (tmp_r.val.x == 0) return E_NAN;
		    if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_NUM) {
			sto.val.x = tmp_l.val.x / tmp_r.val.x;
		    } else if (tmp_r.type == VAL_NUM && tmp_l.type == VAL_MAT) {	
			sto.type = VAL_MAT;
			sto.val.m = new Eigen::MatrixXd(*tmp_r.val.m / tmp_r.val.x);
		    }
		}
	    }
	}
    }
    if (nest_level > 0) {
	printf("expected close paren\n");
	return E_BAD_SYNTAX;
    }

    //last try removing parenthesis
    if (!found_valid_op) {
	//if there isn't a valid parenthetical expression, then we should interpret this as a value string
	if (first_open_ind < 0 || last_close_ind < 0) {
	    str = CGS_trim_whitespace(str, NULL);
	    parse_ercode er = lookup_val(str, sto);
	    if (er != E_SUCCESS) {
		//try interpreting as a number
		sto.val.x = strtod(str, NULL);
		if (errno) return E_BAD_TOKEN;
		sto.type = VAL_NUM;
	    }
	} else if (str[first_open_ind] == '\"' && str[last_close_ind] == '\"') {
	    //this is a string
	    sto.type = VAL_STR;
	    sto.n_els = last_close_ind-first_open_ind;
	    //allocate memory and copy
	    sto.val.s = (char*)malloc(sizeof(char)*sto.n_els);
	    for (size_t i = 0; i < sto.n_els-1; ++i) sto.val.s[i] = str[first_open_ind+i+1];
	    sto.val.s[sto.n_els-1] = 0;
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
		    ret = parse_func(str, first_open_ind, tmp_f, &f_end);
		    if (ret != E_SUCCESS) return ret;
		    //TODO: allow user defined functions
		    if (strcmp(tmp_f.name, "vec") == 0) {
			if (tmp_f.n_args < 3) return E_LACK_TOKENS;
			if (tmp_f.args[0].type != VAL_NUM || tmp_f.args[1].type != VAL_NUM || tmp_f.args[2].type != VAL_NUM) return E_BAD_TOKEN;
			sto.type = VAL_3VEC;
			sto.val.v = new evec3(tmp_f.args[0].val.x, tmp_f.args[1].val.x, tmp_f.args[2].val.x);
			sto.n_els = 3;
			return E_SUCCESS;
		    }
		    str[last_close_ind+1] = term_char;
		    return E_BAD_TOKEN;
		}
	    }
	    //otherwise interpret this as a parenthetical expression
	    str[last_close_ind] = 0;
	    str = CGS_trim_whitespace(str+first_open_ind+1, NULL);
	    ret = parse_value(str, sto);
	}
    }
    return ret;
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
		switch (str[i]) {
		case 'n': str[j++] = '\n';break;
		case 't': str[j++] = '\t';break;
		case '\\': str[j++] = '\\';break;
		case '\"': str[j++] = '\"';break;
		default: er = E_BAD_SYNTAX;
		}
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
	    if (verbatim || (str[i] != ' ' && str[i] != '\t' && str[i] != '\n')) {
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
parse_ercode Scene::parse_list(char* str, Value& sto) const {
    parse_ercode er;
    sto.val.v = new evec3();
    //find the start and the end of the vector
    char* start = strchr(str, '[');
    //make sure that the string is null terminated
    char* end = strchr_block(start+1, ']');
    if (!end) return E_BAD_SYNTAX;
    *end = 0;

    Value tmp_val;

    //read the coordinates separated by spaces
    size_t n_els;
    char** list_els = csv_to_list(start+1, ',', &n_els, er);
    if (er != E_SUCCESS) return er;
    Value* buf = (Value*)malloc(sizeof(Value)*n_els);
    for (size_t i = 0; list_els[i] && i < n_els; ++i) {
	if ((er = parse_value(list_els[i], buf[i])) != E_SUCCESS) return er;
    }
    //cleanup and reset the string
    *end = ']';
    free(list_els);
    //set number of elements and type
    sto.type = VAL_LIST;
    sto.n_els = n_els;
    sto.val.l = buf;
    return E_SUCCESS;
    /*char* save_str;
    char* tok = strtok_r(start+1, ",", &save_str);
    if (!tok) return E_LACK_TOKENS;
    tok = CGS_trim_whitespace(tok, NULL);
    errno = 0;
    if ((er = parse_value(tok, tmp_val)) != E_SUCCESS) return er;
    if (tmp_val.type != VAL_NUM) return E_BAD_VALUE;
    sto.val.v->x() = tmp_val.val.x;
    //y
    tok = strtok_r(NULL, ",", &save_str);
    if (!tok) return E_LACK_TOKENS;
    if ((er = parse_value(tok, tmp_val)) != E_SUCCESS) return er;
    if (tmp_val.type != VAL_NUM) return E_BAD_VALUE;
    sto.val.v->y() = tmp_val.val.x;
    //z
    tok = strtok_r(NULL, ",", &save_str);
    if (!tok) return E_LACK_TOKENS;
    if ((er = parse_value(tok, tmp_val)) != E_SUCCESS) return er;
    if (tmp_val.type != VAL_NUM) return E_BAD_VALUE;
    sto.val.v->z() = tmp_val.val.x;*/

    sto.type = VAL_3VEC;

    return E_SUCCESS;
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
    } else if (strcmp(f.name, "Sphere") == 0) {
	if (f.n_args < 2) return E_LACK_TOKENS;
	//if we have enough tokens make sure we have both elements as vectors
	Value cent = f.args[0].cast_to(VAL_3VEC, er);
	if (er != E_SUCCESS || f.args[1].type != VAL_NUM) return E_BAD_VALUE;
	double rad = f.args[1].val.x;
	if (ptr) *ptr = new Sphere(*(cent.val.v), rad, p_invert);
	if (type) *type = CGS_SPHERE;
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
parse_ercode Scene::parse_func(char* token, long open_par_ind, cgs_func& f, char** end) const {
    parse_ercode er;

    //by default we want to indicate that we didn't get to the end
    if (end) *end = NULL;
    f.n_args = 0;
    //infer the location of the open paren index if the user didn't specify it
    if (open_par_ind < 0 || token[open_par_ind] != '('/*)*/) {
	char* par_char = strchr(token, '('/*)*/);
	//make sure there actually is an open paren
	if (par_char == NULL) return E_BAD_TOKEN;
	open_par_ind = par_char - token;
    }

    //break the string up at the parenthesis and remove surrounding whitespace
    token[open_par_ind] = 0;
    f.name = CGS_trim_whitespace(token, NULL);

    //now remove whitespace from the ends of the string
    char* arg_str = token+open_par_ind+1;
    //make sure that the string is null terminated
    char* term_ptr = strchr_block(arg_str, /*(*/')');
    if (!term_ptr) return E_BAD_SYNTAX;
    *term_ptr = 0;
    if (end) *end = term_ptr+1;

    //read the coordinates separated by spaces
    char** list_els = csv_to_list(arg_str, ',', &(f.n_args), er);
    if (er != E_SUCCESS) return er;
    //make sure that we don't go out of bounds, TODO: let functions accept arbitrarily many arguments?
    if (f.n_args >= ARGS_BUF_SIZE) {
	er = E_NOMEM;
	free(list_els);
	return er;
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
	if ((er = parse_value(list_els[i], f.args[i])) != E_SUCCESS) return er;
    }
    //cleanup and reset the string
    //*term_ptr = /*(*/')';
    free(list_els);
    return E_SUCCESS;

    //keep track of information for individual arguments
    CGS_Stack<type_ind_pair> blk_stk;
    int blk_stk_in = 0;
    int p_block_in = 0;
    size_t j = 0;
    bool tok_start = false;
    char* cur_tok = arg_str;
    char* cur_name = NULL;
    size_t last_non_space = 0;
    size_t i = 0;

    //keep track of the last value
    Value last_val;
    last_val.type = VAL_UNDEF;

    for (; arg_str[i] != 0; ++i) {
	//it is possible that there is an array or function inside one of the arguments, we should ignore field separators in those blocks
	if (arg_str[i] == '[') {
	    blk_stk.push({BLK_SQUARE, i});
	} else if (arg_str[i] == ']') {
	    type_ind_pair start_ind;
	    if (blk_stk.pop(&start_ind) != E_SUCCESS || start_ind.t != BLK_SQUARE) return E_BAD_SYNTAX;
	    //now parse the argument as a vector
	    parse_ercode tmp_er = parse_list(arg_str+start_ind.i, last_val);
	    if (tmp_er != E_SUCCESS) return tmp_er;
	} else if (arg_str[i] == '(') {
	    blk_stk.push({BLK_PAREN, i});
	} else if (arg_str[i] == ')') {
	    type_ind_pair start_ind;
	    //when we reach an end paren without a matching open paren we should stop reading the function
	    if (blk_stk.pop(&start_ind) != E_SUCCESS) break;
	    //if (p_block_in < 0) break;
	} else if (arg_str[i] == '\"') {
	    tok_start = true;
	    arg_str[i] = ' ';
	    size_t str_start_ind = i+1;
	    //skip ahead until we find the end of the quote
	    for (_uint j = i+1; arg_str[j] != 0; ++j) {
		if (arg_str[j] == '\"' && arg_str[j-1] != '\\') {
		    //once we reach the end of the string we need to set the higher level index accordingly
		    last_non_space = j-1;
		    i = j;
		    arg_str[j] = ' ';//assign the character to be ignored
		    break;
		}
	    }
	    //copy the string into the value memory
	    last_val.type = VAL_STR;
	    last_val.n_els = (last_non_space-str_start_ind)+2;
	    last_val.val.s = (char*)malloc(sizeof(char)*last_val.n_els);
	    for (size_t j = 0; j < last_val.n_els; ++j) last_val.val.s[j] = arg_str[str_start_ind+j];
	    last_val.val.s[last_val.n_els-1] = 0;
	    continue;
	} else if (arg_str[i] == '=' && arg_str[i+1] != '=' && arg_str[i+1] != 0) {
	    //this is a named argument, set the name and reset the token to only include the portion after the = sign
	    arg_str[i] = 0;
	    cur_name = CGS_trim_whitespace(cur_tok, NULL);
	    ++i;
	    cur_tok = arg_str+i;
	}
	
	//now if this character is a field separator we add it to the argument list
	if (arg_str[i] == ',' && blk_stk.is_empty()) {
	    if (!tok_start) return E_LACK_TOKENS;
	    arg_str[last_non_space+1] = 0;
	    //check if we need to figure out the value based on context
	    if (last_val.type == VAL_UNDEF) {
		parse_ercode tmp_er = parse_value(cur_tok, last_val);
	    }
	    //set the name and its value
	    if (cur_name) f.arg_names[j] = CGS_trim_whitespace(cur_name, NULL);
	    f.args[j++] = last_val;
	    //now reset the current token
	    cur_name = NULL;
	    cur_tok = arg_str + i + 1;
	    tok_start = false;
	    //make sure we reset the value type
	    last_val.type = VAL_UNDEF;
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
    //if (j > 0) {
        if (!tok_start) return E_LACK_TOKENS;
        arg_str[last_non_space+1] = 0;
	//parse the final value
	if (last_val.type == VAL_UNDEF) {
	    parse_ercode tmp_er = parse_value(cur_tok, last_val);
	}
	if (cur_name) f.arg_names[j] = CGS_trim_whitespace(cur_name, NULL);
        f.args[j++] = last_val;
        arg_str[i] = 0;
    //}
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

parse_ercode Scene::read_file(const char* p_fname) {
    parse_ercode er = E_SUCCESS;
    char buf[BUF_SIZE];
    char func_name[BUF_SIZE];
    size_t n_args = 0;

    //we need to store whatever the current token is
    char* cur_token;

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
	size_t lineno = 1;
	size_t line_len = 0;
	char last_char = 0;
	//iterate over each line in the file
	while (fgets(buf, BUF_SIZE, fp)) {
	    char* red_str = CGS_trim_whitespace(buf, &line_len);
	    cur_token = buf;
	    for (size_t i = 0; i < line_len; ++i) {
		//check the most recent block pushed onto the stack
		block_type cur_type = BLK_UNDEF;
		if (blk_stack.size() > 0) cur_type = blk_stack.peek();
		//only interpret as normal code if we aren't in a comment or literal block
		if (cur_type != BLK_COMMENT && cur_type != BLK_LITERAL) {
		    if (buf[i] == '(' && blk_stack.peek() != BLK_LITERAL) {
			//initialize a new cgs_func with the appropriate arguments
			cgs_func cur_func;
			char* endptr;
			er = parse_func(cur_token, i, cur_func, &endptr);
			switch (er) {
			    case E_BAD_TOKEN:
				printf("Error on line %d: Invalid function name \"%s\"\n", lineno, cur_token);
				return fail_exit(er, fp);
			    case E_BAD_SYNTAX:
				printf("Error on line %d: Invalid syntax\n", lineno);
				return fail_exit(er, fp);
			    default: break;
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
			char* tok = CGS_trim_whitespace(buf, NULL);
			size_t val_len;
			char* val = CGS_trim_whitespace(buf+i+1, &val_len);
			named_items[std::string(tok)] = std::string(val);
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
	    ++lineno;
	}
	fclose(fp);
    } else {
        printf("Error: couldn't open file %s for reading!\n", p_fname);
        return E_NOFILE;
    }

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
    o.named_items.clear();
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
