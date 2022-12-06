#define DEBUG_INFO 1

#ifndef CGS_H
#define CGS_H

#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <unordered_map>
#include <vector>
#include <cstring>
#include <math.h>

//hints for dynamic buffer sizes
#define BUF_SIZE 1024
#define ARGS_BUF_SIZE 256
#define FUNC_BUF_SIZE 8

#define DEF_IM_RES	128
#define IM_DEPTH	255

#define SIDE_END	2
#define SIDE_UNDEF	3

#define FLAG_NORMAL	0
#define FLAG_OPER	1
#define FLAG_COMP	2
#define FLAG_AND	3

//keywords
#define KEY_CLASS_LEN	5
#define KEY_FOR_LEN	3
#define KEY_DEF_LEN	3
#define KEY_IN_LEN	2
#define KEY_IF_LEN	2

typedef enum {OP_EQ, OP_ADD, OP_SUB, OP_MULT, OP_DIV, OP_NOT, OP_OR, OP_AND, OP_GRT, OP_LST, OP_GEQ, OP_LEQ, N_OPTYPES} Optype_e;

typedef unsigned int _uint;
typedef unsigned char _uint8;

typedef enum { E_SUCCESS, E_NOFILE, E_LACK_TOKENS, E_BAD_TOKEN, E_BAD_SYNTAX, E_BAD_VALUE, E_NOMEM, E_EMPTY_STACK, E_NOT_BINARY, E_NAN, E_NOT_DEFINED } parse_ercode;

typedef enum { CGS_UNION, CGS_INTERSECT, CGS_DIFFERENCE, CGS_CMB_NOOP } combine_type;
//note that ROOTS are a special type of COMPOSITES
typedef enum { CGS_UNDEF, CGS_ROOT, CGS_DATA, CGS_COMPOSITE, CGS_SPHERE, CGS_BOX, CGS_PLANE, CGS_CYLINDER } object_type;

typedef Eigen::Matrix3d emat3;
typedef Eigen::Vector3d evec3;
typedef Eigen::Vector4d evec4;

/**
 * A helper class (and enumeration) that can track the context and scope of a curly brace block while parsing a file
 */
typedef enum {BLK_UNDEF, BLK_MISC, BLK_INVERT, BLK_TRANSFORM, BLK_DATA, BLK_ROOT, BLK_COMPOSITE, BLK_FUNC_DEC, BLK_LITERAL, BLK_COMMENT, BLK_SQUARE, BLK_QUOTE, BLK_QUOTE_SING, BLK_PAREN, BLK_CURLY} block_type;
template <typename T>
class CGS_Stack {
protected:
    size_t stack_ptr;
    size_t buf_len;
    T* buf;
    parse_ercode grow(size_t new_size);

public:
    CGS_Stack();
    ~CGS_Stack<T>();
    void swap(CGS_Stack<T>& o);
    CGS_Stack(const CGS_Stack<T>& o);
    CGS_Stack(CGS_Stack<T>&& o);
    CGS_Stack<T>& operator=(CGS_Stack<T> o);

    parse_ercode push(T b);
    parse_ercode pop(T* ptr);
    void reset();
    bool has(const T& key);
    bool is_empty() { return stack_ptr == 0; }
    size_t size() { return stack_ptr; }
    T peek(size_t ind = 1);
};

class CompositeObject;

/**
  * Remove the whitespace surrounding a word
  */
inline char* CGS_trim_whitespace(char* str, size_t* len) {
    if (!str) return NULL;
    size_t start_ind = 0;
    bool started = false;
    _uint last_non = 0;
    for (_uint i = 0; str[i] != 0; ++i) {
        if (str[i] != ' ' && str[i] != '\t' && str[i] != '\n') {
            last_non = i;
            if (!started) {
                start_ind = i;
		started = true;
            }
        }
    }
    str[last_non+1] = 0;
    if (len) *len = last_non - start_ind+1;
    return str+start_ind;
}

/**
  * Remove the whitespace surrounding a word and return a copy which must be deallocated by calling free()
  */
/*inline char* CGS_trim_whitespace(const char* str, size_t* len) {
    if (!str) return NULL;
    size_t start_ind = 0;
    bool started = false;
    _uint last_non = 0;
    for (_uint i = 0; str[i] != 0; ++i) {
        if (str[i] != ' ' && str[i] != '\t' && str[i] != '\n') {
            last_non = i;
            if (!started) {
                start_ind = i;
		started = true;
            }
        }
    }
    _uint ret_len = last_non+2-start_ind;
    char* ret = (char*)malloc(sizeof(char)*ret_len);
    if (ret == NULL) return NULL;
    for (_uint i = start_ind; i <= last_non; ++i) { ret[i-start_ind] = str[i]; }
    ret[ret_len-1] = 0;
    if (len) *len = last_non - start_ind+1;
    return ret;
}*/

/**
 * This acts similar to getline, but stops at a semicolon, newline (unless preceeded by a \), {, or }.
 * bufptr: a pointer to which the buffer is saved. If bufptr is NULL than a new buffer is allocated through malloc()
 * n: a pointer to a size_t with the number of characters in the buffer pointed to by bufptr. The call will return do nothing if n is null but *bufptr is not.
 * fp: file pointer to read from
 * linecount: a pointer to an integer specifying the number of new line characters read.
 * Returns: 0 if the end of the file was reached, 1 otherwise
 */
int read_cgs_line(char** bufptr, size_t* n, FILE* fp, size_t* line);

/*
 * A virtual class which describes a simple shape (such as cubes, cylinders spheres etc)
 */
class Object {
private:
    friend class GeometricObject;
protected:
    emat3 trans_mat;
    int invert;

public:
    Object(int p_invert=0);
    Object(const Eigen::Quaterniond& orientation, int p_invert=0);
    virtual int in(const evec3& r) = 0;
    void set_inversion(int p_invert);
    void rescale(const evec3& components);
    void set_trans_mat(const emat3& new_mat) { trans_mat = new_mat; }
};

class Sphere : public Object {
private:
    evec3 center;
    double rad;

public:
    Sphere(evec3& p_center, double p_rad, int p_invert=0);
    Sphere(evec3& p_center, evec3& p_rad, int p_invert=0);
    int in(const evec3& r);

    evec3 get_center() const { return center; }
    double get_rad() const { return rad; }
};

class Box : public Object {
private:
    evec3 center;
    evec3 offset;

public:
    Box(evec3& p_corner_1, evec3& p_corner_2, int p_invert=0);
    Box(evec3& p_corner_1, evec3& p_corner_2, Eigen::Quaterniond p_orientation, int p_invert=0);
    int in(const evec3& r);

    evec3 get_center() const { return center; }
    evec3 get_offset() const { return offset; }
};

class Plane : public Object {
private:
    evec3 normal;
    double offset;

public:
    Plane(evec3& normal, double offset, int p_invert=0);
    Plane(evec3& point_1, evec3& point_2, evec3& point_3, int p_invert=0);
    int in(const evec3& r);

    evec3 get_normal() const { return normal; }
    double get_offset() const { return offset; }
};

class Cylinder : public Object {
private:
    evec3 center;
    double height;
    double r1_sq;
    double r1_sq_x_h;
    double r2_sq;

public:
    Cylinder(evec3& p_center, double p_height, double p_r1, double p_r2, int p_invert=0);
    int in(const evec3& r);

    evec3 get_center() const { return center; }
    double get_height() const { return height; }
    double get_r1() const { return sqrt(r1_sq); }
    double get_r2() const { return sqrt(r2_sq); }
};

class Scene;
typedef enum {VAL_UNDEF, VAL_STR, VAL_NUM, VAL_LIST, VAL_3VEC, VAL_MAT, VAL_INST} valtype;
class Value;
class cgs_func;
class instance;
union V {
    char* s;
    double x;
    Value* l;
    evec3* v;
    instance* i;
    Eigen::MatrixXd* m;
};
struct Value {
    valtype type;
    V val;
    size_t n_els; //only applicable for string and list types

    Value() { val.x = 0;type=VAL_UNDEF;n_els=0; }
    bool operator==(std::string str);
    bool operator!=(std::string str);
    valtype get_type() { return type; }
    size_t size() { return n_els; }
    V get_val() { return val; }
    char* to_c_str();
    double to_float();
    Value cast_to(valtype type, parse_ercode& er) const;
};
void cleanup_val(Value* o);
Value copy_val(const Value o);
void swap_val(Value* a, Value* b);

struct cgs_func {
    char* name;
    Value args[ARGS_BUF_SIZE];
    char* arg_names[ARGS_BUF_SIZE];
    size_t n_args;
    cgs_func() { name = NULL;n_args = 0; }
};
Value make_val_str(const char* s);
Value make_val_std_str(std::string s);
Value make_val_list(const Value* vs, size_t n_vs);
Value make_val_mat(Eigen::MatrixXd m);
Value make_val_vec3(evec3 vec);
cgs_func parse_func_decl(char* str);
cgs_func copy_func(const cgs_func o);
void cleanup_func(cgs_func* o);
void swap(cgs_func* a, cgs_func* b);

/*
 * A composite object is a node in a binary tree that represents one or more primitives combined via unions and intersections
 */
class CompositeObject : public Object {
    friend class Scene;
protected:
    Object* children[2];
    object_type child_types[2];
    combine_type cmb;
    bool is_leaf;

    //metadata allows users to attach miscelaneous information to an object such as a name or material properties
    std::unordered_map<std::string, Value> metadata;

    //since we need to perform casts to use the appropriate object type, it's helpful to add a layer of abstraction
    int call_child_in(_uint side, const evec3& r);
    //this is a helper function which returns a pointer to a copy of the object pointed to by the child on the specified side
    Object* copy_child(_uint side) const;

public:
    CompositeObject(combine_type p_cmb = CGS_UNION);
    CompositeObject(combine_type p_cmb, const cgs_func& spec, int p_invert);
    CompositeObject(const CompositeObject& o);
    CompositeObject(CompositeObject&& o);
    CompositeObject& operator=(CompositeObject&& o);
    ~CompositeObject();

    //add a child to the composite object by parsing the string description
    void add_child(_uint side, Object* o, object_type p_type);
    int in(const evec3& r);
    const Object* get_child_l() { return children[0]; }
    const Object* get_child_r() { return children[1]; }
    object_type get_child_type_l() const { return child_types[0]; }
    object_type get_child_type_r() const { return child_types[1]; }
    combine_type get_combine_type() const { return cmb; }
    int has_metadata(std::string key) const { return metadata.count(key); }
    Value fetch_metadata(std::string key) { return metadata[key]; }
};

struct type_ind_pair {
public:
    block_type t;
    size_t i;
    type_ind_pair() { t = BLK_UNDEF;i = 0; }
    type_ind_pair(size_t ii) { t = BLK_UNDEF;i = ii; }
    type_ind_pair(block_type tt, size_t ii) { t = tt;i = ii; }
};

/**
 * A helper class for scene which maintains two parallel stacks used when constructing binary trees. The first specifies the integer code for the side occupied and the second stores pointers to objects. Each index specified in the side array is a "tribit" with 0 specifying the left side, 1 the right and 2 the end end of the stack
 * WARNING: The stack only handles pointers to objects. The caller is responsible for managing the lifetime of each object placed on the stack.
 */
struct side_obj_pair {
public:
    _uint side;
    CompositeObject* obj;
    side_obj_pair(_uint p_side, CompositeObject* p_obj) { side = p_side;obj=p_obj; }
    side_obj_pair(const side_obj_pair& o) { side = o.side;obj = o.obj; }
    side_obj_pair(side_obj_pair&& o) { side = o.side;obj = o.obj;o.side = 0;o.obj = NULL; }
    side_obj_pair& operator=(const side_obj_pair& o) { side = o.side;obj = o.obj;return *this; }
};
class ObjectStack : public CGS_Stack<side_obj_pair> {
public:
    //parse_ercode push(_uint side, CompositeObject* obj);
    parse_ercode emplace_obj(Object* obj, object_type p_type);
    parse_ercode pop(_uint* side, CompositeObject** obj);
    _uint look_side();
    CompositeObject* look_obj();
    CompositeObject* get_root();
};
/**
 * A helper class for scene which maintains two parallel stacks used when constructing binary trees. The first specifies the integer code for the side occupied and the second stores pointers to objects. Each index specified in the side array is a "tribit" with 0 specifying the left side, 1 the right and 2 the end end of the stack
 * WARNING: The stack only handles pointers to objects. The caller is responsible for managing the lifetime of each object placed on the stack.
 * NOTE: the context performs a shallow copy of everything pushed onto it and does not perform any cleanup
 */
class name_val_pair {
private:
    char* name;
    Value val;
public:
    void swap(name_val_pair& o);
    void copy(const name_val_pair& o) { name = strdup(o.name);val = copy_val(o.val); }
    name_val_pair(const char* p_name, Value p_val) { name = strdup(p_name);val = copy_val(p_val); }
    name_val_pair(const name_val_pair& o) { copy(o); }
    name_val_pair(name_val_pair&& o) { name = o.name;val = o.val;o.name = NULL;o.val.type = VAL_UNDEF;o.val.val.x = 0; }
    ~name_val_pair() {
	if (name) free(name);
	if (val.type != VAL_UNDEF) cleanup_val(&val);
    }
    name_val_pair& operator=(name_val_pair& o) { swap(o);return *this; }
    name_val_pair& operator=(const name_val_pair& o) { copy(o);return *this; }
    bool name_matches(const char* str) const { if (strcmp(str, name) == 0) return true;return false; }
    Value& get_val() { return val; }
};
class instance {
private:
    std::vector<name_val_pair> fields;
    std::string type;
public:
    instance(cgs_func decl);
};
class context : public CGS_Stack<name_val_pair> {
private:
    Value do_op(char* tok, size_t ind, parse_ercode& er);
public:
    //parse_ercode push(_uint side, CompositeObject* obj);
    void emplace(const char* p_name, Value p_val) { name_val_pair inst(p_name, p_val);push(inst); }
    Value lookup(const char* name) const;
    parse_ercode pop_n(size_t n);
    Value parse_value(char* tok, parse_ercode& er);
    cgs_func parse_func(char* token, long open_par_ind, parse_ercode& f, char** end);
    Value parse_list(char* str, parse_ercode& sto);
    void swap(CGS_Stack<name_val_pair>& o) { CGS_Stack<name_val_pair>::swap(o); }
    parse_ercode set_value(const char* name, Value new_val);
};
/**
 * A class for functions defined by the user along with the implementation code
 */
class user_func {
private:
    cgs_func call_sig;
    //this is a dynamic buffer for the number of lines
    size_t n_lines;
    char** code_lines;
public:
    //read the function with contents stored in the file pointer fp at the current file position
    user_func(cgs_func sig, char** bufptr, size_t* n, FILE* fp);
    ~user_func();
    user_func(const user_func& o);
    user_func(user_func&& o);
    size_t get_n_lines() { return n_lines; }
    Value eval(context& c, cgs_func call, parse_ercode& er);
};
class Scene {
private:
    //std::vector<Object*> objects;
    std::vector<CompositeObject*> roots;
    std::vector<CompositeObject*> data_objs;

    //let users define constants
    context named_items;
    parse_ercode fail_exit(parse_ercode er, FILE* fp);
    parse_ercode read_file(const char* p_fname);

public:
    Scene() {}
    Scene(const char* p_fname, parse_ercode* ercode = NULL);
    Scene(const Scene& o);
    Scene(Scene&& o);
    Scene& operator=(Scene& o);
    ~Scene();

    std::vector<CompositeObject*> get_roots() { return roots; }
    std::vector<CompositeObject*> get_data() { return data_objs; }
    void read();
 
    parse_ercode make_object(const cgs_func& f, Object** ptr, object_type* type, int p_invert) const;
    parse_ercode make_transformation(const cgs_func& f, emat3& res) const;
    context& get_context() { return named_items; }
    //void cleanup_func(cgs_func& f);
};

#endif //CGS_H
