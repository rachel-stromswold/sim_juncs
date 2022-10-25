#define DEBUG_INFO 1

#ifndef CGS_H
#define CGS_H

#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <unordered_map>
#include <vector>
#include <cstring>
#include <math.h>

#define BUF_SIZE 1024
#define ARGS_BUF_SIZE 256

#define SIDE_END	2
#define SIDE_UNDEF	3

#define FLAG_NORMAL	0
#define FLAG_OPER	1
#define FLAG_COMP	2
#define FLAG_AND	3
typedef enum {OP_EQ, OP_ADD, OP_SUB, OP_MULT, OP_DIV, OP_NOT, OP_OR, OP_AND, OP_GRT, OP_LST, OP_GEQ, OP_LEQ, N_OPTYPES} Optype_e;

typedef unsigned int _uint;
typedef unsigned char _uint8;

typedef enum { E_SUCCESS, E_NOFILE, E_LACK_TOKENS, E_BAD_TOKEN, E_BAD_SYNTAX, E_BAD_VALUE, E_NOMEM, E_EMPTY_STACK, E_NOT_BINARY, E_NAN } parse_ercode;

typedef enum { CGS_UNION, CGS_INTERSECT, CGS_DIFFERENCE, CGS_CMB_NOOP } combine_type;
//note that ROOTS are a special type of COMPOSITES
typedef enum { CGS_UNDEF, CGS_ROOT, CGS_DATA, CGS_COMPOSITE, CGS_SPHERE, CGS_BOX, CGS_CYLINDER } object_type;

typedef Eigen::Matrix3d emat3;
typedef Eigen::Vector3d evec3;
typedef Eigen::Vector4d evec4;

class CompositeObject;

/**
  * Remove the whitespace surrounding a word
  */
inline char* CGS_trim_whitespace(char* str, size_t* len) {
    if (!str) return NULL;
    char* start = str;
    bool started = false;
    _uint last_non = 0;
    for (_uint i = 0; str[i] != 0; ++i) {
        if (str[i] != ' ' && str[i] != '\t' && str[i] != '\n') {
            last_non = i;
            if (!started) {
                start = str + i;
		started = true;
            }
        }
    }
    str[last_non+1] = 0;
    if (len) *len = last_non+1;
    return start;
}

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

typedef enum {VAL_UNDEF, VAL_STR, VAL_NUM, VAL_LIST, VAL_3VEC, VAL_MAT} valtype;
class Value;
union V {
    char* s;
    double x;
    Value* l;
    evec3* v;
    Eigen::MatrixXd* m;
};
class Value {
    friend class Scene;
#ifndef DEBUG_INFO
protected:
#else
public:
#endif
    valtype type;
    V val;
    size_t n_els = 1; //only applicable for string and list types
public:
    Value() { type = VAL_UNDEF;val.x = 0;n_els=0; }
    Value(const char* s);
    Value(std::string s);
    Value(double x) { type = VAL_NUM;val.x = x; }
    Value(int x) { type = VAL_NUM;val.x = (double)x; }
    Value(const Value* vs, size_t n_vs);
    Value(Eigen::MatrixXd m);
    Value(evec3 vec);
    void copy(const Value& o);
    void cleanup();
    ~Value();
    Value(const Value& o);//copy 
    Value(Value&& o);//move
    //Value& operator=(const Value& o);
    Value& operator=(Value& o);//assign
    Value& operator=(const Value& o);//assign
    bool operator==(std::string str);
    bool operator!=(std::string str);
    valtype get_type() { return type; }
    size_t size() { return n_els; }
    V get_val() { return val; }
    char* to_c_str();
    double to_float();
    Value cast_to(valtype type, parse_ercode& er) const;
};

typedef struct {
    char* name;
    Value args[ARGS_BUF_SIZE];
    char* arg_names[ARGS_BUF_SIZE];
    size_t n_args = 0;
} cgs_func;

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

/**
 * A helper class (and enumeration) that can track the context and scope of a curly brace block while parsing a file
 */
typedef enum {BLK_UNDEF, BLK_MISC, BLK_INVERT, BLK_TRANSFORM, BLK_DATA, BLK_ROOT, BLK_COMPOSITE, BLK_LITERAL, BLK_COMMENT, BLK_SQUARE, BLK_QUOTE, BLK_QUOTE_SING, BLK_PAREN, BLK_CURLY} block_type;
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

struct type_ind_pair {
public:
    block_type t;
    size_t i;
    type_ind_pair() { t = BLK_UNDEF;i = 0; }
    type_ind_pair(size_t ii) { t = BLK_UNDEF;i = ii; }
    type_ind_pair(block_type tt, size_t ii) { t = tt;i = ii; }
};

struct side_obj_pair {
public:
    _uint side;
    CompositeObject* obj;
    side_obj_pair(_uint p_side, CompositeObject* p_obj) { side = p_side;obj=p_obj; }
    side_obj_pair(const side_obj_pair& o) { side = o.side;obj = o.obj; }
    side_obj_pair(side_obj_pair&& o) { side = o.side;obj = o.obj;o.side = 0;o.obj = NULL; }
    side_obj_pair& operator=(const side_obj_pair& o) { side = o.side;obj = o.obj;return *this; }
};

/**
 * A helper class for scene which maintains two parallel stacks used when constructing binary trees. The first specifies the integer code for the side occupied and the second stores pointers to objects. Each index specified in the side array is a "tribit" with 0 specifying the left side, 1 the right and 2 the end end of the stack
 * WARNING: The stack only handles pointers to objects. The caller is responsible for managing the lifetime of each object placed on the stack.
 */
class ObjectStack : public CGS_Stack<side_obj_pair> {
public:
    //parse_ercode push(_uint side, CompositeObject* obj);
    parse_ercode emplace_obj(Object* obj, object_type p_type);
    parse_ercode pop(_uint* side, CompositeObject** obj);
    _uint look_side();
    CompositeObject* look_obj();
    CompositeObject* get_root();
};

class Scene {
private:
    //std::vector<Object*> objects;
    std::vector<CompositeObject*> roots;
    std::vector<CompositeObject*> data_objs;

    parse_ercode lookup_val(char* tok, Value& sto) const;
    //let users define constants
    std::unordered_map<std::string, Value> named_items;
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

    Value parse_value(char* tok, parse_ercode& er) const;
    parse_ercode parse_list(char* str, Value& sto) const;
    parse_ercode make_object(const cgs_func& f, Object** ptr, object_type* type, int p_invert) const;
    parse_ercode make_transformation(const cgs_func& f, emat3& res) const;
    parse_ercode parse_func(char* token, long open_par_ind, cgs_func& f, char** end) const;
};

#endif //CGS_H
