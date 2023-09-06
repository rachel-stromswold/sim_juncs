#define DEBUG_INFO 1

#ifndef CGS_READ_H
#define CGS_READ_H

#include <stdio.h>
#include <vector>
#include <string>
#include <cstring>
#include <math.h>
#include "geometry.hpp"

//hints for dynamic buffer sizes
#define BUF_SIZE 	1024
#define LINE_SIZE 	512
#define DEF_STACK_SIZE	16
#define ARGS_BUF_SIZE 	256
#define FUNC_BUF_SIZE 	256

//keywords
#define KEY_CLASS_LEN	5
#define KEY_FOR_LEN	3
#define KEY_DEF_LEN	3
#define KEY_IN_LEN	2
#define KEY_IF_LEN	2

//hash params
#define FNV_PRIME	16777619
#define FNV_OFFSET	2166136261
#define TABLE_BITS	8
#define TABLE_MASK(n)	(((u_int32_t)1<<(n))-1)
#define TABLE_SIZE(n)	((u_int32_t)1<<(n))

typedef matrix<3,3> mat3x3;

typedef unsigned int _uint;
typedef unsigned char _uint8;

typedef enum { E_SUCCESS, E_NOFILE, E_LACK_TOKENS, E_BAD_TOKEN, E_BAD_SYNTAX, E_BAD_VALUE, E_BAD_TYPE, E_NOMEM, E_EMPTY_STACK, E_NOT_BINARY, E_NAN, E_NOT_DEFINED, E_OUT_OF_RANGE, E_BAD_FLOAT } parse_ercode;
typedef enum {VAL_UNDEF, VAL_STR, VAL_NUM, VAL_ARRAY, VAL_LIST, VAL_3VEC, VAL_MAT, VAL_FUNC, VAL_INST} valtype;
typedef enum {BLK_UNDEF, BLK_MISC, BLK_INVERT, BLK_TRANSFORM, BLK_DATA, BLK_ROOT, BLK_COMPOSITE, BLK_FUNC_DEC, BLK_LITERAL, BLK_COMMENT, BLK_SQUARE, BLK_QUOTE, BLK_QUOTE_SING, BLK_PAREN, BLK_CURLY} block_type;

inline void* xrealloc(void* p, size_t nsize) {
    void* tmp = realloc(p, nsize);
    if (!tmp) {
	fprintf(stderr, "Insufficient memory to allocate block of size %lu!\n", nsize);
	exit(EXIT_FAILURE);
    }
    return tmp;
}

/** ======================================================== utility functions ======================================================== **/

bool is_char_sep(char c);
bool is_token(const char* str, size_t i, size_t len);
char* find_token_before(char* str, size_t i);
char* strchr_block(char* str, char c);
char* token_block(char* str, const char* comp);
size_t read_cgs_line(char** bufptr, size_t* n, FILE* fp, size_t* lineno);
char* CGS_trim_whitespace(char* str, size_t* len);

/** ============================ line_buffer ============================ **/

/*class cgs_error {
private:
    char* message;
    parse_ercode er;

public:
    cgs_error(cgs_error type, const char* message);
    cgs_error(const cgs_error& o);
    cgs_error(cgs_error&& o);
    cgs_error& operator=(cgs_error& o);
    char* msg() { return message; }
    parse_ercode type() { return er; }
};*/

//store a line number and an offset within that line to describe a position in a line_buffer
struct line_buffer_ind {
    size_t line;
    size_t off;
    line_buffer_ind() { line = 0;off = 0; }
    line_buffer_ind(long pl, long po) { line = pl;off = po; }
    //increase/decrease the line buffer by a specified amount while keeping line number the same
    friend line_buffer_ind operator+(const line_buffer_ind& lhs, const size_t& rhs);
    friend line_buffer_ind operator-(const line_buffer_ind& lhs, const size_t& rhs);
};

class line_buffer {
private:
    char** lines;
    size_t* line_sizes;
    size_t n_lines;
    //helper function for jmp_enclosed and get_enclosed. If the line at index 
#ifndef DEBUG_INFO
    int it_single(char** sto, char start_delim, char end_delim, line_buffer_ind* start, line_buffer_ind* end, int* pdepth, bool include_delimse, bool include_start) const;
public:
#else
public:
    int it_single(char** sto, char start_delim, char end_delim, line_buffer_ind* start, line_buffer_ind* end, int* pdepth, bool include_delimse, bool include_start) const;
#endif
    line_buffer();
    line_buffer(const char* p_fname);
    line_buffer(const char** p_lines, size_t pn_lines);
    line_buffer(const char* line, char sep, const char* ignore_blocks = "\"\"()[]{}");
    line_buffer(const line_buffer& o);
    line_buffer& operator=(line_buffer& o);
    line_buffer& operator=(line_buffer&& o);
    line_buffer(line_buffer&& o);
    ~line_buffer();
    void split(char split_delim);
    line_buffer_ind jmp_enclosed(line_buffer_ind start, char start_delim, char end_delim, bool include_delims=false) const;
    line_buffer get_enclosed(line_buffer_ind start, line_buffer_ind* end, char start_delim, char end_delim, bool include_delims=false, bool include_start=false) const;
    char* get_line(line_buffer_ind p) const;
    char* get_line(size_t i) const { line_buffer_ind p(i,0);return get_line(p); }
    size_t get_line_size(size_t i) const { if (i >= n_lines) { return 0; }return line_sizes[i]; }
    size_t get_n_lines() const { return n_lines; }
    char* flatten(char sep_char = 0) const;
    //increment or decrement the line_buffer_ind p and return whether the operation was successful
    bool inc(line_buffer_ind& p) const;
    bool dec(line_buffer_ind& p) const;
    char get(line_buffer_ind p) const;
};

/** ============================ stack ============================ **/

/**
 * A helper class (and enumeration) that can track the context and scope of a curly brace block while parsing a file
 */
template <typename T>
class stack {
protected:
    size_t stack_ptr;
    size_t buf_len;
    T* buf;
    parse_ercode grow(size_t new_size) {
	size_t old_size = buf_len;
	buf_len = new_size;

	T* tmp_buf = (T*)malloc(sizeof(T)*buf_len);
	//copy memory to new block if allocation was successful or return nomem error
	if (tmp_buf) {
	    for (size_t i = 0; i < stack_ptr; ++i) {
		new(tmp_buf+i) T(std::move(buf[i]));
	    }
	    free(buf);
	    buf = tmp_buf;
	    char* clear_buf = (char*)(tmp_buf+old_size);
	    if (buf_len > old_size) {
		size_t clear_len = sizeof(T)*(buf_len-old_size);
		for (size_t i = 0; i < clear_len; ++i) clear_buf[i] = 0;
	    }
	} else {
	    buf_len = old_size;
	    return E_NOMEM;
	}

	return E_SUCCESS;
    }
public:
    stack() {
	stack_ptr = 0;
	buf_len = DEF_STACK_SIZE;
	buf = (T*)calloc(buf_len, sizeof(T));
	//check that allocation was successful
	if (!buf) {
	    buf = NULL;
	    buf_len = 0;
	    stack_ptr = 0;
	}
    }
    ~stack<T>() {
	if (buf) {
	    for (_uint i = 0; i < stack_ptr; ++i) buf[i].~T();
	    free(buf);
	    buf = NULL;
	}
	stack_ptr = 0;
	buf_len = 0;
    }
    //swap
    void swap(stack<T>& o) {
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
    stack(const stack<T>& o) {
	stack_ptr = o.stack_ptr;
	//to save memory we'll only allocate however many entries the old object had
	buf_len = stack_ptr;
	//make sure that we set aside some space for the buffer
	if (buf_len == 0) buf_len = DEF_STACK_SIZE;
	buf = (T*)calloc(buf_len, sizeof(T));
	if (!buf) {
	    buf = NULL;
	    buf_len = 0;
	    stack_ptr = 0;
	}

	//copy the entries from the old object stack
	for (_uint i = 0; i < stack_ptr; ++i) buf[i] = o.buf[i];
    }
    //move constructor
    stack(stack<T>&& o) {
	stack_ptr = o.stack_ptr;
	buf_len = o.buf_len;
	buf = o.buf;
	//invalidate the object that just got moved
	o.stack_ptr = 0;
	o.buf_len = 0;
	o.buf = NULL;
    }
    //assignment operator
    stack<T>& operator=(stack<T> o) {
	swap(o);
	return *this;
    }
    /**
     * Push a side-object pair onto the stack.
     * returns: An error code if one occurred or E_SUCCESS if the operation was successful. It is possible for this function to fail if there wasn't sufficient space to push the object onto the stack.
     */
    parse_ercode push(T b) {
	//we might need to allocate more memory
	parse_ercode res = E_SUCCESS;
	if (stack_ptr == buf_len)
	    res = grow(2*buf_len);
	if (res == E_SUCCESS)
	    buf[stack_ptr++] = b;

	return res;
    }
    /**
     * Pop a side-object pair into the values pointed to by side and obj respectively. This function is guaranteed to succeed, but if the stack is empty then SIDE_END will be assigned to side and NULL will be assigned to obj. The caller should check that the returned values are valid.
     */
    parse_ercode pop(T* ptr) {
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
    void reset() { stack_ptr = 0; }
    /**
     * Check wheter an entry matching key is already somewhere in the stack
     */
    bool has(const T& key) {
	for (_uint i = 0; i < stack_ptr; ++i) {
	    if (buf[i] == key) return true;
	}
	return false;
    }
    /**
     * Check wheter the stack is empty
     */
    bool is_empty() {
	return stack_ptr == 0;
    }
    /**
     * Get the number of entries on the stack
     */
    size_t size() {
	return stack_ptr;
    }
    /**
     * Returns a copy of the object at the top of the stack. If the stack is empty then the object casted from 0 is returned.
     */
    T peek(size_t ind = 1) {
	if (stack_ptr >= ind && ind > 0)
	    return buf[stack_ptr - ind];
	return (T)0;
    }
};

class CompositeObject;

class value;
class context;
class user_func;
union V {
    char* s;
    double x;
    double* a;
    value* l;
    vec3* v;
    context* c;
    user_func* f;
    mat3x3* m;
};
struct value {
    valtype type;
    V val;
    size_t n_els; //only applicable for string and list types

    value() { val.x = 0;type=VAL_UNDEF;n_els=0; }
    bool operator==(const value& o) const;
    bool operator!=(const value& o) const;
    bool operator==(std::string str) const;
    bool operator!=(std::string str) const;
    valtype get_type() { return type; }
    size_t size() { return n_els; }
    V& get_val() { return val; }
    char* to_c_str();
    double to_float();
    int rep_string(char* sto, size_t n) const;
    value cast_to(valtype type, parse_ercode& er) const;
};
//check whether the value is of the specified type based on the type string
inline bool is_type(value v, valtype t) { return v.type == t; }
bool is_type(value v, const char* type);
void cleanup_val(value* o);
value copy_val(const value o);
void swap_val(value* a, value* b);

struct cgs_func {
    char* name;
    value args[ARGS_BUF_SIZE];
    char* arg_names[ARGS_BUF_SIZE];
    size_t n_args;
    cgs_func() { name = NULL;n_args = 0; }
};
value make_val_undef();
value make_val_num(double x);
value make_val_str(const char* s);
value make_val_std_str(std::string s);
value make_val_array(std::vector<double> a);
value make_val_list(const value* vs, size_t n_vs);
value make_val_mat(mat3x3 m);
value make_val_vec3(vec3 vec);
value make_val_func(const char* name, size_t n_args, value (*p_exec)(context&, cgs_func, parse_ercode&));
cgs_func parse_func_decl(char* str);
cgs_func copy_func(const cgs_func o);
void cleanup_func(cgs_func* o);
void swap(cgs_func* a, cgs_func* b);
value lookup_named(const cgs_func f, const char* name);

//collections (like python dicts) are simply aliases for functions under the hood
typedef context collection;

/**
 * A class which stores a labeled value.
 */
class name_val_pair {
private:
    char* name;
    value val;
public:
    void swap(name_val_pair& o);
    void copy(const name_val_pair& o);
    name_val_pair();
    name_val_pair(int a);
    name_val_pair(const char* p_name, value p_val);
    name_val_pair(const name_val_pair& o);
    name_val_pair(name_val_pair&& o);
    ~name_val_pair();
    name_val_pair& operator=(name_val_pair& o);
    name_val_pair& operator=(const name_val_pair& o);
    bool name_matches(const char* str) const;
    const char* get_name() { return name; }
    value& get_val();
};
struct type_ind_pair {
public:
    block_type t;
    size_t i;
    type_ind_pair() { t = BLK_UNDEF;i = 0; }
    type_ind_pair(size_t ii) { t = BLK_UNDEF;i = ii; }
    type_ind_pair(block_type tt, size_t ii) { t = tt;i = ii; }
};

/** ============================ context ============================ **/
struct name_ind {
    size_t ind;
    name_ind* next;
    name_ind(size_t i, name_ind* ptr=NULL) { ind = i;next=ptr; }
};
class context : public stack<name_val_pair> {
private:
    //members
    context* parent;
    name_ind* table[TABLE_SIZE(TABLE_BITS)];
    //helpers
    struct read_state {
	const line_buffer& b;
	line_buffer_ind pos;
	stack<block_type> blk_stk;
	size_t buf_size;
	char* buf;
	read_state(const line_buffer& pb) : b(pb), pos(0,0) {
	    buf_size = BUF_SIZE;
	    buf = (char*)calloc(buf_size, sizeof(char));
	}
	~read_state() { buf_size=0;free(buf); }
    };
    value do_op(char* tok, size_t ind, parse_ercode& er);
    parse_ercode read_single_line(context::read_state& rs);
    void init() {
	for(size_t i = 0; i < TABLE_SIZE(TABLE_BITS); ++i) { table[i] = NULL; }
    }
    void setup_builtins();

public:
    context() : stack<name_val_pair>() { init();setup_builtins();parent = NULL; }
    context(context* p_parent) : stack<name_val_pair>() { init();parent = p_parent; }
    ~context();
    //void emplace(const char* p_name, value p_val); { name_val_pair inst(p_name, p_val);push(inst); }
    value lookup(const char* name) const;
    parse_ercode pop(name_val_pair* ptr=NULL);
    parse_ercode pop_n(size_t n);
    value parse_value(char* tok, parse_ercode& er);
    value parse_value(const line_buffer& b, line_buffer_ind& pos, parse_ercode& er);
    cgs_func parse_func(char* token, long open_par_ind, parse_ercode& f, char** end, int name_only=0);
    value parse_list(char* str, parse_ercode& sto);
    void swap(stack<name_val_pair>& o) { stack<name_val_pair>::swap(o); }
    parse_ercode set_value(const char* name, value new_val, bool force_push=false);
    parse_ercode read_from_lines(const line_buffer& b);
    void register_func(cgs_func sig, value (*p_exec)(context&, cgs_func, parse_ercode&));
    value peek_val(size_t i=1);
    void error(const char* msg) {
	//TODO: something that isn't dumb
	printf("Error: %s\n", msg);
    }
};

/**
 * A class for functions defined by the user along with the implementation code
 */
class user_func {
private:
    cgs_func call_sig;
    line_buffer code_lines;
    value (*exec)(context&, cgs_func, parse_ercode&);
public:
    //read the function with contents stored in the file pointer fp at the current file position
    user_func(cgs_func sig, line_buffer p_buf);
    user_func(value (*p_exec)(context&, cgs_func, parse_ercode&));
    ~user_func();
    user_func(const user_func& o);
    user_func(user_func&& o);
    line_buffer& get_buffer() { return code_lines; }
    value eval(context& c, cgs_func call, parse_ercode& er);
};

#endif //CGS_READ_H
