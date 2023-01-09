#include "cgs_read.hpp"

/** ======================================================== utility functions ======================================================== **/

/**
 * Helper function for is_token which tests whether the character c is a token terminator
 */
bool is_char_sep(char c) {
    if (c == 0 || c == ' ' || c == '\t' || c == '\n' || c == ';' || c == '+'  || c == '-' || c == '*'  || c == '/')
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
 * Helper function that finds the start of first token before the index i in the string str. Note that this is not null terminated and includes characters including and after str[i] (unless str[i] = 0).
 */
char* find_token_before(char* str, size_t i) {
    while (i > 0) {
	--i;
	if (is_char_sep(str[i])) return str+i+1;
    }
    return str;
}

/**
 * Find the first index of the character c that isn't nested inside a block
 */
char* strchr_block(char* str, char c) {
    stack<size_t> blk_stk;
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
    stack<size_t> blk_stk;
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
	} else if (i > 0 && str[i-1] == '/' && str[i] == '*') {
	    blk_stk.push(i);
	} else if (i > 0 && str[i-1] == '*' && str[i] == '/') {
	    blk_stk.pop(NULL);
	} else if (i > 0 && str[i-1] == '/' && str[i] == '/') {
	    blk_stk.push(i);
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
	size = LINE_SIZE;
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
    stack<type_ind_pair> blk_stk;
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
  * Remove the whitespace surrounding a word
  * Note: this function performs trimming "in place"
  */
char* CGS_trim_whitespace(char* str, size_t* len) {
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

/** ============================ line_buffer ============================ **/

line_buffer::line_buffer(char* p_fname) {
    size_t buf_size = BUF_SIZE;
    lines = (char**)malloc(sizeof(char*)*buf_size);
    line_sizes = (size_t*)malloc(sizeof(size_t)*buf_size);
    n_lines = 0;

    FILE* fp = NULL;
    if (p_fname) {
        fp = fopen(p_fname, "r");
    }
    if (fp) {
	size_t lineno = 1;
	size_t line_len = 0;
	do {
	    //reallocate buffer if necessary
	    if (n_lines >= buf_size) {
		buf_size *= 2;
		lines = (char**)realloc(lines, sizeof(char*)*buf_size);
		line_sizes = (size_t*)realloc(line_sizes, sizeof(size_t)*buf_size);
	    }
	    //these need to be set to zero so that read_cgs_line() will allocate a buffer for us
	    char* this_buf = NULL;
	    line_len = read_cgs_line(&this_buf, &line_len, fp, &lineno);
	    lines[n_lines] = this_buf;
	    line_sizes[n_lines++] = line_len;
	} while (line_len >= 0);
	lines = (char**)realloc(lines, sizeof(char*)*n_lines);
	line_sizes = (size_t*)realloc(line_sizes, sizeof(size_t)*n_lines);
    } else {
        printf("Error: couldn't open file %s for reading!\n", p_fname);
	free(lines);
	lines = NULL;
    }
}
line_buffer::line_buffer(const char** p_lines, size_t pn_lines) {
    n_lines = pn_lines;
    lines = (char**)malloc(sizeof(char*)*n_lines);
    line_sizes = (size_t*)malloc(sizeof(size_t)*n_lines);
    for (size_t i = 0; i < n_lines; ++i) {
	line_sizes[i] = strlen(p_lines[i]);
	lines[i] = strdup(p_lines[i]);
    }
}
/**
 * Create a line buffer from a single line separated by characters of the type sep such that sep are not contained inside any blocks specified by ignore_blocks. Seperation characters are not included.
 * ignore_blocks: if specified then we search for instances of block pairs and make sure we are at the root level. ignore_blocks should have an even number of characters with even index characters specifying open block delimeters and odd index characters specifying close block delimeters.
 */
line_buffer::line_buffer(const char* str, char sep, const char* ignore_blocks) {
    n_lines = 0;
    if (!ignore_blocks) ignore_blocks = "";
    //by doing this we allocate once since we have a guaranteed upper limit on the number of lines that might exist
    size_t line_buf_size = strlen(str);
    lines = (char**)malloc(sizeof(char*)*line_buf_size);
    line_sizes = (size_t*)malloc(sizeof(size_t)*line_buf_size);
    int nest_level = 0;
    size_t last_sep = 0;
    size_t i = 0;
    for (;;++i) {
	//iterate through pairs of open/close block delimeters
	for (size_t j = 0; ignore_blocks[j] && ignore_blocks[j+1]; j += 2) {
	    if (ignore_blocks[j] == ignore_blocks[j+1]) {
		//there's a special case for things like quotations, skip ahead until we find the end
		++i;
		for (;str[i] && str[i] != ignore_blocks[j]; ++i) (void)0;
		break;
	    } else {
		if (str[i] == ignore_blocks[j])
		    ++nest_level;
		else if (str[i] == ignore_blocks[j+1])
		    --nest_level;
	    }
	}
	if (nest_level <= 0 && (str[i] == sep || str[i] == 0)) {
	    line_sizes[n_lines] = i-last_sep;
	    lines[n_lines] = (char*)malloc(sizeof(char)*(line_sizes[n_lines]+1));
	    for (size_t j = last_sep; j < i; ++j)
		lines[n_lines][j] = str[j];
	    lines[n_lines][line_sizes[n_lines]] = 0;
	    last_sep = i+1;
	    ++n_lines;
	    if (str[i] == 0) break;
	}
    }
}
line_buffer::~line_buffer() {
    for (size_t i = 0; i < n_lines; ++i) free(lines[i]);
    free(line_sizes);
    free(lines);
}
line_buffer::line_buffer(const line_buffer& o) {
    n_lines = o.n_lines;
    line_sizes = (size_t*)malloc(sizeof(size_t)*n_lines);
    lines = (char**)malloc(sizeof(char*)*n_lines);
    for (size_t i = 0; i < n_lines; ++i) {
	lines[i] = strdup(o.lines[i]);
	line_sizes[i] = o.line_sizes[i];
    }
}
line_buffer::line_buffer(line_buffer&& o) {
    n_lines = o.n_lines;
    line_sizes = o.line_sizes;
    lines = o.lines;
    o.n_lines = 0;
    o.lines = NULL;
}
/**
 * Find the line buffer starting on line start_line between the first instance of start_delim and the last instance of end_delim respecting nesting (i.e. if lines={"a {", "b {", "}", "} c"} then {"b {", "}"} is returned. Note that the result must be deallocated with a call to free().
 * start_line: the line to start reading from
 * end_line: if this value is not NULL, then the index of the line on which end_delim was found is stored here
 * start_delim: the character to be used as the starting delimiter. This needs to be supplied so that we find a matching end_delim at the root level
 * end_delim: the character to be searched for
 * line_offset: the character on line start_line to start reading from, this defaults to zero. Note that this only applies to the start_line, and not any subsequent lines in the buffer.
 * include_delims: if true then the delimeters are included in the enclosed strings. defualts to false
 * returns: an array of lines that must be deallocated by a call to free(). The number of lines is stored in n_lines.
 */
line_buffer line_buffer::get_enclosed(size_t start_line, long* end_line, char start_delim, char end_delim, size_t line_offset, bool include_delims) {
    //initialization
    line_buffer ret;
    if (end_line) *end_line = -1;
    if (n_lines == 0) {
	ret.n_lines = 0;
	ret.lines = NULL;
	ret.line_sizes = NULL;
	return ret;
    }
    ret.lines = (char**)malloc(sizeof(char*)*n_lines);
    ret.line_sizes = (size_t*)malloc(sizeof(size_t)*n_lines);

    //tracking variables
    int depth = 0;
    bool exit = false;

    //iterate through lines
    size_t k = 0;
    long start_char = -1;
    long j = line_offset;
    for (size_t i = start_line; !exit && i < n_lines; ++i) {
	if (lines[i] == NULL) {
	    ret.n_lines = i;
	    exit = true;
	} else {
	    //this way we ensure that start char is always zero after the first start_delim is read
	    if (start_char >= 0) start_char = 0;
	    //iterate through characters in the line looking for an end_delim without a preceeding start_delim
	    for (; lines[i][j]; ++j) {
		if (lines[i][j] == start_delim) {
		    if (depth == 0) {
			start_char = j;
			if (!include_delims) ++start_char;
		    }
		    ++depth;
		} else if (lines[i][j] == end_delim) {
		    --depth;
		    if (depth <= 0) {
			ret.n_lines = k+1;
			exit = true;
			if (include_delims) ++j;
			if (end_line) *end_line = i;
			break;
		    }
		}
	    }
	    //don't read empty lines
	    if (start_char >= 0) {
		ret.line_sizes[k] = j;
		ret.lines[k++] = strndup(lines[i]+start_char, j-start_char);
	    }
	    j = 0;
	}
    }
    return ret;
}

/**
 * Returns a version of the line buffer which is flattened so that everything fits onto one line.
 * sep_char: if this is not zero, then each newline in the buffer is replaced by a sep_char
 */
char* line_buffer::flatten(char sep_char) {
    //figure out how much memory must be allocated
    size_t tot_size = 1;
    for (size_t i = 0; i < n_lines; ++i) { tot_size += line_sizes[i]; }
    if (sep_char != 0) tot_size += n_lines;
    //allocate the memory
    char* ret = (char*)malloc(sizeof(char)*tot_size);
    size_t k = 0;
    for (size_t i = 0; i < n_lines; ++i) {
	for (size_t j = 0; lines[i][j] && j < line_sizes[i]; ++j) {
	    ret[k++] = lines[i][j];
	}
	if (k != 0 && sep_char != 0) ret[k++] = sep_char;
    }
    ret[k] = 0;
    return ret;
}

/** ======================================================== builtin functions ======================================================== **/
/**
 * Make a vector argument with the x,y, and z coordinates supplied
 */
value make_vec(cgs_func tmp_f, parse_ercode& er) {
    value sto;sto.type = VAL_UNDEF;
    if (tmp_f.n_args < 3) { er = E_LACK_TOKENS;return sto; }
    if (tmp_f.args[0].type != VAL_NUM || tmp_f.args[1].type != VAL_NUM || tmp_f.args[2].type != VAL_NUM) { er = E_BAD_TOKEN;return sto; }
    sto.type = VAL_3VEC;
    sto.val.v = new vec3(tmp_f.args[0].val.x, tmp_f.args[1].val.x, tmp_f.args[2].val.x);
    sto.n_els = 3;
    er = E_SUCCESS;
    return sto;
}
/**
 * Make a range following python syntax. If one argument is supplied then a list with tmp_f.args[0] elements is created starting at index 0 and going up to (but not including) tmp_f.args[0]. If two arguments are supplied then the range is from (tmp_f.args[0], tmp_f.args[1]). If three arguments are supplied then the range (tmp_f.args[0], tmp_f.args[1]) is still returned, but now the spacing between successive elements is tmp_f.args[2].
 */
value make_range(cgs_func tmp_f, parse_ercode& er) {
    value sto;sto.type = VAL_UNDEF;
    if (tmp_f.n_args == 0) {
	er = E_LACK_TOKENS;
	return sto;
    } else if (tmp_f.n_args == 1) {
	if (tmp_f.args[0].type != VAL_NUM) { er = E_BAD_TYPE;return sto; }
	if (tmp_f.args[0].val.x < 0) { er = E_BAD_VALUE;return sto; }
	//interpret the argument as an upper bound starting from 0
	sto.type = VAL_LIST;
	sto.n_els = (size_t)(tmp_f.args[0].val.x);
	sto.val.l = (value*)malloc(sizeof(value)*sto.n_els);
	for (size_t i = 0; i < sto.n_els; ++i) {
	    sto.val.l[i].type = VAL_NUM;
	    sto.val.l[i].val.x = i;
	}
    } else if (tmp_f.n_args == 2) {
	if (tmp_f.args[0].type != VAL_NUM) { er = E_BAD_TYPE;return sto; }
	if (tmp_f.args[1].type != VAL_NUM) { er = E_BAD_TYPE;return sto; }
	//interpret the argument as an upper bound starting from 0
	sto.type = VAL_LIST;
	int list_start = (int)(tmp_f.args[0].val.x);
	int list_end = (int)(tmp_f.args[1].val.x);
	if (list_end < list_start) { er = E_BAD_VALUE;return sto; }
	sto.n_els = list_end - list_start;
	sto.val.l = (value*)malloc(sizeof(value)*sto.n_els);
	size_t j = 0;
	for (int i = list_start; i < list_end; ++i) {
	    sto.val.l[j].type = VAL_NUM;
	    sto.val.l[j].val.x = i;
	    ++j;
	}
    } else {
	if (tmp_f.args[0].type != VAL_NUM) { er = E_BAD_TYPE;return sto; }
	if (tmp_f.args[1].type != VAL_NUM) { er = E_BAD_TYPE;return sto; }
	if (tmp_f.args[2].type != VAL_NUM) { er = E_BAD_TYPE;return sto; }
	//interpret the argument as an upper bound starting from 0
	sto.type = VAL_LIST;
	double list_start = tmp_f.args[0].val.x;
	double list_end = tmp_f.args[1].val.x;
	double inc = tmp_f.args[2].val.x;
	if (list_end < list_start || inc == 0) { er = E_BAD_VALUE;return sto; }
	sto.n_els = (list_end - list_start) / inc;
	sto.val.l = (value*)malloc(sizeof(value)*sto.n_els);
	for (size_t i = 0; i < sto.n_els; ++i) {
	    sto.val.l[i].type = VAL_NUM;
	    sto.val.l[i].val.x = i*inc + list_start;
	}
    }
    er = E_SUCCESS;
    return sto;
}
/**
 * linspace(a, b, n) Create a list of n equally spaced real numbers starting at a and ending at b. This function must be called with three aguments unlike np.linspace. Note that the value b is included in the list
 */
value make_linspace(cgs_func tmp_f, parse_ercode& er) {
    value sto;
    if (tmp_f.n_args < 3) {
        er = E_LACK_TOKENS;
        return sto;
    } else {
        if (tmp_f.args[0].type == VAL_NUM && tmp_f.args[1].type == VAL_NUM && tmp_f.args[2].type == VAL_NUM) {
            sto.type = VAL_LIST;
            sto.n_els = (size_t)(tmp_f.args[2].val.x);
            //prevent divisions by zero
            if (sto.n_els < 2) {
                sto.type = VAL_UNDEF;
                sto.n_els = 0;
                sto.val.l = NULL;
                er = E_BAD_VALUE;
                return sto;
            }
            sto.val.l = (value*)malloc(sizeof(value)*sto.n_els);
            double step = (tmp_f.args[1].val.x - tmp_f.args[0].val.x)/(sto.n_els - 1);
            for (size_t i = 0; i < sto.n_els; ++i) {
                sto.val.l[i].type = VAL_NUM;
                sto.val.l[i].val.x = step*i + tmp_f.args[0].val.x;
            }
            er = E_SUCCESS;
            return sto;
        } else {
            er = E_BAD_TYPE;
            return sto;
        }    
    }
}
/**
 * Take a list value and flatten it so that it has numpy dimensions (n) where n is the sum of the length of each list in the base list. values are copied in order e.g flatten([0,1],[2,3]) -> [0,1,2,3]
 * cgs_func: the function with arguments passed
 */
value flatten_list(cgs_func tmp_f, parse_ercode& er) {
    value sto;sto.type = VAL_UNDEF;
    if (tmp_f.n_args < 1) { er = E_LACK_TOKENS;return sto; }
    if (tmp_f.args[0].type != VAL_LIST) { er = E_BAD_VALUE;return sto; }
    value cur_list = tmp_f.args[0];
    size_t cur_st = 0;
    //this is used for estimating the size of the buffer we need. Take however many elements were needed for this list and assume each sub-list has the same number of elements
    size_t base_n_els = cur_list.n_els;
    //start with the number of elements in the lowest order of the list
    size_t buf_size = cur_list.n_els;
    sto.val.l = (value*)malloc(sizeof(value)*buf_size);
    //there may potentially be nested lists, we need to be able to find our way back to the parent and the index once we're done
    stack<value> lists;
    stack<size_t> inds;
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
		value* tmp_val = (value*)realloc(sto.val.l, sizeof(value)*buf_size);
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
 * print the elements to the console
 */
value print(cgs_func tmp_f, parse_ercode& er) {
    value ret;ret.type = VAL_UNDEF;ret.val.x = 0;ret.n_els = 0;
    for (size_t i = 0; i < tmp_f.n_args; ++i) {
	if (tmp_f.args[i].type == VAL_NUM) {
	    printf("%f", tmp_f.args[i].val.x);
	} else if (tmp_f.args[i].type == VAL_STR) {
	    printf("%s", tmp_f.args[i].val.s);
	} else if (tmp_f.args[i].type == VAL_3VEC) {
	    printf("(%f,%f,%f)", tmp_f.args[i].val.v->x(), tmp_f.args[i].val.v->y(), tmp_f.args[i].val.v->z());
	} else {
	    printf("<object at %p>", tmp_f.args[i].val.l);
	}
    }
    printf("\n");
    return ret;
}
/** ======================================================== value ======================================================== **/
value make_val_undef() {
    value v;
    v.type = VAL_UNDEF;
    v.n_els = 0;
    v.val.x = 0;
    return v;
}
value make_val_num(double x) {
    value v;
    v.type = VAL_NUM;
    v.n_els = 1;
    v.val.x = x;
    return v;
}
value make_val_str(const char* s) {
    value v;
    v.type = VAL_STR;
    v.n_els = strlen(s) + 1;
    v.val.s = (char*)malloc(sizeof(char)*v.n_els);
    for (size_t i = 0; i < v.n_els; ++i) v.val.s[i] = s[i];
    return v;
}
value make_val_std_str(std::string s) {
    value v;
    v.type = VAL_STR;
    v.n_els = s.size();
    v.val.s = (char*)malloc(sizeof(char)*v.n_els);
    for (size_t i = 0; i < v.n_els; ++i) v.val.s[i] = s[i];
    return v;
}
value make_val_list(const value* vs, size_t n_vs) {
    value v;
    v.type = VAL_LIST;
    v.n_els = n_vs;
    v.val.l = (value*)malloc(sizeof(value)*v.n_els);
    for (size_t i = 0; i < v.n_els; ++i) v.val.l[i] = vs[i];
    return v;
}
value make_val_mat(mat3x3 m) {
    value v;
    v.type = VAL_MAT;
    v.n_els = 1;
    v.val.m = new mat3x3(m);
    return v;
}
value make_val_vec3(vec3 vec) {
    value v;
    v.type = VAL_3VEC;
    v.n_els = 1;
    v.val.v = new vec3(vec.el[0], vec.el[1], vec.el[2]);
    return v;
}
/**
 * Add a new callable function with the signature sig and function pointer corresponding to the executed code. This function must accept a function and a pointer to an error code and return a value.
 */
value make_val_func(const char* name, size_t n_args, value (*p_exec)(context&, cgs_func, parse_ercode&)) {
    /*cgs_func sig;
    sig.name = strdup(name);
    sig.n_args = n_args;
    for (size_t i = 0; i < n_args; ++i) {
	sig.args[i] = make_val_num(0);
	sig.arg_names[i] = NULL;
    }*/

    value ret;
    ret.type = VAL_FUNC;
    ret.n_els = n_args;
    ret.val.f = new user_func(p_exec);
    return ret;
}
void cleanup_val(value* v) {
    if (v->type == VAL_STR && v->val.s) {
	free(v->val.s);
    } else if (v->type == VAL_LIST && v->val.l) {
	for (size_t i = 0; i < v->n_els; ++i) cleanup_val(v->val.l + i);
	free(v->val.l);
    } else if (v->type == VAL_MAT && v->val.m) {
	delete v->val.m;
    } else if (v->type == VAL_3VEC && v->val.v) {
	delete v->val.v;
    } else if (v->type == VAL_INST && v->val.c) {
	delete v->val.c;
    } else if (v->type == VAL_FUNC && v->val.f) {
	delete v->val.f;
    }
    v->val.x = 0;
}
value copy_val(const value o) {
    value ret;
    ret.type = o.type;
    ret.n_els = o.n_els;
    //strings or lists must be copied
    if (o.type == VAL_STR) {
	ret.val.s = (char*)malloc(sizeof(char)*o.n_els);
	for (size_t i = 0; i < o.n_els; ++i) ret.val.s[i] = o.val.s[i];
    } else if (o.type == VAL_LIST) {
	ret.val.l = (value*)calloc(o.n_els, sizeof(value));
	for (size_t i = 0; i < o.n_els; ++i) ret.val.l[i] = copy_val(o.val.l[i]);
    } else if (o.type == VAL_MAT) {
	ret.val.m = new mat3x3(*(o.val.m));
    } else if (o.type == VAL_3VEC) {
	ret.val.v = new vec3(o.val.v->el[0], o.val.v->el[1], o.val.v->el[2]);
    } else if (o.type == VAL_INST) {
	ret.val.c = new collection(*(o.val.c));
    } else if (o.type == VAL_FUNC) {
	ret.val.f = new user_func(*(o.val.f));
    } else {
	ret.val.x = o.val.x;
    }
    return ret;
}
void swap_val(value* a, value* b) {
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
bool value::operator==(std::string str) {
    if (type != VAL_STR) return false;
    size_t strlen = str.size();
    if (strlen != n_els-1) return false;
    for (size_t i = 0; i < strlen; ++i) {
	if (str[i] != val.s[i]) return false;
    }
    return true;
}
bool value::operator!=(std::string str) {
    if (type != VAL_STR) return true;
    size_t strlen = str.size();
    if (strlen != n_els-1) return true;
    for (size_t i = 0; i < strlen; ++i) {
	if (str[i] != val.s[i]) return true;
    }
    return false;
}
char* value::to_c_str() {
    if (type == VAL_STR)
	return val.s;
    return NULL;
}
double value::to_float() {
    if (type == VAL_NUM)
	return val.x;
    return 0;
}
/**
 * Convert a value to a string representation. The string is saved to sto.
 * v: the value to convert to a string
 * sto: the string to write to
 * n: the size of the buffer sto
 * returns: the number of characters written excluding the null terminator
 */
int value::rep_string(char* sto, size_t n) const {
	if (type == VAL_STR) {
	    char* end = stpncpy(sto, val.s, n-1);
	    *end = 0;
	    return end-sto;
	} else if (type == VAL_NUM) {
	    return snprintf(sto, n, "%f", val.x);
	} else if (type == VAL_LIST) {
	    //TODO
	} else if (type == VAL_3VEC) {
	    return val.v->to_str(sto, n);
	} else if (type == VAL_MAT) {
	    return val.m->to_str(sto, n);
	}
    return 0;
}
/**
 * Perform an in-place cast of the instance to the type t. An error is returned if a cast is impossible.
 */
value value::cast_to(valtype t, parse_ercode& er) const {
    er = E_SUCCESS;
    if (type == VAL_UNDEF) { er = E_BAD_VALUE;return copy_val(*this); }
    if (type == t) return copy_val(*this);
    value ret;
    ret.type = t;
    ret.n_els = n_els;
    if (t == VAL_3VEC) {
	if (type == VAL_LIST) {
	    value* tmp_lst = val.l;
	    //check that we have at least three numeric values
	    if (n_els < 3) { er = E_LACK_TOKENS;return copy_val(*this); }
	    if (tmp_lst[0].type != VAL_NUM || tmp_lst[1].type != VAL_NUM || tmp_lst[2].type != VAL_NUM) { er = E_BAD_TOKEN;return copy_val(*this); }
	    //actually change data and free the old
	    ret.val.v = new vec3(tmp_lst[0].val.x, tmp_lst[1].val.x, tmp_lst[2].val.x);
	    ret.n_els = 3;
	    return ret;
	}
    } else if (t == VAL_LIST) {
	if (type == VAL_3VEC) {
	    ret.n_els = 3;
	    ret.val.l = (value*)malloc(sizeof(value)*ret.n_els);
	    for (size_t i = 0; i < ret.n_els; ++i) {
		ret.val.l[i].type = VAL_NUM;
		ret.val.l[i].n_els = 1;
		ret.val.l[i].val.x = val.v->el[i];
	    }
	    return ret;
	} else if (type == VAL_MAT) {
	    //TODO
	}
    } else if (t == VAL_STR) {
	ret.val.s = (char*)malloc(sizeof(char)*BUF_SIZE);
	int n_write = rep_string(ret.val.s, BUF_SIZE);
	ret.val.s = (char*)realloc(ret.val.s, sizeof(char)*(n_write+1));
    }
    //if we reach this point in execution then there was an error
    ret.type = VAL_UNDEF;
    ret.n_els = 0;
    return ret;
}

/** ======================================================== cgs_func ======================================================== **/

void cleanup_func(cgs_func* f) {
    if (f) {
	for (size_t i = 0; i < f->n_args; ++i) {
	    cleanup_val(f->args + i);
	    if (f->arg_names[i]) free(f->arg_names[i]);
	}
    }
}
cgs_func copy_func(const cgs_func o) {
    cgs_func f;
    if (o.name) f.name = strdup(o.name);
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
	value tmp_val = sf->args[i];
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
/**
 * Find the named argument with the matching name
 */
value lookup_named(const cgs_func f, const char* name) {
    for (size_t i = 0; i < f.n_args; ++i) {
	if (f.arg_names[i] && strcmp(f.arg_names[i], name) == 0) return f.args[i];
    }
    value ret;
    ret.type = VAL_UNDEF;
    ret.val.x = 0;
    ret.n_els = 0;
    return ret;
}

/** ======================================================== name_val_pair ======================================================== **/

void name_val_pair::swap(name_val_pair& o) {
    char* tmp_name = name;
    name = o.name;
    o.name = tmp_name;
    value tmp_val = val;
    val = o.val;
    o.val = tmp_val;
}

/**
 * Copy the name_val pair into *this
 */
void name_val_pair::copy(const name_val_pair& o) {
    if (o.name)
	name = strdup(o.name);
    else
	name = NULL;
    val = copy_val(o.val);
}
/**
 * constructors
 */
name_val_pair::name_val_pair() {
    name = NULL;
    val.type = VAL_UNDEF;
    val.val.x = 0;
    val.n_els = 0;
}
name_val_pair::name_val_pair(int a) {
    name = NULL;
    val.type = VAL_UNDEF;
    val.val.x = 0;
    val.n_els = 0;
}
name_val_pair::name_val_pair(const char* p_name, value p_val) {
    if (p_name)
	name = strdup(p_name);
    else
	name = NULL;
    val = copy_val(p_val);
}
name_val_pair::name_val_pair(const name_val_pair& o) {
    copy(o);
}
name_val_pair::name_val_pair(name_val_pair&& o) {
    name = o.name;
    val = o.val;
    o.name = NULL;
    o.val.type = VAL_UNDEF;
    o.val.val.x = 0;
}
name_val_pair::~name_val_pair() {
	if (name) free(name);
	if (val.type != VAL_UNDEF) cleanup_val(&val);
    }
name_val_pair& name_val_pair::operator=(name_val_pair& o) {
    swap(o);
    return *this;
}
name_val_pair& name_val_pair::operator=(const name_val_pair& o) {
    copy(o);
    return *this;
}
/**
 * test whether the name of this pair matches the search string
 */
bool name_val_pair::name_matches(const char* str) const {
    if (!name) return false;
    if (strcmp(str, name) == 0) return true;
    return false;
}
/**
 * Fetch a reference to the value in this pair
 */
value& name_val_pair::get_val() {
    return val;
}

/** ======================================================== cgs_func ======================================================== **/

/**
 * Given the string starting at token, and the index of an open paren parse the result into a cgs_func struct.
 * token: a c-string which is modified in place that contains the function
 * open_par_ind: the location of the open parenthesis
 * f: the cgs_func that information should be saved to
 * end: If not NULL, a pointer to the first character after the end of the string is stored here. If an error occurred during parsing end will be set to NULL.
 * name_only: if set to true, then only the names of the function arguments will be set and all values will be set to undefined. This is useful for handling function declarations
 * returns: an errorcode if an invalid string was supplied.
 */
cgs_func context::parse_func(char* token, long open_par_ind, parse_ercode& er, char** end, int name_only) {
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
    if (name_only) {
	for (size_t i = 0; list_els[i] && i < f.n_args; ++i) {
	    f.arg_names[i] = strdup(list_els[i]);
	    f.args[i].type = VAL_UNDEF;
	    f.args[i].n_els = 0;
	    f.args[i].val.x = 0;
	}
    } else {
	for (size_t i = 0; list_els[i] && i < f.n_args; ++i) {
	    //handle named arguments
	    char* eq_loc = strchr_block(list_els[i], '=');
	    if (eq_loc) {
		*eq_loc = 0;
		f.arg_names[i] = strdup(list_els[i]);
		list_els[i] = eq_loc+1;
	    } else {
		f.arg_names[i] = NULL;
	    }
	    f.args[i] = parse_value(list_els[i], er);
	    //if the evaluation of the name failed, then use it as a variable name instead
	    /*if (er != E_SUCCESS) {
		f.args[i].type = VAL_UNDEF;
		f.arg_names[i] = list_els[i];
		er = E_SUCCESS;
	    }*/
	}
    }
    //cleanup and reset the string
    //*term_ptr = /*(*/')';
    free(list_els);
    return f;
}

/** ======================================================== context ======================================================== **/

/**
 * Iterate through the context and return a value to the variable with the matching name.
 * name: the name of the variable to set
 * returns: the matching value, no deep copies are performed
 */
value context::lookup(const char* str) const {
    value ret;
    ret.type = VAL_UNDEF;
    ret.val.x = 0;
    if (stack_ptr == 0) return ret;
    //iterate to find the item highest on the stack with a matching name
    for (long i = stack_ptr-1; i >= 0; --i) {
	if (buf[i].name_matches(str)) return buf[i].get_val();
    }
    //return a default undefined value if nothing was found
    return ret;
}

/**
 * Iterate through the context and assign a new value to the variable with the matching name. NOTE: this function only performs a shallow copy.
 * name: the name of the variable to set
 * new_val: the value to set the variable to
 * returns: E_SUCCESS if the variable with a matching name was found or E_NOT_DEFINED otherwise
 */
parse_ercode context::set_value(const char* str, value new_val) {
    for (long i = stack_ptr-1; i >= 0; --i) {
	if (buf[i].name_matches(str)) {
	    buf[i].get_val() = copy_val(new_val);
	    return E_SUCCESS;
	}
    }
    //if we didn't find the value with the matching name then define it
    emplace(str, new_val);
    return E_SUCCESS;
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
 * Read a string of the format [x, y, z] into an Eigen::Vector3.
 * returns: 0 on success or an error code
 * 	-1: insufficient tokens
 * 	-2: one of the tokens supplied was invalid
 */
value context::parse_list(char* str, parse_ercode& er) {
    value sto;
    //find the start and the end of the vector
    char* start = strchr(str, '[');
    if (!start) { er = E_BAD_SYNTAX;return sto; }
    //make sure that the string is null terminated
    char* end = strchr_block(start+1, ']');
    if (!end) { er = E_BAD_SYNTAX;return sto; }
    *end = 0;

    //read the coordinates separated by spaces
    size_t n_els = 0;
    value* lbuf;
    //check if this is a list interpretation
    char* for_start = token_block(start+1, "for");
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
	value it_list = parse_value(list_expr, er);
	free(list_expr);
	if (er != E_SUCCESS || it_list.type != VAL_LIST) { er = E_BAD_VALUE;return sto; }
	//we now iterate through the list specified, substituting VAL in the expression with the current value
	for_start[0] = 0;
	emplace(var_name, sto);
	n_els = it_list.n_els;
	lbuf = (value*)calloc(n_els, sizeof(value));
	for (size_t i = 0; i < it_list.n_els; ++i) {
	    set_value(var_name, it_list.val.l[i]);
	    char* expr_name = strndup(start+1, expr_len);
	    lbuf[i] = parse_value(expr_name, er);
	    free(expr_name);
	    if (er != E_SUCCESS) { n_els = i;return sto; }
	}
	name_val_pair tmp;
	pop(&tmp);
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
	lbuf = (value*)calloc(n_els, sizeof(value));
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
 * Execute the mathematical operation in the string str at the location op_ind
 */
value context::do_op(char* str, size_t i, parse_ercode& er) {
    value sto;
    sto.type = VAL_UNDEF;
    sto.val.x = 0;
    sto.n_els = 0;
    //some operators (==, >=, <=) take up more than one character, test for these
    int op_width = 1;
    if (str[i+1] == '=') op_width = 2;
    //Store the operation before setting it to zero
    char term_char = str[i];
    str[i] = 0;

    //ternary operators and dereferences are special cases
    if (term_char == '?') {
	//the colon must be present
	char* col_ind = strchr_block(str+i+1, ':');
	if (!col_ind) { er = E_BAD_SYNTAX;return sto; }
	value tmp_l = parse_value(str, er);
	if (er != E_SUCCESS) return sto;
	//false branch
	if (tmp_l.type == VAL_UNDEF || tmp_l.val.x == 0) {
	    sto = parse_value(col_ind+1, er);
	    return sto;
	} else {
	    //true branch
	    *col_ind = 0;
	    sto = parse_value(str+i+op_width, er);
	    *col_ind = ':';
	    return sto;
	}
    } else if (term_char == '.') {
	value inst_val = lookup(find_token_before(str, i));
	if (inst_val.type != VAL_INST) { er = E_BAD_TYPE;return sto; }
	str[i] = term_char;
	return inst_val.val.c->parse_value(str+i+1, er);
    }

    //parse right and left values
    value tmp_l = parse_value(str, er);
    if (er != E_SUCCESS) { cleanup_val(&tmp_l);return sto; }
    value tmp_r = parse_value(str+i+op_width, er);
    if (er != E_SUCCESS) { cleanup_val(&tmp_r);return sto; }
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
	    sto.val.x = (tmp_l.val.x == tmp_r.val.x)? 1: 0;
	}
	++i;
    } else if (term_char == '>') {
	if (tmp_l.type != VAL_NUM || tmp_r.type != VAL_NUM) { cleanup_val(&tmp_l);cleanup_val(&tmp_r);er = E_BAD_VALUE;return sto; }
	if (str[i+1] == '=') {
	    sto.val.x = (tmp_l.val.x >= tmp_r.val.x)? 1: 0;
	    ++i;
	} else {
	    sto.val.x = (tmp_l.val.x > tmp_r.val.x)? 1: 0;
	}
    } else if (term_char == '<') {
	if (tmp_l.type != VAL_NUM || tmp_r.type != VAL_NUM) { cleanup_val(&tmp_l);cleanup_val(&tmp_r);er = E_BAD_VALUE;return sto; }
	if (str[i+1] == '=') {
	    sto.val.x = (tmp_l.val.x <= tmp_r.val.x)? 1: 0;
	    ++i;
	} else {
	    sto.val.x = (tmp_l.val.x < tmp_r.val.x)? 1: 0;
	}
    //handle arithmetic
    } else if (term_char == '+' && (i == 0 || str[i-1] != 'e')) {
        if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_NUM) {	
            sto.val.x = tmp_l.val.x + tmp_r.val.x;
        } else if (tmp_l.type == VAL_MAT && tmp_r.type == VAL_MAT) {	
            sto.type = VAL_MAT;
            sto.val.m = new mat3x3(*tmp_l.val.m + *tmp_r.val.m);
        } else if (tmp_l.type == VAL_STR || tmp_r.type == VAL_STR) {
            sto.type = VAL_STR;
            sto.val.s = (char*)malloc(sizeof(char)*BUF_SIZE);
            int n_write = tmp_l.rep_string(sto.val.s, BUF_SIZE);
            if (n_write > 0 && n_write < BUF_SIZE)
                n_write += tmp_r.rep_string(sto.val.s+n_write, BUF_SIZE-n_write);
	    sto.n_els = n_write+1;
        }
    } else if (term_char == '-') {
	if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_NUM) {
	    sto.val.x = tmp_l.val.x - tmp_r.val.x;
	} else if (tmp_l.type == VAL_MAT && tmp_r.type == VAL_MAT) {	
	    sto.type = VAL_MAT;
	    sto.val.m = new mat3x3(*tmp_l.val.m - *tmp_r.val.m);
	}
    } else if (term_char == '*') {
	if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_NUM) {	
	    sto.val.x = tmp_l.val.x * tmp_r.val.x;
	} else if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_MAT) {	
	    sto.type = VAL_MAT;
	    sto.val.m = new mat3x3((*tmp_r.val.m) * tmp_l.val.x);
	} else if (tmp_l.type == VAL_MAT && tmp_r.type == VAL_MAT) {	
	    sto.type = VAL_MAT;
	    sto.val.m = new mat3x3((*tmp_l.val.m) * (*tmp_r.val.m));
	}
    } else if (term_char == '/') {
	if (tmp_r.val.x == 0) { cleanup_val(&tmp_l);cleanup_val(&tmp_r);er = E_NAN;return sto; }//TODO: return a nan?
	if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_NUM) {
	    sto.val.x = tmp_l.val.x / tmp_r.val.x;
	} else if (tmp_r.type == VAL_NUM && tmp_l.type == VAL_MAT) {	
	    sto.type = VAL_MAT;
	    sto.val.m = new mat3x3(*tmp_r.val.m / tmp_r.val.x);
	}
    } else if (term_char == '^') {
	if (tmp_r.val.x == 0 || tmp_l.type != VAL_NUM || tmp_r.type != VAL_NUM) { er = E_NAN;return sto; }//TODO: return a nan?
	sto.type = VAL_NUM;
	sto.val.x = pow(tmp_l.val.x, tmp_r.val.x);
    }
    str[i] = term_char;
    cleanup_val(&tmp_l);
    cleanup_val(&tmp_r);
    return sto;
}

value context::parse_value(char* str, parse_ercode& er) {
    value sto;
    er = E_SUCCESS;

    //store locations of the first instance of different operators. We do this so we can quickly look up new operators if we didn't find any other operators of a lower precedence (such operators are placed in the tree first).
    int first_open_ind = -1;
    int last_close_ind = -1;

    //keeps track of open and close [], (), {}, and ""
    stack<type_ind_pair> blk_stk;
    type_ind_pair start_ind;

    //iterate until we hit a non whitespace character
    while (str[0] == ' ' || str[0] == '\n' || str[0] == '\t') ++str;
    //variable names are not allowed to start with '+', '-', or a digit and may not contain any '.' symbols. Use this to check whether the value is numeric
    bool is_numeric = (str[0] == '+' || str[0] == '-' || str[0] == '.' || (str[0] >= '0' && str[0] <= '9'));
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
	    if (blk_stk.pop(&start_ind) != E_SUCCESS || start_ind.t != BLK_SQUARE) { printf(/*[*/"unexpected ]\n");er = E_BAD_SYNTAX;return sto; }
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
	    if (blk_stk.pop(&start_ind) != E_SUCCESS || start_ind.t != BLK_PAREN) { printf(/*(*/"unexpected )\n");er = E_BAD_SYNTAX;return sto; }
	    --nest_level;
	    //only set the end paren location if it hasn't been set yet and the stack has no more parenthesis to remove, TODO: make this work with other block types inside a set of parenthesis
	    //if (blk_stk.is_empty() && last_close_ind < 0) last_close_ind = i;
	    last_close_ind = i;
	} else if (str[i] == '{'/*}*/) {
	    //keep track of open and close parenthesis, these will come in handy later
	    blk_stk.push({BLK_CURLY, i});
	    //TODO: decide if this is necessary
	    ++nest_level;
	    //only set the open index if this is the first match
	    if (first_open_ind == -1) { first_open_ind = i; }
	} else if (str[i] == /*{*/'}') {
	    if (blk_stk.pop(&start_ind) != E_SUCCESS || start_ind.t != BLK_CURLY) { printf(/*{*/"unexpected }\n");er = E_BAD_SYNTAX;return sto; }
	    --nest_level;
	    //only set the end paren location if it hasn't been set yet and the stack has no more parenthesis to remove, TODO: make this work with other block types inside a set of parenthesis
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

	    //check if we found a comparison operation symbol
	    if (((str[i] == '=' && str[i+1] == '=') || str[i] == '>' || str[i] == '<') && op_prec < 5) {
		op_prec = 5;
		op_loc = i;
	    } else if (i != 0 && (str[i] == '+' || str[i] == '-') && str[i-1] != 'e' && op_prec < 4) {
		//remember to recurse after we finish looping
		op_prec = 4;
		op_loc = i;
	    } else if (str[i] == '^' && op_prec < 3) {
		op_prec = 3;
		op_loc = i;
	    } else if ((str[i] == '*' || str[i] == '/') && op_prec < 2) {
		op_prec = 2;
		op_loc = i;
	    } else if (op_prec < 1 && (str[i] == '?' || (str[i] == '.' && !is_numeric))) {
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
		//try interpreting as a boolean
		if (strcmp(str, "false") == 0) {
		    sto.type = VAL_NUM;
		    sto.val.x = 0;
		    sto.n_els = 1;
		    return sto;
		} else if (strcmp(str, "true") == 0) {
		    sto.type = VAL_NUM;
		    sto.val.x = 1;
		    sto.n_els = 1;
		    return sto;
		}
		//try interpreting as a number
		errno = 0;
		sto.val.x = strtod(str, NULL);
		if (errno) {
		    printf("undefined token %s\n", str);
		    er = E_NOT_DEFINED;
		    return sto;
		}
		er = E_SUCCESS;
		sto.type = VAL_NUM;
	    }
	    sto = copy_val(sto);
	    str[reset_ind] = ' ';//all whitespace is treated identically so it doesn't matter
	} else if (str[first_open_ind] == '\"' && str[last_close_ind] == '\"') {
	    //this is a string
	    sto.type = VAL_STR;
	    sto.n_els = last_close_ind-first_open_ind;
	    //allocate memory and copy
	    sto.val.s = (char*)malloc(sizeof(char)*sto.n_els);
	    for (size_t j = 0; j < sto.n_els-1; ++j) sto.val.s[j] = str[first_open_ind+j+1];
	    sto.val.s[sto.n_els-1] = 0;
	} else if (str[first_open_ind] == '[' && str[last_close_ind] == ']') {
	    //first check to see if the user is trying to access an element
	    str[first_open_ind] = 0;
	    char* pre_list_name = find_token_before(str, first_open_ind);
	    //if the string is empty then we're creating a new list, otherwise we're accessing an existing list
	    if (pre_list_name[0] == 0) {
		str[first_open_ind] = '[';//]
		sto = parse_list(str+first_open_ind, er);
		sto.type = VAL_LIST;
		if (er != E_SUCCESS) return sto;
	    } else {
		value tmp_lst = lookup(pre_list_name);
		str[first_open_ind] = '[';//]
		if (tmp_lst.type != VAL_LIST) { cleanup_val(&tmp_lst);er = E_BAD_TYPE;return sto; }
		str[last_close_ind] = 0;
		value contents = parse_value(str+first_open_ind+1, er);
		str[last_close_ind] = /*[*/']';
		//check that we found the list and that it was valid
		if (contents.type != VAL_NUM) { cleanup_val(&contents);er = E_BAD_TYPE;return sto; }
		//now figure out the index
		long ind = (long)(contents.val.x);
		if (ind < 0) ind = tmp_lst.n_els+ind;
		if (ind >= tmp_lst.n_els || ind < 0) { er = E_OUT_OF_RANGE;return sto; }
		return copy_val(tmp_lst.val.l[ind]);
	    }
	} else if (str[first_open_ind] == '{' && str[last_close_ind] == '}') {
	    //now parse the argument as a context
	    size_t n_els;
	    str[last_close_ind] = 0;
	    char** list_els = csv_to_list(str+first_open_ind+1, ';', &n_els, er);
	    if (er != E_SUCCESS) { free(list_els);return sto; }
	    //setup the context
	    sto.val.c = new collection(this);
	    sto.n_els = n_els;
	    //insert context members
	    for (size_t j = 0; list_els[j] && j < n_els; ++j) {
		char* cur_name = NULL;
		char* rval = list_els[j];
		char* eq_loc = strchr_block(list_els[j], '=');
		if (eq_loc) {
		    *eq_loc = 0;
		    cur_name = CGS_trim_whitespace(list_els[j], NULL);
		    rval = eq_loc+1;
		}
		value tmp = sto.val.c->parse_value(rval, er);
		if (er != E_SUCCESS) { free(list_els);sto.type=VAL_UNDEF;free(sto.val.c);sto.val.c = NULL;return sto; }
		sto.val.c->emplace(cur_name, tmp);
		cleanup_val(&tmp);
	    }
	    free(list_els);
	    sto.type = VAL_INST;
	    return sto;
	} else if (str[first_open_ind] == '(' && str[last_close_ind] == ')') {
	    //check to see if this is a function call (it is if there are any non-whitespace characters before the open paren
	    int is_func = 0;
	    for (size_t j = 0; j < first_open_ind; ++j) {
		if (str[j] != ' ' && str[j] != '\t' && str[j] != '\n') {
		    //we can't leave this as zero in case the user needs to do some more operations
		    char term_char = str[last_close_ind+1];
		    str[last_close_ind+1] = 0;
		    cgs_func tmp_f;
		    char* f_end;
		    tmp_f = parse_func(str, first_open_ind, er, &f_end);
		    if (er != E_SUCCESS) { return sto; }
		    //TODO: allow user defined functions
		    //er = E_BAD_TOKEN;
		    if (strcmp(tmp_f.name, "vec") == 0) {
			sto = make_vec(tmp_f, er);
		    } else if (strcmp(tmp_f.name, "range") == 0) {
			sto = make_range(tmp_f, er);
		    } else if (strcmp(tmp_f.name, "linspace") == 0) {
			sto = make_linspace(tmp_f, er);
		    } else if (strcmp(tmp_f.name, "flatten") == 0) {
			sto = flatten_list(tmp_f, er);
		    } else if (strcmp(tmp_f.name, "print") == 0) {
			sto = print(tmp_f, er);
		    } else {
			//otherwise lookup the function
			value func_val = lookup(tmp_f.name);
			if (func_val.type == VAL_FUNC) {
			    //make sure that the function was found and that sufficient arguments were provided
			    if (func_val.n_els <= tmp_f.n_args) {
				sto = func_val.val.f->eval(*this, tmp_f, er);
			    }
			} else {
			    er = E_BAD_TYPE;
			}
		    }
		    str[last_close_ind+1] = term_char;
		    cleanup_func(&tmp_f);
		    return sto;
		}
	    }
	    //otherwise interpret this as a parenthetical expression
	    str[last_close_ind] = 0;
	    char* tmp_str = CGS_trim_whitespace(str+first_open_ind+1, &reset_ind);
	    sto = parse_value(tmp_str, er);
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

/**
 * Generate a context from a list of lines. This context will include function declarations, named variables, and subcontexts (instances).
 * lines: the array of lines to read from
 * n_lines: the size of the array
 * returns: an errorcode if one was found or E_SUCCESS on success
 */
parse_ercode context::read_from_lines(line_buffer b) {
    parse_ercode er = E_SUCCESS;
    stack<block_type> blk_stack;

    char last_char = 0;
    size_t buf_size = BUF_SIZE;
    char* buf = (char*)malloc(sizeof(char)*buf_size);
    //iterate over each line in the file
    for (long lineno = 0; lineno < b.get_n_lines(); ++lineno) {
	char* line = b.get_line(lineno);
	//check for class and function declarations
	char* dectype_start = token_block(line, "def");
	if (dectype_start) {
	    char* endptr;
	    cgs_func cur_func = parse_func(dectype_start + KEY_DEF_LEN, -1, er, &endptr, true);
	    size_t func_end = endptr - (dectype_start + KEY_DEF_LEN);
	    //jump ahead until after the end of the function
	    if (er == E_SUCCESS) {
		size_t i = 0;
		for (; endptr[i] && (endptr[i] == ' ' || endptr[i] == '\t' || endptr[i] == '\n'); ++i) (void)0;
		size_t n_enc = b.get_n_lines()-lineno;
		line_buffer lines_enc(b.get_enclosed(lineno, &lineno, '{', '}', func_end));
		if (lineno < 0) { free(buf);return E_BAD_SYNTAX; }
		//now we actually create the function
		user_func tmp(cur_func, lines_enc);
		value v;
		v.type = VAL_FUNC;
		v.val.f = &tmp;
		v.n_els = n_enc;
		emplace(cur_func.name, v);
		//we have to advance to the line after end of the function declaration
	    }
	} else {
	    //char* lval = NULL;
	    size_t rval_start = 0;
	    size_t k = 0;
	    bool started = false; 
	    size_t i = 0;
	    for (; line[i]; ++i) {
		//handle comments
		if (line[i] == '/') {
		    if (line[i+1] == '/') {
			break;
		    } else if (line[i+1] == '*') {
			blk_stack.push(BLK_COMMENT);
		    } else if (i > 0 && line[i-1] == '*') {
			block_type tmp;
			er = blk_stack.pop(&tmp);
			if (tmp != BLK_COMMENT) { return E_BAD_SYNTAX; }
		    }
		}
		block_type cur_blkt = BLK_UNDEF;
		if (blk_stack.size() > 0) cur_blkt = blk_stack.peek();
		if (cur_blkt != BLK_LITERAL && cur_blkt != BLK_COMMENT) {
		    if (k >= buf_size) {
			buf_size *= 2;
			buf = (char*)realloc(buf, sizeof(char)*buf_size);
		    }
		    //ignore preceeding whitespace
		    if (line[i] != ' ' || line[i] != '\t' || line[i] != '\n') started = true;
		    if (started) buf[k++] = line[i];
		    //check for assignment, otherwise handle blocks
		    char match_tok = 0;
		    char sep_tok = 0;
		    if (line[i] == '=') {
			if (blk_stack.is_empty()) {
			    buf[k-1] = 0;
			    //lval = CGS_trim_whitespace(buf, NULL);
			    rval_start = k;
			}
		    } else if (line[i] == '(') {
			match_tok = ')';
		    } else if (line[i] == '[') {
			match_tok = ']';
		    } else if (line[i] == '{') {
			match_tok = '}';
			sep_tok = ';';
		    }
		    if (match_tok != 0 && i > rval_start) {
			line_buffer lines_enc(b.get_enclosed(lineno, &lineno, line[i], match_tok, i, true));
			if (lineno < 0) { free(buf);return E_BAD_SYNTAX; }
			//we have to advance to the end of the function declaration
			char* enc = lines_enc.flatten(sep_tok);
			//concatenate the string using thingies
			size_t dest_len = (k - rval_start) + strlen(enc) + 2;
			buf = (char*)realloc(buf, sizeof(char)*(k+dest_len));
			strcpy(buf+k-1, enc);
			free(enc);
			//parse the value
			value tmp_val = parse_value(buf+rval_start, er);
			if (er != E_SUCCESS) { free(buf);return er; }
			emplace(CGS_trim_whitespace(buf, NULL), tmp_val);
			cleanup_val(&tmp_val);
			break;
		    }
		}
	    }
	    //only set the rval if we haven't done so already
	    if (line[i] == 0) {
		buf[k] = 0;
		value tmp_val = parse_value(buf+rval_start, er);
		if (er != E_SUCCESS) { free(buf);return er; }
		emplace(CGS_trim_whitespace(buf, NULL), tmp_val);
		cleanup_val(&tmp_val);
	    }
	}
    }
    free(buf);

    return er;
}

//

/** ======================================================== user func ======================================================== **/
/**
 * constructor
 * sig: this specifies the signature of the function used when calling it
 * bufptr: a buffer to be used for line reading, see read_cgs_line
 * n: the number of characters currently in the buffer, see read_cgs_line
 * fp: the file pointer to read from
 */
user_func::user_func(cgs_func sig, line_buffer b) : code_lines(b) {
    call_sig = copy_func(sig);
    exec = NULL;
}
/**
 * constructor
 * sig: this specifies the signature of the function used when calling it
 * bufptr: a buffer to be used for line reading, see read_cgs_line
 * n: the number of characters currently in the buffer, see read_cgs_line
 * fp: the file pointer to read from
 */
user_func::user_func(value (*p_exec)(context&, cgs_func, parse_ercode&)) : call_sig() {
    exec = p_exec;
}
//deallocation
user_func::~user_func() {
    cleanup_func(&call_sig);
}
//copy constructor
user_func::user_func(const user_func& o) : code_lines(o.code_lines) {
    call_sig = copy_func(o.call_sig);
    exec = o.exec;
}
//move constructor
user_func::user_func(user_func&& o) : code_lines(o.code_lines) {
    call_sig = o.call_sig;
    exec = o.exec;
    o.call_sig.name = NULL;
    o.exec = NULL;
    o.call_sig.n_args = 0;
}

const char* token_names[] = {"if", "for", "while", "return"};
const size_t n_token_names = sizeof(token_names)/sizeof(char*);
typedef enum {TOK_NONE, TOK_IF, TOK_FOR, TOK_WHILE, TOK_RETURN} token_type;
/**
 * evaluate the function
 */
value user_func::eval(context& c, cgs_func call, parse_ercode& er) {
    value sto;
    if (exec) {
	value ret = (*exec)(c, call, er);
	return ret;
    }
    return sto;
    /*er = E_SUCCESS;
    size_t lineno = 0;
    value bad_val;
    while (lineno < code_lines.get_n_lines()) {
	char* line = code_lines.get_line(lineno);
	//search for tokens
	token_type this_tok = TOK_NONE;
	char* token = NULL;
	for (size_t i = 0; i < n_token_names; ++i) {
	    token = token_block(line, token_names[i]);
	    if (token) {
		this_tok = (token_type)i;
	    }
	}
	switch (this_tok) {
	    //TODO
	    defualt: break;
	}
    }
    //return an undefined value by default
    value ret;ret.type = VAL_UNDEF;ret.val.x = 0;return ret;*/
}
