#include "cgs_read.hpp"

/** ======================================================== utility functions ======================================================== **/

/**
 * check if a character is whitespace
 */
bool is_whitespace(char c) {
    if (c == 0 || c == ' ' || c == '\t' || c == '\n' || c == '\r')
	return true;
    return false;
}

/**
 * Helper function for is_token which tests whether the character c is a token terminator
 */
bool is_char_sep(char c) {
    if (is_whitespace(c) || c == ';' || c == '+'  || c == '-' || c == '*'  || c == '/')
	return true;
    return false;
}

/**
 * Helper function which looks at the string str at index i and tests whether it is a token with the matching name
 */
bool is_token(const char* str, size_t i, size_t len) {
    if (i > 0 && !is_char_sep(str[i-1]))
	return false;
    if (!is_char_sep(str[i+len]))
	return false;
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
 * Returns: the number of characters read (including null termination). On reaching the end of the file, 0 is returned.
 */
size_t read_cgs_line(char** bufptr, size_t* n, FILE* fp, size_t* lineno) {
    //dereference pointers and interpret them correctly
    size_t n_lines = 0;
    if (lineno) n_lines = *lineno;
    size_t size = *n;
    char* buf = *bufptr;
    if (buf) {
	if (size == 0) return 0;
    } else {
	size = LINE_SIZE;
	buf = (char*)malloc(sizeof(char)*size);
    }

    int res = fgetc(fp);
    size_t i = 0;
    for (;; ++i) {
	if (i >= size) {
	    size *= 2;
	    buf = (char*)xrealloc(buf, sizeof(char)*size);
	}
	if (res == EOF) {
	    buf[i] = 0;
	    *n = size;
	    *bufptr = buf;
	    if (lineno) *lineno = n_lines;
	    return 0;
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
	    ret = (char**)xrealloc(ret, sizeof(char*)*(off+1));

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
    ret = (char**)xrealloc(ret, sizeof(char*)*(off+2));
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


line_buffer_ind operator+(const line_buffer_ind& lhs, const size_t& rhs) {
    return line_buffer_ind(lhs.line, lhs.off+rhs);
}
line_buffer_ind operator-(const line_buffer_ind& lhs, const size_t& rhs) {
    line_buffer_ind ret(lhs.line, lhs.off);
    if (rhs > ret.off)
	ret.off = 0;
    else
	ret.off -= rhs;
    return ret;
}

line_buffer::line_buffer() {
    lines = NULL;
    line_sizes = NULL;
    n_lines = 0;
}

line_buffer::line_buffer(const char* p_fname) {
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
	bool go_again = true;
	do {
	    //xreallocate buffer if necessary
	    if (n_lines >= buf_size) {
		buf_size *= 2;
		lines = (char**)xrealloc(lines, sizeof(char*)*buf_size);
		line_sizes = (size_t*)xrealloc(line_sizes, sizeof(size_t)*buf_size);
	    }
	    //read the line until a semicolon, newline or EOF is found
	    size_t this_size = BUF_SIZE;
	    char* this_buf = (char*)malloc(this_size);
	    int res = fgetc(fp);
	    for (line_len = 0; true; ++line_len) {
		//grow the buffer if necessary
		if (line_len >= this_size) {
		    this_size *= 2;
		    this_buf = (char*)xrealloc(this_buf, sizeof(char)*this_size);
		}
		if (res == EOF || (char)res ==';' || (char)res == '\n') {
		    this_buf[line_len] = 0;
		    if ((char)res == '\n')
			++lineno;
		    else if ((char)res == EOF)
			go_again = false;
		    line_len = line_len;
		    break;
		}
		this_buf[line_len] = (char)res;
		res = fgetc(fp);
	    }
	    if (line_len > 0) {
		this_buf = (char*)xrealloc(this_buf, sizeof(char)*(line_len+1));
		lines[n_lines] = this_buf;
		line_sizes[n_lines++] = line_len;
	    } else {
		free(this_buf);
	    }
	} while (go_again);
	lines = (char**)xrealloc(lines, sizeof(char*)*n_lines);
	line_sizes = (size_t*)xrealloc(line_sizes, sizeof(size_t)*n_lines);
	fclose(fp);
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
	    if (str[i] == ignore_blocks[j] && ignore_blocks[j] == ignore_blocks[j+1]) {
		//there's a special case for things like quotations, skip ahead until we find the end
		++i;
		for (;str[i] && str[i] != ignore_blocks[j]; ++i) (void)0;
		break;
	    } else if (str[i] == ignore_blocks[j]) {
		++nest_level;
	    } else if (str[i] == ignore_blocks[j+1]) {
		--nest_level;
	    }
	}
	if (nest_level <= 0 && (str[i] == sep || str[i] == 0)) {
	    line_sizes[n_lines] = i-last_sep;
	    lines[n_lines] = (char*)malloc(sizeof(char)*(line_sizes[n_lines]+1));
	    for (size_t j = last_sep; j < i; ++j)
		lines[n_lines][j-last_sep] = str[j];
	    lines[n_lines][line_sizes[n_lines]] = 0;
	    last_sep = i+1;
	    ++n_lines;
	    if (str[i] == 0) break;
	}
    }
}
line_buffer::~line_buffer() {
    if (lines) {
	for (size_t i = 0; i < n_lines; ++i) {
	    free(lines[i]);
	}
	free(lines);
    }
    if (line_sizes) {
	free(line_sizes);
    }
    n_lines = 0;
    lines = NULL;
    line_sizes = NULL;
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
line_buffer& line_buffer::operator=(line_buffer& o) {
    size_t tmp_n_lines = n_lines;
    n_lines = o.n_lines;
    o.n_lines = tmp_n_lines;
    size_t* tmp_line_sizes = line_sizes;
    line_sizes = o.line_sizes;
    o.line_sizes = tmp_line_sizes;
    char** tmp_lines = lines;
    lines = o.lines;
    o.lines = tmp_lines;
    return *this;
}

line_buffer& line_buffer::operator=(line_buffer&& o) {
    n_lines = o.n_lines;
    lines = o.lines;
    line_sizes = o.line_sizes;
    o.n_lines = 0;
    o.line_sizes = NULL;
    o.lines = NULL;
    return *this;
}

/**
 * Goes through a line buffer and splits into multiple lines at each instance of split_delim.
 * i.e. if the buffer is in the state
 * lines = {"foo; bar", "foobar;"}
 * then split(';') will transform the state to
 * lines = {"foo", " bar", "foobar;"}
 */
void line_buffer::split(char split_delim) {
    for (size_t i = 0; i < n_lines; ++i) {
	for (size_t j = 0; j < line_sizes[i]; ++j) {
	    if (lines[i][j] == split_delim) {
		++n_lines;
		printf("new size = %lu\n", sizeof(char*)*n_lines);
		lines = (char**)xrealloc(lines, sizeof(char*)*n_lines);
		line_sizes = (size_t*)xrealloc(line_sizes, sizeof(size_t)*n_lines);
		//move everything else forward
		for (size_t k = n_lines-1; k > i+1; --k) {
		    lines[k] = lines[k-1];
		    line_sizes[k] = line_sizes[k-1];
		}
		//split this line
		lines[i+1] = strndup(lines[i]+j+1, line_sizes[i]-j);
		lines[i][j] = 0;
		line_sizes[i+1] = line_sizes[i]-j-1;
		line_sizes[i] = j;
	    }
	}
    }
}

/**
 * helper function for get_enclosed and jmp_enclosed. This function reads the line at index k for the line contained between start_delim and end_delim. If only a start_delim is found or start_ind is not NULL, a pointer with a value set to the index of the start of the line is returned.
 * linesto: store the line read into this variable
 * start_delim: the starting delimiter
 * end_delim: the ending delimiter
 * start: the position to start reading. This is updated as soon as the start delimiter is actually found
 * end: the position where reading stopped. This will be at offset 0 on the line after start->line if an end delimiter was not found or the character where reading stopped because of an end delimiter.
 * includ_delims: if true, then the delimiters are included in the string
 * start_ind: a pointer, the value of which is the index in the line k to start reading.
 * depth: keeps track of how many nested pairs of start and end delimeters we've encountered. We only want to exit calling if an end_delim was found. This variable is set to -1 if a zero depth close brace was found to signal that parsing should terminate.
 * end_line: If end delim is set, then the value in end_line is set to k.
 * jp: stores the index where line reading stopped. Either because a NULL terminator was encountered or an end_delim was encountered
 * returns: an index of the position of the index that reading started from or -1 if not found
 */
int line_buffer::it_single(char** linesto, char start_delim, char end_delim, line_buffer_ind* start, line_buffer_ind* end, int* pdepth, bool include_delims, bool include_start) const {
    bool free_after = false;
    if (start == NULL) {
	start = (line_buffer_ind*)malloc(sizeof(line_buffer_ind));
	start->line = 0;
	start->off = 0;
    }
    //setup
    int depth = 0;
    if (pdepth) depth = *pdepth;
    size_t j = start->off;
    size_t init_off = start->off;
    int ret = -1;
    //iterate through characters in the line looking for an end_delim without a preceeding start_delim
    size_t i = start->line;
    end->line = i;
    while (true) {
	if (lines[i][j] == 0) {
	    if(end) {
		end->line = i+1;
		end->off = j;
	    }
	    break;
	} else if (lines[i][j] == start_delim) {
	    //this is a special case, depths may only be toggled
	    if (end_delim == start_delim) {
		if (j == 0 || lines[i][j-1] != '\\') {
		    depth = 1 - depth;
		    if (depth == 1) {
			//start block
			start->line = i;
			start->off = j;
			ret = j;
			if (!include_delims) ++(start->off);
		    } else {
			//end block	
			if (end) {
			    if (include_delims) ++j;
			    end->off = j;
			}
			break;
		    }
		}
	    } else {
		if (depth == 0) {
		    start->line = i;
		    start->off = j;
		    ret = j;
		    if (!include_delims) ++(start->off);
		}
		++depth;
	    }
	} else if (lines[i][j] == end_delim) {
	    --depth;
	    if (depth <= 0) {
		if (include_delims) ++j;
		if (end)
		    end->off = j;
		//--depth;//force escaping from the loop
		break;
	    }
	}
	++j;
    }
    if (pdepth) *pdepth = depth;
    if (linesto) {
	if (include_start)
	    *linesto = lines[i]+init_off;
	else
	    *linesto = lines[i]+(start->off);
    }
    if (include_start)
	start->off = init_off;
    if (free_after) free(start);
    return ret;
}

/**
 * Find the line buffer starting on line start_line between the first instance of start_delim and the last instance of end_delim respecting nesting (i.e. if lines={"a {", "b {", "}", "} c"} then {"b {", "}"} is returned. Note that the result must be deallocated with a call to free().
 * start_line: the line to start reading from
 * end_line: if this value is not NULL, then the index of the line on which end_delim was found is stored here. If the end delimeter is not found, then the line is set to n_lines and the offset is set to zero
 * start_delim: the character to be used as the starting delimiter. This needs to be supplied so that we find a matching end_delim at the root level
 * end_delim: the character to be searched for
 * line_offset: the character on line start_line to start reading from, this defaults to zero. Note that this only applies to the start_line, and not any subsequent lines in the buffer.
 * include_delims: if true, then the delimeters are included in the enclosed strings. defualts to false
 * include_start: if true, then the part preceeding the first instance of start_delim will be included. This value is always false if include_delims is false. If include_delims is true, then this defaults to true.
 * returns: an array of lines that must be deallocated by a call to free(). The number of lines is stored in n_lines.
 */
line_buffer line_buffer::get_enclosed(line_buffer_ind start, line_buffer_ind* pend, char start_delim, char end_delim, bool include_delims, bool include_start) const {
    line_buffer ret;
    //initialization
    if (pend) {
	pend->line = n_lines;
	pend->off = 0;
    }
    if (n_lines == 0) {
	ret.n_lines = 0;
	ret.lines = NULL;
	ret.line_sizes = NULL;
	return ret;
    }
    //set include_start to false if include_delims is false
    include_start &= include_delims;
    ret.lines = (char**)malloc(sizeof(char*)*n_lines);
    ret.line_sizes = (size_t*)malloc(sizeof(size_t)*n_lines);

    //tracking variables
    int depth = 0;

    //iterate through lines
    size_t k = 0;
    size_t start_line = start.line;
    line_buffer_ind end;
    size_t i = start.line;
    bool started = false;
    for (; depth >= 0 && i < n_lines; ++i) {
	if (lines[i] == NULL) {
	    ret.n_lines = i-start_line;
	    break;
	} else {
	    end.line = i;
	    end.off = 0;
	    char* this_line;
	    int start_ind = it_single(&this_line, start_delim, end_delim, &start, &end, &depth, include_delims, include_start);
	    if (start_ind >= 0) started = true;
	    //don't read empty lines
	    if (started) {
		ret.line_sizes[k] = end.off-start.off;
		ret.lines[k++] = strndup(this_line, end.off-start.off);
	    }
	    //This means that an end delimeter was found. In this case, we need to break out of the loop.
	    if (end.line == start.line) break;
	    start.off = 0;
	    ++start.line;
	}
    }
    if (pend) *pend = end;
    ret.n_lines = k;
    return ret;
}

line_buffer_ind line_buffer::jmp_enclosed(line_buffer_ind start, char start_delim, char end_delim, bool include_delims) const {
    int depth = 0;
    for (size_t i = start.line; depth >= 0 && i < n_lines; ++i) {
	line_buffer_ind end(i, 0);
	it_single(NULL, start_delim, end_delim, &start, &end, &depth, include_delims, false);
	if (end.line == start.line) {
	    return end;
	}
	start.off = 0;
	++start.line;
    }
    line_buffer_ind ret(n_lines, 0);
    return ret;
}

/**
 * Return a string with the line contained at index i. This string should be freed with a call to free().
 */
char* line_buffer::get_line(line_buffer_ind p) const {
    if (p.line >= n_lines)
	return NULL;
    size_t line_size = line_sizes[p.line]+1;
    char* ret = (char*)malloc(sizeof(char)*(line_size-p.off));
    for (size_t j = p.off; j < line_size; ++j) ret[j-p.off] = lines[p.line][j];
    ret[line_size-p.off-1] = 0;
    return ret;
}

/**
 * Returns a version of the line buffer which is flattened so that everything fits onto one line.
 * sep_char: if this is not zero, then each newline in the buffer is replaced by a sep_char
 */
char* line_buffer::flatten(char sep_char) const {
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

bool line_buffer::inc(line_buffer_ind& p) const {
    if (p.line >= n_lines) return false;
    if (p.off >= line_sizes[p.line]) {
	if (p.line == n_lines-1) return false;
	p.off = 0;
	p.line += 1;
    } else {
	p.off += 1;
    }
    return true;
}

bool line_buffer::dec(line_buffer_ind& p) const {
    if (p.line > n_lines) return false;
    if (p.off == 0) {
	if (p.line == 0) return false;
	p.line -= 1;
	p.off = line_sizes[p.line];
    } else {
	p.off -= 1;
    }
    return true;
}

/**
 * returns the character at position pos
 */
char line_buffer::get(line_buffer_ind pos) const {
    if (pos.line >= n_lines || pos.off >= line_sizes[pos.line]) return 0;
    return lines[pos.line][pos.off];
}

/** ======================================================== builtin functions ======================================================== **/
/**
 * Get the type of a value
 */
value typeof(context* c, cgs_func tmp_f, parse_ercode& er) {
    value sto;
    sto.type = VAL_STR;
    if (tmp_f.n_args < 1) { er = E_LACK_TOKENS;return sto; }
    switch (tmp_f.args[0].type) {
	case VAL_UNDEF: sto.n_els = strlen("undefined");sto.val.s = strdup("undefined"); break;
	case VAL_STR: sto.n_els = strlen("string");sto.val.s = strdup("string"); break;
	case VAL_NUM: sto.n_els = strlen("number");sto.val.s = strdup("number"); break;
	case VAL_ARRAY: sto.n_els = strlen("array");sto.val.s = strdup("array"); break;
	case VAL_LIST: sto.n_els = strlen("list");sto.val.s = strdup("list"); break;
	case VAL_3VEC: sto.n_els = strlen("vector_3");sto.val.s = strdup("vector_3"); break;
	case VAL_MAT: sto.n_els = strlen("matrix_3x3");sto.val.s = strdup("matrix_3x3"); break;
	case VAL_FUNC: sto.n_els = strlen("function");sto.val.s = strdup("function"); break;
	case VAL_INST:
	   value t = tmp_f.args[0].lookup("__type__");
	   if (t.type == VAL_STR)
	       sto.val.s = strdup(t.val.s);
	   break;
	default: sto.type = VAL_UNDEF;er = E_BAD_VALUE; break;
    }
    return sto;
}
/**
 * Make a vector argument with the x,y, and z coordinates supplied
 */
value make_vec(context* c, cgs_func tmp_f, parse_ercode& er) {
    value sto = make_val_undef();
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
value make_range(context* c, cgs_func tmp_f, parse_ercode& er) {
    value sto = make_val_undef();
    if (tmp_f.n_args == 0) {
	er = E_LACK_TOKENS;
	return sto;
    } else if (tmp_f.n_args == 1) {
	if (tmp_f.args[0].type != VAL_NUM) {
	    printf("Error: ranges can only be specified with numeric types\n");
	    er = E_BAD_TYPE;
	    return sto;
	}
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
	if (tmp_f.args[0].type != VAL_NUM) {
	    printf("Error: ranges can only be specified with numeric types\n");
	    er = E_BAD_TYPE;
	    return sto;
	}
	if (tmp_f.args[1].type != VAL_NUM) {
	    printf("Error: ranges can only be specified with numeric types\n");
	    er = E_BAD_TYPE;
	    return sto;
	}
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
	if (tmp_f.args[0].type != VAL_NUM) {
	    printf("Error: ranges can only be specified with numeric types\n");
	    er = E_BAD_TYPE;return sto; }
	if (tmp_f.args[1].type != VAL_NUM) {
	    printf("Error: ranges can only be specified with numeric types\n");
	    er = E_BAD_TYPE;return sto; }
	if (tmp_f.args[2].type != VAL_NUM) {
	    printf("Error: ranges can only be specified with numeric types\n");
	    er = E_BAD_TYPE;return sto; }
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
value make_linspace(context* c, cgs_func tmp_f, parse_ercode& er) {
    value sto;
    if (tmp_f.n_args < 3) {
        er = E_LACK_TOKENS;
        return sto;
    } else {
        if (tmp_f.args[0].type == VAL_NUM && tmp_f.args[1].type == VAL_NUM && tmp_f.args[2].type == VAL_NUM) {
            sto.type = VAL_ARRAY;
            sto.n_els = (size_t)(tmp_f.args[2].val.x);
            //prevent divisions by zero
            if (sto.n_els < 2) {
                sto.type = VAL_UNDEF;
                sto.n_els = 0;
                sto.val.l = NULL;
                er = E_BAD_VALUE;
                return sto;
            }
            sto.val.a = (double*)malloc(sizeof(double)*sto.n_els);
            double step = (tmp_f.args[1].val.x - tmp_f.args[0].val.x)/(sto.n_els - 1);
            for (size_t i = 0; i < sto.n_els; ++i) {
                sto.val.a[i] = step*i + tmp_f.args[0].val.x;
            }
            er = E_SUCCESS;
            return sto;
        } else {
	    printf("Error: linspaces can only be specified with numeric types");
            er = E_BAD_TYPE;
            return sto;
        }    
    }
}
/**
 * Take a list value and flatten it so that it has numpy dimensions (n) where n is the sum of the length of each list in the base list. values are copied in order e.g flatten([0,1],[2,3]) -> [0,1,2,3]
 * cgs_func: the function with arguments passed
 */
value flatten_list(context* c, cgs_func tmp_f, parse_ercode& er) {
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
		value* tmp_val = (value*)xrealloc(sto.val.l, sizeof(value)*buf_size);
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
value print(context* c, cgs_func tmp_f, parse_ercode& er) {
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

//math functions
value fun_sin(context* c, cgs_func f, parse_ercode& er) {
    if (f.n_args < 1) { er = E_LACK_TOKENS;c->error("expected argument for sin()");return make_val_undef(); }
    if (f.args[0].type != VAL_NUM) { er = E_BAD_TYPE;c->error("sin() only accept numbers");return make_val_undef(); }
    return make_val_num( sin(f.args[0].val.x) );
}
value fun_cos(context* c, cgs_func f, parse_ercode& er) {
    if (f.n_args < 1) { er = E_LACK_TOKENS;c->error("expected argument for cos()");return make_val_undef(); }
    if (f.args[0].type != VAL_NUM) { er = E_BAD_TYPE;c->error("cos() only accept numbers");return make_val_undef(); }
    return make_val_num( cos(f.args[0].val.x) );
}
value fun_tan(context* c, cgs_func f, parse_ercode& er) {
    if (f.n_args < 1) { er = E_LACK_TOKENS;c->error("expected argument for tan()");return make_val_undef(); }
    if (f.args[0].type != VAL_NUM) { er = E_BAD_TYPE;c->error("tan() only accept numbers");return make_val_undef(); }
    return make_val_num( tan(f.args[0].val.x) );
}
value fun_exp(context* c, cgs_func f, parse_ercode& er) {
    if (f.n_args < 1) { er = E_LACK_TOKENS;c->error("expected argument for tan()");return make_val_undef(); }
    if (f.args[0].type != VAL_NUM) { er = E_BAD_TYPE;c->error("tan() only accept numbers");return make_val_undef(); }
    return make_val_num( exp(f.args[0].val.x) );
}
value fun_sqrt(context* c, cgs_func f, parse_ercode& er) {
    if (f.n_args < 1) { er = E_LACK_TOKENS;c->error("expected argument for tan()");return make_val_undef(); }
    if (f.args[0].type != VAL_NUM) { er = E_BAD_TYPE;c->error("tan() only accept numbers");return make_val_undef(); }
    return make_val_num( sqrt(f.args[0].val.x) );
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
    v.n_els = s.size() + 1;
    v.val.s = (char*)malloc(sizeof(char)*v.n_els);
    for (size_t i = 0; i < v.n_els; ++i) v.val.s[i] = s[i];
    return v;
}
value make_val_array(std::vector<double> a) {
    value v;
    v.type = VAL_ARRAY;
    v.n_els = a.size();
    v.val.a = (double*)malloc(sizeof(double)*v.n_els);
    for (size_t i = 0; i < v.n_els; ++i) v.val.a[i] = a[i];
    return v;
}
value make_val_list(const value* vs, size_t n_vs) {
    value v;
    v.type = VAL_LIST;
    v.n_els = n_vs;
    v.val.l = (value*)malloc(sizeof(value)*v.n_els);
    for (size_t i = 0; i < v.n_els; ++i) v.val.l[i] = copy_val(vs[i]);
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
value make_val_func(const char* name, size_t n_args, value (*p_exec)(context*, cgs_func, parse_ercode&)) {
    value ret;
    ret.type = VAL_FUNC;
    ret.n_els = n_args;
    ret.val.f = new user_func(p_exec);
    return ret;
}
/**
 * make an instance object with the given type
 */
value make_val_inst(context* parent, const char* s) {
    value v;
    v.type = VAL_INST;
    v.val.c = new context(parent);
    if (s && s[0] != 0) {
	value tmp = make_val_str(s);
	v.val.c->place_value("__type__", tmp);
    }
    return v;
}

bool is_type(value v, const char* str) {
    if (v.type == VAL_INST) {
	value type_str = v.val.c->lookup("__type__");
	if (type_str.type == VAL_STR)
	    return (strcmp(str, type_str.val.s) == 0);
    } else if (strcmp(str, "undefined") == 0) {
	return (v.type == VAL_UNDEF);
    } else if (strcmp(str, "string") == 0) {
	return (v.type == VAL_STR);
    } else if (strcmp(str, "numeric") == 0) {
	return (v.type == VAL_NUM);
    } else if (strcmp(str, "array") == 0) {
	return (v.type == VAL_ARRAY);
    } else if (strcmp(str, "list") == 0) {
	return (v.type == VAL_LIST);
    } else if (strcmp(str, "vec3") == 0) {
	return (v.type == VAL_3VEC);
    } else if (strcmp(str, "matrix") == 0) {
	return (v.type == VAL_MAT);
    } else if (strcmp(str, "function") == 0) {
	return (v.type == VAL_FUNC);
    } else if (strcmp(str, "instance") == 0) {
	return (v.type == VAL_INST);
    }
    return false;
}

void cleanup_val(value* v) {
    if ((v->type == VAL_STR && v->val.s) || (v->type == VAL_ARRAY && v->val.a)) {
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
    v->type = VAL_UNDEF;
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
    } else if (o.type == VAL_ARRAY) {
	ret.val.a = (double*)malloc(sizeof(double)*o.n_els);
	for (size_t i = 0; i < o.n_els; ++i) ret.val.a[i] = o.val.a[i];
    } else if (o.type == VAL_LIST) {
	if (o.val.l == NULL) {
	    ret.n_els = 0;
	    ret.val.l = NULL;
	} else {
	    ret.val.l = (value*)calloc(o.n_els, sizeof(value));
	    for (size_t i = 0; i < o.n_els; ++i) ret.val.l[i] = copy_val(o.val.l[i]);
	}
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
bool value::operator==(const value& o) const {
    if (type != o.type) return false;
    if (type == VAL_NUM) {
	return (val.x == o.val.x);
    } else if (type == VAL_STR) {
	//first make sure that both strings are not null
	if ((val.s == NULL || o.val.s == NULL) && val.s != o.val.s) return false;
	return (strcmp(val.s, o.val.s) == 0);
    } else if (type == VAL_LIST) {
	if ((val.l == NULL || o.val.l == NULL) && val.l != o.val.l) return false;
	if (n_els != o.n_els) return false;
	for (size_t i = 0; i < n_els; ++i) {
	    bool res = (val.l[i] == o.val.l[i]);
	    if (!res) return false;
	}
	return true;
    } else if (type == VAL_3VEC) {
	if (val.v == o.val.v)
	    return true;
	if (val.v != NULL && o.val.v != NULL)
	    return *(val.v) == *(o.val.v);
	return false;
    } else if  (type == VAL_MAT) {
	if (val.v == o.val.v)
	    return true;
	if (val.v != NULL && o.val.v != NULL)
	    return *(val.m) == *(o.val.m);
	return false;
    }
    //TODO: implement other comparisons
    return false;
}
bool value::operator!=(const value& o) const {
    return !(*this == o);
}
bool value::operator==(const std::string str) const {
    if (type != VAL_STR) return false;
    size_t strlen = str.size();
    if (strlen != n_els-1) return false;
    for (size_t i = 0; i < strlen; ++i) {
	if (str[i] != val.s[i]) return false;
    }
    return true;
}
bool value::operator!=(const std::string str) const {
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
    //exit if there isn't enough space to write the null terminator
    if (n <= 1) return 0;
    if (type == VAL_STR) {
	char* end = stpncpy(sto, val.s, n-1);
	*end = 0;
	return end-sto;
    } else if (type == VAL_ARRAY) {
	int ret = 1;
	sto[0] = '{';//}
	for (size_t i = 0; i < n_els; ++i) {
	    size_t rem = n-ret;
	    int tmp = write_numeric(sto+(size_t)ret, rem, val.a[i]);
	    if (tmp < 0) {
		sto[ret] = 0;
		return ret;
	    }
	    if (tmp >= rem) {
		sto[n-1] = 0;
		return n-1;
	    }
	    ret += (size_t)tmp;
	    if (i+1 < n_els)
		sto[ret++] = ',';
	}
	if (ret < 0) return ret;
	if ((size_t)ret < n)
	    sto[ret++] = '}';
	sto[ret] = 0;
	return ret;
    } else if (type == VAL_NUM) {
	return write_numeric(sto, n, val.x);
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
	    //list -> vector3
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
	    //vector3 -> list
	    ret.n_els = 3;
	    ret.val.l = (value*)malloc(sizeof(value)*ret.n_els);
	    if (!ret.val.l) { ret.type = VAL_UNDEF;ret.n_els = 0;ret.val.x = 0;return ret; }
	    for (size_t i = 0; i < ret.n_els; ++i) {
		ret.val.l[i].type = VAL_NUM;
		ret.val.l[i].n_els = 1;
		ret.val.l[i].val.x = val.v->el[i];
	    }
	    return ret;
	} else if (type == VAL_MAT) {
	    //matrix -> list
	    //TODO
	} else if (type == VAL_INST) {
	    //instance -> list
	    ret.n_els = val.c->size();
	    ret.val.l = (value*)calloc(ret.n_els, sizeof(value));
	    if (!ret.val.l) { ret.type = VAL_UNDEF;ret.n_els = 0;ret.val.x = 0;return ret; }
	    for (size_t i = 1; i < ret.n_els; ++i) {
		ret.val.l[i] = copy_val(val.c->peek_val(i));
	    }
	    ret.n_els = n_els;
	    return ret;
	}
    } else if (t == VAL_ARRAY) {
	if (type == VAL_3VEC) {
	    //vector3 -> array
	    ret.n_els = 3;
	    ret.val.a = (double*)malloc(sizeof(double)*ret.n_els);
	    if (!ret.val.a) { ret.type = VAL_UNDEF;ret.n_els = 0;ret.val.x = 0;return ret; }
	    for (size_t i = 0; i < ret.n_els; ++i) {
		ret.val.a[i] = val.v->el[i];
	    }
	    return ret;
	} else if (type == VAL_LIST) {
	    //list -> array
	    ret.n_els = n_els;
	    ret.val.a = (double*)malloc(sizeof(double)*ret.n_els);
	    if (!ret.val.a) { ret.type = VAL_UNDEF;ret.n_els = 0;ret.val.x = 0;return ret; }
	    parse_ercode tmp_er;
	    for (size_t i = 0; i < ret.n_els; ++i) {
		ret.val.a[i] = val.l[i].cast_to(VAL_NUM, tmp_er).val.x;
		if (tmp_er != E_SUCCESS) { er = tmp_er; return ret; }
	    }
	    return ret;
	}
    } else if (t == VAL_STR) {
	//anything -> string
	ret.val.s = (char*)malloc(sizeof(char)*BUF_SIZE);
	int n_write = rep_string(ret.val.s, BUF_SIZE);
	ret.val.s = (char*)xrealloc(ret.val.s, sizeof(char)*(n_write+1));
    }
    //if we reach this point in execution then there was an error
    ret.type = VAL_UNDEF;
    ret.n_els = 0;
    return ret;
}

void print_spaces(FILE* f, size_t n) {
    for (size_t i = 0; i < n; ++i)
	fprintf(f, " |");
}
/**
 * recursively print out a value and the values it contains
 */
void value::print_hierarchy(FILE* f, size_t depth) const {
    if (f == NULL)
	f = stdout;
    //print the left tree view thingie
    print_spaces(f, depth);
    //now we handle all of the simple (non-recursive) prints
    if (type == VAL_UNDEF) {
	fprintf(f, "UNDEFINED\n");
    } else if (type == VAL_UNINIT) {
	fprintf(f, "UNINITIALIZED\n");
    } else if (type == VAL_NUM) {
	fprintf(f, "%f\n", val.x);
    } else if (type == VAL_STR) {
	fprintf(f, "\"%s\"\n", val.s);
    } else if (type == VAL_3VEC || type == VAL_MAT) {
	char* vec_str = (char*)malloc(sizeof(char)*BUF_SIZE);
	val.v->to_str(vec_str, BUF_SIZE);
	fprintf(f, "%s\n", vec_str);
	free(vec_str);
    } else if (type == VAL_FUNC) {
	fprintf(f, "%p\n", val.f);
    //now handle all the recursive prints
    } else if (type == VAL_LIST) {
	fprintf(f, "[+\n"/*]*/);
	for (size_t i = 0; i < n_els; ++i)
	    val.l[i].print_hierarchy(f, depth+1);
	print_spaces(f, depth);
	fprintf(f, /*[*/"]\n");
    } else if (type == VAL_ARRAY) {
	fprintf(f, "([+\n"/*])*/);
	for (size_t i = 0; i < n_els; ++i) {
	    print_spaces(f, depth+2);
	    fprintf(f, "|%f\n", val.a[i]);
	}
	print_spaces(f, depth+1);
	fprintf(f, /*([*/"])\n");
    } else if (type == VAL_INST) {
	fprintf(f, "{+\n"/*}*/);
	for (size_t i = val.c->size(); i > 0; --i) {
	    name_val_pair pair = val.c->inspect(i);
	    value tmp = pair.get_val();
	    print_spaces(f, depth+1);
	    fprintf(f, "%s:", pair.get_name());
	    //use space economically for simple types
	    if (tmp.type == VAL_LIST || tmp.type == VAL_ARRAY || tmp.type == VAL_INST) {
		fprintf(f, " -V\n");
		tmp.print_hierarchy(f, depth+1);
	    } else {
		fprintf(f, " ");
		tmp.print_hierarchy(f, 0);
	    }
	}
	print_spaces(f, depth);
	fprintf(f, /*{*/"}\n");
    }
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

/** ======================================================== context ======================================================== **/

/**
 * Given the string starting at token, and the index of an open paren parse the result into a cgs_func struct.
 * token: a c-string which is modified in place that contains the function
 * open_par_ind: the location of the open parenthesis
 * f: the cgs_func that information should be saved to
 * end: If not NULL, a pointer to the first character after the end of the string is stored here. If an error occurred during parsing end will be set to NULL.
 * name_only: if set to true, then only the names of the function arguments will be set and all values will be set to undefined. This is useful for handling function declarations
 * returns: an errorcode if an invalid string was supplied.
 */
cgs_func context::parse_func(char* token, long open_par_ind, parse_ercode& er, const char* const * end, int name_only) {
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
 
    //read the arguments separated by spaces
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
	    f.arg_names[i] = list_els[i];
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
	}
    }
    //cleanup and reset the string
    free(list_els);
    return f;
}

void context::copy_context(const context& o) {
    //we can't copy parents, so use the same
    parent = o.parent;
    //allocate the buffer
    grow(o.stack_ptr+1);
    stack_ptr = o.stack_ptr;
    for (size_t i = 0; i < stack_ptr; ++i) {
	buf[i].i = o.buf[i].i;
	buf[i].v = copy_val(o.buf[i].v);
    }
    //now we need to copy over the table
    for (size_t i = 0; i < TABLE_SIZE(TABLE_BITS); ++i) {
	if (o.table[i]) {
	    name_ind* o_cur = o.table[i];
	    table[i] = new name_ind();
	    name_ind* t_cur = table[i];
	    while (true) {
		if (o_cur->name)
		    t_cur->name = strdup(o_cur->name);
		t_cur->ind = o_cur->ind;
		//break if this is the end of the list, otherwise repeat using the next element
		if (o_cur->next == NULL)
		    break;
		t_cur->next = new name_ind();
		o_cur = o_cur->next;
		t_cur = t_cur->next;
	    }
	} else {
	    table[i] = NULL;
	}
    }
}
context::context(context&& o) {
    parent = o.parent;
    stack_ptr = o.stack_ptr;
    for (size_t i = 0; i < stack_ptr; ++i) {
	buf[i] = o.buf[i];
	o.buf[i].v = make_val_undef();
    }
    //moving linked lists is actually trivial
    for (size_t i = 0; i < TABLE_SIZE(TABLE_BITS); ++i) {
	table[i] = o.table[i];
	o.table[i] = NULL;
    }
}

/**
 * cleanup
 */
context::~context() {
    //erase the hash table
    for (size_t i = 0; i < TABLE_SIZE(TABLE_BITS); ++i) {
	if (table[i])
	    delete table[i];
    }
    for (size_t i = 0; i < stack_ptr; ++i)
	cleanup_val(&(buf[i].v));
}

/**
 * setup pointers to all the builtin functions
 */
void context::setup_builtins() {
    value tmp_f = make_val_func("vec", 3, &make_vec);
    place_value("vec", tmp_f);
    tmp_f = make_val_func("range", 1, &make_range);
    place_value("range", tmp_f);
    tmp_f = make_val_func("linspace", 3, &make_linspace);
    place_value("linspace", tmp_f);
    tmp_f = make_val_func("flatten", 1, &flatten_list);
    place_value("flatten", tmp_f);
    tmp_f = make_val_func("print", 1, &print);
    place_value("print", tmp_f);
    tmp_f = make_val_func("sin", 1, &fun_sin);
    place_value("sin", tmp_f);
    tmp_f = make_val_func("cos", 1, &fun_cos);
    place_value("cos", tmp_f);
    tmp_f = make_val_func("tan", 1, &fun_tan);
    place_value("tan", tmp_f);
    tmp_f = make_val_func("exp", 1, &fun_exp);
    place_value("exp", tmp_f);
    tmp_f = make_val_func("sqrt", 1, &fun_sqrt);
    place_value("sqrt", tmp_f);
}

//non-cryptographically hash the string str
size_t fnv_1(const char* str) {
    if (str == NULL)
	return 0;
    size_t ret = FNV_OFFSET;
    for (size_t i = 0; str[i]; ++i) {
	ret = ret^str[i];
	ret = ret*FNV_PRIME;
    }
#if TABLE_BITS > 15
    return (ret >> TABLE_BITS) ^ (ret & TABLE_MASK(TABLE_BITS));
#else
    return ((ret >> TABLE_BITS) ^ ret) & TABLE_MASK(TABLE_BITS);
#endif
}

/**
 * Iterate through the context and assign a new value to the variable with the matching name.
 * name: the name of the variable to set
 * new_val: the value to set the variable to
 * force_push: If set to true, then set_value is guaranteed to increase the stack size by one, even if there is already an element named p_name. This element is guaranteed to receive priority over the existing element, so this may be used to simulate scoped variables.
 * move_assign: If set to true, then the value is directly moved into the context. This can save some time.
 * returns: E_SUCCESS if the variable with a matching name was found or E_NOT_DEFINED otherwise
 */
parse_ercode context::set_value(const char* p_name, value p_val, bool force_push, bool move_assign) {
    //generate a fake name if none was provided
    if (!p_name || p_name[0] == 0) {
	char tmp[BUF_SIZE];
	snprintf(tmp, BUF_SIZE, "\e_%lu", stack_ptr);
	set_value(tmp, p_val, force_push, move_assign);
	return E_SUCCESS;
    }
    size_t ti = fnv_1(p_name);
    //now add it to the hash table
    name_ind* cur = table[ti];
    //iterate through the linked list until we find a matching name
    while(cur && !force_push) {
	if (strcmp(p_name, cur->name) == 0) {
	    cleanup_val(&(buf[cur->ind].v));
	    buf[cur->ind].v = (move_assign)? p_val : copy_val(p_val);
	    return E_SUCCESS;
	}
	cur = cur->next;
    }
    //if we reach this point then no entry with the matching name was found, add it
    table[ti] = new name_ind(p_name, stack_ptr, table[ti]);
    if (force_push)
	table[ti]->next = cur;
    value tmp_v = (move_assign)? p_val : copy_val(p_val);
    val_ind inst(ti, tmp_v);
    push(inst);
    return E_SUCCESS;
}

/**
 * Iterate through the context and return a value to the variable with the matching name.
 * name: the name of the variable to set
 * returns: the matching value, no deep copies are performed
 */
value context::lookup(const char* str) const {
    value ret;
    //lookup
    name_ind* cur = table[fnv_1(str)];
    while(cur) {
	if (strcmp(str, cur->name) == 0)
	    return buf[cur->ind].v;
	cur = cur->next;
    }
    //try searching through the parent if that didn't work
    if (parent)
	return parent->lookup(str);
    //reaching this point in execution means the matching entry wasn't found
    ret.type = VAL_UNDEF;
    ret.val.x = 0;
    return ret;
}

value context::peek_val(size_t i) {
    if (i < stack_ptr)
	return buf[stack_ptr - i].v;
    value ret;ret.type = VAL_UNDEF;ret.n_els = 0;ret.val.x = 0;
    return ret;
}

name_val_pair context::inspect(size_t i) {
    if (i > 0 && i <= stack_ptr) {
	//get the appropriate stack element and search through its corresponding hash table bucket
	size_t j = stack_ptr-i;
	char* name = NULL;
	name_ind* cur = table[buf[j].i];
	while (cur) {
	    if (cur->ind == j) {
		return name_val_pair(cur->name, buf[j].v);
	    }
	    cur = cur->next;
	}
    }
    return name_val_pair(0);
}

parse_ercode context::pop(val_ind* ptr) {
    if (stack_ptr == 0)
	return E_EMPTY_STACK;
    //find the name of the item on top of the stack and look up the current index
    size_t ti = buf[stack_ptr-1].i;
    name_ind* cur = table[ti];
    name_ind* prev = NULL;
    while (cur) {
	if (cur->ind == stack_ptr-1) {
	    if (prev)
		prev->next = cur->next;
	    else
		table[ti] = cur->next;
	    if (ptr)
		*ptr = buf[--stack_ptr];
	    delete cur;
	    break;
	}
    }
    return E_SUCCESS;
}

/**
 * remove the top n items from the stack. If n is too large then no action is taken and an E_EMPTY_STACK error is returned
 */
parse_ercode context::pop_n(size_t n) {
    if (n > stack_ptr)
	return E_EMPTY_STACK;
    for (size_t i = 0; i < n; ++i) {
	pop(NULL);
    }
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
	if (er != E_SUCCESS) return sto;
	if (it_list.type != VAL_ARRAY && it_list.type != VAL_LIST) { er = E_BAD_TYPE;return sto; }
	//the prototype expression needs to be null terminated
	for_start[0] = 0;
	//we need to add a variable with the appropriate name to loop over, we force push because we pop later
	set_value(var_name, sto, true);
	n_els = it_list.n_els;
	lbuf = (value*)calloc(n_els, sizeof(value));
	//we now iterate through the list specified, substituting VAL in the expression with the current value
	for (size_t i = 0; i < it_list.n_els; ++i) {
	    if (it_list.type == VAL_LIST)
		set_value(var_name, it_list.val.l[i]);
	    else
		set_value(var_name, make_val_num(it_list.val.a[i]));
	    char* expr_name = strndup(start+1, expr_len);
	    lbuf[i] = parse_value(expr_name, er);
	    free(expr_name);
	    if (er != E_SUCCESS) { n_els = i;return sto; }
	}
	//we need to remove the the variable we loop over
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
	if (inst_val.type != VAL_INST) {
	    printf("Error: tried to lookup from non instance type\n");
	    er = E_BAD_TYPE;return sto; }
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
	sto.type = VAL_NUM;
	sto.n_els = 1;
	sto.val.x = (tmp_l == tmp_r);
    } else if (term_char == '>') {
	if (tmp_l.type != VAL_NUM || tmp_r.type != VAL_NUM) { cleanup_val(&tmp_l);cleanup_val(&tmp_r);er = E_BAD_VALUE;return sto; }
	if (str[i+1] == '=') {
	    sto.val.x = (tmp_l.val.x >= tmp_r.val.x)? 1: 0;
	    //++i;
	} else {
	    sto.val.x = (tmp_l.val.x > tmp_r.val.x)? 1: 0;
	}
    } else if (term_char == '<') {
	if (tmp_l.type != VAL_NUM || tmp_r.type != VAL_NUM) { cleanup_val(&tmp_l);cleanup_val(&tmp_r);er = E_BAD_VALUE;return sto; }
	if (str[i+1] == '=') {
	    sto.val.x = (tmp_l.val.x <= tmp_r.val.x)? 1: 0;
	    //++i;
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
	if (tmp_r.val.x == 0) {
	    printf("Error: division by zero.\n");
	    cleanup_val(&tmp_l);
	    cleanup_val(&tmp_r);
	    er = E_NAN;
	    return sto;
	}//TODO: return a nan?
	if (tmp_l.type == VAL_NUM && tmp_r.type == VAL_NUM) {
	    sto.val.x = tmp_l.val.x / tmp_r.val.x;
	} else if (tmp_r.type == VAL_NUM && tmp_l.type == VAL_MAT) {	
	    sto.type = VAL_MAT;
	    sto.val.m = new mat3x3(*tmp_r.val.m / tmp_r.val.x);
	}
    } else if (term_char == '^') {
	if (tmp_l.type != VAL_NUM || tmp_r.type != VAL_NUM) { er = E_NAN;return sto; }//TODO: return a nan?
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
	} else if (str[i] == '/') {
	    //ignore everything after a comment
	    if (str[i+1] == '/') break;
	    //TODO: add support for multiline comments
	}

	if (nest_level == 0) {
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
	    for (size_t j = 0; j < sto.n_els-1; ++j)
		sto.val.s[j] = str[first_open_ind+j+1];
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
		if (tmp_lst.type != VAL_LIST) {
		    printf("Error: tried to index from non list type\n");
		    cleanup_val(&tmp_lst);er = E_BAD_TYPE;return sto; }
		str[last_close_ind] = 0;
		value contents = parse_value(str+first_open_ind+1, er);
		str[last_close_ind] = /*[*/']';
		//check that we found the list and that it was valid
		if (contents.type != VAL_NUM) {
		    printf("Error: only integers are valid indices\n");
		    cleanup_val(&contents);er = E_BAD_TYPE;return sto; }
		//now figure out the index
		long tmp = (long)(contents.val.x);
		if (tmp < 0) tmp = tmp_lst.n_els+tmp;
		if (tmp < 0) {
		    printf("Error: index %d is out of range for list of size %lu.\n", tmp, tmp_lst.n_els);
		    er = E_OUT_OF_RANGE;
		    return sto;
		}
		size_t ind = (size_t)tmp;
		if (ind >= tmp_lst.n_els) {
		    printf("Error: index %d is out of range for list of size %lu.\n", tmp, tmp_lst.n_els);
		    er = E_OUT_OF_RANGE;
		    return sto;
		}
		return copy_val(tmp_lst.val.l[ind]);
	    }
	} else if (str[first_open_ind] == '{' && str[last_close_ind] == '}') {
	    //now parse the argument as a context
	    size_t n_els;
	    str[last_close_ind] = 0;
	    char** list_els = csv_to_list(str+first_open_ind+1, ',', &n_els, er);
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
		if (er != E_SUCCESS) { free(list_els);sto.type=VAL_UNDEF;cleanup_val(&sto);return sto; }
		sto.val.c->set_value(cur_name, tmp);
		cleanup_val(&tmp);
	    }
	    free(list_els);
	    sto.type = VAL_INST;
	    return sto;
	} else if (str[first_open_ind] == '(' && str[last_close_ind] == ')') {
	    //check to see if this is a function call (it is if there are any non-whitespace characters before the open paren
	    for (size_t j = 0; j < first_open_ind; ++j) {
		if (str[j] != ' ' && str[j] != '\t' && str[j] != '\n') {
		    cgs_func tmp_f;
		    //check to see if this is a function declaration, this must be handled differently
		    if (strncmp(str+j, "fun", KEY_DEF_LEN) == 0 &&
		    (is_whitespace(str[j+KEY_DEF_LEN]) || str[j+KEY_DEF_LEN] == '('/*)*/)) {
			const char* f_end;
			tmp_f = parse_func(str, first_open_ind, er, &f_end, true);
			//find the contents in the curly brace and separate by semicolons
			line_buffer b(f_end, ';');
			sto.type = VAL_FUNC;
			sto.n_els = tmp_f.n_args;
			sto.val.f = new user_func(tmp_f, b.get_enclosed(rs.pos, &end, '{', '}', func_end));
			er = E_SUCCESS;
			cleanup_func(&tmp_f);
			return sto;
		    } else {
			//we can't leave this as zero in case the user needs to do some more operations
			char term_char = str[last_close_ind+1];
			str[last_close_ind+1] = 0;
			tmp_f = parse_func(str, first_open_ind, er);
			if (er != E_SUCCESS) { return sto; }
			//otherwise lookup the function
			value func_val = lookup(tmp_f.name);
			if (func_val.type == VAL_FUNC) {
			    //make sure that the function was found and that sufficient arguments were provided
			    if (func_val.n_els <= tmp_f.n_args) {
				sto = func_val.val.f->eval(this, tmp_f, er);
			    } else {
				printf("Error: unrecognized function name %s\n", tmp_f.name);
				er = E_LACK_TOKENS;
			    }
			} else {
			    printf("Error: unrecognized function name %s\n", tmp_f.name);
			    er = E_BAD_TYPE;
			}
			cleanup_func(&tmp_f);
			return sto;
		    }
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
 * Easily inspect values using const strings
 */
value context::parse_value(const char* tok) {
    parse_ercode er;
    char* tmp = strdup(tok);
    value ret = parse_value(tmp, er);
    free(tmp);
    if (er != E_SUCCESS) {
	cleanup_val(&ret);
	ret = make_val_undef();
    }
    return ret;
}

/*
 * helper function for parse_value which appends at most n characters from the string str to line while dynamically resizing the buffer if necessary
 * line: the line to save to
 * line_off: the current end of line
 * line_size: the size in memory allocated for line
 * n: the maximum number of characters to write
 */
char* append_to_line(char* line, size_t* line_off, size_t* line_size, const char* str, size_t n) {
    if (!line_off || !line_size || !str || n == 0) return line;
    size_t ls = *line_size;
    size_t i = *line_off;
    if (!line || i + n >= ls) {
	ls = 2*i + n + 1;
	line = (char*)xrealloc(line, ls);
    }
    size_t j = 0;
    for (; j < n; ++j) {
	line[i+j] = str[j];
	//terminate on reaching the string end
	if (!str[j]) {
	    ++j;
	    break;
	}
    }
    *line_off = i+j;
    *line_size = ls;
    return line;
}

/**helper function for read_from lines that reads a single line
 */
parse_ercode context::read_single_line(context::read_state& rs) {
    parse_ercode er = E_SUCCESS;
    //tracking variables
    size_t buf_off = 0;
    size_t len = rs.b.get_line_size(rs.pos.line);
    size_t k = 0;
    bool started = false;
    char* lval = NULL;
    size_t rval_ind = 0;
    line_buffer_ind init = rs.pos;
    for (rs.pos.off = 0;; ++rs.pos.off) {
	//make sure the buffer is large enough
	if (k >= rs.buf_size) {
	    rs.buf_size *= 2;
	    rs.buf = (char*)xrealloc(rs.buf, sizeof(char)*rs.buf_size);
	}
	//exit the loop when we reach the end, but make sure to include whatever parts haven't already been included
	if (rs.pos.off >= len || rs.b.get(rs.pos) == 0) {
	    rs.buf[k] = 0;
	    break;
	}
	//ignore preceeding whitespace
	if (rs.b.get(rs.pos) != ' ' || rs.b.get(rs.pos) != '\t' || rs.b.get(rs.pos) != '\n')
	    started = true;
	if (started) {
	    //handle comments
	    if (rs.b.get(rs.pos) == '/' && rs.pos.off > 0 && rs.b.get(rs.pos-1) == '/') {
		//we don't care about lines that only contain comments, so we should skip over them, but in the other event we need to skip to the end of the line
		if (k == 1)
		    started = false;
		else
		    rs.pos.off = rs.b.get_line_size(rs.pos.line);
		//terminate the expression and move to the next line
		rs.buf[--k] = 0;
		break;
	    } else if (rs.b.get(rs.pos) == '*' && rs.pos.off > 0 && rs.b.get(rs.pos-1) == '/') {
		if (k == 1)
		    started = false;
		rs.buf[--k] = 0;//set the slash to be a null terminator
		while (rs.b.inc(rs.pos)) {
		    if (rs.b.get(rs.pos) == '*' && rs.b.get(rs.pos+1) == '/') {
			rs.pos = rs.pos+2;
			break;
		    }
		}
	    }
	    //handle assignments
	    if (rs.b.get(rs.pos) == '=') {
		rs.buf[k++] = 0;
		lval = CGS_trim_whitespace(rs.buf, NULL);
		init.off = k;
		//buf_off = k;
		rval_ind = k;
		continue;//don't copy the value into rs.buf
	    }
	    //if we encounter a block, then we need to make sure we include all of its contents, even if that block ends on another line
	    char match_tok = 0;
	    switch(rs.b.get(rs.pos)) {
		case '(': match_tok = ')';break;
		case '[': match_tok = ']';break;
		case '{': match_tok = '}';break;
		case '\"': match_tok = '\"';break;
		case '\'': match_tok = '\'';break;
		default: break;
	    }
	    if (match_tok) {
		line_buffer_ind end;
		line_buffer enc = rs.b.get_enclosed(rs.pos, &end, rs.b.get(rs.pos), match_tok, true, true);
		char* tmp = enc.flatten(' ');
		//if we're still on the same line, then we need to continue until we reach the end. Otherwise, save everything inclosed and terminate
		if (end.line == init.line && end.off < len) {
		    rs.buf = append_to_line(rs.buf, &k, &rs.buf_size, tmp, end.off - rs.pos.off);
		    rs.pos = end;
		    init = rs.pos;
		    free(tmp);
		    //continue;//don't copy the value into rs.buf
		} else {
		    rs.buf = append_to_line(rs.buf, &k, &rs.buf_size, tmp, strlen(tmp)+1);
		    free(tmp);
		    rs.pos = end;
		    break;//we're done reading this line since we jumped across a line
		}
	    }
	    //handle everything else
	    rs.buf[k++] = rs.b.get(rs.pos);
	}
    }
    //only set the rval if we haven't done so already
    if (started) {
	value tmp_val = parse_value(rs.buf+rval_ind, er);
	if (er != E_SUCCESS) return er;
	place_value(lval, tmp_val);
    } else {
	rs.pos.line += 1;
	rs.pos.off = 0;
    }
    return E_SUCCESS;
}

/**
 * Generate a context from a list of lines. This context will include function declarations, named variables, and subcontexts (instances).
 * lines: the array of lines to read from
 * n_lines: the size of the array
 * returns: an errorcode if one was found or E_SUCCESS on success
 */
parse_ercode context::read_from_lines(const line_buffer& b) {
    parse_ercode er = E_SUCCESS;

    context::read_state rs(b);
    //iterate over each line in the file
    while (true) {
	char* line = b.get_line(rs.pos.line);
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
		//now we actually create the function
		rs.pos.off = 0;
		line_buffer_ind end;
		user_func tmp(cur_func, b.get_enclosed(rs.pos, &end, '{', '}', func_end));
		if (end.line == b.get_n_lines()) { free(line);return E_BAD_SYNTAX; }
		value v;
		v.type = VAL_FUNC;
		v.val.f = &tmp;
		v.n_els = b.get_n_lines()-rs.pos.line;
		set_value(cur_func.name, v);
		//we have to advance to the line after end of the function declaration
	    }
	} else {
	    er = read_single_line(rs);
	    if (er != E_SUCCESS) { free(line);return er; }
	}
	free(line);
	//if we're at the end of a line, try incrementing. If that doesn't work, then we've reached the end of the file.
	if (rs.pos.off >= rs.b.get_line_size(rs.pos.line)) {
	    if (!b.inc(rs.pos))
		break;
	}
    }

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
user_func::user_func(value (*p_exec)(context*, cgs_func, parse_ercode&)) : call_sig() {
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
value user_func::eval(context* c, cgs_func call, parse_ercode& er) {
    value sto;
    if (exec) {
	value ret = (*exec)(c, call, er);
	return ret;
    } else if (call.n_args == call_sig.n_args) {
	//setup a new scope with function arguments defined
	context func_scope(c);
	for (size_t i = 0; i < call_sig.n_args; ++i) {
	    func_scope.set_value(call_sig.arg_names[i], call.args[i]);
	}
	size_t lineno = 0;
    } else {
	er = E_LACK_TOKENS;
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
	free(line);
    }
    //return an undefined value by default
    value ret;ret.type = VAL_UNDEF;ret.val.x = 0;return ret;*/
}
