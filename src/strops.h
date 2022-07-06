#include <cstring>



/**
  * Remove the whitespace surrounding a word
  */
inline char* trim_whitespace(char* str, size_t* len) {
    char* start = NULL;
    _uint last_non = 0;
    for (_uint i = 0; str[i] != 0; ++i) {
	if (str[i] != ' ' && str[i] != '\t' && str[i] != '\n') {
	    if (start)
		last_non = i;
	    else
		start = str + i;
	}
    }
    str[last_non+1] = 0;
    if (len) *len = last_non+1;
    return start;
}


