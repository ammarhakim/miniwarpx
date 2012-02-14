#ifndef __inpparse__
#define __inpparse__

#include "section.h"

Chapter* parse_input(char *file);
void cleanup_input(Chapter *c);

#endif // __inpparse__
