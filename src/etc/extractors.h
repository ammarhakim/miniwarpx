#ifndef __extractors__2004__
#define __extractors__2004__

/* Given string representation of int (double) returns int (double) */
int extract_int(char *val);
double extract_double(char *val);

/* Given comma seperate sting representation of int (double) returns
 * int (double) array
 */
double* extract_double_p(char *val);
int *extract_int_p(char *val);

#endif // __extractors__2004__
