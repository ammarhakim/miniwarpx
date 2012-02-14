#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "extractors.h"

int
extract_int(char *val) { return atoi(val); }

double
extract_double(char *val) { return atof(val); }

int*
extract_int_p(char *val)
{
    char *cp1, *cp2;
    char *token, delimiter[] = ",";
    int count;
    int *a;

    cp1 = strdup(val);
    cp2 = strdup(val);
    // count how many entries are in array
    token = strtok(cp1, delimiter);
    if(token == NULL)
        return NULL;
    count = 1;
    while(1)
    {
        token = strtok(NULL, delimiter);
        if(token == NULL)
            break;
        count++;
    }
    a = (int*) malloc(sizeof(int)*count);
    // copy these into a
    token = strtok(cp2, delimiter);
    count = 0;
    a[count] = atoi(token);
    while(1)
    {
        token = strtok(NULL, delimiter);
        if(token == NULL)
            break;
        count++;
        a[count] = atoi(token);
    }
    free(cp1);
    free(cp2);
    return a;
}


double*
extract_double_p(char *val)
{
    char *cp1, *cp2;
    char *token, delimiter[] = ",";
    int count;
    double *a;

    cp1 = strdup(val);
    cp2 = strdup(val);
    // count how many entries are in array
    token = strtok(cp1, delimiter);
    if(token == NULL)
        return NULL;
    count = 1;
    while(1)
    {
        token = strtok(NULL, delimiter);
        if(token == NULL)
            break;
        count++;
    }
    a = (double*) malloc(sizeof(double)*count);
    // copy these into a
    token = strtok(cp2, delimiter);
    count = 0;
    a[count] = atof(token);
    while(1)
    {
        token = strtok(NULL, delimiter);
        if(token == NULL)
            break;
        count++;
        a[count] = atof(token);
    }
    free(cp1);
    free(cp2);
    return a;
}
