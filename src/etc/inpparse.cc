#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include "inpparse.h"

extern char * xgetline(FILE *fp, char *buff, int *len);

static
void
copy_substring(char *dest, char *src, int s, int e)
{
    strncpy(dest, src+s, e-s);
    dest[e-s] = '\0';
}

static
regex_t*
compile_re(char *pattern)
{
    int err_no;
    regex_t *p;

    // allocate memory for compiled regular expressions
    p = (regex_t*)malloc(sizeof(regex_t));
    memset(p, 0, sizeof(regex_t));
    // compile regular expression
    err_no = regcomp(p, pattern, REG_EXTENDED);
    // if compile failed print error message and return NULL
    if(err_no != 0)
    {
        size_t length; 
        char *buffer;
        length = regerror (err_no, p, NULL, 0);
        buffer = (char*) malloc(length);
        regerror (err_no, p, buffer, length);
        fprintf(stderr, "%s\n", buffer);
        free(buffer);
        regfree(p);
        return NULL;
    }
    return p;
}

Chapter*
parse_input(char *file)
{
    // pattern1 matches "key = value" lines
    char *pattern1 = "([[:alnum:]_-]+)[[:blank:]]*=[[:blank:]]*(.+)";
    // pattern2 matches section headers
    char *pattern2 = "\[[[:alnum:]_-]+]";
    char *line, *sect, *key, *value;
    int res, len, no_sub1, no_sub2;
    regex_t *p1, *p2;
    regmatch_t *m1, *m2;
    Chapter *c;
    Section *s;
    FILE *fp;

    // compile patterns
    p1 = compile_re(pattern1);
    p2 = compile_re(pattern2);
    if((p1 == NULL) || (p2 == NULL))
    {
        printf("Unable to compile regular expression\n");
        return NULL;
    }
    // allocate memory for match object
    no_sub1 = p1->re_nsub+1;
    no_sub2 = p2->re_nsub+1;
    m1 = (regmatch_t*)malloc(sizeof(regmatch_t)*no_sub1);
    m2 = (regmatch_t*)malloc(sizeof(regmatch_t)*no_sub2);

    // open file for parsing
    fp = fopen(file, "r");
    if(fp == NULL)
    {
        printf("Unable to open input file %s for reading\n", file);
        regfree(p1);  free(p1);
        regfree(p2); free(p2);
        free(m1); free(m2);
        return NULL;
    }
    
    c = new_Chapter(file);
    s = NULL;
    // loop over file reading lines and parsing them
    while((line = xgetline(fp, 0, &len)) != NULL)
    {
        // match line against regular expressions
        res = regexec(p1, line, no_sub1, m1, 0);
        if( res == 0) // matched a data line
        {
            // get key
            len = m1[1].rm_eo - m1[1].rm_so;
            key = (char*) malloc(len+1);
            copy_substring(key, line, m1[1].rm_so, m1[1].rm_eo);
            // get value
            len = m1[2].rm_eo - m1[2].rm_so;
            value = (char*) malloc(len+1);
            copy_substring(value, line, m1[2].rm_so, m1[2].rm_eo);
            // insert data into current section
            insert_data_in_section(s, key, value);
            free(key);
            free(value);
        }
        else
        {
            res = regexec(p2, line, no_sub2, m2, 0);
            if( res == 0) // matched a section header
            {
                // add current section into chapter
                if(s)
                {
                    insert_section_in_chapter(c, s);
                    delete_Section(s);
                }
                // get section name
                len = m2[0].rm_eo - m2[0].rm_so - 2;
                sect = (char*) malloc(len+1);
                copy_substring(sect, line, m2[0].rm_so+1, m2[0].rm_eo-1);
                s = new_Section(sect);
                free(sect);
            }
        }// if line is neither section or data it is ignored
    }
    insert_section_in_chapter(c, s);
    // deallocate memory and close file handle
    regfree(p1);
    free(p1);
    regfree(p2);
    free(p2);
    free(m1);
    free(m2);

    fclose(fp);
    return c;
}
