#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "section.h"

#define BASE_BLOCK_SIZE 4

/* Checks if key exists in Section. If it does returns location of entry 
 * otherwise returns -1
 */
static
int
data_location(Section *s, char *key)
{
    int i;
    for(i = 0; i < s->nentries; i++)
        if(strcmp(key, s->entries[i]->key) == 0)
            return i;
    return -1;
}

Section*
new_Section(char *name)
{
    Section *s;
    s = (Section*) malloc(sizeof(Section));
    s->section_name = strdup(name);
    s->nentries = 0;
    s->entries = (Section_Entry**)malloc(sizeof(Section_Entry)*BASE_BLOCK_SIZE);
    s->curr_alloc = BASE_BLOCK_SIZE;
    return s;
}

void
delete_Section(Section *s)
{
    int i;
    for(i = 0; i < s->nentries; i++)
    {
        free(s->entries[i]->key);
        free(s->entries[i]->value);
        free(s->entries[i]);
    }
    free(s->entries);
    free(s->section_name);
    free(s);
}

/* Inserts data (key, value) in Section. If data already exists nothing is done.
 */
void 
insert_data_in_section(Section *s, char *key, char *value)
{
    Section_Entry *se;
    int loc;
    loc = data_location(s, key);
    if (loc == -1)
    {
        se = (Section_Entry*) malloc(sizeof(Section_Entry));
        se->key = strdup(key);
        se->value = strdup(value);
        if (s->nentries == s->curr_alloc)
        {
            s->entries = (Section_Entry**) realloc(s->entries, 2*s->curr_alloc*sizeof(Section_Entry*));
            s->curr_alloc *= 2;
        }
        s->entries[s->nentries] = (Section_Entry*)malloc(sizeof(Section_Entry));
        s->entries[s->nentries++] = se;
    }
}
/* Returns data associated with key in section if it exists. If key does not
 * exist then NULL is returned
 */
char* 
lookup_in_section(Section *s, char *key)
{
    int loc;
    loc = data_location(s, key);
    if (loc != -1)
        return s->entries[loc]->value;
    else
        return NULL;
}

Chapter*
new_Chapter(char *name)
{
    Chapter *c;
    c = (Chapter*)malloc(sizeof(Chapter));
    c->chapter_name = strdup(name);
    c->nsections = 0;
    c->sections = (Section**)malloc(sizeof(Section*)*BASE_BLOCK_SIZE);
    c->curr_alloc = BASE_BLOCK_SIZE;
    return c;
}

void
delete_Chapter(Chapter *c)
{
    int i;
    for(i = 0; i < c->nsections; i++)
        delete_Section(c->sections[i]);
    free(c->sections);
    free(c->chapter_name);
    free(c);
}

Section* 
lookup_in_chapter(Chapter *c, char *section_name)
{
    int i;
    for(i = 0; i < c->nsections; i++)
        if(strcmp(c->sections[i]->section_name, section_name) == 0)
            return c->sections[i];
    return NULL;
}

static
void
copy_section(Section *src, Section *dest)
{
    int i;
    dest->section_name = strdup(src->section_name);
    dest->nentries = src->nentries;
    dest->curr_alloc = src->curr_alloc;
    dest->entries = (Section_Entry**)malloc(sizeof(Section_Entry*)*dest->curr_alloc);
    for(i = 0; i < dest->nentries; i++)
    {
        dest->entries[i] = (Section_Entry*)malloc(sizeof(Section_Entry));
        dest->entries[i]->key = strdup(src->entries[i]->key);
        dest->entries[i]->value = strdup(src->entries[i]->value);
    }
}

void 
insert_section_in_chapter(Chapter *c, Section *s)
{

    if (c->nsections == c->curr_alloc)
    {
        c->sections = (Section**) realloc(c->sections, 
                                          (c->curr_alloc+BASE_BLOCK_SIZE)*sizeof(Section*));
        c->curr_alloc += BASE_BLOCK_SIZE;
    }
    c->sections[c->nsections] = (Section*)malloc(sizeof(Section));
    copy_section(s, c->sections[c->nsections]);
    c->nsections++;
}
