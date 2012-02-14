#ifndef __section__
#define __section__

struct Section_Entry
{
    char *key;
    char *value;
};
typedef struct Section_Entry Section_Entry;

struct Section
{
    char *section_name;
    int nentries;
    Section_Entry **entries;
    int curr_alloc;
};
typedef struct Section Section;

struct Chapter
{
    char *chapter_name;
    int nsections;
    Section **sections;
    int curr_alloc;
};
typedef struct Chapter Chapter;

Section* new_Section(char *name);
void delete_Section(Section *s);
void insert_data_in_section(Section *s, char *key, char *value);
char* lookup_in_section(Section *s, char *key);

Chapter* new_Chapter(char *name);
void delete_Chapter(Chapter *c);
void insert_section_in_chapter(Chapter *c, Section *s);
Section* lookup_in_chapter(Chapter *c, char *section_name);

#endif // __section__
