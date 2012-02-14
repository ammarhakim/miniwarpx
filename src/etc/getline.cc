#include<ctype.h>
#include<errno.h>
#include<stdlib.h>
#include<stdio.h>

#define	CHUNK	128

static 
int
mygetc(FILE* fp)
{
    int	c;
    static int last = -1;
    
    if (last != -1)
    {
        c = last;
        last = -1;
        return c;
    }
		
    if ((c = getc(fp)) == '\r')
    {
        if ((c = getc(fp)) != '\n')
            last = c;
		
        return('\n');
    }

    return(c);
}

char *
xgetline(FILE* fp, char *buf, int *len)
{
    size_t sz = CHUNK;		/* this keeps track of the current size of buffer */
    size_t i = 0;		/* index into string tracking current position */
    char *ptr;			/* since we may set buf to NULL before returning */
    int	c;			/* to store getc return */

    /* no real good idea for handling this */
    if (!fp)
        return(buf);

    /* start out with buf set to CHUNK + 2 bytes */
    if ((buf = (char*)realloc(buf, CHUNK + 2)) == NULL)
        return(NULL);

    while ((c = mygetc(fp)) != EOF && c != '\n')
    {
        /* the following needed in case we are in cbreak or raw mode */
        if (c != '\b')
            buf[i++] = c;
        else if (i)
            i--;
        
            /* check for buffer overflow */
        if (i >= sz)
            if ((buf = (char*)realloc(buf, (sz += CHUNK) + 2)) == NULL)
                return(NULL);
    }
    
    /* is there anything to return? */
    if (c == EOF && !i)
    {
        free(buf);
        return(NULL);
    }
    *len = i; /* lenght of line excluding newline character */
    buf[i++] = 0;	/* yes I want the ++ */

    /* test for error but don't bother explaining if it fails */
    if ((ptr = (char*)realloc(buf, i)) == NULL)
        ptr = buf;

    return(ptr);
}
// #include <string.h>
// int
// main(void)
// {
//     char *p;
//     FILE *fp;
//     int len;
//     char a[80], b[80];

//     fp = fopen("test.inp", "r");
//     while((p = getline(fp, 0, &len)) != NULL)
//     {
        
//     }

//     return(0);
// }
