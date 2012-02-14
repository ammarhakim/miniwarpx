#include <stdio.h>

extern char * xgetline(FILE *fp, char *buff, int *len);

void
copy_file(char *src, char *dest)
{
    int len;
    char *line;
    FILE *fps, *fpd;

    // open source file for reading
    fps = fopen(src,"r");
    // open desination file for writing
    fpd = fopen(dest,"w");

    // loop over source file, getting each line and writing it to the
    // destination file
    while((line = xgetline(fps, 0, &len)) != NULL)
        fprintf(fpd, "%s\n", line);
    
    // close the files
    fclose(fps); fclose(fpd);
}
