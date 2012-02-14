#include <stdio.h>
#include <unistd.h>
#include <fnmatch.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include <errno.h>

#include "miniwarpx.h"

/* Creates output directory to write data. Cleans up directory 
   if it already exists. Returns 1 if initialization  succeeded
   or -1 otherwise
 */
int 
init_output(Run_Data& ri)
{
    char *dir1;
    int len, res, curr_err;
    char *pattern, *temp;
    DIR *dir;
    struct dirent *entry;
    
    len = strlen(ri.run_name) + 10;
    dir1 = (char*) malloc(len);
    // create directory to store output, if it does not already exist
    sprintf(dir1, "./%s", ri.run_name);
    res = mkdir(dir1, S_IRWXU | S_IRWXG | S_IRWXO);
    curr_err = errno; // copy this so that unexpected things don't occur
    if(res == 0)
    {
        free(dir1);
        return 1;
    }
    else
    {
        // directory creation failed. See what really happened
        if(curr_err == ENOENT)
        {
            printf("Directory path in input file incorrect or output directory does not exist\n");
            return -1;
        }
        else if(curr_err == EACCES)
        {
            printf("Permision denied to create directory %s\n", dir1);
            return -1;
        }
        else if(curr_err == EEXIST) // directory exists
        {
            printf("Directory %s exists\n", dir1);
            printf("Deleting contents...\n");
            // delete all frame.* files
            pattern = (char*)malloc(strlen("frame")+5); // for globbing
            temp = (char*) malloc(len+20);
            sprintf(pattern, "frame.*");
            dir = opendir(dir1);
            if(dir)
            {
                while((entry = readdir(dir)))
                {
                    if(fnmatch(pattern, entry->d_name, FNM_NOESCAPE) == 0)
                    {
                        sprintf(temp, "%s/%s", dir1, entry->d_name);
                        res = remove(temp);
                        if(res == -1)
                            printf("Unable to delete %s\n", entry->d_name);
                    }
                }
            }
            else
            {
                printf("Unable to clean up %s\n", dir1);
                return -1;
            }
            free(temp);
            free(pattern);
            free(dir1);
        }
        else
        {
            free(dir1);
            return -1;
        }
    }
    return 1;
}
