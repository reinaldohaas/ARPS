#include <sys/types.h>
#include <sys/stat.h>

#ifdef UNDERSCORE
#define c_mkdir c_mkdir_
#endif

void c_mkdir(const char *path_name, int *status )
{
    *status = mkdir(path_name, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

/*    printf("Status: %d\n", status); */

    return; 
}
