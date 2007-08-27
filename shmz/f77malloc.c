#include <bigmalloc.h>

void *malloc_(int *sz)
{
    void *ptr;
    ptr = Malloc(*sz);
    return ptr;
}

void free_(void **ptr)
{
    Free(*ptr);
}
