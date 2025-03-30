#include "cu_util.h"

#ifndef C_CU_FQ_READ
#define C_CU_FQ_READ

typedef struct {
    char* id;
    char* seq;
    char* qual;
    char* metadata;
    size_t size;
} cu_fq_read;

#endif