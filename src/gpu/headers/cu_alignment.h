#include "cu_util.h"

#ifndef C_CU_ALIGNMENT
#define C_CU_ALIGNMENT

typedef struct {
    char* cigar;
    int end_ref;
    int end_read;
    int score;
} cu_alignment;

#endif