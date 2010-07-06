#ifndef __OUT_UTILS_H
#define __OUT_UTILS_H

#include "basic.h"

/* pairwise blocks */
void print_align(FILE *fw);
void print_align_mcl(FILE *fw);

/* multiple blocks */
void print_geneSet(FILE *fw, const geneSet &g);
void print_params(FILE *fw);
void print_POG_memory(FILE *fw, const POG_order &ref, int block);
void print_POG_block(FILE *fw, const POG_order &ref, int block, int col);

#endif
