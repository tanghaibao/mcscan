#ifndef __POG_H
#define __POG_H

#include "basic.h"

void POG_main(FILE *fw);

// out_utils
extern void print_params(FILE *fw);
extern void print_geneSet(FILE *fw, const geneSet &g);
extern void print_POG_memory(FILE *fw, const POG_order &ref, int block);
extern void print_POG_block(FILE *fw, const POG_order &ref, int block, int col);

#endif
