/*
 * Author: Haibao Tang <bao@uga.edu> May 17, 2007
 *
 * Simple output subroutines, controls what to print in the final .blocks file.
 * Also contains some helper functions for other data structures.
 */

#include "out_utils.h"

void print_params(FILE *fw)
/* print parameters */
{
    fprintf( fw, "############### Parameters ###############\n");
    fprintf( fw, "# MATCH_SCORE: %d\n", MATCH_SCORE );
    fprintf( fw, "# MATCH_SIZE: %d\n", MATCH_SIZE );
    fprintf( fw, "# UNIT_DIST: %d\n", UNIT_DIST );
    fprintf( fw, "# GAP_SCORE: %d\n", GAP_SCORE );
    fprintf( fw, "# OVERLAP_WINDOW: %d\n", OVERLAP_WINDOW );
    fprintf( fw, "# EXTENSION_DIST: %d\n", EXTENSION_DIST );
    fprintf( fw, "# E_VALUE: %g\n", E_VALUE );
    fprintf( fw, "# PIVOT: %s\n", PIVOT.c_str() );
    fprintf( fw, "##########################################\n\n");
}
void print_align(FILE* fw)
/* print alignment */
{
    int i, j, pid;
    int nseg = seg_list.size(), nanchor;
    Seg_feat *s;

    print_params(fw);

    for (i=0; i<nseg; i++)
    {
        s = &seg_list[i];
        nanchor = s->pids.size();
        fprintf(fw, "## Alignment %d: score=%.1f e_value=%.2g N=%d %s %s\n",
                i, s->score, s->e_value, nanchor, s->mol_pair.c_str(),
                s->sameStrand?"plus":"minus");
        for (j=0; j<nanchor; j++)
        {
            pid = s->pids[j];
            fprintf(fw, "%3d-%3d:\t%s\t%s\t%7.1g\n",
                    i, j, match_list[pid].gene1.c_str(),
                    match_list[pid].gene2.c_str(), match_list[pid].score);
        }
    }
}

void print_align_mcl(FILE* fw)
/* sometimes we wish to print to a simple three-column file for mcl clustering */
{
    int i, j, pid;
    int nseg = seg_list.size(), nanchor;
    Seg_feat *s;
    for (i=0; i<nseg; i++)
    {
        s = &seg_list[i];
        nanchor = s->pids.size();
        for (j=0; j<nanchor; j++)
        {
            pid = s->pids[j];
            fprintf(fw, "%s\t%s\t%.1g\n",
                    match_list[pid].gene1.c_str(),
                    match_list[pid].gene2.c_str(), match_list[pid].score);
        }
    }
}

void print_geneSet(FILE *fw, const geneSet &g)
/* helper function to print out geneSet(*/
{
    geneSet::const_iterator i=g.begin();
    if (g.empty()) fprintf(fw, ".");
    for (; i!=g.end(); i++)
    {
        if (i!=g.begin()) fprintf(fw, ";");
        fprintf(fw, "%s", (*i)->name.c_str());
    }
}

void print_POG_memory(FILE *fw, const POG_order &ref, int block)
/* print verbose info about ref for debugging */
{
    POG_order::const_iterator it;
    set<POG_node *>::const_iterator p;
    int j = 0;
    for (it=ref.begin(); it!=ref.end(); it++)
    {
        fprintf(fw, "%3d-%4d:\t", block, j++);
        print_geneSet(fw, (*it)->master_genes);
        fprintf(fw, "\t");
        if ((*it)->fusion.empty()) fprintf(fw, ".");
        for (p=(*it)->fusion.begin(); p!=(*it)->fusion.end(); p++)
        {
            print_geneSet(fw, (*p)->genes);
            fprintf(fw, "|");
        }
        fprintf(fw, "\t[%p]\t", (void*)*it);
        for (p=(*it)->next.begin(); p!=(*it)->next.end(); p++)
            fprintf(fw, "`%p", (void*)*p);
        fprintf(fw, "\n");
    }
}

void print_POG_block(FILE *fw, const POG_order &ref, int block, int cols)
/* multiple blocks output */
{
    POG_order::const_iterator it;
    set<POG_node *>::const_iterator p;
    vector<POG_node *> v(cols);
    int j = 0, k;
    for (it=ref.begin(); it!=ref.end(); it++)
    {
        fprintf(fw, "%3d-%4d:\t", block, j++);
        print_geneSet(fw, (*it)->master_genes);
        /* by now the columns for the syntenic region has been assigned */
        for (k=0; k<cols; k++) v[k] = NULL;
        for (p=(*it)->fusion.begin(); p!=(*it)->fusion.end(); p++)
            v[(*p)->r->col] = *p;
        for (k=0; k<cols; k++)
        {
            fprintf(fw, "\t");
            if (v[k] == NULL) fprintf(fw, ".");
            else print_geneSet(fw, v[k]->genes);
        }
        fprintf(fw, "\n");
    }
}

