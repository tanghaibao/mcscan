#ifndef __BASIC_H
#define __BASIC_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <cstdio>
#include <cstring>
#include <cstdarg>
#include <cerrno>
#include <climits>
#include <cmath>
#include <ctime>
#include <sys/time.h>
#include <argp.h>
#include <string>
#include <vector>
#include <list>
#include <queue>
#include <set>
#include <map>
#include <algorithm>

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

/***** Useful macros *****/

#define sameString(a, b) (strcmp((a), (b))==0)
#define MAX(a, b) ((a)>(b)?(a):(b))
#define MIN(a, b) ((a)<(b)?(a):(b))
#define LABEL_LEN 256

#define all(c) (c).begin(),(c).end() 
#define tr(c,i) for(i=(c).begin();i!=(c).end();i++)

/* Compatibility of __attribute__ with non-GNU */
#ifndef __GNUC__
#define __attribute__(x)
#endif

/***** Data structures *****/
struct Blast_record
{
    string gene1, gene2;
    string mol_pair;
    int pair_id;
    int node;
    double score;
};
struct Gene_feat
{
    string name;
    string mol;
    int mid;
    int node;
    bool operator < (const Gene_feat &g) const
    {
        return (mol == g.mol && mid < g.mid) || mol < g.mol;
    }
};
struct geneCmp
{
    bool operator() (const Gene_feat *a, const Gene_feat *b) const
    {
        return (a->mol == b->mol && a->mid < b->mid) || a->mol < b->mol;
    }
};
typedef vector<Gene_feat *> geneVec;
typedef set<Gene_feat *, geneCmp> geneSet;

struct Seg_feat
{
    vector<int> pids;
    Gene_feat *s1, *t1, *s2, *t2;
    double score, e_value;
    string mol_pair;
    bool sameStrand;
};

// dagchainer
struct Cell_t
{
    float raw;
int score :
    30;
unsigned from :
    2;
};
struct Score_t
{
    int pairID;  // identifier of match pair
    int x, y;  // x,y coordinates
    float score;
    bool operator< (const Score_t & node) const
    {
        return  (x < node.x || (x == node.x && y < node.y));
    }
};
struct Path_t
{
    float score;
    int rc;  // sum of row and column of last entry
    int sub;
};

// pog (partial order graph)
struct Syn_region
{
    Seg_feat *s;
    int col;
    int score;
    bool match1;
};

struct POG_node
{
    geneSet master_genes, genes;
    set<POG_node *> fusion, next;
    int node;
    bool visited;
    Syn_region *r;
};
typedef list<POG_node *> POG_order;

struct DP
{
    POG_node *s, *t;
    DP *from;
    int score;
};
typedef vector<DP > DPVec;

struct End_point
{
    Syn_region *s;
    bool start;
    POG_node *a;
    int ref_index;
    bool operator < (const End_point &g) const
    {
        return ref_index < g.ref_index ||
               (ref_index==g.ref_index && s->score > g.s->score);
    }
};

/***** All data *****/
// use gene name to search its node, mol, mid
extern map<string, Gene_feat> gene_map;
extern vector<Blast_record> match_list;
extern vector<Seg_feat> seg_list;
extern map<string, int> mol_pairs;
extern map<string, geneSet > chr_map;

/***** CONSTANTS *****/
// match bonus
extern int MATCH_SCORE;
extern int MATCH_SIZE;
// gap extension penalty
extern int GAP_SCORE;
// length for a single gap in basepairs
extern int GAP_SIZE;
// The filter window for linking locally repetitive hits
extern int OVERLAP_WINDOW;
// significance cut-off
extern double E_VALUE;
// reference genome
extern string PIVOT;
// intergenic distance
extern int UNIT_DIST;
// segment extension limit
extern int EXTENSION_DIST;
// alignment significance score
extern int CUTOFF_SCORE;
extern bool IN_SYNTENY;
// use base pair distance rather than gene ranks
extern bool USE_BP;

// direction in the 2d dynamic matrix
enum { DIAG, UP, LEFT, DEL };

/***** Helper functions (Some from James Kent library) *****/
void progress(const char *format, ...)
/* Print progress message */
__attribute__((format(printf, 1, 2)));

void err(const char *format, ...)
/* Print error message but do not exit */
__attribute__((format(printf, 1, 2)));

void warn(const char *format, ...)
/* Print error message but do not exit */
__attribute__((format(printf, 1, 2)));

void errAbort(const char *format, ...)
/* Print error message to stderr and exit */
__attribute__((noreturn, format(printf, 1, 2)));

long clock1000();
/* A millisecond clock. */

void uglyTime(const char *label, ...)
/* Print label and how long it's been since last call.  Call with
 * a NULL label to initialize. */
__attribute__((format(printf, 1, 2)));

FILE *mustOpen(const char *fileName, const char *mode);
/* Open a file or die */

#endif
