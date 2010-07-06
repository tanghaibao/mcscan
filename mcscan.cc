/*
 * Author: Haibao Tang <bao@uga.edu> May 10, 2007
 * Main entry point for the executable mcscan
*/

#include "mcscan.h"

static bool IS_PAIRWISE;
static bool BUILD_MCL;
static char prefix_fn[LABEL_LEN];


const char *argp_program_version = "MCSCAN 0.8";
const char *argp_program_bug_address = "<bao@uga.edu>";

/* Program documentation */
static char doc[] = "MCSCAN -- multiple collinearity scan"
                    " (compiled "__DATE__" "__TIME__")\n\n"
                    "Reference:\n"
                    " Tang,H., Wang,X., Bowers,J.E., Ming,R., Alam,M., Paterson,A.H.\n"
                    " Unraveling ancient hexaploidy through"
                    " multiply-aligned angiosperm gene maps\n"
                    " Genome Research (2008) 18, 1944-1954.\n";

/*  A description of the arguments we accept. */
const unsigned int nargs = 1;
static char args_doc[] = "prefix_fn";
static char *args[nargs];

/* The options we understand. */
static struct argp_option options[] =
{
    {"match_score", 'k', "MATCH_SCORE", 0,
        "final score=MATCH_SCORE+NUM_GAPS*GAP_SCORE" },
    {"gap_score", 'g', "GAP_SCORE", 0, "gap penalty" },
    {"e_value", 'e', "E_VALUE", 0, "alignment significance" },
    {"pivot", 'p', "PIVOT", 0,
     "PIVOT is the reference genome, make it two letter prefix in"\
     "your .bed file, everything else will be aligned to the reference" },
    {"unit_dist", 'u', "UNIT_DIST", 0, "average intergenic distance" },
    {0, 'A', 0, 0, "use base pair dist instead of gene ranks" },
    {0, 'a', 0, 0, "only builds the pairwise blocks (.aligns file)" },
    {0, 'b', 0, 0, "limit within genome synteny (e.g. Vv-Vv) mapping" },
    { 0 }
};

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
    /* Get the input argument from argp_parse, which we
      know is a pointer to our arguments structure. */
    switch (key)
    {
    case 'k':
        MATCH_SCORE = atoi(arg);
        break;
    case 'g':
        GAP_SCORE = atoi(arg);
        break;
    case 's':
        MATCH_SIZE = atoi(arg);
        break;
    case 'e':
        E_VALUE = atof(arg);
        break;
    case 'p':
        PIVOT = string(arg);
        break;
    case 'u':
        UNIT_DIST = atoi(arg);
        break;
    case 'a':
        IS_PAIRWISE = true;
        break;
    case 'b':
        IN_SYNTENY = true;
        break;
    case 'A':
        USE_BP = true;
        break;

    case ARGP_KEY_ARG:
        if (state->arg_num >= nargs)
            /* Too many arguments. */
            argp_usage (state);

        args[state->arg_num] = arg;

        break;
    case ARGP_KEY_END:
        if (state->arg_num < nargs)
            /* Not enough arguments. */
            argp_usage (state);

        break;

    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}


/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

static int read_opt (int argc, char **argv)
{
    /* Default values. */
    // match bonus, final score=MATCH_SCORE+GAPS*GAP_SCORE
    MATCH_SCORE = 40;
    // the number of genes required to call synteny, sometimes more
    MATCH_SIZE = 5;
    // gap extension penalty
    GAP_SCORE = -2;
    // alignment significance
    E_VALUE = 1e-5;
    // align with a reference genome (occurs as first column in .blocks file)
    PIVOT = "ALL";
    UNIT_DIST = 0;

    IS_PAIRWISE = false;
    BUILD_MCL = false;
    IN_SYNTENY = false;
    USE_BP = false;

    /* Parse our arguments; every option seen by parse_opt will
      be reflected in arguments. */
    argp_parse (&argp, argc, argv, 0, 0, 0);
    strcpy(prefix_fn, args[0]);


    // default unit values for the distance calculation
    if (USE_BP) 
    {
        if (UNIT_DIST==0) UNIT_DIST = 10000;
    }
    else
    {
        if (UNIT_DIST==0) UNIT_DIST = 2;
    }

    OVERLAP_WINDOW = MATCH_SCORE*UNIT_DIST/10;
    EXTENSION_DIST = MATCH_SCORE*UNIT_DIST/2;
    CUTOFF_SCORE = MATCH_SCORE*MATCH_SIZE;

    return 0;
}


int main(int argc, char *argv[])
{
    /* Start the timer */
    uglyTime(NULL);

    char align_fn[LABEL_LEN], block_fn[LABEL_LEN];
    FILE *fw;

    read_opt(argc, argv);
    print_params(stdout);

    read_bed(prefix_fn);
    if (!IS_PAIRWISE) read_mcl(prefix_fn);
    read_blast(prefix_fn);

    sprintf(align_fn, "%s.aligns", prefix_fn);
    fw = mustOpen(align_fn, "w");

    progress("%d pairwise comparisons", (int) mol_pairs.size());

    map<string, int>::const_iterator ip;
    for (ip=mol_pairs.begin(); ip!=mol_pairs.end(); ip++)
    {
        if (ip->second >= MATCH_SIZE) feed_dag(string(ip->first));
    }

    progress("%d alignments generated", (int) seg_list.size());
    if (BUILD_MCL) print_align_mcl(fw);
    else print_align(fw);

    fclose(fw);
    uglyTime("Pairwise synteny written to %s", align_fn);

    if (IS_PAIRWISE) return 0;

    sprintf(block_fn, "%s.blocks", prefix_fn);
    fw = mustOpen(block_fn, "w");

    POG_main(fw);

    fclose(fw);
    uglyTime("Multiple synteny written to %s", block_fn);

    return 0;
}

