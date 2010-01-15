/*
 * Author: Haibao Tang <bao@uga.edu> May 10, 2007
 * Main entry point for the executable mcscan
*/

#include "mcscan.h"

static bool IS_PAIRWISE;
static bool BUILD_MCL;
static char prefix_fn[LABEL_LEN];

static void print_banner()
{
    progress("\nMCSCAN %.1f: multiple collinearity scan\n"
             " (compiled "__DATE__" "__TIME__")\n\n"
             "Reference:\n"
             " Tang,H., Wang,X., Bowers,J.E., Ming,R., Alam,M., Paterson,A.H.\n"
             " Unraveling ancient hexaploidy through"
             " multiply-aligned angiosperm gene maps\n"
             " Genome Research (2008) 18, 1944-1954.\n", VER);
}

static void init_opt()
{
    // match bonus, final score=MATCH_SCORE+GAPS*GAP_SCORE
    MATCH_SCORE = 50;
    // the number of genes required to call synteny, sometimes more
    MATCH_SIZE = 6;
    // gap extension penalty
    GAP_SCORE = -3;
    // alignment significance
    E_VALUE = 1e-5;
    // align with a reference genome (occurs as first column in .blocks file)
    PIVOT = "ALL";
    // this variable is dependent on gene density
    UNIT_DIST = 10000;

    IS_PAIRWISE = false;
    BUILD_MCL = false;
    IN_SYNTENY = false;
}

static void print_help(const char *prg)
{
    progress("[Usage] %s prefix_fn [options]\n"
             " -k  MATCH_SCORE, final score=MATCH_SCORE+NUM_GAPS*GAP_SCORE\n"
             "     (default: %d)\n"
             " -g  GAP_SCORE, gap penalty (default: %d)\n"
             " -s  MATCH_SIZE, number of genes required to call synteny\n"
             "     (default: %d)\n"
             " -e  E_VALUE, alignment significance (default: %g)\n"
             " -p  PIVOT, pivot is the reference genome, make it two letter prefix\n"
             "     your .gff file, everything else will be aligned to the reference\n"
             "     (default: %s)\n"
             " -u  UNIT_DIST, average intergenic distance (default: %d)\n"
             " -a  only builds the pairwise blocks (.aligns file)\n"
             " -b  limit within genome synteny (e.g. Vv-Vv) mapping\n"
             " -h  print this help page\n",
             prg, MATCH_SCORE, GAP_SCORE, MATCH_SIZE, E_VALUE, PIVOT.c_str(), UNIT_DIST);
    exit(1);
}

static void read_opt(int argc, char *argv[])
{
    int c;
    opterr = 0;

    if (argc < 2) print_help(argv[0]);

    while ((c = getopt(argc, argv, "k:g:s:e:p:u:abch")) != -1)
        switch (c)
        {
        case 'k':
            MATCH_SCORE = atoi(optarg);
            break;
        case 'g':
            GAP_SCORE = atoi(optarg);
            break;
        case 's':
            MATCH_SIZE = atoi(optarg);
            break;
        case 'e':
            E_VALUE = atof(optarg);
            break;
        case 'p':
            PIVOT = string(optarg);
            break;
        case 'u':
            UNIT_DIST = atoi(optarg);
            break;
        case 'a':
            IS_PAIRWISE = true;
            break;
        case 'b':
            IN_SYNTENY = true;
            break;
        case 'c':
            BUILD_MCL = IS_PAIRWISE = true;
            break;
        case '?':
            if (optopt=='k' || optopt=='s' || optopt=='g' || optopt=='e' || optopt=='p')
                errAbort("Option -%c requires an argument.", optopt);
            else if (isprint (optopt))
                errAbort("Unknown option `-%c'.", optopt);
            else
                errAbort("Unknown option character `\\x%x'.", optopt);
        default:
            print_help(argv[0]);
            break;
        }

    if (optind==argc) errAbort("Please enter your input file");
    else strcpy(prefix_fn, argv[optind]);

    OVERLAP_WINDOW = MATCH_SCORE*UNIT_DIST/10;
    EXTENSION_DIST = MATCH_SCORE*UNIT_DIST/2;
    CUTOFF_SCORE = MATCH_SCORE*MATCH_SIZE;
}

int main(int argc, char *argv[])
{
    /* Start the timer */
    uglyTime(NULL);

    print_banner();
    char align_fn[LABEL_LEN], block_fn[LABEL_LEN];
    FILE *fw;

    init_opt();
    read_opt(argc, argv);

    read_gff(prefix_fn);
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

