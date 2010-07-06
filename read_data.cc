/*
 * Author: Haibao Tang <bao@uga.edu> May 10, 2007
 *
 * Data input module for mcscan, contains several procedures
 * Read blast output file, formatted by -m8, is a bunch of hits
 * Read MCL cluster, which is retrieved by clustering the above blast file
 * Read GFF file, which includes the chromosome, position information
 */


#include "read_data.h"

// incremental sorting y coord
static bool cmp_y (const Score_t& t1, const Score_t& t2)
{
    return t1.y < t2.y ||
           (t1.y == t2.y && t1.x < t2.x);
}

// incremental sorting e-value
static bool cmp_ev (const Score_t& t1, const Score_t& t2)
{
    return t1.score < t2.score;
}

// filter the blast -m8 output by the following threshold:
// lexically sorted, gene #1 < gene #2
// non-self blast match
// both be present in the mcl output file and in the same group
void read_blast(const char *prefix_fn, bool gff_flag=true)
{
    char fn[LABEL_LEN], g1[LABEL_LEN], g2[LABEL_LEN];
    double score;
    Blast_record br;
    int i;

    sprintf(fn, "%s.blast", prefix_fn);
    FILE *fp = mustOpen(fn, "r");

    int pair_id = 0;
    int total_num = 0;
    map<string, Gene_feat>::iterator it1, it2;
    Gene_feat *gf1, *gf2;
    while ( fscanf(fp, "%s%s%lg",
                   &g1[0], &g2[0], &score)==3 )
    {
        total_num++;
        // swap lexically and ignore self match
        i = strcmp(g1, g2);
        if (i < 0)
        {
            br.gene1.assign(g1);
            br.gene2.assign(g2);
        }
        else if (i > 0)
        {
            // interchange
            br.gene1.assign(g2);
            br.gene2.assign(g1);
        }
        else continue;  // bug fixed by bao, May 22nd 2009
        it1 = gene_map.find(br.gene1);
        it2 = gene_map.find(br.gene2);
        if (it1==gene_map.end() || it2==gene_map.end()) continue;
        gf1 = &(it1->second), gf2 = &(it2->second);

        // assert both has the same MCL node id
        br.node = gf1->node;
        if (gff_flag && br.node != gf2->node) continue;

        if (gf1->mol.empty() || gf2->mol.empty()) continue;
        br.mol_pair = gf1->mol+"&"+gf2->mol;
        mol_pairs[br.mol_pair]++;

        br.pair_id = pair_id++;
        br.score = score;
        match_list.push_back(br);
    }

    int selected_num = match_list.size();
    progress("%d matches imported (%d discarded)",
             selected_num, total_num - selected_num);

    fclose(fp);
}

void read_mcl(const char *prefix_fn)
{
    char delims[] = " \t\r\n";
    char *atom = NULL;
    char fn[LABEL_LEN];

    sprintf(fn, "%s.mcl", prefix_fn);
    FILE *fp = mustOpen(fn, "r");

    int node_num = 0;
    size_t n = 0;
    char *line;
    map<string, Gene_feat>::iterator it;
    Gene_feat *gf;
    while (getline(&line, &n, fp)>=0)
    {
        atom = strtok(line, delims);
        while (atom != NULL)
        {
            if ((it=gene_map.find(string(atom))) != gene_map.end())
            {
                gf = &(it->second);
                gf->node = node_num;
                if (!gf->mol.empty()) chr_map[gf->mol].insert(gf);
            }
            atom = strtok(NULL, delims);
        }
        node_num++;
    }

    fclose(fp);
}

void read_bed(const char *prefix_fn)
{
    char fn[LABEL_LEN], gn[LABEL_LEN], mol[LABEL_LEN];
    int end5, end3;
    Gene_feat gf;
    // default position for genes are based on gene ranks
    vector<Gene_feat> bed;

    sprintf(fn, "%s.bed", prefix_fn);
    FILE *fp = mustOpen(fn, "r");

    while (fscanf(fp, "%s%d%d%s",
                  &mol[0], &end5, &end3, &gn[0]) == 4)
    {
        gf.mol = string(mol);
        gf.name = string(gn);
        gf.mid = end5;
        bed.push_back(gf);
    }

    fclose(fp);
    
    // sort bed with respect to chromosome and position
    sort(all(bed));

    vector<Gene_feat>::iterator bi;
    unsigned int i = 0;
    tr(bed, bi) 
    { 
        if (! USE_BP) bi->mid = i++; 
        gene_map[bi->name] = *bi;
        //printf("%s\n", bi->name.c_str());
    }
}

static void filter_matches_x ()
{
    // match_bin is a list of records that are potentially repetitive
    vector<Score_t> match_bin, score_cpy;
    vector<Score_t>::const_iterator it, prev_rec;

    sort(score.begin(), score.end());
    prev_rec = it = score.begin();
    it++;
    match_bin.push_back(*(prev_rec));
    for (; it != score.end(); it++)
    {
        // scan whether it has a linking window with previous one
        if ((prev_rec->x != it->x) ||
                (it->y - prev_rec->y) > OVERLAP_WINDOW)
        {
            // record last match_bin, take only least e-value
            score_cpy.push_back(*min_element(match_bin.begin(),
                                             match_bin.end(), cmp_ev));
            // start a new match_bin
            match_bin.clear();
        }
        match_bin.push_back(*it);
        prev_rec = it;
    }
    // don't forget the last match_bin
    score_cpy.push_back(*min_element(match_bin.begin(),
                                     match_bin.end(), cmp_ev));
    match_bin.clear();

    // copy into score
    score.clear();
    score = score_cpy;
    score_cpy.clear();
}

static void filter_matches_y ()
{
    // match_bin is a list of records that are potentially repetitive
    vector<Score_t> match_bin, score_cpy;
    vector<Score_t>::const_iterator it, prev_rec;

    sort(score.begin(), score.end(), cmp_y);
    prev_rec = it = score.begin();
    it++;
    match_bin.push_back(*(prev_rec));
    for (; it != score.end(); it++)
    {
        // scan whether it has a linking window with previous one
        if ((prev_rec->y != it->y) ||
                (it->x - prev_rec->x) > OVERLAP_WINDOW)
        {
            // record last match_bin, take only least e-value
            score_cpy.push_back(*min_element(match_bin.begin(),
                                             match_bin.end(), cmp_ev));
            // start a new match_bin
            match_bin.clear();
        }
        match_bin.push_back(*it);
        prev_rec = it;
    }
    // don't forget the last match_bin
    score_cpy.push_back(*min_element(match_bin.begin(),
                                     match_bin.end(), cmp_ev));
    match_bin.clear();

    // copy into score
    score.clear();
    score = score_cpy;
    score_cpy.clear();
}

// feed into dagchainer
void feed_dag(const string &mol_pair)
{
    // two additional filters will be applied here
    // best hsp (least e-value)
    // non-repetitive in a window of 50kb region
    vector<Blast_record>::const_iterator it;
    Score_t cur_score;

    for (it = match_list.begin(); it < match_list.end(); it++)
    {
        if (it->mol_pair != mol_pair) continue;

        cur_score.pairID = it->pair_id;
        cur_score.x = gene_map[it->gene1].mid;
        cur_score.y = gene_map[it->gene2].mid;
        cur_score.score = MATCH_SCORE;

        score.push_back(cur_score);
    }

    // sort by both axis and remove redundant matches within
    // a given window length (default 50kb)
    filter_matches_x();
    filter_matches_y();

    dag_main(score, mol_pair);
}

