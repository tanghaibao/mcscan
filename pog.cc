/*
 * Author: Haibao Tang <bao@uga.edu>, Jan.15, 2008
 *
 * Performs partial order graph alignment of gene orders
 * this offers an improvement over the consensus method in versions <0.8
 *
 * First, consecutive tandems on each chromosome are merged
 * then for each reference chromosome, the syntenic region is added
 * sequentially, beginning from the best-scoring, this is considered
 * as re-alignment of gene orders, but utilizes the partial order graph
 * data structure
 */

#include "pog.h"

static POG_order ref, ref_slave, master, slave;
static vector<POG_node *> memory_pool;
static DPVec v, track_v;
static vector<Syn_region> Q;
static vector<End_point> endpoints;
static Syn_region *syn;
static int cols;

/* Comparator for sorting */
bool synCmp (const Syn_region &a, const Syn_region &b)
{
    return a.s->score > b.s->score;
}

static bool check_self_genome(const string &s)
/* whether a mol_pair is self comparison, e.g. "Vv1&Vv14" */
{
    int pos = s.find('&');
    return s.substr(0, 2) == s.substr(pos+1, 2);
}

static void init_POG(POG_order &g, const geneSet &s)
/* convert geneSet to POG_order by merging consecutive tandems*/
{
    g.clear();
    geneSet::const_iterator i = s.begin();
    POG_node *t = new POG_node;
    memory_pool.push_back(t);
    t->node = (*i)->node;
    t->master_genes.insert(*i);
    g.push_back(t);

    for (i++; i!=s.end(); i++)
    {
        if ((*i)->node != t->node)
        {
            t = new POG_node;
            memory_pool.push_back(t);
            t->node = (*i)->node;
            g.push_back(t);
        }
        t->master_genes.insert(*i);
    }
}

static void link_POG(POG_order &g)
/* populate the directed edge in the graph */
{
    POG_order::const_iterator i = g.begin();
    POG_node *t = *i;
    for (i++; i!=g.end(); i++)
    {
        t->next.insert(*i);
        t = *i;
    }
}

static void init_synteny(POG_order &g, POG_order &t,
                         Gene_feat *a, Gene_feat *b)
/* collects all the genes in range [*a, *b] */
{
    POG_node *p;
    POG_order::const_iterator it=t.begin();
    g.clear();
    for (; it!=t.end(); it++)
    {
        p = *it;
        if (p->node == a->node &&
                p->master_genes.find(a) != p->master_genes.end()) break;
    }
    for (; it!=t.end(); it++)
    {
        p = *it;
        g.push_back(p);
        if (p->node == b->node &&
                p->master_genes.find(b) != p->master_genes.end()) break;
    }
}
static void init_master(Gene_feat *a, Gene_feat *b)
/* master version - collects all the genes in range [*a, *b] */
{
    POG_order &g = master;
    printf(" search between %s - %s\n", a->name.c_str(), b->name.c_str());

    init_synteny(g, ref, a, b);

    printf(" master contains %d elements.\n", (int)master.size());
}

static void init_slave(Gene_feat *a, Gene_feat *b, bool sameStrand)
/* slave version - collects all the genes in range [*a, *b] */
{
    POG_order &g = slave;
    init_POG(ref_slave, chr_map[a->mol]);
    printf(" search between %s - %s\n", a->name.c_str(), b->name.c_str());

    init_synteny(g, ref_slave, a, b);
    if (!sameStrand) reverse(g.begin(), g.end());
    link_POG(slave);

    /* slave regions do not require master_genes */
    POG_order::iterator j = g.begin();
    for (; j!=g.end(); j++)
    {
        (*j)->genes = (*j)->master_genes;
        (*j)->master_genes.clear();
        (*j)->fusion.insert(*j);
        (*j)->r = syn;
    }

    printf(" slave contains %d elements.\n", (int)slave.size());
}

static void refresh_POG(POG_order &g)
/* helper function to use with DFS */
{
    POG_order::iterator i = g.begin();
    for (; i!=g.end(); i++) (*i)->visited = false;
}

static void DFS(POG_node *src, POG_node *target,
                int score, int &max_score)
/* performs a depth-first-search to get distance between src and target */
{
    /* prune when it is impossible to get more better score  */
    if (score <= max_score) return;
    else if (src == target)
    {
        max_score = score;
        return;
    }
    set<POG_node *>::const_iterator p;
    for (p=src->next.begin(); p!=src->next.end(); p++)
    {
        if ((*p)->visited) continue;
        (*p)->visited = true;
        DFS(*p, target, score+GAP_SCORE, max_score);
    }
}

static void fuse_POG_node(POG_node *t, POG_node *g)
/* dump the genes in fused POG node for gene retrieval */
{
    t->fusion.insert(g);
}

static void fuse_POG()
/* backtracking through the best path and fuse aligned POGs */
{
    int n=v.size(), max_score = 0, max_i = -1, i;
    for (i=0; i<n; i++)
    {
        if (v[i].score > max_score)
        {
            max_score = v[i].score;
            max_i = i;
        }
    }
    //for (i=0; i<n; i++) printf("%d ", v[i].score); puts("");
    printf(" best pog path score %d\n", max_score);
    if (max_score < CUTOFF_SCORE) return;
    syn->score = max_score;

    track_v.clear();
    DP *a = &v[max_i], *b;
    while (a != NULL)
    {
        track_v.push_back(*a);
        a = a->from;
    }
    reverse(track_v.begin(), track_v.end());

    n = track_v.size();
    POG_order::iterator is, it, ix, iy;
    for (i=0; i<n-1; i++)
    {
        a = &track_v[i];
        b = &track_v[i+1];
        is = find(ref.begin(), ref.end(), a->s);
        is++;
        it = find(ref.begin(), ref.end(), b->s);
        ix = find(slave.begin(), slave.end(), a->t);
        ix++;
        iy = find(slave.begin(), slave.end(), b->t);

        /* insert the interleaved slave nodes before master nodes */
        ref.insert(is, ix, iy);

        /* fix links and merge nodes */
        fuse_POG_node(a->s, a->t);
        fuse_POG_node(b->s, b->t);

        /* make sure the interleaved slave nodes are non-empty */
        if (ix != iy)
        {
            a->s->next.insert(*ix);
            iy--;
            (*iy)->next.clear();
            (*iy)->next.insert(b->s);
        }
        else
            a->s->next.insert(b->s);
    }

    /* for .blocks layout */
    End_point ep;
    ep.s = syn;
    ep.a = track_v.begin()->s;
    ep.start = true;
    endpoints.push_back(ep);
    ep.s = syn;
    ep.a = track_v.rbegin()->s;
    ep.start = false;
    endpoints.push_back(ep);
}

static void align_POG()
/* core algorithm, one dimensional dynamic programming */
{
    POG_order::const_iterator i=master.begin(), j=slave.begin();
    //for (; i!=master.end(); i++) printf("%d ", (*i)->node); puts("");
    //for (; j!=slave.end(); j++) printf("%d ", (*j)->node); puts("");

    DP p;
    /* collect matching nodes for sparse dynamic programming */
    for (i=master.begin(); i!=master.end(); i++)
    {
        int node = (*i)->node;
        for (j=slave.begin(); j!=slave.end(); j++)
        {
            if (node == (*j)->node)
            {
                p.from = NULL, p.s = *i, p.t = *j, p.score = MATCH_SCORE;
                v.push_back(p);
            }
        }
    }
    /* distances between matches are computed and plugged in formula*/
    int n=v.size(), aa, bb, del_x, del_y, del;
    DP *a, *b;
    for (aa=0; aa<n; aa++)
    {
        a = &v[aa];
        for (bb=aa+1; bb<n; bb++)
        {
            b = &v[bb];
            if (a->s == b->s || a->t == b->t ) continue;
            refresh_POG(master), refresh_POG(slave);

            del_x = del_y = -MATCH_SCORE;
            DFS(a->s, b->s, MATCH_SCORE, del_x);
            if (del_x == -MATCH_SCORE) break;
            DFS(a->t, b->t, MATCH_SCORE, del_y);
            if (del_y == -MATCH_SCORE) continue;
            del = a->score + MIN(del_x, del_y);

            if (del > b->score)
            {
                b->score = del;
                b->from = a;
            }
        }
    }
    fuse_POG();

    v.clear();
}

static void cluster_POG(const string &mol)
/* collect threaded alignments from dagchainer and re-align */
{
    bool match1, match2;
    int n = seg_list.size(), i;

    Syn_region r;
    Seg_feat *s;

    for (i=0; i<n; i++)
    {
        s = &seg_list[i];
        match1 = mol==s->s1->mol;
        match2 = mol==s->s2->mol;
        r.s = s;
        if (match1)
        {
            r.match1 = true;
            Q.push_back(r);
        }
        if (match2)
        {
            r.match1 = false;
            Q.push_back(r);
        }
    }
    n = Q.size();
    sort(Q.begin(), Q.end(), synCmp);

    for (i=0; i<n; i++)
    {
        syn = &Q[i];
        s = syn->s;
        if (IN_SYNTENY && check_self_genome(s->mol_pair)) continue;
        printf(" original dagchainer score %.1f\n", s->score);
        if (syn->match1)
        {
            init_master(s->s1, s->t1);
            init_slave(s->s2, s->t2, s->sameStrand);
            align_POG();
        }
        else
        {
            init_master(s->s2, s->t2);
            init_slave(s->s1, s->t1, s->sameStrand);
            align_POG();
        }
    }
}

static void layout_POG()
{
    POG_order::const_iterator it;
    vector<End_point>::iterator ip;
    int i = 0;

    for (ip=endpoints.begin(); ip!=endpoints.end(); ip++)
    {
        for (i=0, it=ref.begin(); it!=ref.end(); it++, i++)
        {
            if (*it == ip->a)
            {
                ip->ref_index = i;
                break;
            }
        }
    }
    sort(endpoints.begin(), endpoints.end());

    ip = endpoints.begin();
    /* assign columns for each syntenic region */
    priority_queue<int, vector<int>, greater<int> > pool;
    cols = 0;
    for (it=ref.begin(); it!=ref.end(); it++)
    {
        if (ip == endpoints.end()) break;
        while (*it == ip->a)
        {
            if (ip->start)
            {
                if (pool.empty())
                {
                    ip->s->col = cols++;
                }
                else
                {
                    ip->s->col = pool.top();
                    pool.pop();
                }
            }
            else
            {
                pool.push(ip->s->col);
            }
            ip++;
        }
    }
}

void POG_main(FILE *fw)
{
    map<string, geneSet >::const_iterator it;
    vector<POG_node *>::iterator iq;
    int i=0;

    print_params(fw);

    for (it=chr_map.begin(); it!=chr_map.end(); it++)
    {
        if (PIVOT!="ALL" && it->first.find(PIVOT)==string::npos) continue;
        if ((int)it->second.size() < MATCH_SIZE) continue;

        string query = it->first;
        init_POG(ref, it->second);
        link_POG(ref);

        printf("## pivot %s contains %d tandem clusters\n",
               query.c_str(), (int)ref.size());
        cluster_POG(query);
        fprintf(fw, "## View %d: pivot %s\n", i, query.c_str());

        //print_POG_memory(fw, ref, i);
        layout_POG();
        print_POG_block(fw, ref, i, cols);

        fprintf(fw, "\n");

        ref.clear(), ref_slave.clear(), master.clear(), slave.clear();
        /* release the memory held by partial order graph */
        for (iq=memory_pool.begin(); iq!=memory_pool.end(); iq++)
            delete *iq;
        memory_pool.clear();
        Q.clear(), endpoints.clear();
        i++;
    }
}

