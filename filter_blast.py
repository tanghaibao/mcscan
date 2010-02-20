#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Filter the 3-column BLAST output according to following rules
1. Remove self matches, multiple matches
2. Re-order gene pairs lexicographically

"""

import os
import sys

if __name__ == '__main__':
        
    from optparse import OptionParser

    usage = "%prog infile outfile"
    parser = OptionParser(usage)
    (options, args) = parser.parse_args()

    try:
        infile = args[0]
        outfile = args[1]
    except:
        sys.exit(parser.print_help())
    
    assert os.path.exists(infile)

    fp = file(infile)
    pairs = {}
    j = 0
    for row in fp:
        j += 1
        a, b, e = row.split()
        e = float(e)
        if a==b: continue
        if a > b: a, b = b, a
        pair_name = (a, b)
        if pair_name not in pairs or (pair_name in pairs and e<pairs[pair_name]): 
            pairs[pair_name] = e

    print >>sys.stderr, j, "records read"
    print >>sys.stderr, len(pairs), "records after filtering"

    fw = file(outfile,"w")
    for k, e in sorted(pairs.items()):
        a, b = k
        fw.write("%s\t%s\t%g\n" % (a, b, e))

    fp.close()
    fw.close()
