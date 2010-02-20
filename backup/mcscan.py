#!/usr/bin/env python

import sys
import os
import os.path as op
from time import asctime
from subprocess import call # py2.4

def sh(cmd, info):
	print cmd
	print "==> %s [%s]"%(info, asctime(localtime()))
	call(cmd, shell=True)


if __name__ == '__main__':
    
    from optparse import OptionParser
    usage = "%prog prefix_name \n" \
            "- <prefix_name> should agree with your prefix.blast and prefix.gff"

    parser = OptionParser(usage)
    (options, args) = parser.parse_args()

    try: 
        prefix = sys.argv[1]
    except:
        print usage; sys.exit(1)

    if not (op.exists("%s.blast" % prefix) and op.exists("%s.gff" % prefix)):
        print >>sys.stderr, \
                "Please make sure that you have both %s.blast and %s.gff in the same folder." % (prefix, prefix)
        sys.exit(1)

    # markov clustering, make sure you have mcl in your current directory
    cur_dir = op.dirname(__file__)
    if not os.path.exists("mcl"):
        print >>sys.stderr, \
                "Please make sure that you have mcl in the current directory"
        sys.exit(1)

    if not op.exists("%s.mcl" % prefix):
        sh("more %s.blast | %s/mcl - --abc --abc-neg-log -abc-tf 'mul(0.4343), ceil(200)' -o %s.mcl" % (prefix, cur_dir, prefix), "MCL clustering")

    # multiple alignment
    sh("./mcscan %s"%prefix, "Multiple alignment")
