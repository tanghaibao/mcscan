#!/usr/bin/env python

from time import *
import sys
import os
try:
	from subprocess import call # py2.4
except:
	def call(cmd, **kwargs):
		os.system(cmd)

def log(cmd, info):
	print cmd
	print "==> %s [%s]"%(info, asctime(localtime()))
	call(cmd, shell=True)

usage = """
[USAGE] python mcscan.py prefix_name
- <prefix_name> should agree with your prefix.blast and prefix.gff
"""

try: 
	prefix = sys.argv[1]
except:
	print usage; sys.exit(1)

if not (os.path.exists("%s.blast"%prefix) and os.path.exists("%s.gff"%prefix)):
	print "Please make sure that you have both %s.blast and %s.gff in the same folder."%\
            (prefix, prefix); sys.exit(1)

# markov clustering, make sure you have mcl in your current directory
cur_dir = os.getcwd()
if not os.path.exists("mcl"):
    print "Please make sure that you have mcl in the current directory"; sys.exit(1)

if not os.path.exists("%s.mcl"%prefix):
    log("more %s.blast | %s/mcl - --abc --abc-neg-log -abc-tf 'mul(0.4343), ceil(200)' -o %s.mcl"%\
            (prefix, cur_dir, prefix), "MCL clustering")

# multiple alignment
log("./mcscan %s"%prefix, "Multiple alignment")
