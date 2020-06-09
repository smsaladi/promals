#!/usr1/local/bin/python

import os, sys, re

# change working directory here
work_dir = ''
if not work_dir: work_dir = os.getcwd()

# options
args = sys.argv

input_file = args[1]
outfile = re.sub("\.fa$", "", input_file)
outfile += ".mummals.aln"

ss = 3
solv = 1
unaligned = 1

id_thr = 0.6

for i in range(2,len(args)):
    if args[i] == '-ss':
	ss = int(args[i+1])
	continue
    if args[i] == '-solv':
	solv = int(args[i+1])
	continue
    
    if args[i] == '-unaligned':
	unaligned = int(args[i+1])
	continue

    if args[i] == '-outfile':
	outfile = args[i+1]
	continue

    if args[i] == '-id_thr':
	id_thr = args[i+1]
	continue

if ss not in [1, 3]:
    print "Error in option -ss:"
    print "secondary structure option (ss) must be 1 or 3"
    sys.exit(0)

if solv not in [1, 2, 3]:
    print "Error in option -solv:"
    print "solvent accessibility option (solv) must be 1, 2 or 3"
    sys.exit(0)

if unaligned not in [0, 1]:
    print "Error in option -unaligned:"
    print "unaligned option must be 0 or 1"
    sys.exit(0)

param_dir = work_dir + "/hmm_parameters"

if unaligned:
    param_file = "%s/dataset_0.20_0.40_0.60_abcd.dali.solv%d_ss%d.mat" %(param_dir, solv, ss)
else:
    param_file = "%s/dataset_0.20_0.40_0.60_abcd.dali0.solv%d_ss%d.mat" %(param_dir, solv, ss)

command = "%s/mummals %s -ss %d -solv %d -unaligned %d -param %s -outfile %s -id_thr %s" %(work_dir, input_file, ss, solv, unaligned, param_file, outfile, id_thr)

print command
print

os.system(command)

