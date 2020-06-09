#!/usr1/local/bin/python

Document = "\
\n\n\
This python script runs mummals with the same options as the c++\n\
executable of mummals, except that the parameter file name does\n\
not need to be specified.\n\
\n\
usage:\n\
      python run_mummals.py input_fasta_file [options]\n\
\n\
options:\n\
         -ss: secondary structure types [1,3]\n\
         -solv: solvent accessibility categories [1,2,3]\n\
         -unaligned: number of distinct unaligned match state [0,1]\n\
         -outfile: output file name\n\
\n\
example:\n\
      python run_mummals.py tmp1.fa -ss 3 -solv 1 -unaligned 1 -outfile tmp1.mummals.aln\n\
\n\
";


import os, sys, re

# change working directory here
work_dir = ''
if not work_dir: work_dir = os.getcwd()

# options
args = sys.argv

if len(args)<2: 
	print Document
	sys.exit(0)

input_file = args[1]
if input_file[0] == '-': 
	print Document
	sys.exit(0)

if not os.path.isfile(input_file):
	print "input fasta file does not exist:", input_file
	sys.exit(0)

outfile = re.sub("\.fa$", "", input_file)
outfile += ".mummals.aln"

ss = 3
solv = 1
unaligned = 1

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

if not os.path.isfile(param_file):
	print "Error: parameter file does not exist:", 
	print param_file
	exit(0)

command = "%s/mummals %s -ss %d -solv %d -unaligned %d -param %s -outfile %s" %(work_dir, input_file, ss, solv, unaligned, param_file, outfile)

print command
print

os.system(command)

