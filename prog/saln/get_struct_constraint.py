#!/usr1/local/bin/python

# this script calculate the alignments using several structure
# comparison programs

# and output the alignments into a constraint file
# fasta format pairwise alignments are separated by "@"

# arguments:
# DAT1 DAT2 [options]

import os, sys, glob, re, shutil

#sys.path.append( "/home/jpei/evolvable_database/bin/" )

import database_util

# pdb1: pdb file name; start1: starting number; end1: ending number (sequentially from 1..naa)
def get_pdb_from_range(pdb1, start1, end1):
        fp = open(pdb1)
        ofp = open(pdb1+".tmp", "w")
        current_resnum = ""
        count=0
        for line in fp:
                if line[0:4]!='ATOM': continue
                if current_resnum != line[22:27]:
                        current_resnum = line[22:27]
                        count+=1
                if count<start1: continue
                if count > end1: break
                ofp.write(line)
        ofp.close()
        fp.close()
        os.system("mv %s.tmp %s" %(pdb1, pdb1))

# from two aligned sequences seq1 and seq2, specified ranges [start1, end1] and [start2, end2]
# get the aligned residue pairs in pos1 and pos2
# if checkseqs is nonzero, check if the sequences matches originalseq1 and originalseq2
def get_aligned_pairs(seq1, seq2, originalseq1, originalseq2, checkseqs, start1, end1, start2, end2):

        # do a few checks
        if len(seq1)!=len(seq2):
                print("The two sequences in the alignment do not have the same length")
                sys.exit(0)
        if start1 == -1:  # special case, the whole sequence
                start1 = 1
                end1 = 1000000
        if start2 == -1:  # special case, the whole sequence
                start2 = 1
                end2 = 1000000
        if checkseqs:
                tmpseq = seq1.replace("-", "").upper()
                if tmpseq != originalseq1[start1-1:end1]:
                        print("Error: sequences do not match for seq1")
                        print("    sequence in alignment: ", tmpseq)
                        print("    sequence in original:  ", originalseq1[start1:end1])
                        sys.exit()
                tmpseq = seq2.replace("-", "").upper()
                if tmpseq != originalseq2[start2-1:end2]:
                        print("Error: sequences do not match for seq2")
                        print("    sequence in alignment: ", tmpseq)
                        print("    sequence in original:  ", originalseq1[start2:end2])
                        sys.exit()
        count1 = 0
        count2 = 0
        pos1 = []
        pos2 = []
        for i, j in enumerate(seq1):
                k = seq2[i]
                if j!='-': count1+=1
                if k!='-': count2+=1
                if j.isupper() and k.isupper():
                        pos1.append(count1+start1-1)
                        pos2.append(count2+start2-1)
        return pos1, pos2

# pdb and sequence file directory
pdb_fa_dir = "/home/jpei/promals/src_structure/structure_db/promals_new/blastresults"

if len(sys.argv)<3:
        print("get_struct_constraint.py - get the structural comparison constraint by dali, fast or tmalign")
        print("Usage:")
        print("get_struct_constraint.py id1 id2 -[dft] -dirname dirname -start1 start1 -start2 start2 -end1 end1 -end2 end2 -clear [01] -pdbdir pdb_fa_dir")
        print()
        print("pdb and other records are available at: pdb_fa_dir")
        sys.exit(0)


# four digit id plus chain id
id1 = sys.argv[1]
id2 = sys.argv[2]

# more options
use_dali = 0
use_fast = 0
use_tmalign = 0
filecreate = 0
fileappend = 0
dirname = ""
resultfile = ""
start1 = start2 = -1
end1 = end2 = 1000000
cleardir = 0
program_dir = "/usr1/bin"
for i in range(len(sys.argv)):
        if sys.argv[i] == '-d' or sys.argv[i] == '-dali': use_dali = 1
        if sys.argv[i] == '-f' or sys.argv[i] == '-fast': use_fast = 1
        if sys.argv[i] == '-t' or sys.argv[i] == '-tmalign': use_tmalign = 1
        if sys.argv[i] == '-resultfile': resultfile = sys.argv[i+1]
        if sys.argv[i] == '-dirname': dirname = sys.argv[i+1]
        if sys.argv[i] == '-append': fileappend = 1
        if sys.argv[i] == '-create': filecreate = 1
        if sys.argv[i] == '-start1': start1 = int(sys.argv[i+1])
        if sys.argv[i] == '-start2': start2 = int(sys.argv[i+1])
        if sys.argv[i] == '-end1': end1 = int(sys.argv[i+1])
        if sys.argv[i] == '-end2': end2 = int(sys.argv[i+1])
        if sys.argv[i] == '-clear': cleardir = int(sys.argv[i+1])
        if sys.argv[i] == '-pdbdir': pdb_fa_dir = sys.argv[i+1]
        if sys.argv[i] == '-progdir': program_dir = sys.argv[i+1]


# go to the directory made up of id names
if os.path.isdir(dirname): os.chdir(dirname)
targetdir = "%s_%s_%s" %(id1, id2, os.getpid())
if not os.path.isdir(targetdir): os.mkdir(targetdir)
os.chdir(targetdir)

# open the result file
ofp = ""
if resultfile:
        if fileappend:
                ofp = open(resultfile, "a")
        elif filecreate:
                ofp = open(resultfile, "w")

# pdbids
pdbid1 = id1[:4]
pdbid2 = id2[:4]

# dat files
dat1 = id1 + ".dat"
dat2 = id2 + ".dat"

# sequence files
seqf1 = pdbid1 + ".fa"
seqf2 = pdbid2 + ".fa"

# pdb files
pdb1 = pdbid1 + ".pdb"
pdb2 = pdbid2 + ".pdb"

# get the files
command = "cp %s/%s/%s.fa ." %(pdb_fa_dir, pdbid1[0:3], pdbid1)
os.system(command)
command = "cp %s/%s/%s.fa ." %(pdb_fa_dir, pdbid2[0:3], pdbid2)
os.system(command)
command = "cp %s/%s/%s.pdb ." %(pdb_fa_dir, pdbid1[0:3], pdbid1)
os.system(command)
command = "cp %s/%s/%s.pdb ." %(pdb_fa_dir, pdbid2[0:3], pdbid2)
os.system(command)

# sequences
# first check if the id corresponds to a scop domain or not,
# if is a scop domain, then use the whole pdb, setting start to -1
fp1 =  open(seqf1)
tmpline = fp1.readline().strip()
assert(len(tmpline.split()) == 2)
structid = tmpline.split()[1]
print(structid, end=' ')
if structid[0] == 'D': start1 = -1; # do not consider scop domains separately, otherwise would be 'd'
fp1.close()
fp2 =  open(seqf2)
tmpline = fp2.readline().strip()
assert(len(tmpline.split()) == 2)
structid = tmpline.split()[1]
print(structid)
if structid[0] == 'D': start2 = -1; # do not consider scop domains separately, otherwise would be 'd'
fp2.close()
# now read the sequences
if start1==-1:
        seq1 = open(seqf1).readlines()[1].replace("\n", "")
        originalseq1 = seq1
else:
        originalseq1 = open(seqf1).readlines()[1].replace("\n", "")
        seq1 = originalseq1[start1-1:end1]
if start2==-1:
        seq2 = open(seqf2).readlines()[1].replace("\n", "")
        originalseq2 = seq2
else:
        originalseq2 = open(seqf2).readlines()[1].replace("\n", "")
        seq2 = originalseq2[start2-1:end2]

# structures
if start1 > 0:
        get_pdb_from_range(pdb1, start1, end1)
if start2 > 0:
        get_pdb_from_range(pdb2, start2, end2)

#print start1, start2

# calculate dali alignment
if use_dali:
        # get the .dat file
        if start1 > 0:
                command = "%s/DaliLite -r %s %s 2>/dev/null 1>/dev/null" %(program_dir, pdb1, pdbid1)
                #print command
                os.system(command)
        else:
                if not os.path.isdir("DAT"): os.mkdir("DAT")
                datlist = glob.glob("%s/%s/%s*.dat"  %(pdb_fa_dir, pdbid1[0:3], pdbid1) )
                if len(datlist):
                        command = "cp %s/%s/%s*.dat DAT" %(pdb_fa_dir, pdbid1[0:3], pdbid1)
                        print(command)
                        os.system(command)
                else:
                        command = "%s/DaliLite -r %s %s 2>/dev/null 1>/dev/null" %(program_dir, pdb1, pdbid1)
                        print(command)
                        os.system(command)
        if start2 > 0:
                command = "%s/DaliLite -r %s %s 2>/dev/null 1>/dev/null" %(program_dir, pdb2, pdbid2)
                #print command
                os.system(command)
        else:
                if not os.path.isdir("DAT"): os.mkdir("DAT")
                datlist = glob.glob("%s/%s/%s*.dat"  %(pdb_fa_dir, pdbid2[0:3], pdbid2) )
                if len(datlist):
                        command = "cp %s/%s/%s*.dat DAT" %(pdb_fa_dir, pdbid2[0:3], pdbid2)
                        #print command
                        os.system(command)
                else:
                        command = "%s/DaliLite -r %s %s 2>/dev/null 1>/dev/null" %(program_dir, pdb2, pdbid2)
                        #print command
                        os.system(command)
        datlist = glob.glob("DAT/%s*.dat" %id1)
        if len(datlist)==0:
                print("Error: cannot get the dat file for", id1)
                sys.exit()
        datlist = glob.glob("DAT/%s*.dat" %id2)
        if len(datlist)==0:
                print("Error: cannot get the dat file for", id2)
                sys.exit()
        id1 = glob.glob("DAT/%s*.dat" %id1)[0].replace("DAT/", "").replace(".dat", "")
        id2 = glob.glob("DAT/%s*.dat" %id2)[0].replace("DAT/", "").replace(".dat", "")
        dali_results = database_util.run_dalilite(id1, id2, seq1, seq2)
        #print dali_results
        if len(dali_results)==6:
                if ofp:
                        ofp.write("@\n")
                        ofp.write(">" + pdbid1 + "\n");
                        ofp.write(dali_results[4] + "\n");
                        ofp.write(">" + pdbid2 + "\n");
                        ofp.write(dali_results[5] + "\n");
                else:
                        print(">seq1:")
                        print(dali_results[4])
                        print(">seq2:")
                        print(dali_results[5])
                        pos1, pos2 = get_aligned_pairs(dali_results[4], dali_results[5], originalseq1, originalseq2, 1, start1, end1, start2, end2)
                        print("aligned positions", len(pos1))
                        for i, j in enumerate(pos1):
                                print(j, pos2[i])


if use_fast:
        fast_prog = "%s/fast" %program_dir
        fast_results = database_util.run_fast(pdb1, pdb2, seq1, seq2, fast_prog)
        if len(fast_results)==6:
                if ofp:
                        ofp.write("@\n")
                        ofp.write(">" + pdbid1 + "\n");
                        ofp.write(fast_results[4] + "\n");
                        ofp.write(">" + pdbid2 + "\n");
                        ofp.write(fast_results[5] + "\n");
                else:
                        print(">seq1:")
                        print(fast_results[4])
                        print(">seq2:")
                        print(fast_results[5])
                        pos1, pos2 = get_aligned_pairs(fast_results[4], fast_results[5], originalseq1, originalseq2, 1, start1, end1, start2, end2)
                        print("aligned positions", len(pos1))
                        for i, j in enumerate(pos1):
                                print(j, pos2[i])
if use_tmalign:
        tmalign_prog = "%s/TMalign" %program_dir
        tmalign_results = database_util.run_tmalign(pdb1, pdb2, seq1, seq2, tmalign_prog)
        if len(tmalign_results)==6:
                if ofp:
                        ofp.write("@\n")
                        ofp.write(">" + pdbid1 + "\n");
                        ofp.write(tmalign_results[4] + "\n");
                        ofp.write(">" + pdbid2 + "\n");
                        ofp.write(tmalign_results[5] + "\n");
                else:
                        print(">seq1:")
                        print(tmalign_results[4])
                        print(">seq2:")
                        print(tmalign_results[5])
                        pos1, pos2 = get_aligned_pairs(tmalign_results[4], tmalign_results[5], originalseq1, originalseq2, 1, start1, end1, start2, end2)
                        print("aligned positions", len(pos1))
                        for i, j in enumerate(pos1):
                                print(j, pos2[i])
if ofp:
        ofp.close()

if cleardir:
        os.chdir("..")
        shutil.rmtree(targetdir)
