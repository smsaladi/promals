
import os, sys, re

clean = 1
install = 1
for i in sys.argv:
    if i == '-unclean': clean = 0
    if i == '-uninstall': install = 0

python_prog = sys.executable 

# 1. set up blastpgp 
print "Determine system and processor types ..."

uname_s = os.popen("uname -s").readlines()[0].replace("\n", "").lower()
print "Your operating system is:", uname_s

pgp = ""
makemat = ""
#fast = ""
if uname_s == 'sunos':
    isainfok = os.popen("isainfo -k").readlines()[0].replace("\n", "").lower()
    #if 'sparc' in isainfok:
    if re.search('sparc', isainfok):
        uname_r = os.popen("uname -r").readlines()[0].replace("\n", "").lower()
        if uname_r == '5.8':
            pgp = "prog/blastpgp/blastpgp_sparc64-solaris8"
            makemat = "prog/blastpgp/makemat_sparc64-solaris8"
            #fast = "prog/fast/fast-sun-solaris/fast"
        elif uname_r == '5.10':
            pgp = "prog/blastpgp/blastpgp_sparc64-solaris10"
            makemat = "prog/blastpgp/makemat_sparc64-solaris10"
            #fast = "prog/fast/fast-sun-solaris/fast"
    else:
        isainfob = os.popen("isainfo -b").readlines()[0].replace("\n", "").lower()
        if isainfob == '32':
            pgp = "prog/blastpgp/blastpgp_ia32-solaris9"
            makemat = "prog/blastpgp/makemat_ia32-solaris9"
            #fast = "prog/fast/fast-sun-solaris/fast"
        if isainfob == '64':
            pgp = "prog/blastpgp/blastpgp_x64-solaris10"    
            makemat = "prog/blastpgp/makemat_x64-solaris10"
            #fast = "prog/fast/fast-sun-solaris/fast" # not available
#elif 'irix' in uname_s:
elif re.search('irix', uname_s):
    pgp = "prog/blastpgp/blastpgp_mips64-irix"
    makemat = "prog/blastpgp/makemat_mips64-irix"
    #fast = "prog/fast/fast-sun-solaris/fast" # not available
#elif 'linux' in uname_s:
elif re.search('linux', uname_s):
    uname_i = os.popen("uname -i").readlines()[0].replace("\n", "").lower()
    if uname_i == 'i386':
        pgp = "prog/blastpgp/blastpgp_ia32-linux"
        makemat = "prog/blastpgp/makemat_ia32-linux"
        #fast = "prog/fast/fast-linux-32/fast"
    elif uname_i == "ia64":
        pgp = "prog/blastpgp/blastpgp_ia64-linux"
        makemat = "prog/blastpgp/makemat_ia64-linux"
        #fast = "prog/fast/fast-linux-64/fast"
    elif uname_i == "x86_64":
        pgp = "prog/blastpgp/blastpgp_x64-linux"
        makemat = "prog/blastpgp/makemat_x64-linux"
        #fast = "prog/fast/fast-linux-64/fast"
    else:
        pgp = "prog/blastpgp/blastpgp_ia32-linux"
        makemat = "prog/blastpgp/makemat_ia32-linux"
        #fast = "prog/fast/fast-linux-32/fast"
#elif 'freebsd' in uname_s:
elif re.search('freebsd', uname_s):
    pgp = "prog/blastpgp/blastpgp_ia32-freebsd"
    makemat = "prog/blastpgp/makemat_ia32-freebsd"
    #fast = "prog/fast/fast-sun-solaris/fast" # not available
#elif 'darwin' in uname_s:
elif re.search('darwin', uname_s):
    pgp = "prog/blastpgp/blastpgp_universal-macosx"
    makemat = "prog/blastpgp/makemat_universal-macosx"
    #fast = "prog/fast/fast-mac-osx/fast"
#elif 'osf' in uname_s:
elif re.search('osf', uname_s):
    pgp = "prog/blastpgp/blastpgp_axp64-tru64"
    makemat = "prog/blastpgp/makemat_axp64-tru64"
    #fast = "prog/fast/fast-alpha64/fast"

if not pgp:
    print "Cannot determine system/processor types"
    print "please install a blastpgp executable file in bin/ directory"
else: 
    print "The blastpgp version is: " + pgp
    print "The makemat version is: " + makemat

# clean up the bin directory
os.chdir("bin/")
os.system("rm -f promals_c promals TMalign progress_for_web.py cd-hit blastpgp makemat runpsipred1 al2co_consensus mafft 2>/dev/null")
os.chdir("../")

# clean up example directory
os.system("rm -r -f example/*")
os.system("cp src/Yfp.fa example/yfp.fa")

# copy the proper blastpgp to bin/ directory; and BLOSUM62
os.system("cp " + pgp + " bin/blastpgp") 
os.chmod("bin/blastpgp", 0755)
os.system("cp " + makemat + " bin/makemat") 
os.chmod("bin/makemat", 0755)
#os.system("cp prog/blastpgp/BLOSUM62 bin")
#os.system("cp " + fast + " bin/fast")  # for program fast too
#os.chmod("bin/fast", 0755)

# set up the .ncbirc file
'''
homedir = os.environ["HOME"]
ncbircfile = homedir + "/.ncbirc"
if not os.path.isfile(ncbircfile):
    tmpfp = open(ncbircfile, "w")
    tmpfp.write("[NCBI]\n\tDATA=")
    data_dir = os.getcwd() + "/prog/blastpgp"
    tmpfp.write(data_dir + "\n")
    tmpfp.close()
'''
ncbircfile = "bin/.ncbirc"
tmpfp = open(ncbircfile, "w")
tmpfp.write("[NCBI]\n\tDATA=")
data_dir = os.getcwd() + "/prog/blastpgp"
tmpfp.write(data_dir + "\n")
tmpfp.close()

# 2. set up psipred
print "\nInstalling psipred ..."
os.chdir("prog/psipred/bin")
os.system("rm -f * 2>/dev/null")
os.chdir("../src")
if clean: os.system("make clean")
if install: os.system("make")
if install: os.system("make install")
os.chdir("../")

# check if psipred is in place
if install:
    if not os.path.isfile("bin/psipred"):
        print "Error: psipred package failed to install correctly"
        sys.exit(0)

os.chdir("../../")
print

# 3. set up uniref90

# check if uniref90 is in place
#if not os.path.isfile("db/uniref90/uniref90_filt"):
    #print "Error: uniref90 not found"
    #sys.exit(0)

# 4. change program directory name in param.c
current_dir = os.path.abspath(".")
fp = open("src/param.c")
if not fp: 
    print "Error: cannot find param.c in current directory"
    sys.exit(0)
fp1 = open("src/tmp_param.c", "w")
if not fp1:
    print "Error: cannot create tmp_param.c in current directory"
    sys.exit(0)
myline = ""
for line in fp.readlines():
    #if "strcpy(program_dir" in line:
    if re.search("strcpy\(program_dir", line):
        myline = "\tstrcpy(program_dir, \"" + current_dir + '");\n'
        fp1.write(myline)
    else:
        fp1.write(line)
fp1.close()
fp.close()
os.system("mv src/tmp_param.c src/param.c")

#sys.exit()

# 5. make the program
os.chdir("src")
print "Installing promals program ..."
if clean: os.system("make clean")
makestatus = 0
if install: 
    makestatus = os.system("make")
    if makestatus:
        print "-Error with compiling PROMALS."
        sys.exit(1)
    os.system("cp promals ../bin")
    os.system("cp progress_for_web.py ../bin")
    os.system("cp promals_c ../bin")
    os.system("cp al2co_consensus ../bin")
os.chdir("..")

# 6. set up the promals script
if install:
    pfp = os.popen("which python")
    python_des = pfp.readline().strip()
    python_des = os.path.abspath(python_des)
    os.chdir("bin")
    fp = open("promals")
    fp1 = open("tmp_promals", "w")
    if not fp:
        print "Error: promals not find in current directory"
        sys.exit(0)
    if not fp1:
        print "Error: cannot create tmp_promals"
        sys.exit(0)
    
    for line in fp.readlines():
        if re.match("#!.*python", line):
            line = "#!%s\n" %python_des
        if re.match("prog_dir", line):
            line = "prog_dir = '%s/bin/'\n" %(current_dir)
        if re.match("#!", line):
            line = "#!%s\n" %(python_prog)
        fp1.write(line)
    fp.close()
    fp1.close()
    os.system("mv tmp_promals promals")
    os.chmod("promals", 0755);
    os.chdir("..")

# 7. set up progress_for_web.py
if install:
    os.chdir("bin")
    fp = open("progress_for_web.py")
    fp1 = open("tmp_progress_for_web.py", "w")
    if not fp:
        print "Error: promals not find in current directory"
        sys.exit(0)
    if not fp1:
        print "Error: cannot create tmp_promals"
        sys.exit(0)
    for line in fp.readlines():
        #if "al2co" in line:
        if re.match("#!.*python", line):
            line = "#!%s\n" %python_des
        if re.search("al2co", line):
            line = re.sub("\/.*al2co", current_dir+"/bin/al2co", line)
            #line = line.replace("/home/jpei/bin", current_dir+"/bin")
        fp1.write(line)
    fp.close()
    fp1.close()
    os.system("mv tmp_progress_for_web.py progress_for_web.py") 
    os.chdir("..")

# 8. set up runpsipred1
# determine tcsh
command = "which tcsh"
pfp = os.popen(command)
tcsh = pfp.readline().strip()
pfp.close()
print "tcsh is available (or not) at: %s" %tcsh

makemat_dir = current_dir + "/bin"
psipred_bin_dir = current_dir + "/prog/psipred/bin"
psipred_data_dir = current_dir + "/prog/psipred/data"
os.chdir("prog/psipred")
fp = open("runpsipred1")
fp1 = open("tmp_runpsipred1", "w")
if not fp:
    print "Error: runpsipred1 not find in current directory"
    sys.exit(0)
if not fp1:
    print "Error: tmp_runpsipred1 cannot be opened"
    sys.exit(0)

mylist = []
for line in fp.readlines():
    # the first line contains tcsh location
    if '#!/bin/tcsh' in line:
        line = '#!%s' %tcsh
        line +='\n'
    #if "set ncbidir = " in line:
    if re.search("set ncbidir = ", line):
        mylist = line.split()
        mylist[3] = makemat_dir
        line = ' '.join(mylist)
        line += '\n'
    if re.search("set execdir = ", line):
        mylist = line.split()
        mylist[3] = psipred_bin_dir
        line = ' '.join(mylist)
        line += '\n'
    if re.search("set datadir = ", line):
        mylist = line.split()
        mylist[3] = psipred_data_dir
        line = ' '.join(mylist)
        line += '\n'
    fp1.write(line)
os.system("mv tmp_runpsipred1 runpsipred1")
os.chmod("runpsipred1", 0755)
fp.close()
fp1.close()
os.system("cp runpsipred1 ../../bin")
os.chdir("../..")

# now test makemat
print "Testing makemat ..."
os.chdir("test")
os.system("rm -f psitmp*")
if not os.path.isfile("seq7.fa"):
    print "Warning: seq7.fa is not present in test/"
    print "makemat is not tested"
elif not os.path.isfile("seq7.chk"):
    print "Warning: seq7.chk is not present in test/"
    print "makemat is not tested"
else:
    print "testing makemat on files seq7.fa and seq7.chk"
    os.system("cp seq7.chk psitmp.chk")
    os.system("echo psitmp.chk > psitmp.pn")
    os.system("echo seq7.fa > psitmp.sn")
    if not os.path.isfile("../bin/makemat"):
        print "Error - ../bin/makemat does not exist"
    print "Execute makemat..."
    print "../bin/makemat -P psitmp >& makemat.log"
    os.system("../bin/makemat -P psitmp >& makemat.log")
    if not os.path.isfile("psitmp.mtx"):
        print "Error - psitmp.mtx not generated by makemat"
    else:
        print "psitmp.mtx generated"
        print "makemat seems OK"
    print "Below is the makemat.log file:"
    if os.path.isfile("makemat.log"):
        fp = open("makemat.log")
        for line in fp: print line.strip('\n')
        os.system("rm -f makemat.log")
        fp.close()
    else:
        print "Error - makemat.log not generated"
    # clean up
    os.system("rm -f psitmp*")
print
os.chdir("../")

# now test runpsipred1
print "Testing runpsipred1 ..."
os.chdir("test")
if os.path.isfile("seq7.ss2"): os.system("rm -f seq7.ss2")
if os.path.isfile("seq7.horiz"): os.system("rm -f seq7.horiz")
if not os.path.isfile("seq7.fa"):
    print "Warning: seq7.fa is not present in test/"
    print "runpsipred1 is not tested"
elif not os.path.isfile("seq7.chk"):
    print "Warning: seq7.chk is not present in test/"
    print "runpsipred1 is not tested"
else:
    print "testing on files seq7.fa and seq7.chk"
    os.system("cp seq7.chk psitmp.chk")
    if not os.path.isfile("../bin/runpsipred1"):
        print "Error - ../bin/runpsipred1 does not exist"
    print "run psipred..."
    print "../bin/runpsipred1 seq7.fa >& seq7.log"
    os.system("../bin/runpsipred1 seq7.fa >& seq7.log")
    if not os.path.isfile("seq7.ss2"):
        print "Error - seq7.ss2 not generated by runpsipred1"
    elif not os.path.isfile("seq7.horiz"):
        print "Error - seq7.horiz not generated by runpsipred1"
    else:
        print "seq7.ss2 and seq7.horiz generated"
        print "runpsipred1 seems OK"
    print "Below is the seq7.log file:"
    if os.path.isfile("seq7.log"):
        fp = open("seq7.log")
        for line in fp: print line.strip('\n')
        fp.close()
        os.system("rm -f seq7.log")
        #os.system("rm -f seq7.ss2"); os.system("rm -f seq7.horiz")
    else:
        print "Error - seq7.log not generated"
print
os.chdir("../")

# 9. set up cd-hit
print "\nInstalling cd-hit ..."
os.chdir("prog/cd-hit")
if clean: os.system("make clean")
if install: 
    os.system("make")
    os.system("cp cd-hit ../../bin")
os.chdir("../../")
print

#10. set up mafft
print "\nInstalling mafft ..."
os.chdir("prog/mafft/src")
if clean: os.system("make clean")
if install: os.system("make")
os.chdir("../")
fp = open("mafft_template")
fp1 = open("tmp_mafft", "w")
if not fp:
    print "Error: mafft not find in current directory"
    sys.exit(0)
if not fp1:
    print "Error: tmp_mafft cannot be opened"
    sys.exit(0)

mylist = []
for line in fp.readlines():
    #if "set ncbidir = " in line:
    if re.search("prefix=", line):
        if 'prefix="$MAFFT_BINARIES"' in line:
            pass
        else:
            line = "\tprefix=%s/prog/mafft/src\n" %current_dir
    fp1.write(line)
fp.close()
fp1.close()
os.system("mv tmp_mafft mafft_template")
os.chmod("mafft_template", 0755)
os.system("cp mafft_template ../../bin/mafft")
os.chdir("../..")

# 11. set up tmalign
print "\nInstalling TMalign ..."
os.chdir("prog/TMalign")
command = "rm -f TMalign 2>/dev/null"
os.system(command)
if install: 
    print "install by g77 ..."
    command = "which g77"
    os.system(command)
    #command = "g77 TMalign.f -o TMalign"
    command = "g77 -static -O3 -lm -o TMalign TMalign.f"
    os.system(command)
    if not os.path.isfile("TMalign"):
        print "install by g77 failed"
        print "install by f77 ..."
        command = "which f77"
        os.system(command)
        #command = "f77 TMalign.f -o TMalign"
        command = "f77 -static -O3 -lm -o TMalign TMalign.f" 
        os.system(command)
        if not os.path.isfile("TMalign"):
            print "install by f77 failed"
            print "install by gfortran ..."
            command = "which gfortran"
            os.system(command)
            command = "gfortran -static -O3 -ffast-math -lm -o TMalign TMalign.f"
            os.system(command)
            if not os.path.isfile("TMalign"):
                print "install by gfortran failed"
    command = "cp TMalign ../../bin"
    os.system(command)
os.chdir("../..")
print


# 12. check if everything is installed or uninstalled
if not install:
    os.chdir("src/example")
    os.system("rm -f yfp.fa.* 2>/dev/null")
    os.system("rm -f yfp.promals.* 2>/dev/null")
    os.chdir("../../")

os.chdir("bin")
if not install:
    os.system("rm -f promals_c promals TMalign progress_for_web.py cd-hit blastpgp makemat runpsipred1 al2co_consensus mafft 2>/dev/null")
    sys.exit()
print "Checking executables in bin/ ...\n"
if not os.path.isfile("promals"):
    print "Error: promals is not installed"
    print
else:
    print "promals is present"
if not os.path.isfile("promals_c"):
    print "Error: promals_c is not installed"
    print
else: print "promals_c is present"
if not os.path.isfile("progress_for_web.py"):
    print "Error: progress_for_web.py is not installed"
    print
else:
    print "progress_for_web.py is present"
#if not os.path.isfile("al2co"):
#    print "Error: al2co is not installed"
#    print
if not os.path.isfile("al2co_consensus"):
    print "Error: al2co_consensus is not installed"
    print
else: print "al2co_consensus is present"
if not os.path.isfile("cd-hit"):
    print "Error: cd-hit is not installed"
    print
else: print "cd-hit is present"
if not os.path.isfile("blastpgp"):
    print "Error: blastpgp is not installed"
    print "'blastpgp' executable needs to be in the bin/ directory"
    print
else: print "blastpgp is present"
if not os.path.isfile("makemat"):
    print "Error: makemat is not installed"
    print "'makemat' executable needs to be in the bin/ directory"
    print
else: print "makemat is present"
if not os.path.isfile("runpsipred1"):
    print "Error: runpsipred1 is not installed"
    print "'runpsipred1' executable needs to be in the bin/ directory"
    print
else: print "runpsipred1 is present"
if not os.path.isfile("mafft"):
    print "Error: mafft is not installed"
    print "'mafft' executable needs to be in the bin/ directory"
    print
else: print "mafft is present"
print
print "Below check structural alignment programs: TMalign, 'fast', DaliLite"
print " -- default setting of promals uses TMalign and 'fast' to obtain available"
print "    3D structural constraints (-fast 1 -tmalign 1 -dali 0)"
print 
if not os.path.isfile("TMalign"):
    print "Warning: TMalign is not installed"
    print "'TMalign' executable needs to be in the bin/ directory for using TMalign structural alignment"
    print
else: 
    print "TMalign is present"
    print "Checking TMalign..."
    mark = 0
    command = "./TMalign ../test/test.pdb ../test/test.pdb"
    pfp = os.popen(command)
    for line in pfp:
        if "VATCRPDEFQCSDGNCIHGSRQCDREYDCKDMSDEVGCVN" in line:
            mark = 1
    pfp.close()
    if mark==0:
        print "Error with test of 'TMalign' executable"
    else:
        print "OK"
if not os.path.isfile("fast"):
    print "Warning: fast is not present"
    print "'fast' executable needs to be in the bin/ directory for using 'fast' structural alignment"
    print "'fast' is available at: http://biowulf.bu.edu/FAST/download.htm"
    print
else: 
    print "'fast' is present"
    print "Checking 'fast'..."
    mark = 0
    command = "./fast ../test/test.pdb ../test/test.pdb"
    pfp = os.popen(command)
    for line in pfp:
        if "VATCRPDEFQCSDGNCIHGSRQCDREYDCKDMSDEVGCVN" in line:
            mark = 1
    pfp.close()
    if mark==0:
        print "Error with test of 'fast' executable"
    print "OK"
if not os.path.isfile("DaliLite"):
    print "Warning: DaliLite is not present"
    print "'DaliLite' executable needs to be in the bin/ directory for using 'DaliLite' structural alignment"
    print "'DaliLite' is available at: http://ekhidna.biocenter.helsinki.fi/dali_lite/downloads"
    print
os.chdir("..")


print
print "\nFinished installing PROMALS\n"

