
import os, sys, re, popen2

python_prog = sys.executable 

# 1. set up blastpgp 
print "Determine system and processor types ..."

uname_s = popen2.popen2("uname -s")[0].readlines()[0].replace("\n", "").lower()
print "Your operating system is:", uname_s

pgp = ""
makemat = ""
if uname_s == 'sunos':
	isainfok = popen2.popen2("isainfo -k")[0].readlines()[0].replace("\n", "").lower()
	#if 'sparc' in isainfok:
	if re.search('sparc', isainfok):
		uname_r = popen2.popen2("uname -r")[0].readlines()[0].replace("\n", "").lower()
		if uname_r == '5.8':
			pgp = "blastpgp/blastpgp_sparc64-solaris8"
			makemat = "blastpgp/makemat_sparc64-solaris8"
		elif uname_r == '5.10':
			pgp = "blastpgp/blastpgp_sparc64-solaris10"
			makemat = "blastpgp/makemat_sparc64-solaris10"
	else:
		isainfob = popen2.popen2("isainfo -b")[0].readlines()[0].replace("\n", "").lower()
		if isainfob == '32':
			pgp = "blastpgp/blastpgp_ia32-solaris9"
			makemat = "blastpgp/makemat_ia32-solaris9"
		if isainfob == '64':
			pgp = "blastpgp/blastpgp_x64-solaris10"		
			makemat = "blastpgp/makemat_x64-solaris10"
#elif 'irix' in uname_s:
elif re.search('irix', uname_s):
	pgp = "blastpgp/blastpgp_mips64-irix"
	makemat = "blastpgp/makemat_mips64-irix"
#elif 'linux' in uname_s:
elif re.search('linux', uname_s):
	uname_i = popen2.popen2("uname -i")[0].readlines()[0].replace("\n", "").lower()
	if uname_i == 'i386':
		pgp = "blastpgp/blastpgp_ia32-linux"
		makemat = "blastpgp/makemat_ia32-linux"
	elif uname_i == "ia64":
		pgp = "blastpgp/blastpgp_ia64-linux"
		makemat = "blastpgp/makemat_ia64-linux"
	elif uname_i == "x86_64":
		pgp = "blastpgp/blastpgp_x64-linux"
		makemat = "blastpgp/makemat_x64-linux"
#elif 'freebsd' in uname_s:
elif re.search('freebsd', uname_s):
	pgp = "blastpgp/blastpgp_ia32-freebsd"
	makemat = "blastpgp/makemat_ia32-freebsd"
#elif 'darwin' in uname_s:
elif re.search('darwin', uname_s):
	pgp = "blastpgp/blastpgp_universal-macosx"
	makemat = "blastpgp/makemat_universal-macosx"
#elif 'osf' in uname_s:
elif re.search('osf', uname_s):
	pgp = "blastpgp/blastpgp_axp64-tru64"
	makemat = "blastpgp/makemat_axp64-tru64"

if not pgp:
	print "Cannot determine system/processor types"
	print "please install a blastpgp executable file in bin/ directory"
else: print "The blastpgp version is: " + pgp

# copy the proper blastpgp to bin/ directory; and BLOSUM62
os.system("cp " + pgp + " bin/blastpgp") 
os.chmod("bin/blastpgp", 0755)
os.system("cp " + makemat + " bin/makemat") 
os.chmod("bin/makemat", 0755)
os.system("cp blastpgp/BLOSUM62 bin")

# set up the .ncbirc file if it is not available in the home directory
homedir = os.environ["HOME"]
ncbircfile = homedir + "/.ncbirc"
if not os.path.isfile(ncbircfile):
	tmpfp = open(ncbircfile, "w")
	tmpfp.write("[NCBI]\n\tDATA=")
	data_dir = os.getcwd() + "/data"
	tmpfp.write(data_dir + "\n")
	tmpfp.close()

# 2. set up psipred
print "\nInstalling psipred ..."
os.chdir("psipred/bin")
os.system("rm *")
os.chdir("../src")
os.system("make clean")
os.system("make")
os.system("make install")
os.chdir("../")

# check if psipred is in place
if not os.path.isfile("bin/psipred"):
	print "Error: psipred package failed to install correctly"
	sys.exit(0)

os.chdir("..")
print

# 3. set up uniref90

# check if uniref90 is in place
if not os.path.isfile("uniref90/uniref90_filt"):
	print "Error: uniref90 not found"
	sys.exit(0)

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

# 5. make the program
os.chdir("src")
print "Compiling the PROMALS program ..."
os.system("make clean")
os.system("make")
os.system("cp promals ../bin")
os.system("cp progress_for_web.py ../bin")
os.system("cp promals_c ../bin")
os.system("cp al2co ../bin")
os.chdir("..")

# 6. set up the promals script
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
	if re.search("al2co", line):
		line = re.sub("\/.*al2co", current_dir+"/bin/al2co", line)
		#line = line.replace("/home/jpei/bin", current_dir+"/bin")
	fp1.write(line)
fp.close()
fp1.close()
os.system("mv tmp_progress_for_web.py progress_for_web.py") 
os.chdir("..")

# 8. set up runpsipred1
makemat_dir = current_dir + "/bin"
psipred_bin_dir = current_dir + "/psipred/bin"
psipred_data_dir = current_dir + "/psipred/data"
os.chdir("psipred")
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
os.chdir("..")

print "\nFinished installing PROMALS\n"

