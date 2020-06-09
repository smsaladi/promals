
from __future__ import print_function

import os
import sys
import subprocess
import shutil

# Compile PROMALS
os.chdir("src")
subprocess.check_call(["make", "-j2"])
os.system("cp promals ../bin")
os.system("cp progress_for_web.py ../bin")
os.system("cp promals_c ../bin")
os.system("cp al2co_consensus ../bin")
os.chdir("..")

# Check for promals/bin in PATH
deps = ['progress_for_web.py', 'al2co_consensus']
for d in deps:
    if not shutil.which(d):
        print("`{}` not found in PATH. Make sure promals/bin is added to PATH".format(d))


# Check for dependencies in PATH
deps = ['blastpgp', 'makemat', 'runpsipredplus', 'cd-hit', 'mafft', 'TMalign']
for d in deps:
    if not shutil.which(d):
        raise ValueError("`{}` not found in PATH. Is it installed?".format(d))

# Check optional structural alignment programs
# -- default setting of promals uses TMalign and 'fast' to obtain available"
#    3D structural constraints (-fast 1 -tmalign 1 -dali 0)"
deps = ['fast', 'DaliLite']
for d in deps:
    if not shutil.which(d):
        print("Optional dependency: `{}` not found in path. is it installed?".format(d),
              file=sys.stderr)

os.chdir("..")

