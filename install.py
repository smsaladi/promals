
from __future__ import print_function

import os
import re
import shutil

# 1. Check for legacy blast
if not shutil.which("blastpgp"):
    raise ValueError("`blastpgp` not found in PATH. Make sure legacy blast is installed")

os.system("cp data/BLOSUM62 bin")
# Skip setting up .ncbirc file

# 2. Check for psipred
if not shutil.which("runpsipred"):
    raise ValueError("`runpsipred` not found in PATH")

# 3. Check if uniref90 is in place
if not os.path.isfile("uniref90/uniref90_filt"):
	raise ValueError("Error: uniref90 not found")

# 4. make the program
os.chdir("src")

print("Compiling the PROMALS program ...")

os.system("make clean")
os.system("make")
os.system("make install")

print("Finished installing PROMALS")

