#
#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#
# trusty makefile template from:
# https://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html
#

CC = gcc
CXX = g++
CFLAGS = -g -O3
CXXFLAGS = -Wall -g -O3

INCLUDES = 
LFLAGS =
LIBS = -lm

SRCS = progress.c kmer_dist.c util.c amino.c mathfunc.c subalign.c sequences.c tnode.c btree_template.h smallnumber.c progressiveAlignHMM.h sparsematrix.c time.c mm.c hmm_profpair.c param.c hmm_local.c hmm_multim.c refinegap.c consv1.c multiple.c ss_prof.c regularizer.c hmm_psipred.c blastpgp.c gap_refining.c refine.c seq_str_aln.c constraint.c structutil.c filter_similar.c

OBJS = $(SRCS:.c=.o)

MAIN = promals_c

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:    $(MAIN)

refinegap: refinegap.c consv1.o mathfunc.o subalign.o util.o amino.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LFLAGS) $(LIBS) refinegap.c consv1.o mathfunc.o subalign.o util.o amino.o -o refinegap


$(MAIN): $(OBJS) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

.c.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN) al2co_consensus

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
