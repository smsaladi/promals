install: mummals meta_align

CC = g++
OTHERFLAGS = -DNumInsertStates=1 -DVERSION="1.09"
#CXXFLAGS = -O3 -W -Wall -pedantic -DNDEBUG $(OTHERFLAGS) -funroll-loops
#CXXFLAGS = -O3 -W -pedantic -DNDEBUG $(OTHERFLAGS) -funroll-loops
CXXFLAGS = -O3
CFLAGS = -c


mummals: mummals.o kmer_dist.o util.o amino.o mathfunc.o subalign.o sequences.o tnode.o btree_template.h hmm_profpair.o smallnumber.o progressiveAlignHMM.h sparsematrix.o time.o mm.o hmm_profpair1.o param.o hmm_local.o hmm_multim.o refinegap.o consv1.o
	g++ mummals.o sequences.o kmer_dist.o util.o amino.o mathfunc.o subalign.o tnode.o hmm_profpair.o smallnumber.o sparsematrix.o time.o mm.o hmm_profpair1.o param.o hmm_local.o hmm_multim.o refinegap.o consv1.o -o mummals -lm $(CXXFLAGS)

meta_align: meta_align.o kmer_dist.o util.o amino.o mathfunc.o subalign.o sequences.o tnode.o btree_template.h hmm_profpair.o smallnumber.o progressiveAlignHMM.h sparsematrix.o time.o mm.o hmm_profpair1.o param.o hmm_local.o hmm_multim.o refinegap.o consv1.o
	g++ meta_align.o sequences.o kmer_dist.o util.o amino.o mathfunc.o subalign.o tnode.o hmm_profpair.o smallnumber.o sparsematrix.o time.o mm.o hmm_profpair1.o param.o hmm_local.o hmm_multim.o refinegap.o consv1.o -o meta_align -lm $(CXXFLAGS)

mummals.o: mummals.c btree_template.h progressiveAlignHMM.h
	$(CC) $(CFLAGS) mummals.c $(CXXFLAGS)

meta_align.o: meta_align.c
	$(CC) $(CFLAGS) meta_align.c $(CXXFLAGS)

consv1.o: consv1.c consv1.h
	$(CC) $(CFLAGS)  consv1.c $(CXXFLAGS)

smallnumber.o: smallnumber.c smallnumber.h
	$(CC) $(CFLAGS) smallnumber.c $(CXXFLAGS)

mathfunc.o: mathfunc.c mathfunc.h
	$(CC) $(CFLAGS) mathfunc.c $(CXXFLAGS)

util.o: util.h util.c
	$(CC) $(CFLAGS) util.c $(CXXFLAGS)

kmer_dist.o: kmer_dist.c kmer_dist.h
	$(CC) $(CFLAGS) kmer_dist.c $(CXXFLAGS)

sequences.o: sequences.c sequences.h
	$(CC) $(CFLAGS) sequences.c $(CXXFLAGS)

tnode.o: tnode.c tnode.h
	$(CC) $(CFLAGS) tnode.c $(CXXFLAGS)

#hmm_profpair.o: hmm_profpair.h hmm_profpair.c
#	$(CC) $(CFLAGS) hmm_profpair.c $(CXXFLAGS)

hmm_profpair1.o: hmm_profpair1.h hmm_profpair1.c
	$(CC) $(CFLAGS) hmm_profpair1.c $(CXXFLAGS)

hmm_local.o: hmm_local.h hmm_local.c
	$(CC) $(CFLAGS) hmm_local.c $(CXXFLAGS)

hmm_multim.o: hmm_multim.h hmm_multim.c
	$(CC) $(CFLAGS) hmm_multim.c $(CXXFLAGS)

sparsematrix.o: sparsematrix.c sparsematrix.h
	$(CC) $(CFLAGS) sparsematrix.c $(CXXFLAGS)

time.o: time.c time.h
	$(CC) $(CFLAGS) time.c $(CXXFLAGS)

param.o: param.c param.h
	$(CC) $(CFLAGS) param.c  $(CXXFLAGS)

mm.o: mm.h mm.c
	$(CC) $(CFLAGS) mm.c $(CXXFLAGS)

refinegap.o: refinegap.c refinegap.h
	$(CC) $(CFLAGS) refinegap.c $(CXXFLAGS)

.c.o :
	$(CC) $(CFLAGS) $(CXXFLAGS) $?

clean: 
	rm *.o
