#install: testPurge
#install: merge_align
#install: mummals
install: promals_c al2co

CC = g++
OTHERFLAGS = -DNumInsertStates=1 -DVERSION="1.09"
#CXXFLAGS = -O3 -W -Wall -pedantic -DNDEBUG $(OTHERFLAGS) -funroll-loops
#CXXFLAGS = -O3 -W -pedantic -DNDEBUG $(OTHERFLAGS) -funroll-loops
CXXFLAGS = -O3
CFLAGS = -c

promals_c: progress.o kmer_dist.o util.o amino.o mathfunc.o subalign.o sequences.o tnode.o btree_template.h smallnumber.o progressiveAlignHMM.h sparsematrix.o time.o mm.o hmm_profpair1.o param.o hmm_local.o hmm_multim.o refinegap.o consv1.o multiple.o ss_prof.o regularizer.o hmm_psipred.o blastpgp.o gap_refining.o refine.o 
	g++ progress.o  sequences.o kmer_dist.o util.o amino.o mathfunc.o subalign.o tnode.o smallnumber.o sparsematrix.o time.o mm.o hmm_profpair1.o param.o hmm_local.o hmm_multim.o refinegap.o consv1.o multiple.o ss_prof.o regularizer.o hmm_psipred.o blastpgp.o gap_refining.o refine.o -o promals_c -lm $(CXXFLAGS)

refine.o: refine.c refine.h
	$(CC) $(CFLAGS)  refine.c $(CXXFLAGS)
refinegap: refinegap.c consv1.o mathfunc.o subalign.o util.o amino.o
	g++ refinegap.c consv1.o mathfunc.o subalign.o util.o amino.o -o refinegap -lm

consv1.o: consv1.c consv1.h
	$(CC) $(CFLAGS)  consv1.c $(CXXFLAGS)

gap_refining.o: gap_refining.c gap_refining.h
	$(CC) $(CFLAGS) gap_refining.c $(CXXFLAGS)

multiple.o: multiple.c multiple.h progressiveAlignHMM.h btree_template.h
	$(CC) $(CFLAGS)  multiple.c $(CXXFLAGS)

blastpgp.o: blastpgp.c blastpgp.h
	$(CC) $(CFLAGS) blastpgp.c $(CXXFLAGS)

test_pairwise_psipred_aln.o: test_pairwise_psipred_aln.c hmm_psipred.h
	$(CC) $(CFLAGS)  test_pairwise_psipred_aln.c $(CXXFLAGS)

hmm_psipred.o: hmm_psipred.c hmm_psipred.h ScoreType.h
	$(CC) $(CFLAGS)  hmm_psipred.c $(CXXFLAGS)

regularizer.o: regularizer.c 
	$(CC) $(CFLAGS)  regularizer.c $(CXXFLAGS)

ss_prof.o: ss_prof.c  ss_prof.h
	$(CC) $(CFLAGS)  ss_prof.c $(CXXFLAGS)

profilehmm.o: profilehmm.c  profilehmm.h
	$(CC) $(CFLAGS)  profilehmm.c $(CXXFLAGS)

progress.o: progress.c 
	$(CC) $(CFLAGS)  progress.c $(CXXFLAGS)

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

subalign.o: subalign.h subalign.c
	$(CC) $(CFLAGS) subalign.c  $(CXXFLAGS)

al2co: al2co.c
	gcc al2co.c -o al2co $(CXXFLAGS) -lm

.c.o :
	$(CC) $(CFLAGS) $(CXXFLAGS) $?

clean: 
	rm *.o al2co promals_c
