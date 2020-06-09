#!/usr1/local/bin/python

#This is simple module for alignment database.

import sys, os, time, math, tempfile, shutil, glob, cPickle, re, math
import struct

############################
# Setting
############################
#home dir
user_home = '/home/jpei/'
db_home = user_home + 'projects/evolvable_database/'

#USER id!!
user = os.getenv( "USER" ) #supposed to be kim

#log file
log_file = db_home + 'log'

#update_job_directory
update_job_dir = db_home + 'update_jobs/'

#binary base dir 
binary_base_dir = user_home

#program to Build backbone atoms (maxsprout)

#command to make DaliLite DAT file in the current directory under DAT directory
#dali_dat_dalilite_command = binary_base_dir + 'local/DaliLite_2.4.4/DaliLite'
dali_dat_dalilite_command = "/home/jpei/local/DaliLite_2.4.4/DaliLite"
dali_dat_command = '%s -r %%s %%s' % dali_dat_dalilite_command #DaliLite -r pdb_fn uaid

seq_dir = '/home/jpei/evolvable_database/sequences/'
ca_str_dir = '/home/jpei/evolvable_database/ca_structures/'
dat_dir = '/home/jpei/evolvable_database/DAT/'
hhm_dir = '/home/jpei/evolvable_database/hhm/'
a3m_dir = '/home/jpei/evolvable_database/a3m/'
compass_profile_dir = '/home/jpei/evolvable_database/compass_profile/'
numerical_profile_dir = '/home/jpei/evolvable_database/numerical_profile/'

#command to run DaliLite
#dalilite command suppose to have DAT directory to run
dalilite_command = dali_dat_dalilite_command 
dali_command = '%s -a %%s %%s' % dalilite_command

verbose_flag = 2

#DaliLite 2.4.4 wrapper

#FAST wrapper

#TMalign wrapper

########################

def unid2uaid( unid ) :
	'''convert unique number id into unique alphabet id'''
	n = int(unid)

	if n/(26**4) :
		return '????'

	if n <= 0 :
		return '????'

	n = n-1 #adjusting starting from 1
	# 1 -> aaaa, 456976 (or 26*26*26*26) -> zzzz

	a3 = n/(26**3)
	r = n%(26**3)
	a2 = r/(26**2)
	r = r%(26**2)
	a1 = r/26
	r = r%26
	a0 = r

	acount = 'abcdefghijklmnopqrstuvwxyz'
	return acount[a3]+acount[a2]+acount[a1]+acount[a0]

def uaid2unid( uaid ) :
	if len(uaid) != 4 :
		return 0
	
	acount = 'abcdefghijklmnopqrstuvwxyz'
	n3 = acount.index( uaid[0] )
	n2 = acount.index( uaid[1] )
	n1 = acount.index( uaid[2] )
	n0 = acount.index( uaid[3] )
	return n3*(26**3) + n2*(26**2) + n1*(26) + n0 + 1

###################################################
###################################################
##
## Alignment Preparation Utility Functions
##
###################################################
###################################################

def remove_gapped_columns( content, label_length=32) :
	#get columns to be deleted
	columns = [1]*len(content[0])
	for i, c in enumerate( content[0] ) :
		if i >= label_length and c == '-' :
			columns[ i ] = 0

	new_content = []
	for l in content :
		new_line = ''
		for i, c in enumerate(l) :
			if columns[i] :
				new_line += c
		new_content.append( new_line )

	return new_content


#################################
# simple pdb checking routines
#################################

def check_multiple_models( content ) :
	for l in content :
		if l[:5] == 'MODEL' :
			return 1
	return 0

def extract_first_model( content ) :
	new_content = []
	inmodel = 0
	for l in content :
		if l[:5] == 'MODEL' :
			inmodel = 1
			continue
		if inmodel and l[:6] == 'ENDMDL' :
			break
		else :
			new_content.append(l)
	return new_content

def find_first_chainid( content ) :
	#chains = set()
        first_chain_id = ""
	for l in content :
		if l[12:16] == ' CA ' and ( l[0:6] == 'ATOM  ' or l[0:6] == 'HETATM' ) :
                        if first_chain_id == "": 
                                first_chain_id = l[21]
                                break
			#chains.add( l[21] )
	#if len(chains) > 1 : return first_chain_id
	#else : return 0
        return first_chain_id

def check_multiple_chains( content ) :
	chains = set()
	for l in content :
		if l[12:16] == ' CA ' and ( l[0:6] == 'ATOM  ' or l[0:6] == 'HETATM' ) :
			chains.add( l[21] )
	if len(chains) > 1 :
		return 1
	else :
		return 0

def merge_multiple_chains( content ) :
	atmnum = 0
        resnum = 0
	new_content = []
        prn = '    ' #previous residue number
        for l in content :
                if l[:6] == 'HETATM' or l[:6] == 'ATOM  ' or l[:3] == 'TER' :
                        if prn == l[21:26] :
                                pass
                        else :
                                prn = l[21:26]
                                resnum += 1
                        atmnum += 1
                        if l[21].isalnum() :
                                new_content.append( l[:6]+'%5d'%(atmnum)+l[11:21] + 'A' + '%4d'%(resnum)+l[26:] )
                        else :
                                print >>sys.stderr, 'Error: no chain_id is detected!', l

                else :
			new_content.append( l )

	return new_content

def select_chain( content, chain_id ) :
	atmnum = 0
        resnum = 0
	new_content = []
        prn = '    ' #previous residue number
        for l in content :
                if l[:6] == 'HETATM' or l[:6] == 'ATOM  ' :
                        if l[21]!=chain_id: continue
                        if prn == l[21:26] :
                                pass
                        else :
                                prn = l[21:26]
                                resnum += 1
                        atmnum += 1
                        if l[21].isalnum() or l[21]==' ':
                                #new_content.append( l[:6]+'%5d'%(atmnum)+l[11:21] + 'A' + '%4d'%(resnum)+l[26:] )
                                new_content.append( l[:6]+'%5d'%(atmnum)+l[11:21] + chain_id + '%4d'%(resnum)+l[26:] )
                        else :
                                print >>sys.stderr, 'Error: no chain_id is detected!', l

                else :
			new_content.append( l )

	return new_content

def convert_selenomet_into_met( content ) :
	new_content = []
	for l in content :
		if l[17:20] == 'MSE' and l[:6] == 'ATOM  ':
			new_content.append( l[:17] + 'MET' + l[20:] )
		elif l[17:20] == 'MSE' and l[:6] == 'HETATM' :
			new_content.append( 'ATOM  '+l[6:17] + 'MET' + l[20:] )
		else :
			new_content.append( l )

	return new_content

def extract_ca_only_record( content ) :
	new_content = []
	for l in content :
		if (l[:4] == 'ATOM' or l[:6] == 'HETATM') and l[12:16] == ' CA ' and (l[11] == ' ' or l[11]=='A' or l[11] =='1') :
			new_content.append( l ) 
	return new_content

def save_sequence( seq, uaid ) :
	fp = open( seq_dir + '/' + uaid + '.fa', 'w' )
	print >>fp, '>' + uaid + '\n' + seq 
	fp.close()

# JP modification: add directory name to parameters
def save_sequence_JP( seq, uaid, seq_dir_JP ) :
	fp = open( seq_dir_JP + '/' + uaid + '.fa', 'w' )
	print >>fp, '>' + uaid + '\n' + seq 
	fp.close()

def save_ca_str_pdb_file( ca_only_content, uaid ) :
	pdb_fn = ca_str_dir + '/' + uaid + '.pdb'
	fp = open( pdb_fn, 'w' )
	new_content = []
	for l in ca_only_content :
		if l[:6] == 'ATOM  ' and l[17:20] in three2one_common :
			new_content.append( l )
		elif l[:6] == 'HETATM' and l[17:20] in three2one_common :
			new_content.append( 'ATOM  ' + l[6:] )
		else :
			new_content.append( 'ATOM  ' + l[6:17] + 'ALA' + l[20:] )
	
	fp.writelines( new_content ) #save only ca records only ATOM and not non-common amino acids...
				     #All non-common aa changed into Alanine!!

# JP modification: add directory name in parameters
def save_ca_str_pdb_file_JP( ca_only_content, uaid, ca_dir) :
	pdb_fn = ca_dir + '/' + uaid + '.pdb'
	fp = open( pdb_fn, 'w' )
	new_content = []
	for l in ca_only_content :
		if l[:6] == 'ATOM  ' and l[17:20] in three2one_common :
			new_content.append( l )
		elif l[:6] == 'HETATM' and l[17:20] in three2one_common :
			new_content.append( 'ATOM  ' + l[6:] )
		else :
			new_content.append( 'ATOM  ' + l[6:17] + 'ALA' + l[20:] )
	
	fp.writelines( new_content ) #save only ca records only ATOM and not non-common amino acids...
				     #All non-common aa changed into Alanine!!

def check_ca_only_simple( content ) :
	ca_number = None
	for i, l in enumerate(content) :
		if l[:4] == 'ATOM' and l[12:16] == ' CA ' :
			if ca_number == None :
				ca_number = i
			else :
				if i - ca_number == 1 :
					return 1
				else :
					continue
	return 0

def build_backbone( pdb_content, id='junk' ) :
	temp_dir = tempfile.mkdtemp()
	cur_dir = os.getcwd()

	os.chdir( temp_dir )
        print temp_dir;

	pdb_fn_base = id
	pdb_fn = id + '.pdb'
	open( pdb_fn, 'w' ).writelines( pdb_content ) #save pdb_content as junk.pdb
	
	maxsprout_read_cmd = "/home/jpei/local/maxsprout_new/readbrk -pdb %s -rd ./ -wd ./ >& /dev/null"
	maxsprout_buildbackbone_cmd = '/home/jpei/local/maxsprout_new/buildbackbone -pdb %s -pl /home/jpei/local/maxsprout_new/dglp.list -d Y >& /dev/null'

	os.system( maxsprout_read_cmd % pdb_fn )
	os.system( maxsprout_buildbackbone_cmd % pdb_fn_base )

	content = open( pdb_fn_base + ".brk_mod" ).readlines()
	#sys.stdout.writelines( content )
	new_content = []
	anum = 1
	for l in content :
		if l[:6] == 'ATOM  ' :
			nl = l[:6] + "%5d"%anum + l[11:]
			anum += 1
			new_content.append( nl )

	os.chdir( cur_dir )
	shutil.rmtree( temp_dir,1 )

	return new_content

def change_chainid_as_A( content ) :
	new_content = []
	for l in content :
		if l[:6] == "ATOM  " or l[:6] == 'HETATM' :
			new_content.append( l[:21]+ 'A' + l[22:] )
		else :
			new_content.append( l )
	return new_content

def change_chainid_as_X( content, chainX ) :
	new_content = []
	for l in content :
		if l[:6] == "ATOM  " or l[:6] == 'HETATM' :
			new_content.append( l[:21]+ chainX + l[22:] )
		else :
			new_content.append( l )
	return new_content

def build_dat_file( pdb_content, uaid, cmd_template=dali_dat_command ) :
	cwd = os.getcwd()
	tempdir = tempfile.mkdtemp()
	os.chdir( tempdir )
	pdb_fn = uaid + ".pdb"
	fp = open( pdb_fn, 'w' )
	pdb_content = change_chainid_as_A( pdb_content )
	fp.writelines( pdb_content )
	fp.close()
	os.system( cmd_template% (pdb_fn, uaid) )

	if verbose_flag :
		print "### building dat file for DaliLite", cmd_template%(pdb_fn, uaid)

	if os.path.exists( './DAT/' ) :
		dat_fn = './DAT/' + uaid + 'A.dat'
	else :
		print >>sys.stderr, "Error! No DAT file were generated!", uaid

	shutil.copyfile( dat_fn, dat_dir + '/' + uaid+'A.dat' )
	os.chdir( cwd )
	shutil.rmtree( tempdir, 1 )

# JP modification: add dat_dir in parameters
def build_dat_file_JP( pdb_content, uaid, dat_dir_JP, chain_id, cmd_template=dali_dat_command ) :
	cwd = os.getcwd()

        # judge if the path is absolute or not
        if dat_dir_JP[0] != '/':
                dat_dir_JP = "%s/%s" %(cwd, dat_dir_JP)
        
	tempdir = tempfile.mkdtemp()
        os.system("chmod 755 " + tempdir)
	os.chdir( tempdir )
	pdb_fn = uaid + ".pdb"
	fp = open( pdb_fn, 'w' )
        chain_id = chain_id.upper()
        if chain_id == '-': chain_id = " "
	pdb_content = change_chainid_as_X( pdb_content, chain_id )
	#pdb_content = change_chainid_as_A( pdb_content)
	fp.writelines( pdb_content )
	fp.close()
	os.system( (cmd_template% (pdb_fn, uaid)) + ">& tmp.err" )
        print cmd_template% (pdb_fn, uaid)
        #sys.exit(0)

	if verbose_flag :
		print "### building dat file for DaliLite", cmd_template%(pdb_fn, uaid)

        dat_fn = ""
	if os.path.exists( './DAT/' ) :
		#dat_fn = './DAT/' + uaid + 'A.dat'
                if chain_id == ' ':
        		dat_fn = './DAT/' + uaid + '_.dat'
                else:
		        dat_fn = './DAT/' + uaid + chain_id + '.dat'
	else :
		print >>sys.stderr, "Error! No DAT file were generated!", uaid

        #print dat_fn, dat_dir_JP
        #print dat_dir_JP + '/' + uaid+chain_id+'.dat'
        #os.system("pwd")
        print "|%s|" %dat_dir_JP
        if chain_id != ' ':
        	shutil.copyfile( dat_fn, dat_dir_JP + '/' + uaid+chain_id+'.dat' )
                #print "cp %s %s/%s%s.dat" %(dat_fn, dat_dir_JP, uaid, chain_id)
                #os.system("cp %s %s/%s%s.dat" %(dat_fn, dat_dir_JP, uaid, chain_id) )
        else:
                shutil.copyfile( dat_fn, dat_dir_JP + '/' + uaid+'_.dat' )
	os.chdir( cwd )
	shutil.rmtree( tempdir, 1 )

one2three_map = {'Z':'GLX', 'X':'UNK','G':'GLY', 'P':'PRO', 'A':'ALA', 'V':'VAL', 'L':'LEU', 'I':'ILE', 'M':'MET', 'C':'CYS', 'F':'PHE', 'Y':'TYR', 'W':'TRP', 'H':'HIS', 'K':'LYS', 'R':'ARG', 'Q':'GLN', 'N':'ASN', 'E':'GLU', 'D':'ASP', 'S':'SER', 'T':'THR'}

three2one_map = { 'GLX':'Z', 'UNK':'X', 'GLY':'G', 'PRO':'P', 'ALA':'A', 'VAL':'V', 'LEU':'L', 'ILE':'I', 'MET':'M', 'CYS':'C', 'PHE':'F', 'TYR':'Y', 'TRP':'W', 'HIS':'H', 'LYS':'K', 'ARG':'R', 'GLN':'Q', 'ASN':'N', 'GLU':'E', 'ASP':'D', 'SER':'S', 'THR':'T'}
three2one_common = { 'GLY':'G', 'PRO':'P', 'ALA':'A', 'VAL':'V', 'LEU':'L', 'ILE':'I', 'MET':'M', 'CYS':'C', 'PHE':'F', 'TYR':'Y', 'TRP':'W', 'HIS':'H', 'LYS':'K', 'ARG':'R', 'GLN':'Q', 'ASN':'N', 'GLU':'E', 'ASP':'D', 'SER':'S', 'THR':'T'}

def one2three( one ) :
        if len( one ) != 1 :
                return "UNK"
        else :
                return one2three_map[one]

def extract_sequence( ca_only_content ) :
	if not check_ca_only_simple( ca_only_content) :
		ca_only_content = extract_ca_only_record( ca_only_content )

	seq = []
	for l in ca_only_content :
		if l[17:20] in three2one_map :
			seq.append( three2one_map[ l[17:20] ] )
		else :
			seq.append( 'X' )
	return ''.join( seq )

def extract_sequence_from_dat( uaid ) :
	fn = dat_dir + '/' + uaid + "A.dat"
	content = open( fn ).readlines()
	for l in content :
		if l[:9] == '-sequence' :
			return l.split()[-1][1:]
	return None

def extract_equence_from_dat_file( fn ) :
	content = open( fn ).readlines()
	for l in content :
		if l[:9] == '-sequence' :
			return l.split()[-1][1:]
	return None
	

def extract_sequence_and_coordinates_from_dat( uaid ) :
	fn = dat_dir + '/' + uaid + 'A.dat'
	content = open( fn ).readlines()
	seq = ''
	coord = []
	for l in content :
		if l[:9] == '-sequence' :
			seq = l.split()[-1][1:]
	for l in content :
		if l[:4] == '-acc':
			#read dat coordinates by columns
			coord.append( (float(l[22:34]), float(l[34:46]), float(l[46:58])) ) 
	return seq, coord

def extract_sequence_and_coordinates_from_dat_JP( uaid, dat_dir_JP, chain_id ) :
        if chain_id == ' ':
	        fn = dat_dir_JP + '/' + uaid + '_.dat'
        else: 
                fn = dat_dir_JP + '/' + uaid + chain_id + '.dat'
	content = open( fn ).readlines()
	seq = ''
	coord = []
	for l in content :
		if l[:9] == '-sequence' :
			seq = l.split()[-1][1:]
	for l in content :
		if l[:4] == '-acc':
			#read dat coordinates by columns
			coord.append( (float(l[22:34]), float(l[34:46]), float(l[46:58])) ) 
	return seq, coord

def parse_dat( fn ) :
	content = open( fn ).readlines()
	seq = ''
	coord = []
	for l in content :
		if l[:9] == '-sequence' :
			seq = l.split()[-1][1:]
	for l in content :
		if l[:4] == '-acc':
			#read dat coordinates by columns
			coord.append( (float(l[22:34]), float(l[34:46]), float(l[46:58])) ) 
	return seq, coord

def make_atom_line( anum, aname, rname, rnum, coord ):
        rname = one2three( rname )
        altLoc = ' '
        chain = 'A'
        iCode = ' '
        occu = 1.0
        bfac = 10.0
        segID = '    '
        charge = '  '
        elem = ' C'
        line = "%6s%5d %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s" % ( "ATOM  ", anum, aname, altLoc, rname, chain, rnum, iCode, coord[0], coord[1], coord[2],occu, bfac, segID, elem, charge)
        #line = "%s%5d%s%s%s %s%s%s%8.3f%8.3f%8.3f%6.2f%6.2f%s%s%s" % ( "ATOM  ", anum, aname, altLoc, rname, chain, rnum, iCode, coord[0], coord[1], coord[2],occu, bfac, segID, elem, charge)
        return line

def make_atom_line_JP( anum, aname, rname, rnum, coord , chain_id):
        rname = one2three( rname )
        altLoc = ' '
        chain = chain_id
        iCode = ' '
        occu = 1.0
        bfac = 10.0
        segID = '    '
        charge = '  '
        elem = ' C'
        line = "%6s%5d %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s" % ( "ATOM  ", anum, aname, altLoc, rname, chain, rnum, iCode, coord[0], coord[1], coord[2],occu, bfac, segID, elem, charge)
        #line = "%s%5d%s%s%s %s%s%s%8.3f%8.3f%8.3f%6.2f%6.2f%s%s%s" % ( "ATOM  ", anum, aname, altLoc, rname, chain, rnum, iCode, coord[0], coord[1], coord[2],occu, bfac, segID, elem, charge)
        return line

def make_ca_only_content_from_dat_sequence_and_coordinates( ca_coord, sequence ):

	if len( ca_coord) == len(sequence) :
		pass
	else :
		print "Error! ca_only coordinates and sequence have different lengths!!"
		sys.exit()

	content = []
	for i in range( len(sequence ) ) :
		res = sequence[i]
		coord = ca_coord[i]

		if res == 'X' or res == 'Z' :
			#sys.stderr.write( 'WARNING: contains %s at sequence position %d and it is changed into Ala.\n' %(res,i) )
			res = 'A'
                if res.islower(): res = 'C'

		line = make_atom_line( i+1, ' CA ', res, i+1, coord )
		content.append( line + '\n' )
	return content

def make_ca_only_content_from_dat_sequence_and_coordinates_JP( ca_coord, sequence , chain_id):

	if len( ca_coord) == len(sequence) :
		pass
	else :
		print "Error! ca_only coordinates and sequence have different lengths!!"
		sys.exit()

	content = []
	for i in range( len(sequence ) ) :
		res = sequence[i]
		coord = ca_coord[i]

		if res == 'X' or res == 'Z' :
			#sys.stderr.write( 'WARNING: contains %s at sequence position %d and it is changed into Ala.\n' %(res,i) )
			res = 'A'
                if res.islower(): res = 'C'

		line = make_atom_line_JP( i+1, ' CA ', res, i+1, coord, chain_id )
		content.append( line + '\n' )
	return content

#returns list of a3m files
#also save the a3m file into the a3m directory!!!
def build_a3m( seq, uaid ) :
	#write temp sequence file
	uaid_fn = uaid + '.fa'
	fp = open( uaid_fn, 'w' )
	print >>fp, '>' + uaid + '\n' + seq
	fp.close()

	hhsearch_buildali_cmd = '/home/jpei/local/hhsearch1.5_2/buildali.pl -cn %s'
	os.system( hhsearch_buildali_cmd % (uaid_fn) )
	
	a3m_list = glob.glob( uaid + '*.a3m' )
	a3m_list.sort()
	save_dir = a3m_dir + '/' + a3m_list[0][:-6] + '/' #assuming aaaa.1.a3m
	if os.path.exists( save_dir ) :
		pass
	else :
		os.makedirs( save_dir )
	
	for fn in a3m_list :
		shutil.copy( fn, save_dir )

	return a3m_list

def remove_query_gaps( read_fn, save_fn ) :
	content = open( read_fn ).readlines()
	label_length = 32
	#get columns to be deleted
	columns = [1]*len(content[0])
	for i, c in enumerate( content[0] ) :
		if i >= label_length and c == '-' :
			columns[ i ] = 0

	new_content = []
	for l in content :
		new_line = ''
		for i, c in enumerate(l) :
			if columns[i] :
				new_line += c
		new_content.append( new_line )

	open( save_fn, 'w' ).writelines( new_content )


def convert_a3m_into_psi_and_remove_query_gaps( a3m_list ) :
	temp_fn = 'junk'
	aln_list = []
	for a3m_fn in a3m_list :
		save_fn = a3m_fn[:-4] + '.aln'
		reformat_cmd = '/home/jpei/local/hhsearch1.5_2/reformat.pl a3m psi %s %s'
		os.system( reformat_cmd % (a3m_fn, temp_fn) )
		remove_query_gaps( temp_fn, save_fn )
		aln_list.append( save_fn )

	return aln_list
	
def write_dummy_list_file( fn, dummy='dummy_list' ) :
	open( dummy, 'w' ).write( fn + '\n' )

def convert_aln_into_compass_profile( aln_list ) :
	mk_compass_db_cmd = '/home/jpei/local/compass/mk_compass_db -i dummy_list -o %s '
	for aln in aln_list :
		write_dummy_list_file( aln )
		save_fn = aln[:-4] + '.cnp'
		os.system( mk_compass_db_cmd % save_fn )
	save_dir = compass_profile_dir + '/' + aln_list[0][:4] + '/'
	if os.path.exists( save_dir ) :
		pass
	else :
		os.makedirs( save_dir )
	
	os.system( "cp *.cnp *.cnp.len " + save_dir )

def convert_aln_into_numerical_profile( aln_list ) :
	import profile_score_module
	for aln in aln_list :
		msa_seqs = [ l.split()[-1] for l in open( aln ).readlines() ]

	        neffs, Qs, pssm = profile_score_module.get_neffs_Qs_and_pssm( msa_seqs )
		
		save_fp = open( aln[:-4] + '.pnp', 'w' )
		cPickle.dump( neffs, save_fp, -1 )
		cPickle.dump( Qs, save_fp, -1 )
		cPickle.dump( pssm, save_fp, -1 )
		save_fp.close()

	save_dir = numerical_profile_dir + '/'+ aln_list[0][:4] + '/' 
	if os.path.exists( save_dir  ) :
		pass
	else :
		os.makedirs( save_dir )

	os.system( "cp *.pnp " + save_dir ) 

def read_fas( fas_fn ) :
	content = open( fas_fn ).readlines()
	headers = []
	sequences = []
	seq = ''
	for l in content :
		if l[0] == '>' :
			if seq :
				sequences.append( seq )
			headers.append( l[1:-1] ) 
			seq = ''
		else :
			seq += l[:-1]
	else :
		sequences.append( seq )

	return headers, sequences

def remove_query_gaps_from_simple_msa( content ) :
	#get columns to be deleted
	columns = [1]*len(content[0])
	for i, c in enumerate( content[0] ) :
		if c == '-' :
			columns[ i ] = 0

	new_content = []
	for l in content :
		new_line = ''
		for i, c in enumerate(l) :
			if columns[i] :
				new_line += c
		new_content.append( new_line )

	return new_content

def save_fas( headers, sequences, save_fn ) :
	fp = open( save_fn, 'w' )
	for header, sequence in zip( headers, sequences ) :
		print >> fp, '>' + header + '\n' + sequence


def convert_aln_into_hhm( a3m_list ) :
	reformat_cmd = '/home/jpei/local/hhsearch1.5_2/reformat.pl a3m fas %s %s'
	hhmake_and_calibration_cmd = '/home/jpei/local/hhsearch1.5_2/hhmake -i %s ;/home/jpei/local/hhsearch1.5_2/hhsearch -cal -d /home/jpei/local/hhsearch/cal.hhm -i %s'
	for a3m in a3m_list :
		fas_fn = a3m[:-4] + ".fas"
		os.system( reformat_cmd % ( a3m, 'junk.fas' ) )
		headers, sequences = read_fas( 'junk.fas' )
		sequences = remove_query_gaps_from_simple_msa( sequences )
		save_fas( headers, sequences, fas_fn )
	
		hhm_fn = a3m[:-4] + '.hhm' 
		os.system( hhmake_and_calibration_cmd % (fas_fn, hhm_fn) )
	
	save_dir = hhm_dir + '/' + a3m_list[0][:4] + '/'
	if os.path.exists( save_dir ) :
		pass
	else :
		os.makedirs( save_dir )

	os.system( "cp *.hhm " + save_dir )
	

#######################################################
#######################################################
##
## Program Running Utilities
##
#######################################################
#######################################################

#Robot class
# writes script to run in cluster
# submit script to queue
# check the queue to find errors
# if a single job does not end within 12hrs then report the error to the
# designated email address

class alignment_runner :

	header='''#!/bin/bash
#$ -cwd
#$ -S /bin/tcsh
#$ -j y

mkdir /local_scratch/%s/; 
mkdir /local_scratch/%s/%%s; 
cd /local_scratch/%s/%%s;'''%(user, user, user)

	footer='''cp /local_scratch/%s/%%s/*.result %%s
rm -rf /local_scratch/%s/%%s/''' %(user, user)

	queue_submit_cmd = 'qsub'
	queue_check_cmd = "qstat"
	check_interval = 60 #seconds
	
	def __init__( self, name_base, cmd_template, data_dirs, save_dir, pair_list, nruns_per_job ) :
		self.jobs = []
		self.name_base = name_base
		self.cmd_template = cmd_template

		self.job_dir = update_job_dir + '/%d.%02d.%02d.%02d.%02d'%time.localtime()[:5]
		os.makedirs( self.job_dir ) 
			
		self.data_dirs = data_dirs
		self.save_dir = save_dir
		self.pair_list = pair_list #list of strings "id1 id2"
		self.number_of_runs_per_job = nruns_per_job
		self.order = '_%0'+str( int(math.log( len(pair_list) / nruns_per_job, 10 )) ) + 'd'

	def write_jobs( self ):
		jcount = 0
		current_job = None
		job_fp = None
		job_name = None
		
		for rcount in xrange( len(self.pair_list) ) :
			if not rcount%self.number_of_runs_per_job :

				#write footer 
				if jcount :
					print >>job_fp, self.footer % (job_name, self.save_dir, job_name)
				jcount += 1

				job_name = self.name_base + self.order%jcount + '.job'
				result_fn = self.name_base + self.order%jcount + '.result'

				current_job = [ self.job_dir + "/" + job_name ]
				self.jobs.append( current_job )
				job_fp = open( current_job[0], 'w' )
				
				#write header
				print >>job_fp, self.header % (job_name, job_name)

				#write data directory copying
				for data_dir in self.data_dirs :
					print >>job_fp, "cp -r ", data_dir, "./"
			
			print >>job_fp, self.cmd_template % (self.pair_list[rcount], result_fn)
		else :
			print >>job_fp, self.footer % (job_name, self.save_dir, job_name)

	def submit_job( self, job_fn ) :
		pp = os.popen( self.queue_submit_cmd + ' ' + job_fn )
		job_line = pp.readlines()[0].split()
		return int(job_line[2])
		
	def submit_jobs( self ) :
		pwd = os.getcwd()
		os.chdir( self.job_dir )

		for job in self.jobs :
			job.append( self.submit_job( job[0] ) )
			job.append( 0 ) #submitting time!!

		os.chdir( pwd )

	##resubmitting!!??
			
	def get_current_job_list( self ) :
		pp = os.popen( self.queue_check_cmd + ' -u ' + user )
		ids = [ int(i.split()[0]) for i in pp.readlines()[2:] ]
		return ids

	def wait_for_jobs( self ) :
		remaining_job_ids = [ i[1] for i in self.jobs ]
		while remaining_job_ids :
			time.sleep( self.check_interval )
			remaining_job_ids = self.get_current_job_list()
			
 
#######################################################
#######################################################
##
## Program Result Processing and Update Database
##
#######################################################

'''
def get_x_postion_list_from_sequence( seq ) :
	xpositions = []
	for i, a in enumerate( seq ) :
		if a == 'X' or a == 'x' :
			xpositions.append( i )
	return xpositions
'''

'''
def remove_xmasks_from_aln( seq, aln ) :
	acount = 0
	new_aln = ''
	for a in aln :
		if a.isalpha() :
			if (a == 'A' or a == 'a') and (seq[acount] == 'X' or seq[acount] == 'x')  :
				new_aln += 'X' 
			else :
				new_aln += a
			acount += 1
		else :
			new_aln += a
	return new_aln
'''			

def remove_xmasks_from_aln( seq, aln ) :
        acount = -1
        new_aln = list( aln )
        xpos = [ (a == 'X' or a=='x') for a in seq ]

        '''#debugging
        print len(seq), seq
        seq_from_aln = []
        for a in aln :
                if a.isalpha() :
                        acount += 1
                        seq_from_aln.append(a)
        print acount, aln
        print acout, ''.join( seq_from_aln )
        acount = -1
        '''

        for i, a in enumerate(aln) :
                if a.isalpha() :
                        acount += 1
                        if xpos[acount] :
                                if ( a == 'A' ) :
                                        new_aln[acount] = 'X'
                                elif (a == 'x' ) :
                                        new_aln[acount] = 'x'
                                else :
                                        print >>sys.stderr, "Error! this message should not be seen!", seq, aln

        return ''.join(new_aln)
		

#######################################################
#DaliLite realted processing functions
######################################################

############################################
# parse alignment definition in dccp file
# simple l.split() can have serious bug
# since the alignment definition is written 
# fixed length 
############################################
def parse_alignment_line( line ) :
        lengthrec = 10 # the alignment definition '????  ????' 4+2+4
        if line[-1] == '\n' :
                line = line[:-1]

        nrec = len( line ) / lengthrec #
        if len(line) != lengthrec*nrec :
                print >>sys.stderr, "Error the alignment line has some strange bug!!"
                print >>sys.stderr, line
                sys.exit()

        parsed = []
        for i in xrange( nrec ) :
                temp = line[i*10:(i+1)*10]
                l1, l2 = temp.split()
                parsed.append( l1 )
                parsed.append( l2 )

        return parsed

##############################################################
# parse dccp file
##############################################################


#bug found!!!
def parse_dccp( fp, code2 ) :
        dccp_array = [] #temporary array to save dccp information
        raw_score_array = []
        switch_flag = 0

        dccp_index = 0  #index for max_z ---- 1 based!!! careful!!!
        max_index = 0

        l = fp.readline()
        while l :
                temp = []
                if l[1:5] == 'DCCP':
                        dccp_index += 1
                        raw_score = float( l[9:18] )
                        Z_score = float( l[27:34] )
                        rmsd = float( l[19:23] )
                        raw_score_array.append( raw_score )
                        if not max_index :
                                max_index = 1
                        elif Z_score > dccp_array[max_index-1][2] :
                                max_index = dccp_index

                        #previously
                        #rmsd < rmsd_array[max_index-1] :
                        #was the second condition

                        elif Z_score == dccp_array[max_index-1][2] and raw_score > raw_score_array[max_index-1]  :
                                max_index = dccp_index

                        if  l[69:73] == code2[:4] : #'mol2' :
                                child_id = l[69:74]
                                switch_flag = 1
                        else :
                                child_id = l[75:80]
                                switch_flag = 0

                        temp.append( child_id )
                        temp.append( raw_score )
                        temp.append( Z_score )

                #read alignment
                l = fp.readline()
                l = fp.readline()
                parent = []
                child = []
                alignments = []
                while not l[1:5] == 'DCCP' :
                        #Bug!!! fixed length should be used
                        #l = l.split()
                        l = parse_alignment_line( l ) #fixed!!

                        for i in range( len(l) ) :
                                l[i] = int( l[i] )
                        alignments.append( l )
                        l = fp.readline()
                        #end of loop when the file is read
                        if not l :
                                break

                if not switch_flag :
                        for i in range( len(alignments) / 2 ) :
                                for j in range( len( alignments[i] ) / 2 ) :
                                        parent.append( [alignments[i][j*2],alignments[i][j*2+1] ] )
                        for i in range( len(alignments) / 2, len(alignments) ) :
                                for j in range( len( alignments[i] ) / 2 ) :
                                        child.append( [alignments[i][j*2], alignments[i][j*2+1]] )
                else :
                        for i in range( len(alignments) / 2, len( alignments ) ) :
                                for j in range( len( alignments[i] ) / 2 ) :
                                        parent.append( [alignments[i][j*2],alignments[i][j*2+1] ] )
                        for i in range( len(alignments) / 2 ) :
                                for j in range( len( alignments[i] ) / 2 ) :
                                        child.append( [alignments[i][j*2], alignments[i][j*2+1]] )

                if len( parent ) == len( child ) :
                        pass
                else :
                        print >>sys.stderr, "Error, while parsing dccp file!"

                temp.append( parent )
                temp.append( child )
                dccp_array.append( temp )

        if not max_index :
                return None

        #debug
        #print dccp_array[max_index-1]
        return dccp_array[max_index-1]

        #returns only maximum Z score alignment information
###########################################################
# Alignment building function
###########################################################
def build_alignment( seq1, seq2, def1, def2 ) :
        alignment = ['','']

        prep = 0
        prec = 0
        for (ps,pe),(cs,ce) in zip( def1, def2 ):
                #adding gaps
                if (ps-cs-prep+prec) > 0 :
                        alignment[1] += '-'*(ps-cs-prep+prec)
                elif (cs-ps+prep-prec) > 0 :
                        alignment[0] += '-'*(cs-ps+prep-prec)
                #adding non-equivalent parts
                alignment[0] += seq1[prep:ps-1].lower()
                alignment[1] += seq2[prec:cs-1].lower()
                #equivalent parts
                alignment[0] += seq1[ps-1:pe].upper()
                alignment[1] += seq2[cs-1:ce].upper()
                #saving position
                prep = pe
                prec = ce

        ##
        ##Adding last part beyond definition
        ##

        ps = len(seq1)+1
        pe = len(seq1)+1
        cs = len(seq2)+1
        ce = len(seq2)+1
        #adding gaps
        if (ps-cs-prep+prec) > 0 :
                alignment[1] += '-'*(ps-cs-prep+prec)
        elif (cs-ps+prep-prec) > 0 :
                alignment[0] += '-'*(cs-ps+prep-prec)
        #adding non-equivalent parts
        alignment[0] += seq1[prep:ps-1].lower()
        alignment[1] += seq2[prec:cs-1].lower()
        #equivalent parts
        alignment[0] += seq1[ps-1:pe].upper()
        alignment[1] += seq2[cs-1:ce].upper()
        #saving position
        prep = pe
        prec = ce

        if len( alignment[0] ) != len( alignment[1] ) :
                print >>sys.stderr, "Error! alignment lengths are not equal.", len( alignment[0]), len( alignment[1])
                print >>sys.stderr, alignment[0]
                print >>sys.stderr, alignment[1]

        return alignment

def check_log_file( fn ) :
	pp = os.popen( "egrep -i 'fatal|error' " + fn )
	err_msg = pp.readlines()

	if err_msg :
		print >>sys.stderr, 'WARNING:!', fn, 'has error message in it!'
		return 1 #ERROR!!
	else :
		return 0 #NO ERROR!!
		

def run_dalilite( code1, code2, seq1, seq2 ) :
	id1 = code1[:4]
	id2 = code2[:4]
	#make and move into temporary directory
	cwd = os.getcwd()
	tempdir = tempfile.mkdtemp()
	os.symlink( cwd+ '/DAT', tempdir+'/DAT' )
	os.chdir(tempdir)

	###############
	#RUN DALILITE COMMAND
	os.system( dali_command %(code1, code2) + " >> %s 2>&1" % (code1+code2+'.log') )
	###############
	
	#seq1, coord1 = parse_dat( './DAT/'+code1+'.dat' )
	#seq2, coord2 = parse_dat( './DAT/'+code2+'.dat' )
	
	dccp = parse_dccp( open( code1+'.dccp' ), code2 )
	if not dccp :
		return None
		
	if dccp[0] == code2 :
		pass
        elif code2[4] == '_' and dccp[0][0:4] == code2[0:4]:
                pass
	else :
		print >>sys.stderr, "Error in running DaliLite", code1, code2, tempdir
		print >>sys.stderr, "code2 supposed to be, ", code2," and from file,", dccp[0], "is different!"
		sys.exit()

	'''
	#building alignment and calculate the dali score
	#get equivalence mapping between parent and child
	eq_map = [[],[]]
	for s, e in dccp[3] :
		#s means start offset 1 based
		#e means end offset 1 based
		s = s-1
		e = e

		for i in range( s, e ):
			eq_map[0].append( i )

	for s, e in dccp[4] :
		#s means start offset 1 based
		#e means end offset 1 based
		s = s-1
		e = e

		for i in range( s, e ) :
			eq_map[1].append( i )

	#verifying the equivalences
	if len( eq_map[0] ) != len( eq_map[1] ) :
		print >>sys.stderr, "Error! Equivalence mapping is wrong!"
		#print >>sys.stderr, eq_map[0]
		#print >>sys.stderr, eq_map[1]


	if len(seq[pid]) < eq_map[0][-1] :
		print pid, 'mapping is wrong'
		for i in range(len(eq_map[0]) ) :
			print i, eq_map[0][i]
		for i in range(len(coord[pid]) ) :
			print i, coord[pid][i]
		print dccp[pid]
	if len( seq[cid]) < eq_map[1][-1] :
		print cid, 'mapping is wroing'
		for i in range(len(eq_map[1]) ) :
			print i, eq_map[1][i]
		for i in range(len(coord[cid]) ) :
			print i, coord[cid][i]
		print dccp[cid]
	'''

	aln1, aln2 = build_alignment( seq1, seq2 , dccp[3],dccp[4] )

	check_log_file( code1+code2+'.log' )

	#check log file!!
	#move back into the current directory and cleanup!!
	os.chdir( cwd )
	shutil.rmtree( tempdir )

	#id1, id2, raw_score, Zscore, start position1, start pos2, aln1, aln2
	return (dccp[1], dccp[2], 0, 0, aln1, aln2) 

def get_position( id2, index_fp ) :
	long_byte = 4
	nrecords = 2
	index_fp.seek( (id2-1)*long_byte*nrecords )
	return struct.unpack( 'l', index_fp.read( long_byte ) )
	

def get_old_dali( uaid1, uaid2 ) :
	id1 = uaid2unid( uaid1 )
	id2 = uaid2unid( uaid2 )

	dali_old_dir = '/usr2/db/dalilite_aln_1.69/'
	pid = '%04d'%id1
	cid = '%04d'%id2
	
	#prepare file
	if os.path.exists( pid+ '.result' ) :
		pass
	else :
		shutil.copy(  dali_old_dir + pid + '.result', './' )
	
	if os.path.exists( pid+ '.result.idx' ) :
		pass
	else :
		shutil.copy( dali_old_dir + pid + '.result.idx', './' )
	
	#read index!!
	pos = get_position( id2, open( pid + '.result.idx' ) ) 
	fp = open( pid + '.result' )
	fp.seek( pos[0] ) #from tuple to integer

	header = fp.readline()
	dccp_line = fp.readline()
	if dccp_line == '\n' :
		return None

        raw_score = float( dccp_line[9:18] )
        Z_score = float( dccp_line[27:34] )
        rmsd = float( dccp_line[19:23] )
	
	aln1 = fp.readline()[:-1]
	aln2 = fp.readline()[:-1]

	fp.close()

	return ( raw_score, Z_score, 0, 0, aln1, aln2 )
	

######################################################
# END DaliLite result parsing code
######################################################
######################################################


######################################################
######################################################
# FAST result parsing code
######################################################
######################################################
head_pattern = re.compile( 'FAST ALIGNMENT: (.{4}).\S+\s(.{4}).\S+' )
score_pattern = re.compile( 'L=(\S+)\sSX=(\S+)\sSN=(\S+)\sL1=\S+\sL2=\S+\sRMSD=(\S+)' )
alignment_pattern = re.compile( ' [12]:\t(\S*)\**' )
def parse_fast_record( fp ) :
        #first line
        l = fp.readline()
        #case of no more record
        if not l:
                return None
        #print l
        match = head_pattern.match( l )
        #checking for starting point
        if not match :
                print >>sys.stderr, 'WARNING: wrong starting point!!' ,l
                while( not match ) :
                        l = fp.readline()
                        if not l :
                                print >>sys.stderr, 'End of file reached!'
                                return None
                        print >>sys.stderr, 'following lines:', l
                        match = head_pattern.match( l )

        parent = match.group(1)
        child = match.group(2)
        #second line
        l = fp.readline()
        match = score_pattern.match( l )
        aligned_length = int(match.group(1))
        fast_score = match.group(2)
        norm_score = match.group(3)
        rmsd = match.group(4)

        if not aligned_length :
                l1 = fp.readline()
                l2 = fp.readline()
                if l1 == l2 == '\n' :
                        return parent, child, fast_score, rmsd, aligned_length, norm_score, None, None
                else:
                        print >> sys.stderr, "Error! while parsing fast record. Aborted!", parent, child

        l = fp.readline()
        if l != '\n':
                print >>sys.stderr, "Error, unexpected pattern at line3:", parent, child, l

        #alignment process
        paseq = ''
        caseq = ''
        while ( True ) :
                l = fp.readline()
                if l == '\n' :
                        break
                match = alignment_pattern.match( l )
                if not match :
                        break
                paseq += match.group(1)

                l = fp.readline()
                match = alignment_pattern.match( l )
                caseq += match.group(1)

                l = fp.readline()
                if l == '\n' :
                        pass
                else :
                        print >>sys.stderr, "Error, unexpected pattern in alignment region!", parent, child, paseq, caseq


        #need to consume the trailing equvalence parts
        while ( True ) :
                l = fp.readline()
                if l == '\n' :
                        l = fp.readline()
                        if l == '\n' :
                                break
                        else :
                                print >>sys.stderr, 'Error, unexpected pattern at the end!', parent, child, l

	'''
        eq = ''
        for p,c in zip( paseq, caseq ) :
                if p.isalpha() and c.isalpha() :
                        eq += ':'
                else :
                        eq += ' '
	'''

        return (parent, child, fast_score, rmsd, aligned_length, norm_score, paseq[:-1], caseq[:-1] )

def run_fast( id1, id2, seq1, seq2, fast_program ) :
	pp = os.popen( '%s %s %s'%(fast_program, id1, id2) )
	result = parse_fast_record( pp )
	
	if not result :
		return None
	if 'X' in seq1  or 'x' in seq1 :
		aln1 = remove_xmasks_from_aln( seq1, result[6] )
	else :
		aln1 = result[6]

	if 'X' in seq2  or 'x' in seq2:
		aln2 = remove_xmasks_from_aln( seq2, result[7] )
	else :
		aln2 = result[7]

	return ( float(result[2]), float(result[5]), 0, 0, aln1, aln2 )

#END of FAST result parsing code
######################################################


######################################################
######################################################
#TMalign result parsing code
######################################################

def parse_tm_record( fp ) :
        for i in range( 18 ) :
                l = fp.readline()
                #case of no more record
                if not l and i == 0:
                        return None
                if not l :
                        return None, None
                if i == 8 :
                        child =  l[13:17]
                elif i == 9 :
                        parent = l[13:17]
                elif i == 11 :
                        tm_score =  float( l[43:50] )  #tm_align score
                        rmsd = float( l[26:32] )  #RMSD
                        aligned_length = int( l[15:19] )  # aligned elgnth
                        identity = float( l[55:60] )  #identity
                elif i == 14 :
                        length = len( l )
                        c_alignment = l[:-1]
                elif i == 15 :
                        equivalence = l[:-1]
                        if len(l) == length :
                                pass
                        else :
                                print >>sys.stderr, 'WARNING: aligned length is not correct\n' + 'parsed_len:', length, 'real:', len(l[:-1]) ,'\n'+l[:-1]
                elif i==16 :
                        p_alignment = l[:-1]
                        if len(l) == length :
                                pass
                        else :
                                print >>sys.stderr, 'WARNING: aligned length is not correct\n' + 'parsed_len:', length, 'real:', len(l[:-1]) ,'\n'+l[:-1]

        return (parent, child, tm_score, rmsd, aligned_length, identity, p_alignment, c_alignment )

def run_tmalign( pdb1, pdb2, seq1, seq2, TMalign_prog ) :
	tmalign_cmd = "%s %s %s" 
	pp = os.popen( tmalign_cmd%(TMalign_prog, pdb2, pdb1 ) ) #be careful the tmalign has reversed input
	
	result = parse_tm_record( pp )

	if not result :
		return None

	if 'X' in seq1 or 'x' in seq1 :
		aln1 = remove_xmasks_from_aln( seq1, result[6] )
	else :
		aln1 = result[6]

	if 'X' in seq2 or 'x' in seq2 :
		aln2 = remove_xmasks_from_aln( seq2, result[7] )
	else :
		aln2 =  result[7]

	return ( result[2]*len(seq1), result[2], 0, 0, aln1, aln2 )

#END of TMalign parsing code
######################################################

######################################################
# HHsearch parsing code
def parse_hhsearch( content ) :
	for i, l in enumerate(content) :
		if l[:7] == 'Probab=' :
			l = l.split()
			probability = float(l[0][7:])/100
			evalue = float( l[1][8:] )
			content = content[i+1:]
			break


	length = len(content)
	if length%11 == 2 :
		pass
	else :
		print >>sys.stderr, "Error!! The format is not right!!"
		sys.exit()

	nlines = length/11
	aln1 = ''
	aln2 = ''
	start1 = 0
	start2 = 0
	for i in xrange( nlines ) :
		#process query line
		l = content[i*11+3]
		l = l.split()
		if l[0] == 'Q' :
			pass
		else :
			print >>sys.stderr, "Error!! The format is not right!!"
			sys.exit()

		if not aln1 :
			start1 = int(l[2]) -1 #correction for 0 based numbering
		aln1 += l[3]
		#process subject line
		l = content[i*11+3+4]
		l = l.split()
		if l[0] == 'T' :
			pass
		else :
			print >>sys.stderr, "Error!! The format is not right!!"
			sys.stderr.writelines( content )
			sys.exit()

		if not aln2 :
			start2 = int(l[2]) - 1 #correction for 0 based numbering
		aln2 += l[3]
	return ( evalue, probability, start1, start2, aln1, aln2 )

def run_hhsearch( id1_fn, id2_fn ) :
	hhsearch_cmd = "/home/jpei/local/hhsearch1.5_2/hhsearch -i %s -d %s -o stdout -b 1 -B 1 -z 1 -Z 1"
	pp = os.popen( hhsearch_cmd % (id1_fn, id2_fn ) )
	stdout = pp.readlines()
	#sys.stdout.writelines( stdout )
	#content = open( 'temp.hhr' ).readlines()
	return parse_hhsearch( stdout )

#End of Hhsearch parse
#############################################################################

#############################################################################
#COMPASS
#############################################################################
def parse_compass_result( content ) :
	new_content = []
	smith_waterman_score = 0
	negloge = 0

	in_flag = 0
	for l in content :
		if l[:23] == "Smith-Waterman score = " :
			in_flag = 1
			score_line = l.split() 
			try :
				smith_waterman_score = float(score_line[3])
				negloge = -math.log(float( score_line[6] ))
			except :
				print >>sys.stderr, 'WARNING! possible problem in compass score line', score_line
				negloge = -350.0

			continue
		if l[:11] == "Parameters:" :
			in_flag = 0
			break
		if in_flag and l != '\n' and l[0] != ' ' :
			new_content.append( l )

	l = new_content[0].split()
	aln1 = l[2]
	start1 = int( l[1] ) -1 #correction for 0 based numbering

	l = new_content[1].split()
	aln2 = l[2]
	start2 = int( l[1] ) -1 #correction for 0 based numbering

	for i in xrange( 1, len(new_content)/2 ) :
		l = new_content[i*2].split()
		aln1 += l[1]

		l = new_content[i*2+1].split()
		aln2 += l[1]
	return ( smith_waterman_score, negloge, start1, start2, aln1, aln2 )


def run_compass( id1_fn, id2_fn ) :
	compass_cmd = '/home/jpei/local/compass/compass_241_db1Xdb2 -i %s -j %s'
	pp = os.popen( compass_cmd % (id1_fn, id2_fn ) )
	content = pp.readlines()
	return parse_compass_result( content )

def main() :
	name_base = "test"
	cmd_template = "test_run %s > %s"
	data_dirs = [ 'data1', 'data2' ]
	save_dir = 'save_dir'
	pair_list = [ 'abcd abc2', 'abc2 abbb', 'abde adsa' ]
	nruns_per_job = 1
		
	a = alignment_runner( name_base, cmd_template, data_dirs, save_dir, pair_list, nruns_per_job )
	a.write_jobs()
	a.submit_jobs()
	a.wait_for_jobs()

if __name__== '__main__' :
	main()
