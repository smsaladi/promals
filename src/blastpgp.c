#include "blastpgp.h"
#include "param.h"

//char uniref90_file[500];

char blastpgp_options[500];

// given a file name, get the blastpgp alignment
// 1. read the alignment file blast_dir/base_name.aln
// 2. if unsuccessful, run blastpgp in blast_dir/
//        then read blastpgp output br file, output it to output_file_name(has a random number)
//        then read output_file_name to a subalign, get rid of positions with gaps in the first sequence
//             the missing ends of blastpgp alignment are filled
//        then purge the subalign by deleting divergent, highly similar, gappy sequences
//        then save the final alignment to a file: query_name.aln
//        then return the final alignment

string tab = "\t";
	  
subalign *get_blastpgp_alignment(char *basename, char *query_name, char *query_seq){

	int i, j, k;

	subalign *x;

	//cout << blast_dir << tab << string(blast_dir) << "\t" << string(basename) << tab << endl;
	k = strlen(blast_dir);
	if(blast_dir[k-1]!='/') {
		assert(k<500-1);
		blast_dir[k] = '/';
		blast_dir[k+1] = '\0';
	}

	//cout << blast_dir << tab << string(blast_dir) << "\t" << string(basename) << tab << endl;
	//cout << ( string(blast_dir) + string(basename)+string(".aln")).c_str() << endl;
	//exit(0);
	x = read_blastpgp_alignment(( string(blast_dir) + string(basename)+string(".aln")).c_str(), query_name, query_seq);
	if(x!=NULL) return x;

	//get_blastpgp_cmd();
	strcpy(blastpgp_cmd, blastpgp_command);
	cout <<  "blastpgp_cmd: " << blastpgp_cmd << endl;
	//set_blastpgp_options("");
	char optionStr[500];
	sprintf(optionStr, " -j %d -e %e -h %e -v 1000 -b 1000 -a 2 -m 6 ", iteration_number, evalue, evalue); 
	set_blastpgp_options(optionStr);
	//set_blastpgp_options("");
	//set_uniref90(uniref90_file);

	x = run_blastpgp(query_name, query_seq);
	return x;

}

// clean blastpgp directory for files/directories beginning with a certain sequene name
void clean_blast_psipred(char *seqname) {

	char cmdname[500];
	
	// sprintf(cmdname, "rm -fr %s%s*", blast_dir, seqname);
	
	//cout << cmdname << endl;
	// system(cmdname);

}

// return 0 if successfully get the blastpgp command in blastpgp_cmd
// otherwise return 1
int get_blastpgp_cmd() {

	int i;

	if(blastpgp_cmd[0]!='\0') {
		return 0;
	}

	int random_int = rand();

	char command[500], blastpgp_file[500];
	sprintf(command, "which blastpgp > blastpgp%d.tmpfile", random_int);
	sprintf(blastpgp_file, "blastpgp%d.tmpfile", random_int);

	i = system(command);

	if(i!=0) {
		cout << "blastpgp is not found in your paths\n" << endl;
		system((string("rm ")+string(blastpgp_file)).c_str());
		return 1;
	}

	ifstream ifp(blastpgp_file, ios::in);
	ifp >> blastpgp_cmd;
	ifp.close();
	
	system((string("rm ")+string(blastpgp_file)).c_str());

	return 0;
}
	
void set_uniref90(char *file1) {

	int i, j;

	strcpy(uniref90_file, file1);

}

void set_blastpgp_options(char *options) {

	int i, j;

	char default_options[500];
	if(strlen(options) == 0 ) {
		strcpy(default_options, " -j 3 -e 0.001 -h 0.001 -v 1000 -b 1000 -a 2 -m 6 ");
		strcpy(blastpgp_options, default_options);
	}
	else {
		strcpy(blastpgp_options, options);
	}
	cout << "blastpgp options: " << blastpgp_options << endl;
}

subalign * run_blastpgp(char *query_name, char *query_seq) {

	int i, j, k;

	char query_OJoj[strlen(query_seq)+1];

	// get the fasta file
	int rand_number = rand();	
	char fasta_file_name[500];
	int check_file_status = 0;
	// check if the fasta file name exists as a file
	while(!check_file_status) {
		sprintf(fasta_file_name, "%s%s.%d.fa", blast_dir, query_name, rand_number);
		FILE *fp = fopen(fasta_file_name, "r");
		if(fp) {
			fclose(fp);
			rand_number = rand();
			continue;
		}
		//fclose(fp);
		check_file_status = 1;
	}

	FILE *fp1 = fopen(fasta_file_name, "w");
	if(!fp1) {
		cout << "cannot write to fasta file " << fasta_file_name << endl;
		exit(0);
	}
	fputs(">", fp1); fputs(query_name, fp1); fputs("\n", fp1);
	for(i=0;i<strlen(query_seq);i++) {
		if( (query_seq[i]=='O') || (query_seq[i]=='o') ||(query_seq[i]=='J') ||(query_seq[i]=='j') ) { 
			fputs("A", fp1); 
			query_OJoj[i] = 'A';
		}
                else {
			fprintf(fp1, "%c", query_seq[i]);
			query_OJoj[i] = query_seq[i];
		}
	}
	//old verison: fputs(query_seq, fp1);
	fputs("\n", fp1);
	fclose(fp1);

	// run blastpgp
	char command[500];
	//sprintf(command, "%s -i %s -o %s.br %s -d %s -C %s%s.chk -M %s", blastpgp_cmd, fasta_file_name, fasta_file_name, blastpgp_options, uniref90_file, blast_dir, query_name, blosum62_file);
	sprintf(command, "%s -i %s -o %s.br %s -d %s -C %s%s.chk", blastpgp_cmd, fasta_file_name, fasta_file_name, blastpgp_options, uniref90_file, blast_dir, query_name);
	//cout << "blastpgp command: " << endl;
	//cout << "   " << command << endl;
	system(command);

	// get blastpgp results
	char br_file_name[500];
	strcpy(br_file_name, (string(fasta_file_name)+string(".br")).c_str());
	subalign *x1 = read_blastpgp_result(br_file_name, query_name, query_seq);
	// purge alignment, removing divergent sequences, higly similar sequences, gappy sequences
	// and keep the maximum sequence number to be 50
	//subalign *x = x1->purge_align(0.25, 0.9, 50, 0.5);
	//subalign *x = x1->purge_align(0.2, 0.9, 300, 0.9);
	subalign *x = x1->purge_align(low_id_thr, 0.9, max_num_sequences, 0.9);
	//cout << "low_id_thr: " << low_id_thr << " max_num_sequences: " << max_num_sequences << endl;

	// output x in the blast_dir directory
	ifstream fp2( (string(blast_dir)+string(query_name)+string(".aln") ).c_str(), ios::in);
	if(!fp2.good()) { x->printali((string(blast_dir)+string(query_name)+string(".aln") ).c_str(), 60); }
	fp2.close();	
	
	// move the checkpoint file to blast_dir
	//sprintf(command, "mv %s.chk %s", query_name, blast_dir);
	//system(command);

	// remove the blast output file
	strcpy(command, "rm -f ");
	strcat(command, br_file_name); 
    // system(command);

	// system( (string("rm -f ") + string(br_file_name) + string(".aln")).c_str());
	// system( (string("rm -f ") + string(fasta_file_name) + string(".aln")).c_str());
	// system( (string("rm -f ") + string(fasta_file_name)).c_str());

	delete x1;

	return x;
}

/*
subalign *get_blastpgp_result(const char *blastpgp_result, char *query_name, char *query_seq) {

	subalign *x = read_blastpgp_result(blastpgp_result, query_name, query_seq);
	
	if(x!=NULL) { return x; }

	x = run_blastpgp(query_name, query_seq);

	if(x==NULL) {
		cout << "Something is wrong with blast file reading"<< endl;
		cout << "returning a null pointer" << endl;
	}

	return x;
}
*/

subalign *read_blastpgp_result(const char *blastpgp_result, char *query_name, char *query_seq) {


   int i,j,k;
   int round = 0;
   char output_aln_file[500];
   char str[300];

   FILE *fp;
   if((fp=fopen(blastpgp_result, "r"))==NULL) {
        //fprintf(stderr, "Error message: Cannot open the blast output file\n\n");
        //fprintf(stderr, "This code reads the blast output file in m6 format\n");
        //fprintf(stderr, "and output the alignment from the last round of iteration\n\n");
        //fprintf(stderr, "Usage: \n");
        fprintf(stderr, "\t getaln brfilename\n\n");
	return NULL;
        // exit(0);
   }
   while(fgets(str, 200, fp)!=NULL) {
      if(strncmp(str, "Results from", 12)==0) {
        round++;
      }
   }
   //cout << "round: " << round << endl;
   rewind(fp);
   while(fgets(str, 200, fp)!=NULL) {
     if(strncmp(str, "Results from", 12)==0) {
        round--;
        if(round==0) break;
     }
   }

   strcpy(output_aln_file, (string(blastpgp_result) + string(".aln")).c_str());

   FILE *fp1 = fopen(output_aln_file, "w");
	
   if(!fp1) {
	cout << "cannot open the temporary alignment file for write " << output_aln_file << endl;
	exit(0);
   }

   int flag = 0;
   int flag1 = 0;
   while(fgets(str, 200, fp)!=NULL) {
     if(strncmp(str, "QUERY ", 6)==0) {
        flag = 1;
	flag1 = 1;
     }
     if(strncmp(str, " ***** No hits found ******", 27)==0) {flag = 0; }
     if(strncmp(str, "  Database: ",12)==0) {
        flag = 0;
     }
     if(flag==1) {
        //fprintf(stdout, "%s", str);
	if(strncmp(str, "Searching...", 12 )==0) continue;
	fprintf(fp1, "%s", str);
     }
   }
   fclose(fp);
   fclose(fp1);
  
   subalign *x;

   if(flag1 == 0) {
	cout << "nothing is read from blast output file " << output_aln_file << endl;
	cout << "use the query sequene as the alignment - " << query_name<< endl;

        x = oneSeq2subalign(query_seq, query_name);
        return x;
   }

   x = new subalign(output_aln_file);

   char *tmpstring = new char [x->alilen+1];
   j = 0;
   for(i=0;i<x->alilen;i++) {
	if(x->aseq[0][i]!='-') {
		tmpstring[j] = x->aseq[0][i];
		j++;
	}
   }
   tmpstring[j] = '\0';

   int right_numbers, left_numbers;

   right_numbers = string(query_seq).find(tmpstring, 0);

   if(right_numbers == string::npos) {
        cout << "Warning: cannot find the sequence in " << query_name << endl;
	cout << query_seq << endl;
	cout << tmpstring << endl;
        return NULL;
   }

   left_numbers =strlen(query_seq) - right_numbers - strlen(tmpstring);

   assert(right_numbers>=0);
   assert(left_numbers >=0);

   // the first sequence of blast matches the query sequence exactly
   if( (right_numbers+left_numbers) ==  0) {
	purge_alignment_by_first_sequence(x);
	return x;  
   }

   cout << "|" << query_seq << "|" << endl;

   string Nterm, Cterm, Nterm_gap, Cterm_gap;
   Nterm = string(query_seq).substr(0, right_numbers);
   Cterm = string(query_seq).substr(right_numbers+strlen(tmpstring), left_numbers);
   Nterm_gap = string("");
   for(i=0;i<Nterm.size();i++) Nterm_gap += "-";
   Cterm_gap = string("");
   for(i=0;i<Cterm.size();i++) Cterm_gap += "-";

   cout << "|" << Nterm << "|" << endl;
   cout << "|" << Cterm << "|" << endl;
   //cout << "N-cap: |" << Nterm_gap << "|" << endl;
   //cout << "C-cap: |" << Cterm_gap << "|" << endl;

   char **tmpseq = cmatrix(x->nal, x->alilen+right_numbers+left_numbers+1);
   for(i=0;i<x->nal;i++) {
	if(i==0) {
		strcpy(tmpseq[i], Nterm.c_str());
		strcat(tmpseq[i], x->aseq[i]);
	}
	else {
		strcpy(tmpseq[i], Nterm_gap.c_str());
		strcat(tmpseq[i], x->aseq[i]);
	}
   }
   //cout << "nstr: " << Nterm_gap.c_str() << endl;
   //cout << "nstr: " << Nterm.c_str() << endl;
   if( left_numbers!=0 ) {
	//cout << "cstr: " <<  Cterm_gap.c_str() << endl;
	//cout << "cstr: " <<  Cterm.c_str() << endl;
	for(i=0;i<x->nal;i++) {
		if(i==0) {
			strcat(tmpseq[i], Cterm.c_str());
			tmpseq[i][x->alilen+right_numbers+left_numbers] = '\0';
		}
		else {
			strcat(tmpseq[i], Cterm_gap.c_str());
			tmpseq[i][x->alilen+right_numbers+left_numbers] = '\0';
		}
	}
   }
   for(i=0;i<x->nal;i++) {
	delete [] x->aseq[i];
	x->aseq[i] = tmpseq[i];
	//cout << x->aseq[i] << endl;
   }

   //cout << "x->alilen: " << x->alilen << endl;
   x->alilen = x->alilen + right_numbers + left_numbers;
   //cout << "x->alilen: " << x->alilen << endl;
   for(i=1;i<=x->nal;i++) {
	delete [] x->alignment[i];
	x->alignment[i] = new int [x->alilen+1];
	for(j=1;j<=x->alilen;j++) {
                x->alignment[i][j] = am2num(x->aseq[i-1][j-1]);
           }
   }

   //x->printali(60);

   purge_alignment_by_first_sequence(x);

   // delete the blastpgp output alignment file
   char command[500];
   strcpy(command, "rm -f ");
   strcat(command, output_aln_file);
   system(command);

   return x;
}

int check_blastpgp_alignment(const char *blastpgp_alignment) {

	int i, j, k;

	FILE *fp;
	char tmpstring[10000];

	if((fp = fopen(blastpgp_alignment, "r"))==NULL) {
		cout << blastpgp_alignment << " does not exist" << endl;
		return 0;
	}
	i = 0;
	fclose(fp);

	// check if the file contains only empty lines
	ifstream f1(blastpgp_alignment, ios::in);
	char *ss;
	while(f1.good()) {
		f1 >> tmpstring;
		for(ss=tmpstring;isspace(*ss);ss++);
		if(tmpstring) return 1;
	}
	f1.close();

	return 0;

}

subalign *read_blastpgp_alignment(const char *blastpgp_alignment, char *query_name, char *query_seq) {


   int i,j,k;
  
   subalign *x;

   // file does not exist; file contain only empty lines
   if(check_blastpgp_alignment(blastpgp_alignment) == 0) {
	return NULL;
   }

   x = new subalign(blastpgp_alignment);

   char *tmpstring = new char [x->alilen+1];
   j = 0;
   for(i=0;i<x->alilen;i++) {
	if(x->aseq[0][i]!='-') {
		tmpstring[j] = x->aseq[0][i];
		j++;
	}
   }
   tmpstring[j] = '\0';

   int right_numbers, left_numbers;

	//cout << query_seq << endl;
	//cout << tmpstring << endl;

   right_numbers = string(query_seq).find(tmpstring, 0);

   if(right_numbers == string::npos) {
        cout << "Warning: cannot find the sequence in " << query_name << endl;
	cout << query_seq << endl;
	cout << tmpstring << endl;
        return NULL;
   }


   left_numbers =strlen(query_seq) - right_numbers - strlen(tmpstring);

   //cout << "right_numbers: " << right_numbers << endl;
   //cout << "left_numbers: " << left_numbers << endl;

   assert(right_numbers>=0);
   assert(left_numbers >=0);

   // the first sequence of blast matches the query sequence exactly
   if( (right_numbers+left_numbers) ==  0) return x;  

   string Nterm, Cterm, Nterm_gap, Cterm_gap;
   Nterm = string(query_seq).substr(0, right_numbers);
   Cterm = string(query_seq).substr(right_numbers+strlen(tmpstring), left_numbers);
   Nterm_gap = string("");
   for(i=0;i<Nterm.size();i++) Nterm_gap += "-";
   Cterm_gap = string("");
   for(i=0;i<Cterm.size();i++) Cterm_gap += "-";

   char **tmpseq = cmatrix(x->nal, x->alilen+right_numbers+left_numbers+1);
   for(i=0;i<x->nal;i++) {
	if(i==0) {
		strcpy(tmpseq[i], Nterm.c_str());
		strcat(tmpseq[i], x->aseq[i]);
	}
	else {
		strcpy(tmpseq[i], Nterm_gap.c_str());
		strcat(tmpseq[i], x->aseq[i]);
	}
   }
   //cout << "nstr: " << Nterm_gap.c_str() << endl;
   //cout << "nstr: " << Nterm.c_str() << endl;
   if( left_numbers!=0 ) {
	//cout << "cstr: " <<  Cterm_gap.c_str() << endl;
	//cout << "cstr: " <<  Cterm.c_str() << endl;
	for(i=0;i<x->nal;i++) {
		if(i==0) {
			strcat(tmpseq[i], Cterm.c_str());
			tmpseq[i][x->alilen+right_numbers+left_numbers] = '\0';
		}
		else {
			strcat(tmpseq[i], Cterm_gap.c_str());
			tmpseq[i][x->alilen+right_numbers+left_numbers] = '\0';
		}
	}
   }
   for(i=0;i<x->nal;i++) {
	delete [] x->aseq[i];
	x->aseq[i] = tmpseq[i];
	//cout << x->aseq[i] << endl;
   }

   //cout << "x->alilen: " << x->alilen << endl;
   x->alilen = x->alilen + right_numbers + left_numbers;
   //cout << "x->alilen: " << x->alilen << endl;
   for(i=1;i<=x->nal;i++) {
	delete [] x->alignment[i];
	x->alignment[i] = new int [x->alilen+1];
	for(j=1;j<=x->alilen;j++) {
                x->alignment[i][j] = am2num(x->aseq[i-1][j-1]);
           }
   }

   //x->printali(60);

   purge_alignment_by_first_sequence(x);

   return x;
}

void purge_alignment_by_first_sequence(subalign *x) {

	int i, j, k;

	int count = 0;
	for(i=0;i<x->alilen;i++) {
		if(x->aseq[0][i]!='-') {count++;}
	}

	if(count == x->alilen) return;


        char **tmpseq = cmatrix(x->nal, count);

	count = 0;
	for(i=0;i<x->alilen;i++) {
		if(x->aseq[0][i]!='-') {
			for(j=0;j<x->nal;j++) {
				tmpseq[j][count] = x->aseq[j][i];
			}
			count++;
		}
	}
	//cout << "count: " << count << endl;
	for(i=0;i<x->nal;i++) {
		delete [] x->aseq[i];
		x->aseq[i] = tmpseq[i];
		x->aseq[i][count] = '\0';
	}

	x->alilen = count;

	for(i=1;i<=x->nal;i++) {
        	delete [] x->alignment[i];
        	x->alignment[i] = new int [x->alilen+1];
        	for(j=1;j<=x->alilen;j++) {
                	x->alignment[i][j] = am2num(x->aseq[i-1][j-1]);
           	}
   	}
}
