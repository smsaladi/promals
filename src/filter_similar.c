#include "multiple.h"
#include "sequences.h"
#include "btree_template.h"
#include "progressiveAlignHMM.h"
#include "multiple.h"
#include "constraint.h"
#include "param.h"

void run_cdhit(char *filefilter, float thr);

// filename is the name of the input file
// allseqs will be changed
void do_filter_similar(sequences &allseqs, char *filename, float thr, vector<subalign *> &similaraln) {

        int i, j, k;
        char command[500];

        // 1. output the sequences to a new fasta file with seq names potentially changed
        char filefilter[500];
        sprintf(filefilter, "%s_fiLtEr", filename);
        allseqs.output_fasta(filefilter);
        
        // 2. use cdhit to get the clusters
        print_section_info("Below filter sequences by cd-hit");
        cout << "\tnumber of input sequences: " << allseqs.nseqs << endl;
        run_cdhit(filefilter, thr);
        //sprintf(command, "%s -i %s -o %s.clus -c %f -d 30 -aL 0.95 -p 1 1>/dev/null 2>/dev/null", cdhit, filefilter, filefilter, thr);
        //printinfo(command, 0);
        //system(command);
        //exit(0);

        // 3. read clstr file and aln similar sequences within each cluster
        char clstrfile[500];
        sprintf(clstrfile, "%s.clus.clstr", filefilter);
        ifstream fp(clstrfile, ios::in);
        if(!fp) {
                cout << "cannot open cd-hit clstr file " << clstrfile << endl;
                exit(0);
        }
        int clustercount = 0;
        char tmpstr[1000];
        char *s, *ss;
        char *seqname;
        vector<char *> seqnames;
        seqnames.clear();
        //vector<subalign *> similaraln;
        while( fp.getline(tmpstr, 1000) ) {
                // cluster line
                if(tmpstr[0] == '>') {
                        clustercount++;
                        if(clustercount==1) continue;
                        if(seqnames.size()>=2) if(exclude_similar==0)  
                                similaraln.push_back( get_similar_aln(seqnames, allseqs, filefilter) );
                        //for(i=0;i<seqnames.size();i++) cout << i << " " << seqnames[i] << endl;
                        //cout << "------------" << endl;
                        for(i=0;i<seqnames.size();i++) {
                                //cout << seqnames[i] << endl;
                                delete [] seqnames[i];
                        }
                        //cout << "=========" << endl;
                        seqnames.clear();
                        continue;
                }
                // below are for sequence lines
                for(s=tmpstr; *s!='>'; s++);
                s++;
                for(ss=s;*ss!=' ';ss++);
                ss--;ss--;ss--;
                *ss = '\0';
                seqname = new char [26];
                //strcpy(seqname, s);
                strncpy(seqname, s, 26); // pay attention to strncpy ending, '\0' is not automatically added at the end
                //cout << s << endl;
                //cout << "|" << seqname << "|" << endl;
                seqnames.push_back(seqname);
        }
        fp.close();
        if(seqnames.size()>=2) if(exclude_similar==0) 
                similaraln.push_back( get_similar_aln(seqnames, allseqs, filefilter) );
        for(i=0;i<seqnames.size();i++) delete [] seqnames[i];
        seqnames.clear();

        // 4. modify the allseqs object, so that only cluster representatives are kept
        char clusfile[500];
        sprintf(clusfile, "%s.clus", filefilter);
        allseqs.seq.clear();
        allseqs.name.clear();
        allseqs.readFasta(clusfile, 1, 0);
        //allseqs.printSeqs();
 
        // 5. clear up the temporary files
        sprintf(command, "rm -f %s*", filefilter);
        //system(command);

        // 6. return the similaraln vector
        // debug
        //cout << "similaraln size: " << similaraln.size() <<endl;
        print_section_info("Below are simiaraln");
        for(i=0;i<similaraln.size();i++) {
                cout << "aln " << i << endl;
                similaraln[i]->printali(100);
        }
        //return similaraln;
}

// run cdhit to get the cluster files, dynamically adjust the similarity threshold util
// the number of clusters is smaller than a certain threshold (10000)
void run_cdhit(char *filefilter, float thr) {
        
        int i, j, k;
        float coption=0.95, aLoption=0.95;
        int noption=5;

        // first use the defined thr
        if(thr >= 0.7) noption =5;
        else if(thr>=0.6) noption =4;
        else if(thr>=0.5) noption = 3;
        else if(thr>=0.4) noption = 2;
        else {
                cout << "cd-hit threshold should be no less than 0.4" << endl;
                exit(0);
        }
        char command[500];
        sprintf(command, "%s -i %s -o %s.clus -c %f -d 30 -aL 0.95 -p 1 -n %d 1>%s.info 2>/dev/null", cdhit, filefilter, filefilter, thr, noption, filefilter);
        system(command);
        // check the number of clusters
        char info[500];
        char tmpstr[500];
        int clusternum;
        sprintf(info, "%s.info", filefilter);
        ifstream fp(info, ios::in);
        while(fp.good() ) {
                fp >> tmpstr;
                if(strcmp(tmpstr, "finished")==0) {
                        fp >> clusternum;
                }
        }
        fp.close(); fp.clear();
        if(clusternum <= 10000) {
                cout << "\tc: " << thr << " aL: " << aLoption << " n: " << noption << endl;
                cout << "\tnumber of clusters: " << clusternum << endl;
                return;
        }

        for(i=1;i<100;i++) {
                coption = 1 - 0.05 * i;
                if(coption < 0.4) break;
                thr = coption;
                if(thr >= 0.7) noption =5;
                else if(thr>=0.6) noption =4;
                else if(thr>=0.5) noption = 3;
                else if(thr>=0.4) noption = 2;

                if(coption>=0.85) {
                        aLoption = coption;
                }
                else {
                        aLoption = 0.8;
                }
                char command[500];
                sprintf(command, "%s -i %s -o %s.clus -c %f -d 30 -aL %f -p 1 -n %d 1>%s.info 2>/dev/null", cdhit, filefilter, filefilter, coption, aLoption, noption, filefilter);
                system(command);
                printinfo(command, 0);
                fp.open(info, ios::in);
                while(fp.good() ) {
                        fp >> tmpstr;
                        if(strcmp(tmpstr, "finished")==0) {
                                fp >> clusternum;
                        }
                }
                fp.close(); fp.clear();
                cout << "\tc: " << coption << " aL: " << aLoption << " n: " << noption << endl;
                cout << "\tnumber of clusters: " << clusternum << endl;
                if(clusternum <= 10000) {
                        //cout << "number of clusters: " << clusternum << endl;
                        return;
                }
        }

        if(clusternum>10000) {
                cout << "c: " << coption << " aL: " << aLoption << " n: " << noption << endl;
                cout << "number of clusters: " << clusternum << " is too large" << endl;
                exit(0);
        }
        
}

subalign *get_similar_aln(vector<char *> seqnames, sequences &allseqs, char *filefilter) {

        int i, j, k;

        // write the sequences to a fasta file
        char individualclstrfilename[500];
        sprintf(individualclstrfilename, "%s.fa", filefilter);
        ofstream ofp(individualclstrfilename, ios::out);
        if(!ofp) {
                cout << "cannot open the individualclstrfilename " << individualclstrfilename << endl;
                exit(0);
        }
        int mark;
        for(i=0;i<seqnames.size();i++) {
                // find the sequence with the name seqnames[i] in allseqs
                mark = 0;
                for(j=1;j<allseqs.seq.size();j++) {
                        if(strcmp(seqnames[i], allseqs.name[j].c_str())==0) {
                                mark = 1;
                                break;
                        }
                }
                if(mark==0) {
                        cout <<"cannot find " << seqnames[i] << " in allseqs" << endl;
                        exit(0);
                }
                
                // now write the sequence to the individualclstrfilename file
                ofp << ">" << seqnames[i] << endl;
                ofp << allseqs.seq[j] << endl;
        }
        ofp.close();

        // run mafft on the fasta file
        char command[500];
        sprintf(command, "%s --auto %s 1>%s.aln 2>/dev/null", mafft, individualclstrfilename, individualclstrfilename);
        printinfo(command, 0);
        system(command);
        char individualaln[500];
        sprintf(individualaln, "%s.aln", individualclstrfilename);

        // obtain the subalign object of the alignment file
        subalign *tmpaln = new subalign(individualaln, "fasta", 25);
        return tmpaln;
}

void get_rep_names(sequences allseqs, vector<subalign *> &similaraln, vector<char *> &repnames) {

        int i, j, k;

        int found;
        for(i=0;i<similaraln.size();i++) {
                found=0;
                //for(k=1;k<=similaraln[i]->nal;k++) { cout << similaraln[i]->aname[k-1] << endl; }
                for(j=1;j<allseqs.seq.size();j++) {
                        for(k=1;k<=similaraln[i]->nal;k++) {
                                if(strcmp(allseqs.name[j].c_str(), similaraln[i]->aname[k-1])==0) {
                                        repnames.push_back(similaraln[i]->aname[k-1]);
                                        found=1;
                                        break;
                                }
                        }
                        if(found) break;
                }
                if(!found) {
                        cout << "cannot find the rep name for similaraln number " << i << endl;
                        exit(0);
                }
        }
}
                                                
                

