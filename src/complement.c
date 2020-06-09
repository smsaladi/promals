#include "all.h"

int isp(int i, int include_histidine) {
       if(include_histidine) {
               if(i>=18) return 1;
               else return 0;
       }
       else {
               if(i>=19) return 1;
               else return 0;
       }
} 

int isn(int i) {
       if(i==16) return 1;
       else if(i==17) return 1;
       else return 0;
} 

int isc(int i, int include_histidine) {
       if(include_histidine) {
               if(i>=16) return 1;
               else return 0;
       }
       else {
               if(i>=19) return 1;
               else if(i==16) return 1;
               else if(i==17) return 1;
               else return 0;
       }
} 

int main(int argc, char **argv) {

        int i, j, k, l;

        subalign *taln = new subalign(argv[1]);

        taln->prof(1);
        taln->prof_freq = dmatrix(taln->prof_len, 20);
        taln->pseudoCounts(taln->prof_effn, taln->prof_nef, taln->prof_len, taln->prof_freq);

        //for(i=1;i<=taln->alilen;i++) { for(j=1;j<=20;j++) { fprintf(stdout, "%4.3f ",  taln->prof_freq[i][j]); } fprintf(stdout, "\n"); }

        // now calculate frequencies of positive and negative residues in each position
        double **profreq = taln->prof_freq;
        float positivefreq[taln->alilen+1];
        float negativefreq[taln->alilen+1];
        for(i=1;i<=taln->alilen;i++) {
                positivefreq[i] = profreq[i][19]+profreq[i][20];
                negativefreq[i] = profreq[i][16]+profreq[i][17];
        }
        int nal = taln->nal;
        int mark[taln->nal+1];
        int **alignment = taln->alignment;
        float totalweight, pw1, nw1, pw2, nw2, ow1, ow2; // weights of total marked sequences, positive/negative amino acids in two positions
        float p1n2w, n1p2w; // positive/negative pairs
        double *prof_hwt_all = taln->prof_hwt_all;
        float logodds;
        float othercombo, complement, background_complement, background_othercombo;
        int restricted_to_both_charge = 1;
        int include_histidine = 0;
        for(i=1;i<=taln->alilen;i++) {
                for(j=i+1;j<=taln->alilen;j++) {
                        // mark sequences where both positions are amino acids
                        totalweight = pw1 = nw1 = pw2 = nw2 = ow1 = ow2 = 0;
                        p1n2w = n1p2w = 0;
                        for(k=1;k<=nal;k++) {
                                if( (alignment[k][i]) && (alignment[k][j]) ) {
                                        if(restricted_to_both_charge) {
                                                if(!isc(alignment[k][i], include_histidine) ) continue;
                                                if(!isc(alignment[k][j], include_histidine) ) continue;
                                        }
                                        mark[k] = 1;
                                        totalweight += prof_hwt_all[k];
                                        if( isn(alignment[k][i]) ) {
                                                nw1 += prof_hwt_all[k];
                                        }
                                        else if( isp(alignment[k][i], include_histidine) ) {
                                                pw1 += prof_hwt_all[k];
                                        }
                                        //else ow1 += prof_hwt_all[k];
                                        if( isn(alignment[k][j]) ) {
                                                nw2 += prof_hwt_all[k];
                                        }
                                        else if( isp(alignment[k][j], include_histidine) ) {
                                                pw2 += prof_hwt_all[k];
                                        }
                                        //else ow2 += prof_hwt_all[k];
                                        if(isn(alignment[k][i]) && isp(alignment[k][j], include_histidine) ) n1p2w += prof_hwt_all[k];
                                        if( isp(alignment[k][i], include_histidine)&&isn(alignment[k][j]) ) p1n2w += prof_hwt_all[k];

                                }
                                else {mark[k] = 0;}
                        }
                        if(totalweight==0) continue;
                        pw1 /= totalweight;
                        nw1 /= totalweight;
                        ow1 /= totalweight;
                        pw2 /= totalweight;
                        nw2 /= totalweight;
                        ow2 /= totalweight;
                        p1n2w /= totalweight;
                        n1p2w /= totalweight;

                        complement = p1n2w + n1p2w;
                        othercombo = 1 - p1n2w - n1p2w;
                        background_complement = pw1*nw2 + nw1*pw2;
                        background_othercombo = 1 - pw1*nw2 - nw1*pw2;

                        if(background_complement==0) continue;
                        if(pw1<0.1) continue;
                        if(pw2<0.1) continue;
                        if(nw1<0.1) continue;
                        if(nw2<0.1) continue;

                        if(totalweight<0.2) continue;
                        
                        fprintf(stdout, "%-5d %-5d %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %7.4f\n", i, j, totalweight, pw1, nw1, pw2, nw2, n1p2w, p1n2w, n1p2w + p1n2w, pw1*nw2 + nw1*pw2, log( complement/background_complement/(othercombo/background_othercombo) ) );
                        continue;

                        
                        if(pw1<0.1) continue;
                        if(pw2<0.1) continue;
                        if(nw1<0.1) continue;
                        if(nw2<0.1) continue;
                        if(pw1+nw1 < 0.3) continue;
                        if(pw2+nw2 < 0.3) continue;
                        if(pw1*nw2 + nw1*pw2<0.2) continue;
                        if(pw1*nw2 + nw1*pw2==0) continue;
                        if(n1p2w + p1n2w==0) continue;
                        fprintf(stdout, "%d %d %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f\n", i, j, pw1, nw1, pw2, nw2, n1p2w, p1n2w, n1p2w + p1n2w, pw1*nw2 + nw1*pw2, log( (n1p2w + p1n2w)/(pw1*nw2 + nw1*pw2) ) );
                }
        }
                                        


}
