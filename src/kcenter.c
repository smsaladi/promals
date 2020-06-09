#include "header_cpp.c"

// retricted kcenter method
// simmat: simlarity matrix of dimXdim, nc: number of clusters; 
// maxelem: maximum number of elements in a cluster
int **kcenter(int **simmat, int dim, int nc, int maxelem, int maxiter) {

        int i, j, k, iter, maxclusterindex, maxsim;
        int sum_maxsumsim, maxsumsim, sumsim, tmpvalue;

        assert(nc < dim);

        // allocate memory
        int **clusters = imatrix(nc, maxelem);
        int *clustersize = ivector(nc);

        // initialize the cluster centers
        // initial cluster centers are first nc elements
        for(i=1;i<=nc;i++) {
                clusters[i][1] = i;
        }

    for(iter=1;iter<=maxiter;iter++) {
    
        // assign elements to centers
        for(i=1;i<=nc;i++) clustersize[i] = 1;
        for(i=1;i<=dim;i++) {
                maxclusterindex = -1;
                maxsim = -100000;
                for(j=1;j<=nc;j++) {
                        if(clustersize[i] == maxelem) continue;
                        if(simmat[i][j]>maxsim) {
                                maxsim = simmat[i][j];
                                maxclusterindex = j;
                        }
                }
                clustersize[maxclusterindex]++;
                clusters[maxclusterindex][clustersize[maxclusterindex]] = j;
        }

        // update the cluster centers
        sum_maxsumsim = 0;
        for(i=1;i<=nc;i++) {
                maxsumsim = -1000000;
                sumsim = 0;
                for(j=1;j<=clustersize[i];j++) {
                        sumsim = 0;
                        for(k=1;k<=clustersize[i];k++) {
                                sumsim += simmat[clusters[i][j]][clusters[i][k]];
                        }
                        if(sumsim>maxsumsim) {
                                maxsumsim = sumsim;
                                centerindex = j;
                        }
                }
                sum_maxsumsim += maxsumsim;
                tmpvalue = clusters[1];
                clusters[1] = clusters[centerindex];
                clusters[centerindex] = tmpvalue;
        }

        // print information
        cout << "after round " << iter << " : " << sum_maxsumsim << endl;
        for(i=1;i<=nc;i++) {
                cout << "cluster " << i << endl;
                for(j=1;j<=clusterisze[i];j++) {
                        cout << "\t" << clusters[i][j] << endl;
                }
        }
                        
    }
}

