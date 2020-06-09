float **map_structure_fast(seq_str_aln *a1, seq_str_aln *a2) {

	int i, j, k;

	int **salnpos; // structure alignment positions; two arrays of integer values
	// [0..1][0..numalignedp-1]
	// Example:
	// aCKE-EKaa
 	// -CEEKEKa-
	// array 1: 2 3 4 5 6
	// array 2: 1 2 3 5 6
	// if no structural alignment with z-score > 5, salnpos = 0
	int numalignedp;

	float **tmpMat = gmatrix<float>(a1->len, a2->len);
	for(i=1;i<=a1->len;i++) for(j=1;j<=a2->len;j++) tmpMat[i][j] = 0;
	// p1 and p2 are the correponding position numbers of the queries
	int p1, p2;
	//cout << "here " << endl;
	//cout << a1->query << endl;
	//cout << a2->query << endl;
	//cout << a1->nhits << " " << a2->nhits << endl;
	////cout << a1->id[1] << "  ----  " <<  a2->id[1] << endl;
	int tmpMat_count = 0;
	for(i=1;i<=a1->nhits;i++) {
		for(j=1;j<=a2->nhits;j++) {
			//cout << "HHHHH" << endl;
			//cout << a1->id[i] << "  ----  " <<  a2->id[j] << endl;
			salnpos = get_salnpos_fast(a1->id[i], a2->id[j], numalignedp);
			if(!salnpos) continue;
			//cout << numalignedp << endl;
			for(k=0;k<numalignedp;k++) {
				//cout << "1: " << salnpos[0][k] <<  " " << a1->aln[i][salnpos[0][k]] << endl;
				p1 = a1->aln[i][salnpos[0][k]];
				if(!p1) {
                                        //cout << "not p1" << endl;
                                        continue;
                                }
				//cout << "2: " << endl; cout << salnpos[1][k] <<  endl; cout <<  " " << a2->aln[j][salnpos[1][k]] << endl;
				p2 = a2->aln[j][salnpos[1][k]];
				if(!p2) continue;
				if(tmpMat[p1][p2]==0) {
					tmpMat[p1][p2] = 1;
					tmpMat_count++;
				}
			}
			free_imatrix(salnpos, 1, numalignedp);
		}
	}
	//cout << "tmpMat_count: " << tmpMat_count << endl;
	if(tmpMat_count) return tmpMat;
	else {
		free_gmatrix<float>(tmpMat, a1->len, a2->len);
		return NULL;
	}
        //cout << "This place" << endl;

