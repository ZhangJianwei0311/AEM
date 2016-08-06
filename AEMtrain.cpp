#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "topiclib.h"
#define MAXLEN 50000 // maximum length of the line is 50000 characters
#define LEN 50 // the constant memory size

/* Online Belief Propagation (OBP) algorithm */
void OBP(double ALPHA, double BETA, int W, int J, int D, double *pr, int *ir, int *jc, 
	double *phi, double *phitot, int NNZ, double TD, double TK,double *theta) 
{
	int wi, di, i, j, k, topic, iter, ii, J2 = (int)J, D2 = (int) (TD*D),thredJ=(int)(J*0.1);
	int *ind_rk, *ind_rd,wordNum,idx=1;
	double mutot, totprob, xi, rsum, xtot = 0.0,*thetatot,entropy,temp_entropy,perp;
	double JALPHA = (double) (J*ALPHA), WBETA = (double) (W*BETA);
	double *mu, *rk, *rd, *munew;
        double tempPerplex,totProb;
        double oldPerp=0;

	mu = dvec(J*NNZ);
        thetatot=dvec(D);
	ind_rk = ivec(J*D);
	ind_rd = ivec(D);
	munew = dvec(J);
	rk = dvec(J*D);
	rd = dvec(D);
	/* initialize ind_rk and ind_rd */
	for (di=0; di<D; di++) {
		ind_rd[di] = di;
		for (j=0; j<J; j++) {
			ind_rk[di*J + j] = j;
		}
	}
   
	/* random initialization */
	for (di=0; di<D; di++) {
		for (i=jc[di]; i<jc[di + 1]; i++) {
			wi = (int) ir[i];
			xi = pr[i];
			xtot += xi;
			// pick a random topic 0..J-1
			topic = (int) (J*drand());
			mu[i*J + topic] = 1.0; // assign this word token to this topic
			phi[wi*J + topic] += xi; // increment phi count matrix
			theta[di*J + topic] += xi; // increment theta count matrix
                        thetatot[di]+=xi;
			phitot[topic] += xi; // increment phitot matrix
		}
	}
        
	for (iter=0; iter<500; iter++) {
		/* passing message mu */
		/* iteration 0 */
		if (iter == 0) {//useTime=etime();
			for (di=0; di<D; di++) {
				for (i=jc[di]; i<jc[di + 1]; i++) {
					wi = (int) ir[i];
					xi = pr[i];
					mutot = 0;
					for (j=0; j<J; j++) {
						phi[wi*J + j] -= xi*mu[i*J + j];
						phitot[j] -= xi*mu[i*J + j];
						theta[di*J + j] -= xi*mu[i*J + j];
                                                
						munew[j] = (phi[wi*J + j] + BETA)/(phitot[j] + WBETA)*(theta[di*J + j] + ALPHA);
						mutot += munew[j];
					}
					for (j=0; j<J; j++) {
						munew[j] /= mutot;
						rk[di*J + j] += xi*fabs(munew[j] - mu[i*J + j]);
						rd[di] += xi*fabs(munew[j] - mu[i*J + j]);
						mu[i*J + j] = munew[j];
						phi[wi*J + j] += xi*mu[i*J + j];
						theta[di*J + j] += xi*mu[i*J + j];
                                                phitot[j] += xi*mu[i*J + j];
					}
				}
				psort(rk + di*J, ind_rk + di*J, 0, J - 1, J2);
			}
			psort(rd, ind_rd, 0, D - 1, D2);
                       

		} else { /* iteration > 0 */
		
			if((iter+1)%3==0&&J2>thredJ){
				
				for (ii=0; ii<D2; ii++) {
					di = (int) ind_rd[ii];
					insertionsort(J, rk + di*J, ind_rk + di*J);
					rd[di]=0;
					for (j=0; j<J2; j++) {
						k = (int) ind_rk[di*J + j];
						rk[di*J + k] = 0.0;
                                 	}
				}
				if(J2>thredJ){
					J2=(int)(J2*TK);
				}
				//printf("%d\n",J2);
			}
			for (ii=0; ii<D2; ii++) {
				di = (int) ind_rd[ii];
				for (j=0; j<J2; j++) {				
					k = (int) ind_rk[di*J + j];
					rd[di] -= rk[di*J + k];
					rk[di*J + k] = 0.0;
				}
				for (i=jc[di]; i<jc[di + 1]; i++) {
					wi = (int) ir[i];
					xi = pr[i];
					for (j=0, mutot=0.0, totprob=0.0; j<J2; j++) {
						k = (int) ind_rk[di*J + j];
						phi[wi*J + k] -= xi*mu[i*J + k];
						phitot[k] -= xi*mu[i*J + k];
						theta[di*J + k] -= xi*mu[i*J + k];
                                                	
						totprob += mu[i*J + k];
						munew[k] = (phi[wi*J + k] + BETA)/(phitot[k] + WBETA)*(theta[di*J + k] + ALPHA);
						mutot += munew[k];
					}
					for (j=0; j<J2; j++) {
						k = (int) ind_rk[di*J + j];
						munew[k] /= mutot;
						munew[k] *= totprob;
						rk[di*J + k] += xi*fabs(munew[k] - mu[i*J + k]);
						mu[i*J + k] = munew[k];
						phi[wi*J + k] += xi*mu[i*J + k];
						phitot[k] += xi*mu[i*J + k];
						theta[di*J + k] += xi*mu[i*J + k];
                                                
					}
				}
				for (j=0; j<J2; j++) {
					k = (int) ind_rk[di*J + j];
					rd[di] += rk[di*J + k];
				}
                         //       insertionsort(J, rk + di*J, ind_rk + di*J);
			//	psort(rk + di*J, ind_rk + di*J, 0, J - 1, J2);
			}
                        insertionsort(D, rd, ind_rd);
	//		psort(rd, ind_rd, 0, D - 1, D2);
		}
		
		/* check the convergence condition */
		if((iter+1)%10==0){
			tempPerplex=0;
			wordNum=0;
			for(di=0;di<D;di++){
				for(i=jc[di];i<jc[di+1];i++){
					wi=(int)ir[i];
					xi=pr[i];
					totProb=0;
					for(j=0;j<J;j++){
						totProb+=((theta[di*J+j]+ALPHA)/(thetatot[di]+J*ALPHA)*(phi[wi*J+j]+BETA)/(phitot[j]+WBETA));
					}
					tempPerplex-=(log(totProb)*xi);
					wordNum+=(int)xi;
                          	}
                      	}
			tempPerplex=exp(tempPerplex/(double)wordNum);
			//printf("%f\n",fabs(tempPerplex));
			if(fabs(tempPerplex-oldPerp)<10.0){
                        	printf("break earily,use iteration:%d,perplexity:%f\n",iter,tempPerplex);
                        	break;
			}
			
			oldPerp=tempPerplex;
		}
		//printf("asd\n");
	}

	free(mu);
	free(ind_rk);
	free(ind_rd);
	free(munew);
	
        free(thetatot);
	free(rk);
	free(rd);
}

/*==========================================
* main
*========================================== */
int main(int argc, char* argv[])
{
	int W, K, D, DS, NNZ, numMinibatch, SEED, m, s,i,j;
	int *ir, *jc;
	double TD, ALPHA, BETA;
	double TK;
	//double TH;
	double *phi, *phitot, *pr;
	double *theta,allTime=0,returnTime;
	double WBETA,KALPHA,*returnTheta;
	char *pch, mystring[MAXLEN];
	char *PhiName="topic_word_matrix.txt";
	//char *ThetaName="document_topic_matrix.txt";
	int JD;
	FILE *fp;
    	printf("%d\n",argc);
	// Usage
	if (argc != 9) {
		fprintf(stderr, "usage: %s IN OUT K (DS TD TK) ALPHA BETA SEED\n", argv[0]);
		//exit(-1);
		return -1;
	}

	// Read parameters
	K       = atoi(argv[2]);
	DS      = atoi(argv[3]);
	TD      = atof(argv[4]);
	TK      = atof(argv[5]);
	//TH      = atof(argv[6]);
	ALPHA   = atof(argv[6]);
	BETA    = atof(argv[7]);
	SEED    = 2*atoi(argv[8]) + 1; // seeding only works on uneven numbers

	assert(K > 0);
	assert(DS > 0);
	assert(TD > 0);
	assert(TK > 0);
	//assert(TH > 0);
	assert(ALPHA >= 0);
	assert(BETA >= 0);
	assert(SEED > 0);

	seedMT(SEED);

	// open document file
	fp = fopen(argv[1], "r"); 
	if (fp == NULL) {fprintf(stderr, "Text file cannot be found.\n");return -1;}

	// read the file header
	fgets(mystring, sizeof(mystring), fp);

	// Unix file contains a special space at each line of the file
	pch = strtok(mystring, " ");
	if ((pch != NULL) && (atoi(pch) != 0)){
		D = atoi(pch);
		printf("#Document: %d \n", D);
	} else {
		fprintf(stderr, "File header error.\n");
		//exit(-1);
		return -1;
	}

	pch = strtok(NULL, " ");
	if ((pch != NULL) && (atoi(pch) != 0)){
		W = atoi(pch);
		printf("#Vocabulary: %d \n", W);
	} else {
		fprintf(stderr, "File header error.\n");
		//exit(-1);
		return -1;
	}

	pch = strtok(NULL, " ");
	if ((pch != NULL) && (atoi(pch) != 0)){
		NNZ = atoi(pch);
		printf("#NNZ: %d \n", NNZ);
	} else {
		fprintf(stderr, "File header error.\n");
		//exit(-1);
		return -1;
	}
	if (NNZ == 0){ 
		fprintf(stderr, "Empty file.\n");
		//exit(-1);
		return -1;
	}
    
	/* initializae global parameters */
	phitot  = dvec(K);
	phi     = dvec(K*W);

	// The number of mini-batches
	numMinibatch = (int) (D/DS); 
    
	// copy each mini-batch in the file to sparse matrix
	// run OBP algorithm
	// allocate enough memory for each mini-batch
	// we assume that each fp contains ten times word indices than that on average
	jc     = ivec(DS+1);
	ir     = ivec((int) (DS*(NNZ/D)*LEN)); 
	pr     = dvec((int) (DS*(NNZ/D)*LEN)); 
    	KALPHA=K*ALPHA;
	
	theta = dvec(K*DS);
	for (m=0; m<=numMinibatch; m++) {
        /* read "DS" docs from the file */
		NNZ = 0;
		
		for (s=0; s<DS; s++) {
			jc[s] = NNZ;
			if (fgets(mystring, sizeof(mystring), fp) == NULL) break;
			//printf("5");
			pch = strtok(mystring, " ");
			//printf("4");
			if ((pch != NULL) && (atoi(pch) != 0)) {
				ir[NNZ] = atoi(pch) - 1;
				pch = strtok(NULL, " ");
				pr[NNZ] = atof(pch);
				NNZ++;
			}
			// Unix file contains a special space at each line of the file
			while ((pch != NULL) && (atoi(pch) != 0)) {
				pch = strtok(NULL, " ");
				if ((pch != NULL) && (atoi(pch) != 0)) {
					ir[NNZ] = atoi(pch) - 1;
					pch = strtok(NULL, " ");
					pr[NNZ] = atof(pch);
					NNZ++;
				}
			}
		}
		jc[s] = NNZ;
		
		/* run the learning algorithm */
		for(i=0;i<(K*DS);i++)
			theta[i]=0;
		printf("Minibatch No.%d start!\n",m+1);
		OBP(ALPHA, BETA, W, K, s, pr, ir, jc, phi, phitot, NNZ, TD, TK,theta);
		printf("Minibatch No.%d end!\n",m+1);
     
	}
	free(theta);
	
	fclose(fp);
        
	/* output */
    	write_dvec(K*W, K, W, phi,PhiName);
	free(phi);
	free(phitot);
	free(ir);
	free(jc);
	free(pr);
	return 0;
}

