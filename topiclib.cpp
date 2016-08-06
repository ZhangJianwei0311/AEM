/*---------------------------------------------------
* file:    topiclib.c
* purpose: various utility routines for c code
* version: 1.0
* author:  newman@uci.edu
* date:    12/20/05
*-------------------------------------------------*/

// grep //$ topiclib.c

#include "topiclib.h"
#include "math.h"
/*------------------------------------------
* private to util.c
*------------------------------------------ */

static int icomp(const void *, const void *); /* comparison for isort */
static int dcomp(const void *, const void *); /* comparison for dsort */
static int    *icomp_vec;                     /*  data used for isort */
static double *dcomp_vec;                     /*  data used for dsort */

/*------------------------------------------
* allocation routines
* imat
* dmat
* ivec
* dvec
*------------------------------------------ */

int **imat(int nr, int nc) //
{
	int ntot = nr*nc;
	int *tmp = (int*) calloc(ntot,sizeof(int));
	int **x  = (int**)calloc(nr,sizeof(int*));
	int r;
	assert(tmp);
	assert(x);
	for (r = 0; r < nr; r++) x[r] = tmp + nc*r;
	return x;
}

void free_imat(int **x) //
{
	free(x[0]);
	free(x);
}

double **dmat(int nr, int nc) //
{
	int ntot = nr*nc;
	double *tmp = (double*) calloc(ntot,sizeof(double));
	double **x  = (double**)calloc(nr,sizeof(double*));
	int r;
	assert(tmp);
	assert(x);
	for (r = 0; r < nr; r++) x[r] = tmp + nc*r;
	return x;
}

void free_dmat(double **x) //
{
	free(x[0]);
	free(x);
}

float **fmat(int nr, int nc) //
{
	int ntot = nr*nc;
	float *tmp = (float*) calloc(ntot,sizeof(float));
	float **x  = (float**)calloc(nr,sizeof(float*));
	int r;
	assert(tmp);
	assert(x);
	for (r = 0; r < nr; r++) x[r] = tmp + nc*r;
	return x;
}

void free_fmat(float **x) //
{
	free(x[0]);
	free(x);
}

int *ivec(int n) //
{
	int *x = (int*)calloc(n,sizeof(int));
	assert(x);
	return x;
}

double *dvec(int n) //
{
	double *x = (double*)calloc(n,sizeof(double));
	assert(x);
	return x;
}

float *fvec(int n) //
{
	float *x = (float*)calloc(n,sizeof(float));
	assert(x);
	return x;
}

/*------------------------------------------
* vector routines
* imax
* dmax
*------------------------------------------ */

int imax(int n, int *x) //
{
	int i, xmax=x[0];
	for (i = 0; i < n; i++) xmax = MAX(xmax,x[i]);
	return xmax;
}

double dmax(int n, double *x) //
{
	int i;
	double  xmax=x[0];
	for (i = 0; i < n; i++) xmax = MAX(xmax,x[i]);
	return xmax;
}

int imin(int n, int *x) //
{
	int i, xmin=x[0];
	for (i = 0; i < n; i++) xmin = MIN(xmin,x[i]);
	return xmin;
}

double dmin(int n, double *x) //
{
	int i;
	double  xmin=x[0];
	for (i = 0; i < n; i++) xmin = MIN(xmin,x[i]);
	return xmin;
}

int isum(int n, int *x) //
{
	int i, xsum=0;
	for (i = 0; i < n; i++) xsum += x[i];
	return xsum;
}

double dsum(int n, double *x) //
{
	int i;
	double xsum=0;
	for (i = 0; i < n; i++) xsum += x[i];
	return xsum;
}

double ddot(int n, double *x, double *y) //
{
	int i;
	double xdot=0;
	for (i = 0; i < n; i++) xdot += x[i]*y[i];
	return xdot;
}

/*------------------------------------------
* countlines
* 
*------------------------------------------ */
int countlines(char *fname) //
{
	int lines = 0;
	char buf[BUFSIZ];
	FILE *fp = fopen(fname ,"r"); assert(fp);
	while (fgets(buf, BUFSIZ, fp)) lines++;
	fclose(fp);
	lines -= 3; // less 3 header lines
	assert(lines>0);
	return lines;
}

int countntot(char *fname) //
{
	int i, count, ntot = 0;
	char buf[BUFSIZ];
	FILE *fp = fopen(fname ,"r"); assert(fp);
	for (i = 0; i < 3; i++) fgets(buf, BUFSIZ, fp); // skip 3 header lines
	while (fscanf(fp, "%*d%*d%d", &count) != EOF) ntot += count;
	fclose(fp);
	assert(ntot>0);
	return ntot;
}

/*------------------------------------------
* sort: call qsort library function
* isort
* dsort
* psort
* insertionsort
*------------------------------------------ */

int partition(double *x, int *indx, int left, int right, int pivotIndx) 
{
	double pivotValue = x[indx[pivotIndx]];
	int i, tmp, storeIndx = left;

	/* move pivot to end */
	tmp = indx[pivotIndx];
	indx[pivotIndx] = indx[right];
	indx[right] = tmp;

	for (i=left; i<right; i++) {
		if (x[indx[i]] > pivotValue) {
			tmp = indx[i];
			indx[i] = indx[storeIndx];
			indx[storeIndx] = tmp;
			storeIndx++;
		}
	}

	/* move pivot to storeIndx */
	tmp = indx[storeIndx];
	indx[storeIndx] = indx[right];
	indx[right] = tmp;

	return storeIndx;
}

/* *index should be in the range [0, length(x)] */
void psort(double *x, int *indx, int left, int right, int k)
{
	int pivotIndx, pivotNewIndx;

	if (right > left) {
		pivotIndx = (int) ((right - left)*drand());
		pivotNewIndx = partition(x, indx, left, right, pivotIndx);
		if (pivotNewIndx > left + k) {
			psort(x, indx, left, pivotNewIndx - 1, k);
		}
		if (pivotNewIndx < left + k) {
			psort(x, indx, pivotNewIndx + 1, right, k + left - pivotNewIndx - 1); 
		}
	}
}

void insertionsort(int n, double *x, int *indx) //
{
	int tmp, i, k;

	for (i=1; i<n; i++) {
		for (k=i; k>0 && x[indx[k]]>x[indx[k-1]]; k--) {
			tmp = indx[k];
			indx[k] = indx[k-1];
			indx[k-1] = tmp;
		}
	}
}

void isort(int n, int *x, int direction, int *indx) //
{
	int i;
	assert(direction*direction==1);
	icomp_vec = ivec(n);
	for (i=0; i<n; i++) {
		icomp_vec[i] = direction*x[i];
		indx[i] = i;
	}
	qsort(indx,n,sizeof(int),icomp);
	free(icomp_vec);
}

static int icomp(const void *pl, const void *p2)
{
	int i = * (int *) pl;
	int j = * (int *) p2;
	return (icomp_vec[i] - icomp_vec[j]);
}

void dsort(int n, double *x, int direction, int *indx) //
{
	int i;
	assert(direction*direction==1);
	dcomp_vec = dvec(n);
	for (i=0; i<n; i++) {
		dcomp_vec[i] = direction*x[i];
		indx[i] = i;
	}
	qsort(indx, n, sizeof(int), dcomp);
	free(dcomp_vec);
}

static int dcomp(const void *p1, const void *p2)
{
	int i = * (int *) p1;
	int j = * (int *) p2;
	return dcomp_vec[i] > dcomp_vec[j] ? 1:-1;
}


void read_docwordcountbin(int nnz, int *w, int *d, int *c, char *fname) //
{
	FILE *fp;
	int chk;
	fp = fopen(fname,"r"); assert(fp);
	fscanf(fp,"%d", &chk);
	fscanf(fp,"%d", &chk);
	fscanf(fp,"%d", &chk);
	chk = fread(w,sizeof(int),nnz,fp); assert(chk==nnz);
	chk = fread(d,sizeof(int),nnz,fp); assert(chk==nnz);
	chk = fread(c,sizeof(int),nnz,fp); assert(chk==nnz);
	fclose(fp);
}


int countnnz(int nr, int nc, int **x) //
{
	int i, j, nnz=0;
	for (i = 0; i < nr; i++) 
		for (j = 0; j < nc; j++) 
			if (x[i][j] > 0) nnz++;
	return nnz;
}

void write_sparse(int nr, int nc, int **x, char *fname) //
{
	FILE *fp = fopen(fname, "w");
	int i, j;
	assert(fp);
	fprintf(fp, "%d\n", nr);
	fprintf(fp, "%d\n", nc);
	fprintf(fp, "%d\n", countnnz(nr,nc,x));
	for (i = 0; i < nr; i++) 
		for (j = 0; j < nc; j++) 
			if (x[i][j] > 0) fprintf(fp, "%d %d %d\n", i+1 , j+1 , x[i][j]);
	fclose(fp);
}

void write_sparsebin(int nr, int nc, int **x, char *fname) //
{
	int i, j, k, chk;
	int nnz   = countnnz(nr,nc,x);
	int *col1 = ivec(nnz);
	int *col2 = ivec(nnz);
	int *col3 = ivec(nnz);
	FILE *fp  = fopen(fname,"w"); assert(fp);
	for (i = 0, k = 0; i < nr; i++) 
		for (j = 0; j < nc; j++) 
			if (x[i][j] > 0) {
				col1[k] = i;
				col2[k] = j;
				col3[k] = x[i][j];
				k++;
			}
			assert(k==nnz);
			fwrite(&nr, sizeof(int),1,fp);
			fwrite(&nc, sizeof(int),1,fp);
			fwrite(&nnz,sizeof(int),1,fp);
			chk = fwrite(col1,sizeof(int),nnz,fp); assert(chk==nnz);
			chk = fwrite(col2,sizeof(int),nnz,fp); assert(chk==nnz);
			chk = fwrite(col3,sizeof(int),nnz,fp); assert(chk==nnz);
			fclose(fp);
			free(col1);
			free(col2);
			free(col3);
}

int **read_sparse(char *fname, int *nr_, int *nc_) //
{
	FILE *fp = fopen(fname,"r");
	int i, j, c, nr, nc, nnz;
	int **x;
	assert(fp);
	fscanf(fp,"%d", &nr);  assert(nr>0);
	fscanf(fp,"%d", &nc);  assert(nc>0);
	fscanf(fp,"%d", &nnz); assert(nnz>0);
	x = imat(nr,nc);
	while (fscanf(fp, "%d%d%d", &i, &j, &c) != EOF) {
		i--;
		j--;
		assert(i<nr);
		assert(j<nc);
		assert(c>0);
		x[i][j] = c;
	}
	fclose(fp);
	*nr_ = nr;
	*nc_ = nc;
	return x;
}

void read_dw(char *fname, int *d, int *w, int *D, int *W) //
{
	int i,wt,dt,ct,count,nnz;
	FILE *fp = fopen(fname ,"r"); assert(fp);
	count = 0;
	fscanf(fp,"%d", D);    assert(*D>0);
	fscanf(fp,"%d", W);    assert(*W>0);
	fscanf(fp,"%d", &nnz); assert(nnz>0);
	while (fscanf(fp, "%d%d%d", &dt, &wt, &ct) != EOF) {
		for (i = count; i < count+ct; i++) {
			w[i] = wt-1;
			d[i] = dt-1;
		}
		count += ct;
	}
	fclose(fp);
}

void fill_dtot(int ntot, int *d, int *dtot) //
{
	int i;
	for (i = 0; i < ntot; i++) dtot[d[i]]++;
}

void read_dwc(char *fname, int *d, int *w, int *c, int *D, int *W) //
{
	FILE *fp = fopen(fname,"r");
	int i=0, dd, ww, cc, nnz;
	assert(fp);
	fscanf(fp,"%d", D);    assert(*D>0);
	fscanf(fp,"%d", W);    assert(*W>0);
	fscanf(fp,"%d", &nnz); assert(nnz>0);
	while (fscanf(fp, "%d%d%d", &dd, &ww, &cc) != EOF) {
		d[i] = --dd;
		w[i] = --ww;
		c[i] = cc;
		i++;
	}
	assert(i==nnz);
	fclose(fp);
}

int read_nnzbin(char *fname) //
{
	int nnz;
	FILE *fp = fopen(fname,"r"); assert(fp);
	assert(fread(&nnz,sizeof(int),1,fp)); // nr
	assert(fread(&nnz,sizeof(int),1,fp)); // nc
	assert(fread(&nnz,sizeof(int),1,fp)); // nnz
	fclose(fp);
	return nnz;
}

int read_nnz(char *fname) //
{
	int nnz;
	FILE *fp = fopen(fname,"r"); assert(fp);
	fscanf(fp,"%d", &nnz); // nr
	fscanf(fp,"%d", &nnz); // nc
	fscanf(fp,"%d", &nnz); // nnz
	fclose(fp);
	return nnz;
}

void read_sparsebin(char *fname, int *col1, int *col2, int *col3) //
{
	int nr, nc, nnz, chk;
	FILE *fp = fopen(fname,"r"); assert(fp);
	assert(fread(&nr, sizeof(int),1,fp)); assert(nr>0);
	assert(fread(&nc, sizeof(int),1,fp)); assert(nc>0);
	assert(fread(&nnz,sizeof(int),1,fp)); assert(nnz>0);
	chk = fread(col1,sizeof(int),nnz,fp); assert(chk==nnz);
	chk = fread(col2,sizeof(int),nnz,fp); assert(chk==nnz);
	chk = fread(col3,sizeof(int),nnz,fp); assert(chk==nnz);
	fclose(fp);
}

void write_ivec (int n, int *x, char *fname) //
{
	FILE *fp = fopen(fname,"w");
	int i;
	assert(fp);
	for (i = 0; i < n; i++)  fprintf(fp, "%d\n", x[i]+1 );
	fclose(fp);
}

void read_ivec (int n, int *x, char *fname) //
{
	FILE *fp = fopen(fname,"r");
	int i;
	assert(fp);
	for (i = 0; i < n; i++)  { fscanf(fp, "%d", x+i ); x[i]--; }
	fclose(fp);
}

void write_dvec (int n, int nrow, int ncol, double *x, char *fname) //
{
	FILE *fp = fopen(fname,"ab");
	int i,j;
	printf("1");
	assert(fp);
	printf("2");
	fprintf(fp, "%d %d\n", nrow, ncol);
	for (i = 0; i < nrow; i++){
		for (j=0;j< ncol;j++)
	       fprintf(fp, "%f ", x[j*nrow+i] );
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void read_dvec (int n, double *x, char *fname) //
{
	FILE *fp = fopen(fname,"r");
	int i;
	assert(fp);
	for (i = 0; i < n; i++)  { fscanf(fp, "%f", x+i ); }
	fclose(fp);
}

/*------------------------------------------
* randperm
*------------------------------------------ */

int *randperm(int n) //
{
	int *order = ivec(n);
	int k, nn, takeanumber, temp;
	for (k=0; k<n; k++) order[ k ] = k;
	nn = n;
	for (k=0; k<n; k++) {
		// take a number between 0 and nn-1
		takeanumber = (int) (nn*drand());
		temp = order[ nn-1 ];
		order[ nn-1 ] = order[ takeanumber ];
		order[ takeanumber ] = temp;
		nn--;
	}
	return order;
}

/*------------------------------------------
* randomassignment
*------------------------------------------ */
void randomassignment(int ntot, int T, int *w, int *d, int *z, int **wp, int **dp, int *ztot) //
{
	int i, t;
	for (i = 0; i < ntot; i++) {
		t = (int)(T*drand());
		z[i] = t;
		wp[w[i]][t]++;
		dp[d[i]][t]++;
		ztot[t]++;
	}
}

void randomassignment_rank(int ntot, int T, int *w, int *d, int *drank, int *z, int **wp, int **dp, int *ztot) //
{
	int i, t;
	for (i = 0; i < ntot; i++) {
		t = (int)(T*drand());
		z[i] = t;
		wp[w[i]][t] += drank[d[i]];
		dp[d[i]][t] += drank[d[i]];
		ztot[t] += drank[d[i]];
	}
}

void assignment(int ntot, int *w, int *d, int *z, int **wp, int **dp, int *ztot) //
{
	int i, t;
	for (i = 0; i < ntot; i++) {
		t = z[i];
		wp[w[i]][t]++;
		dp[d[i]][t]++;
		ztot[t]++;
	}
}

void randomassignment_2layer(int ntot, int T, int S, int *w, int *d, int *z, int *y, int **wp, int **zy, int **dp, int *ztot, int *ytot) //
{
	int i, t, s;
	for (i = 0; i < ntot; i++) {
		t = (int)(T*drand());
		s = (int)(S*drand());
		z[i] = t;
		y[i] = s;
		wp[w[i]][t]++;
		zy[t][s]++;
		dp[d[i]][s]++;
		ztot[t]++;
		ytot[s]++;
	}
}

void randomassignment2(int ntot, int T, int *d, int *z, int **dp) //
{
	int i, t;
	for (i = 0; i < ntot; i++) {
		t = (int)(T*drand());
		z[i] = t;
		dp[d[i]][t]++;
	}
}

/*------------------------------------------
* sample_chain
*------------------------------------------ */
void sample_chain (int ntot, int W, int T, double alpha, double beta, int *w, int *d, int *z, int **wp, int **dp, int *ztot, int *order) //
{
	int ii, i, t;
	double totprob, maxprob, currprob;
	double *probs = dvec(T);
	double wbeta = W*beta;

	for (ii = 0; ii < ntot; ii++) {

		i = order[ ii ];

		t = z[i];      // take the current topic assignment to word token i
		ztot[t]--;     // and substract that from the counts
		wp[w[i]][t]--;
		dp[d[i]][t]--;
		totprob = 0;

		for (t = 0; t < T; t++) {
			probs[t] = (wp[w[i]][t] +  beta)/(ztot[t]+  wbeta)*(dp[d[i]][t]+  alpha);
			totprob += probs[t];
		}

		maxprob  = drand()*totprob;
		currprob = probs[0];
		t = 0;

		// sample a topic t from the distribution
		while (maxprob>currprob) {
			t++;
			currprob += probs[t];
		}

		z[i] = t;      // assign current word token i to topic t
		wp[w[i]][t]++; // and update counts
		dp[d[i]][t]++;
		ztot[t]++;
	}

	free(probs);
}

/*------------------------------------------
* sample_chain_rank
*------------------------------------------ */
void sample_chain_rank (int ntot, int W, int T, double alpha, double beta, int *w, int *d, int *drank, int *z, int **wp, int **dp, int *ztot, int *order) //
{
	int ii, i, t;
	double totprob, maxprob, currprob;
	double *probs = dvec(T);
	double wbeta = W*beta;

	for (ii = 0; ii < ntot; ii++) {

		i = order[ ii ];

		t = z[i];      // take the current topic assignment to word token i
		ztot[t] -= drank[d[i]];
		wp[w[i]][t] -= drank[d[i]];
		dp[d[i]][t] -= drank[d[i]];
		totprob = 0;

		for (t = 0; t < T; t++) {
			probs[t] = (wp[w[i]][t] +  beta)/(ztot[t]+  wbeta)*(dp[d[i]][t]+  alpha);
			totprob += probs[t];
		}

		maxprob  = drand()*totprob;
		currprob = probs[0];
		t = 0;

		// sample a topic t from the distribution
		while (maxprob>currprob) {
			t++;
			currprob += probs[t];
		}

		z[i] = t;      // assign current word token i to topic t
		wp[w[i]][t] += drank[d[i]];
		dp[d[i]][t] += drank[d[i]];
		ztot[t] += drank[d[i]];
	}

	free(probs);
}

void sample_chain0 (int ntot, int W, int T, double alpha, double beta, int *w, int *d, int *z, int **wp, int **dp, int *ztot) //
{
	int i, t;
	double totprob, maxprob, currprob;
	double *probs = dvec(T);
	double wbeta = W*beta;


	for (i = 0; i < ntot; i++) {

		t = z[i];      // take the current topic assignment to word token i
		ztot[t]--;     // and substract that from the counts
		wp[w[i]][t]--;
		dp[d[i]][t]--;

		for (t = 0, totprob = 0.0; t < T; t++) {
			probs[t] = (dp[d[i]][t] + alpha) * (wp[w[i]][t] + beta) / (ztot[t] + wbeta);
			totprob += probs[t];
		}

#if 1
		assert(totprob > 0.0);
#endif

		maxprob  = drand()*totprob;
		currprob = probs[0];
		t = 0;

		// sample a topic t from the distribution
		while (maxprob>currprob) {
			t++;
			currprob += probs[t];
		}

		z[i] = t;      // assign current word token i to topic t
		wp[w[i]][t]++; // and update counts
		dp[d[i]][t]++;
		ztot[t]++;
	}

	free(probs);  
}

void sample_chain_2layer (int ntot, int W, int T, int S, double alpha, double beta, double gamma, int *w, int *d, int *z, int *y, int **wp, int **zy, int **dp, int *ztot, int *ytot) //
{
	int    i, t, s;
	double totprob, maxprob, currprob, term1, term2, term3;
	double wbeta  = W*beta;
	double tgamma = T*gamma;
	double **probs = dmat(T,S);

	for (i = 0; i < ntot; i++) {

		t = z[i];
		s = y[i];
		ztot[t]--;
		ytot[s]--;
		wp[w[i]][t]--;
		zy[t][s]--;
		dp[d[i]][s]--;

		totprob = 0;      
		for (t = 0; t < T; t++) {
			for (s = 0; s < S; s++) {
				term1 = (wp[w[i]][t] + beta) / (ztot[t] + wbeta);
				term2 = (zy[t][s]  + gamma)  / (ytot[s] + tgamma);
				term3 = (dp[d[i]][s] + alpha);
				probs[t][s] = term1*term2*term3;
				totprob += probs[t][s];
			}
		}

		maxprob  = drand()*totprob;
		currprob = probs[0][0];
		t = 0;
		s = 0;
		while (maxprob>currprob) {
			t++;
			if (t >= T) { s++; t=0; }
			currprob += probs[t][s];
		}

		z[i] = t;
		y[i] = s;
		ztot[t]++;
		ytot[s]++;
		wp[w[i]][t]++;
		zy[t][s]++;
		dp[d[s]][t]++;
	}

	free_dmat(probs);

}

void resample_chain (int ntot, int W, int T, double alpha, double beta, int *w, int *d, int *z, int **wp, int **dp, int *ztot) //
{
	int i, t;
	double totprob, maxprob, currprob;
	double *probs = dvec(T);
	double wbeta = W*beta;

	for (i = 0; i < ntot; i++) {

		t = z[i];
		dp[d[i]][t]--;
		totprob = 0;

		for (t = 0; t < T; t++) {
			probs[t] = (wp[w[i]][t] +  beta)/(ztot[t] + wbeta)*(dp[d[i]][t] + alpha);
			totprob += probs[t];
		}

		maxprob  = drand()*totprob;
		currprob = probs[0];
		t = 0;

		// sample a topic t from the distribution
		while (maxprob>currprob) {
			t++;
			currprob += probs[t];
		}

		z[i] = t;
		dp[d[i]][t]++;
	}

	free(probs);
}

void oversample_dp (int ntot, int W, int T, double alpha, double beta, int *w, int *d, int *z, int **wp, int **dp, int *ztot) //
{
	int i, t, k, ntimes=4;
	double totprob, maxprob, currprob;
	double *probs = dvec(T);
	double wbeta = W*beta;

	for (i = 0; i < ntot; i++) {

		totprob = 0;
		for (t = 0; t < T; t++) {
			probs[t] = (wp[w[i]][t] +  beta)/(ztot[t] + wbeta)*(dp[d[i]][t] + alpha);
			totprob += probs[t];
		}

		for (k = 0; k < ntimes; k++) {
			maxprob  = drand()*totprob;
			currprob = probs[0];
			t = 0;
			while (maxprob>currprob) {
				t++;
				currprob += probs[t];
			}
			dp[d[i]][t]++;
		}
	}

	free(probs);
}

void loglike (int ntot, int W, int D, int T, double alpha, double beta, int *w, int *d, int **wp, int **dp, int *ztot, int *dtot) //
{
	int    i, j, t;
	double llike;
	static int init = 0;
	static double **prob_w_given_t;
	static double **prob_t_given_d;
	static double *dtot_;
	double ztot_;

	if (init==0) {
		init = 1;
		prob_w_given_t = dmat(W,T);
		prob_t_given_d = dmat(D,T);
		dtot_ = dvec(D);
		for (j = 0; j < D; j++) dtot_[j] = dtot[j] + T*alpha;
	}

	for (t = 0; t < T; t++) {
		ztot_ = ztot[t] + W*beta;
		for (i = 0; i < W; i++) prob_w_given_t[i][t] = (wp[i][t]+beta) / ztot_;
		for (j = 0; j < D; j++) prob_t_given_d[j][t] = (dp[j][t]+alpha)/ dtot_[j];
	}

	llike = 0;
	for (i = 0; i < ntot; i++)
		llike += log(ddot(T, prob_w_given_t[w[i]], prob_t_given_d[d[i]]));

	printf(">>> llike = %.6e    ", llike);
	printf("pplex = %.4f\n", exp(-llike/ntot));
}

void chksum (int n, int T, int **x, int *sumx) //
{
	int i, t, sum;
	for (t = 0; t < T; t++) {
		sum = 0;
		for (i = 0; i < n; i++) sum += x[i][t];
		assert(sum==sumx[t]);
	}
}

void getztot (int n, int T, int **x, int *ztot) //
{
	int i, t, sum;
	for (t = 0; t < T; t++) {
		sum = 0;
		for (i = 0; i < n; i++) sum += x[i][t];
		ztot[t] = sum;
	}
}

double etime() //
{
	static double last_clock = 0;
	static double now_time = 0;
	last_clock = now_time;
	now_time = (double) clock ();
	return (double) (now_time - last_clock) / CLOCKS_PER_SEC;
}

/* digamma funciton implemented by David Blei */
double digamma(double x)
{
	double p;
	x=x+6;
	p=1/(x*x);
	p=(((0.004166666666667*p-0.003968253986254)*p+
		0.008333333333333)*p-0.083333333333333)*p;
	p=p+log(x)-0.5/x-1/(x-1)-1/(x-2)-1/(x-3)-1/(x-4)-1/(x-5)-1/(x-6);
	return p;
}

typedef unsigned long uint32;

#define N              (624)                 // length of state vector
#define M              (397)                 // a period parameter
#define K              (0x9908B0DFU)         // a magic constant
#define hiBit(u)       ((u) & 0x80000000U)   // mask all but highest   bit of u
#define loBit(u)       ((u) & 0x00000001U)   // mask all but lowest    bit of u
#define loBits(u)      ((u) & 0x7FFFFFFFU)   // mask     the highest   bit of u
#define mixBits(u, v)  (hiBit(u)|loBits(v))  // move hi bit of u to hi bit of v

static uint32   state[N+1];     // state vector + 1 extra to not violate ANSI C
static uint32   *next;          // next random value is computed from here
static int      left = -1;      // can *next++ this many times before reloading


void seedMT(uint32 seed)
 {
    //
    // We initialize state[0..(N-1)] via the generator
    //
    //   x_new = (69069 * x_old) mod 2^32
    //
    // from Line 15 of Table 1, p. 106, Sec. 3.3.4 of Knuth's
    // _The Art of Computer Programming_, Volume 2, 3rd ed.
    //
    // Notes (SJC): I do not know what the initial state requirements
    // of the Mersenne Twister are, but it seems this seeding generator
    // could be better.  It achieves the maximum period for its modulus
    // (2^30) iff x_initial is odd (p. 20-21, Sec. 3.2.1.2, Knuth); if
    // x_initial can be even, you have sequences like 0, 0, 0, ...;
    // 2^31, 2^31, 2^31, ...; 2^30, 2^30, 2^30, ...; 2^29, 2^29 + 2^31,
    // 2^29, 2^29 + 2^31, ..., etc. so I force seed to be odd below.
    //
    // Even if x_initial is odd, if x_initial is 1 mod 4 then
    //
    //   the          lowest bit of x is always 1,
    //   the  next-to-lowest bit of x is always 0,
    //   the 2nd-from-lowest bit of x alternates      ... 0 1 0 1 0 1 0 1 ... ,
    //   the 3rd-from-lowest bit of x 4-cycles        ... 0 1 1 0 0 1 1 0 ... ,
    //   the 4th-from-lowest bit of x has the 8-cycle ... 0 0 0 1 1 1 1 0 ... ,
    //    ...
    //
    // and if x_initial is 3 mod 4 then
    //
    //   the          lowest bit of x is always 1,
    //   the  next-to-lowest bit of x is always 1,
    //   the 2nd-from-lowest bit of x alternates      ... 0 1 0 1 0 1 0 1 ... ,
    //   the 3rd-from-lowest bit of x 4-cycles        ... 0 0 1 1 0 0 1 1 ... ,
    //   the 4th-from-lowest bit of x has the 8-cycle ... 0 0 1 1 1 1 0 0 ... ,
    //    ...
    //
    // The generator's potency (min. s>=0 with (69069-1)^s = 0 mod 2^32) is
    // 16, which seems to be alright by p. 25, Sec. 3.2.1.3 of Knuth.  It
    // also does well in the dimension 2..5 spectral tests, but it could be
    // better in dimension 6 (Line 15, Table 1, p. 106, Sec. 3.3.4, Knuth).
    //
    // Note that the random number user does not see the values generated
    // here directly since reloadMT() will always munge them first, so maybe
    // none of all of this matters.  In fact, the seed values made here could
    // even be extra-special desirable if the Mersenne Twister theory says
    // so-- that's why the only change I made is to restrict to odd seeds.
    //

    register uint32 x = (seed | 1U) & 0xFFFFFFFFU, *s = state;
    register int    j;

    for(left=0, *s++=x, j=N; --j;
        *s++ = (x*=69069U) & 0xFFFFFFFFU);
 }


uint32 reloadMT(void)
 {
    register uint32 *p0=state, *p2=state+2, *pM=state+M, s0, s1;
    register int    j;

    if(left < -1)
        seedMT(4357U);

    left=N-1, next=state+1;

    for(s0=state[0], s1=state[1], j=N-M+1; --j; s0=s1, s1=*p2++)
        *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);

    for(pM=state, j=M; --j; s0=s1, s1=*p2++)
        *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);

    s1=state[0], *p0 = *pM ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);
    s1 ^= (s1 >> 11);
    s1 ^= (s1 <<  7) & 0x9D2C5680U;
    s1 ^= (s1 << 15) & 0xEFC60000U;
    return(s1 ^ (s1 >> 18));
 }


inline uint32 randomMT(void)
 {
    uint32 y;

    if(--left < 0)
        return(reloadMT());

    y  = *next++;
    y ^= (y >> 11);
    y ^= (y <<  7) & 0x9D2C5680U;
    y ^= (y << 15) & 0xEFC60000U;
    y ^= (y >> 18);
    return(y);
 }
