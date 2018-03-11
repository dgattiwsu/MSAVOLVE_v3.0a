/* PSICOV - Protein Sparse Inverse COVariance analysis program */

/* by David T. Jones August 2011 - Copyright (C) 2011 University College London */

/* Version 1.05 - Last Edit 13/2/12 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#define FALSE 0
#define TRUE 1

#define SQR(x) ((x)*(x))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

#define MAXSEQLEN 5000
#define MINSEQS 50
#define MINEFSEQS 50


extern glasso_(int *, double *, double *, int *, int *, int *, int *, double *, int *, double *, double *, int *, double *, int *);


/* Dump a rude message to standard error and exit */
void
                fail(char *errstr)
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
}

/* Convert AA letter to numeric code (0-21) */
int
                aanum(int ch)
{
    const static int aacvs[] =
    {
	999, 0, 3, 4, 3, 6, 13, 7, 8, 9, 21, 11, 10, 12, 2,
	21, 14, 5, 1, 15, 16, 21, 19, 17, 21, 18, 6
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 20);
}

/* Allocate matrix */
void           *allocmat(int rows, int columns, int size)
{
    int             i;
    void          **p, *rp;

    rp = malloc(rows * sizeof(void *) + sizeof(int));

    if (rp == NULL)
	fail("allocmat: malloc [] failed!");

    *((int *)rp) = rows;

    p = rp + sizeof(int);

    for (i = 0; i < rows; i++)
	if ((p[i] = calloc(columns, size)) == NULL)
	    fail("allocmat: malloc [][] failed!");

    return p;
}

/* Free matrix */
void
                freemat(void *rp)
{
    int             rows;
    void **p = rp;

    rows = *((int *)(rp - sizeof(int)));

    while (rows--)
	free(p[rows]);

    free(rp - sizeof(int));
}

/* Allocate vector */
void           *allocvec(int columns, int size)
{
    void          *p;

    p = calloc(columns, size);

    if (p == NULL)
	fail("allocvec: calloc failed!");

    return p;
}

struct sc_entry
{
    float sc;
    int i, j;
} *sclist;

/* Sort descending */
int cmpfn(const void *a, const void *b)
{
    if (((struct sc_entry *)a)->sc == ((struct sc_entry *)b)->sc)
	return 0;

    if (((struct sc_entry *)a)->sc < ((struct sc_entry *)b)->sc)
	return 1;

    return -1;
}
    

int             main(int argc, char **argv)
{
    int             a, b, i, j, k, seqlen, nids, s, nseqs, ncon, opt, ndim, approxflg=0, initflg=0, debugflg=0, diagpenflg=1, apcflg=1, maxit=10000, npair, nnzero, niter, jerr, shrinkflg=1, rawscflg = 1, pseudoc = 1, minseqsep = 5;
    unsigned int *wtcount;
    double thresh=1e-4, del, sum, score, (**pab)[21][21], **pa, wtsum, pc, **pcmat, *pcsum, pcmean, rhodefault = -1.0, lambda, smean, fnzero, lastfnzero, trialrho, rfact, r2, targfnzero = 0.0, scsum, scsumsq, mean, sd, zscore, ppv;    
    float *weight, idthresh = -1.0, maxgapf = 0.9;
    char            buf[4096], seq[MAXSEQLEN], *blockfn = NULL, **aln;
    FILE *ifp;

    while ((opt = getopt(argc, argv, "alnpr:b:i:t:c:g:d:j:")) >= 0)
	switch (opt)
	{
	case 'a':
	    approxflg = 1;
	    break;
	case 'n':
	    shrinkflg = 0;
	    break;
	case 'p':
	    rawscflg = 0;
	    break;
	case 'l':
	    apcflg = 0;
	    break;
	case 'r':
	    rhodefault = atof(optarg);
	    break;
	case 'd':
	    targfnzero = atof(optarg);
	    break;
	case 't':
	    thresh = atof(optarg);
	    break;
	case 'i':
	    idthresh = 1.0 - atof(optarg)/100.0;
	    break;
	case 'c':
	    pseudoc = atoi(optarg);
	    break;
	case 'j':
	    minseqsep = atoi(optarg);
	    break;
	case 'b':
	    blockfn = strdup(optarg);
	    break;
	case 'g':
	    maxgapf = atof(optarg);
	    break;
	case '?':
	    exit(-1);
	}

    if (optind >= argc)
	fail("Usage: psicov [options] alnfile\n\nOptions:\n-a\t: use approximate Lasso algorithm\n-n\t: don't pre-shrink the sample covariance matrix\n-p\t: output PPV estimates rather than raw scores\n-l\t: don't apply APC to Lasso output\n-r nnn\t: set initial rho parameter\n-d nnn\t: set target precision matrix sparsity (default 0 = not specified)\n-t nnn\t: set Lasso convergence threshold (default 1e-4)\n-i nnn\t: select BLOSUM weighting with given identity threshold (default selects threshold automatically)\n-c nnn\t: set pseudocount value (default 1)\n-j nnn\t: set minimum sequence separation (default 5)\n-g nnn\t: maximum fraction of gaps (default 0.9)\n-b file\t: read rho parameter file\n");

    ifp = fopen(argv[optind], "r");
    if (!ifp)
	fail("Unable to open alignment file!");

    for (nseqs=0;; nseqs++)
	if (!fgets(seq, MAXSEQLEN, ifp))
	    break;

    aln = allocvec(nseqs, sizeof(char *));
    
    weight = allocvec(nseqs, sizeof(float));

    wtcount = allocvec(nseqs, sizeof(unsigned int));
    
    rewind(ifp);
    
    if (!fgets(seq, MAXSEQLEN, ifp))
	fail("Bad alignment file!");
    
    seqlen = strlen(seq)-1;

    if (nseqs < MINSEQS)
	fail("Alignment too small - not enough homologous sequences to proceed (or change MINSEQS at your own risk!)");

    if (!(aln[0] = malloc(seqlen)))
	fail("Out of memory!");

    for (j=0; j<seqlen; j++)
	aln[0][j] = aanum(seq[j]);
    
    for (i=1; i<nseqs; i++)
    {
	if (!fgets(seq, MAXSEQLEN, ifp))
	    break;
	
	if (seqlen != strlen(seq)-1)
	    fail("Length mismatch in alignment file!");
	
	if (!(aln[i] = malloc(seqlen)))
	    fail("Out of memory!");
	
	for (j=0; j<seqlen; j++)
	    aln[i][j] = aanum(seq[j]);
    }


    /* Calculate sequence weights */

    if (idthresh < 0.0)
    {
	double meanfracid = 0.0;
	
	for (i=0; i<nseqs; i++)
	    for (j=i+1; j<nseqs; j++)
	    {
		int nids;
		float fracid;
		
		for (nids=k=0; k<seqlen; k++)
		    if (aln[i][k] == aln[j][k])
			nids++;
		
		fracid = (float)nids / seqlen;
		
		meanfracid += fracid;
	    }
	
	meanfracid /= 0.5 * nseqs * (nseqs - 1.0);

	idthresh = 0.38 * 0.32 / meanfracid;
    }

    for (i=0; i<nseqs; i++)
	for (j=i+1; j<nseqs; j++)
	{
	    int nthresh = seqlen * idthresh;

	    for (k=0; nthresh > 0 && k<seqlen; k++)
		if (aln[i][k] != aln[j][k])
		    nthresh--;
	    
	    if (nthresh > 0)
	    {
		wtcount[i]++;
		wtcount[j]++;
	    }
	}

    for (wtsum=i=0; i<nseqs; i++)
	wtsum += (weight[i] = 1.0 / (1 + wtcount[i]));

    if (wtsum < MINEFSEQS)
	puts("\n*** WARNING - not enough sequence variation - or change MINEFSEQS at your own risk! ***\n");

    pa = allocmat(seqlen, 21, sizeof(double));
    pab = allocmat(seqlen, seqlen, 21*21*sizeof(double));

    /* Calculate singlet frequencies with pseudocount */
    for (i=0; i<seqlen; i++)
    {
	for (a=0; a<21; a++)
	    pa[i][a] = pseudoc;
	
	for (k=0; k<nseqs; k++)
	{
	    a = aln[k][i];
	    if (a < 21)
		pa[i][a] += weight[k];
	}
	
	for (a=0; a<21; a++)
	    pa[i][a] /= pseudoc * 21.0 + wtsum;
    }

    /* Calculate pair frequencies with pseudocount */
    for (i=0; i<seqlen; i++)
    {
	for (j=i+1; j<seqlen; j++)
	{
	    for (a=0; a<21; a++)
		for (b=0; b<21; b++)
		    pab[i][j][a][b] = pseudoc / 21.0;

	    for (k=0; k<nseqs; k++)
	    {
		a = aln[k][i];
		b = aln[k][j];
		if (a < 21 && b < 21)
		    pab[i][j][a][b] += weight[k];
	    }
	    
	    for (a=0; a<21; a++)
		for (b=0; b<21; b++)
		{
		    pab[i][j][a][b] /= pseudoc * 21.0 + wtsum;
		    pab[j][i][b][a] = pab[i][j][a][b];

//		    printf("%d/%d %d/%d %f %f %f %f\n", i+1, a, j+1, b, pab[i][j][a][b], pa[i][a] , pa[j][b], pab[i][j][a][b] - pa[i][a] * pa[j][b]);
		}
	}
    }

    for (i=0; i<seqlen; i++)
	for (a=0; a<21; a++)
	    for (b=0; b<21; b++)
		pab[i][i][a][b] = (a == b) ? pa[i][a] : 0.0;

    gsl_matrix *cmat, *rho, *ww, *wwi, *tempmat;

    ndim = seqlen * 21;

    cmat = gsl_matrix_calloc(ndim, ndim);

    /* Form the covariance matrix */
    for (i=0; i<seqlen; i++)
	for (j=0; j<seqlen; j++)
	    for (a=0; a<21; a++)
		for (b=0; b<21; b++)
		    if (i != j)
			gsl_matrix_set(cmat, i*21+a, j*21+b, pab[i][j][a][b] - pa[i][a] * pa[j][b]);
		    else if (a == b)
			gsl_matrix_set(cmat, i*21+a, j*21+b, pab[i][j][a][b] - pa[i][a] * pa[j][b]);

    freemat(pab);

    /* Shrink sample covariance matrix towards shrinkage target F = Diag(1,1,1,...,1) * smean */

    if (shrinkflg)
    {
	for (smean=i=0; i<ndim; i++)
	    smean += gsl_matrix_get(cmat, i, i);
	
	smean /= (float)ndim;
	lambda = 0.1;

//	smean = 1;

	tempmat = gsl_matrix_calloc(ndim, ndim);

	gsl_set_error_handler_off();
	
	for (;;)
	{
	    gsl_matrix_memcpy(tempmat, cmat);
	    
	    /* Test if positive definite using Cholesky decomposition */
	    if (!gsl_linalg_cholesky_decomp(tempmat))
		break;
	    
	    for (i=0; i<seqlen; i++)
		for (j=0; j<seqlen; j++)
		    for (a=0; a<21; a++)
			for (b=0; b<21; b++)
			    if (i != j)
				gsl_matrix_set(cmat, i*21+a, j*21+b, (1.0 - lambda) * gsl_matrix_get(cmat, i*21+a, j*21+b));
			    else if (a == b)
				gsl_matrix_set(cmat, i*21+a, j*21+b, smean * lambda + (1.0 - lambda) * gsl_matrix_get(cmat, i*21+a, j*21+b));
	}

	gsl_matrix_free(tempmat);
    }

    rho = gsl_matrix_alloc(ndim, ndim);
    ww = gsl_matrix_alloc(ndim, ndim);
    wwi = gsl_matrix_alloc(ndim, ndim);

    lastfnzero=0.0;

    /* Guess at a reasonable starting rho value if undefined */
    if (rhodefault < 0.0)
	trialrho = MAX(0.001, 1.0 / wtsum);
    else
	trialrho = rhodefault;

    rfact = 0.0;

    for (;;)
    {
	if (trialrho <= 0.0 || trialrho >= 1.0)
	    fail("Sorry - failed to find suitable value for rho (0 < rho < 1)!");

	gsl_matrix_set_all(rho, trialrho);

	for (i=0; i<seqlen; i++)
	    for (j=0; j<seqlen; j++)
		for (a=0; a<21; a++)
		    for (b=0; b<21; b++)
			if ((a != b && i == j) || pa[i][20] > maxgapf || pa[j][20] > maxgapf)
			    gsl_matrix_set(rho, i*21+a, j*21+b, 1e9);
	
	/* Mask out regions if block-out list provided */
	if (blockfn != NULL)
	{
	    ifp = fopen(blockfn, "r");
	    
	    for (;;)
	    {
		if (fscanf(ifp, "%d %d %lf", &i, &j, &score) != 3)
		    break;
		
		for (a=0; a<21; a++)
		    for (b=0; b<21; b++)
		    {
			gsl_matrix_set(rho, (i-1)*21+a, (j-1)*21+b, score);
			gsl_matrix_set(rho, (j-1)*21+b, (i-1)*21+a, score);
		    }
	    }
	    
	    fclose(ifp);
	}
    
	/* All matrices are symmetric so no need to transpose before/after calling Fortran code */

	glasso_(&ndim, cmat->data, rho->data, &approxflg, &initflg, &debugflg, &diagpenflg, &thresh, &maxit, ww->data, wwi->data, &niter, &del, &jerr);

	if (targfnzero <= 0.0)
	    break;
	
	for (npair=nnzero=i=0; i<ndim; i++)
	    for (j=i+1; j<ndim; j++,npair++)
		if (gsl_matrix_get(wwi, i, j) != 0.0)
		    nnzero++;

	fnzero = (double) nnzero / npair;

//      printf("rho=%f fnzero = %f\n", trialrho, fnzero);

	/* Stop iterating if we have achieved the target sparsity level */
	if (fabs(fnzero - targfnzero)/targfnzero < 0.01)
	    break;
	
	if (fnzero == 0.0)
	{
	    /* As we have guessed far too high, halve rho and try again */
	    trialrho *= 0.5;
	    continue;
	}
	
	if (lastfnzero > 0.0 && fnzero != lastfnzero)
	{
//	    printf("fnzero=%f lastfnzero=%f trialrho=%f oldtrialrho=%f\n", fnzero, lastfnzero, trialrho, trialrho/rfact);
	    
	    rfact = pow(rfact, log(targfnzero / fnzero) / log(fnzero / lastfnzero));

//	    printf("New rfact = %f\n", rfact);
	}

	lastfnzero = fnzero;

	/* Make a small trial step in the appropriate direction */

	if (rfact == 0.0)
	    rfact = (fnzero < targfnzero) ? 0.9 : 1.1;
	
	trialrho *= rfact;
    }

    gsl_matrix_free(rho);
    gsl_matrix_free(ww);

    /* Calculate background corrected scores using average product correction */

    pcmat = allocmat(seqlen, seqlen, sizeof(double));
    pcsum = allocvec(seqlen, sizeof(double));
    
    pcmean = 0.0;
    
    for (i=0; i<seqlen; i++)
	for (j=i+1; j<seqlen; j++)
	{	
	    for (pc=a=0; a<20; a++)
		for (b=0; b<20; b++)
		    pc += fabs(gsl_matrix_get(wwi, i*21+a, j*21+b));

	    pcmat[i][j] = pcmat[j][i] = pc;
	    pcsum[i] += pc;
	    pcsum[j] += pc;

	    pcmean += pc;
	}

    pcmean /= seqlen * (seqlen - 1) * 0.5;

    /* Build final list of predicted contacts */

    sclist = allocvec(seqlen * (seqlen - 1) / 2, sizeof(struct sc_entry));

    for (scsum=scsumsq=ncon=i=0; i<seqlen; i++)
	for (j=i+minseqsep; j<seqlen; j++)
	    if (pcmat[i][j] > 0.0)
	    {
		/* Calculate APC score */
		if (apcflg)
		    sclist[ncon].sc = pcmat[i][j] - pcsum[i] * pcsum[j] / SQR(seqlen - 1.0) / pcmean;
		else
		    sclist[ncon].sc = pcmat[i][j];
		scsum += sclist[ncon].sc;
		scsumsq += SQR(sclist[ncon].sc);
		sclist[ncon].i = i;
		sclist[ncon++].j = j;
	    }

    qsort(sclist, ncon, sizeof(struct sc_entry), cmpfn);

    mean = scsum / ncon;
    sd = sqrt(scsumsq / ncon - SQR(mean));

    /* Print output in CASP RR format with optional PPV estimated from final Z-score */
    if (rawscflg)
	for (i=0; i<ncon; i++)
	    printf("%d %d 0 8 %f\n", sclist[i].i+1, sclist[i].j+1, sclist[i].sc);
    else
	for (i=0; i<ncon; i++)
	{
	    zscore = (sclist[i].sc - mean) / sd;
	    ppv = 0.904 / (1.0 + 16.61 * exp(-0.8105 * zscore));
	    printf("%d %d 0 8 %f\n", sclist[i].i+1, sclist[i].j+1, ppv);
	}
    
    return 0;
}
