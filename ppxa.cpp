#include "ppxa.h"

int sign(double v)
{
    return v > 0 ? 1 : (v < 0 ? -1 : 0);
}

void set_zero(double*v, int size)
{
    for (int i = 0; i< size; i++) v[i] = 0;
}

//iD is a lower trianguler matrix, with diagonal elements and above diagonal equal 1
void *iDxVector (double* vector, int vector_size, double *result)
 {
     result[0] = vector[0];
     for (int i=1; i<vector_size; i++) {result[i] = vector[i] + result[i-1];}
 }

 
void *iDtxVector (double* vector, int vector_size, double *result)
 {
     result[vector_size-1] = vector[vector_size-1];
     for (int i=vector_size-2; i>=0; i--) {result[i] = vector[i] + result[i+1];}
 }

/*
	Function implements codebook update using PPXA+ algorithm
	
	PARAMETERS
	Input: q, data, id, nCodevectors, sData, de, typ
	Output: q
	
	data - pointer to the data
    q - pointer to codevector
	sData - data's size
	nCodevectors - number of codevectors 
	de - delta of projector onto [de, infinity]^{nCodevectors -1} (see [1] for details)
	typ - 1 - quadratic norm, 2 - mean absolute value criterion

*/

void codebook_ppxa_update(
            double *Q_val,
            double * data, // orginal data
            int *id,
            int Q,
            int ps,
            double de,
            int typ
            )
{
    int i = 0, j =0; // counter
    double I = 0;
    int stop = 0; // loop control
    double lambda = 1.5;
    int gamma = 1;

    double om1 = 0.5; 
    double om2 = 0.5;


    double const1 = 1/om1;
    double const2 = 1/(const1+1);
    double prec = 0.001;  //Precision parameter to stop the iterations
    int nit = 1; // number of iterations
    int signum ; // for gsl_linalg_LU_decomp

    
	double *NI = new double [Q];
	
    gsl_vector *suI = gsl_vector_alloc(Q); //suI = zeros(Q,1); % vector with sum of pixels values with label n
    gsl_vector *y = gsl_vector_alloc(Q); //suI = zeros(Q,1); % vector with sum of pixels values with label n
    gsl_vector *yold = gsl_vector_alloc(Q);
    gsl_vector *hv = gsl_vector_alloc(Q); // help vector
    gsl_matrix *iLQ = gsl_matrix_alloc(Q,Q);
    gsl_vector *c = gsl_vector_alloc(Q);

	
	double *t1 = new double [ps]; 
    double *t2 = new double [Q]; 
    double *p1 = new double [ps]; 
    double *p2 = new double [Q]; 
    double *cy = new double [Q]; 
    double *cy1 = new double [Q]; 

    /*-------------- Init --------------*/
    om1 = om1/gamma;
    om2 = om2/gamma;

    set_zero(NI,Q);
    gsl_vector_set_zero(suI);
    double *dp2 = gsl_vector_ptr (suI,0);

    for (i=0; i<ps; i++)
    {
        j = *(id+i);
        *(NI+j) = *(NI+j)+1; // ind = find(il == n); NI(n) = length(ind);
        *(dp2+j) = *(dp2+j) + (double)*(data+i); //muI(n) = (f(ind));
    }

    for (i=0; i<Q; i++) if (*(NI+j)==0) *(NI+j)=1; // male oszustwo

    //init LQ matrix

    gsl_matrix *LQ = gsl_matrix_alloc(Q,Q);
    dp2 = gsl_matrix_ptr (LQ,0,0); //LQ = iD' diag(NI) * iD
    *(dp2+Q-1)=*(NI+Q-1);
    *(dp2+(Q-1)*Q)=*(NI+Q-1);
    for (i=Q-2; i>=0; i--)
    {
        *(dp2+i) = *(NI+i) + *(dp2+i+1); // first raw
        *(dp2+Q*i)= *(dp2+i); // first column
    }
    for (j=1; j<Q; j++)
    {
        for (i=1; i<j; i++)  *(dp2+j*Q+i)=*(dp2+j);
        for (i=j; i<Q; i++)  *(dp2+j*Q+i)=*(dp2+(j-1)*Q+i);
    }
    gsl_matrix_scale(LQ,om1); // LQ = om1*LQ+om2*eye(Q);
    gsl_matrix_add_diagonal(LQ, om2);

    //init iLQ matrix

    gsl_permutation * perm = gsl_permutation_alloc (Q);
    gsl_linalg_LU_decomp (LQ, perm, &signum);
    gsl_linalg_LU_invert (LQ, perm, iLQ);
    gsl_permutation_free (perm);
    gsl_matrix_free (LQ); // we will not use LQ anymore

    // init t1 vector

    for (i=0; i<ps; i++) *(t1+i) = *(data+i);

    // init t2 vector

    for (i=0; i<Q; i++) t2[i]=0;

    // init y vector

    double *dp1 = gsl_vector_ptr (hv,0); //y = iLQ*(om1*iD'*suI+om2*t2);
    dp2 = gsl_vector_ptr (suI,0);
    iDtxVector(dp2,Q,dp1); // hv = iD'*suI
    gsl_vector_scale(hv,om1); // hv = om1*hv;
    //gsl_vector_scale(t2,om2);for non zero vector t2
    //gsl_vector_add(y,t2) // for non zero vector t2
    gsl_blas_dgemv( CblasNoTrans, 1, iLQ, hv, 0, y ); // y = iLQ *hv

    /*-------------- Main loop --------------*/
    stop = 0;
    while (stop == 0)
    {
        if (typ==2)
        {
            for (i=0; i<ps; i++) //p1 = (f/om1+t1)/(1/om1+1);
                *(p1+i) = (*(data+i)/om1+ *(t1+i))* const2;
        }
        else if (typ==1)
        {
            for (i=0; i<ps; i++) // p1 = sign(t1-f).*max(abs(t1-f)-1/om1,0)+f;
            {
                I = *(t1+i) - *(data+i);
                *(p1+i) = sign(I) * MAX(ABS(I)-const1,0) + *(data+i);
            }

        }

        for (i=0; i<Q; i++) p2[i] = t2[i];
        for (i=1; i<Q; i++) {if (t2[i] < de ) p2[i] = de;}

        dp2 = gsl_vector_ptr (suI,0); // suI update
        gsl_vector_set_zero(suI);

        for (i=0; i<ps; i++)
        {
            j = (uint8_t)*(id+i);
            *(dp2+j) = *(dp2+j) + *(p1+i); //suI(n) = (p1(ind));
        }

        dp1 = gsl_vector_ptr (hv,0); //c = iLQ*(om1*iD'*suI+om2*p2);

        iDtxVector(dp2,Q,dp1); //hv = iD'*suI

        for (i=0; i<Q; i++) // hv = om1* hv+om2*p2
             dp1[i] = om1 * dp1[i] + om2 * p2[i];


        gsl_blas_dgemv( CblasNoTrans, 1, iLQ, hv, 0, c ); // c = iLQ *hv

        dp1 = gsl_vector_ptr (c,0); //cy = 2*c-y;
        dp2 = gsl_vector_ptr (y,0);
        for (i=0; i<Q; i++) cy[i] = 2*dp1[i]-dp2[i];

        iDxVector(cy,Q,cy1); //cy1 = iD*cy;

        for (i=0; i<ps; i++) //t1 = t1+lambda*(cy1(il)-p1);
        {
            j = (uint8_t)id[i];
            t1[i] = t1[i] + lambda * (cy1[j] - p1[i]);
        }

        for (i=0; i<Q; i++) t2[i] = t2[i] +lambda*(cy[i]-p2[i]); //t2 = t2+lambda*(cy-p2);

        gsl_blas_dcopy (y,yold);  //yold = y;

        dp1 = gsl_vector_ptr (c,0); //y = y+lambda*(c-y);
        dp2 = gsl_vector_ptr (y,0);
        for (i=0; i<Q; i++) dp2[i] = dp2[i] + lambda * (dp1[i] - dp2[i]);

        gsl_vector_sub(yold,y); //-(y-yold)
        I = gsl_blas_dnrm2(yold); //norm(y-yold)

        if ( nit > 1 && I < prec) stop = 1;
        nit++;
        if ( nit > 5000) stop = 1;

    }

    dp2 = gsl_vector_ptr (y,0);

    iDxVector(dp2,Q,Q_val);

    gsl_vector_free(c);
	
	delete [] cy;
	delete [] cy1;
	delete [] t1;
	delete [] t2;
	delete [] p1;
	delete [] p2;
	delete [] NI;
 
    gsl_vector_free (hv);
    gsl_vector_free (y);
    gsl_vector_free (yold);
    gsl_vector_free (suI);
    gsl_matrix_free (iLQ);

}
