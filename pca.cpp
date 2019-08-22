
#include "pca.h"
/*
Copyright 2012 Anna Jezierska
Function finds a linear operator T using Principal Component Analysis of data [ch1, ch2,ch3]
INPUT: pointer data to channels : ch1, ch2, ch3 (vectors of length given by "size"), size
OUTPUT: T (pointer to the matrix 3x3)

Algorithm:
1. Find the mean of ch1, ch2 and ch3.
2. Substruct the mean from the data
3. Find the covariance matrix A
4. Find the proincipal components using SVD of A
*/


void getPCA (
			  double ** data, // pointer to data : ch1, ch2, ch3
              int size,   // channel size
              double * T // PCA linear operator
              )
{
    double *ch1 = data[0];
    double *ch2 = data[1];
    double *ch3 = data[2];
	//  Find the mean of ch1, ch2 and ch3.
	double m1 =0, m2 = 0, m3 = 0;
	for (int i = 0; i<size; i++)
	{
		m1+= ch1[i]/size;
		m2+= ch2[i]/size;
		m3+= ch3[i]/size;
	}
    
	// Substruct the mean from the data
	gsl_matrix *I = gsl_matrix_alloc(3,size);
	
	for (int i = 0; i<size; i++)
	{
		gsl_matrix_set (I, 0, i, ch1[i]-m1);
		gsl_matrix_set (I, 1, i, ch2[i]-m2);
		gsl_matrix_set (I, 2, i, ch3[i]-m3);
	}
	
	// Find the covariance matrix A
	gsl_matrix *A = gsl_matrix_alloc(3,3);
	gsl_blas_dgemm (CblasNoTrans,CblasTrans, 1,I, I, 0, A); //A = I * I^T 
	
	
	// Find the proincipal components using SVD of A
	gsl_matrix *V = gsl_matrix_alloc(3,3);
	gsl_vector *work = gsl_vector_alloc(3);
	gsl_vector *S = gsl_vector_alloc(3);
	
	gsl_matrix_set_all (V, 1.0);
	
	gsl_linalg_SV_decomp (A, V, S, work);
	
    // Set result to T
	for (int i = 0; i <3; i++)
	{
		for (int j = 0; j <3; j++)
		{
			T[i*3+j] = gsl_matrix_get(V,i,j);
		}
	}
	
	// ELEMENTS OF VECTOR S CORRESPOND TO THE VARIANCE
	
	gsl_matrix_free (I);
	gsl_matrix_free (A);
	gsl_matrix_free (V);
	gsl_vector_free (S);
    gsl_vector_free (work);
}

/*
Convert data [ch1, ch2, ch3] using linear operator T
Operation is done in place.
Input: pointer data to channels : ch1, ch2, ch3, size T
Output: ch1, ch2, ch3
*/
void D2T (
			  double ** data,
              int size,
              double * T
			)
{
    
	double *ch1 = data[0];
    double *ch2 = data[1];
    double *ch3 = data[2];
    double T1 = T[0];
	double T2 = T[1];
	double T3 = T[2];
	double T4 = T[3];
	double T5 = T[4];
	double T6 = T[5];
	double T7 = T[6];
	double T8 = T[7];
	double T9 = T[8];
	
	for (int i = 0; i< size; ++i)
    {
		double d1 = ch1[i];
		double d2 = ch2[i];
		double d3 = ch3[i];
		
        ch1[i] = d1 * T1   + d2 * T4 + d3 * T7 ;
		ch2[i] = d1 * T2   + d2 * T5 + d3 * T8;
		ch3[i] = d1 * T3   + d2 * T6 + d3 * T9;
    }
}	
/*
Find inverse matrix T^-1 
The T matrix is assumed to be square [size (n x n)]
The operation is done in place.
Input: T, n
Output: T
*/

void invT (
           double *T,
		   int n
		   )
{
    //read input
	gsl_matrix *M = gsl_matrix_alloc(n,n);
	for (int i = 0; i <n; i++)
	{
		for (int j = 0; j <n; j++)
		{
			gsl_matrix_set (M, i, j, T[i*n+j]);
		}
	}
	gsl_permutation * p = gsl_permutation_alloc (n);
	
	// invert
	int s;
	gsl_matrix *iM = gsl_matrix_alloc(n,n);
	gsl_linalg_LU_decomp (M, p, &s);    
	gsl_linalg_LU_invert (M, p, iM);
	
	//set output
	for (int i = 0; i <n; i++)
	{
		for (int j = 0; j <n; j++)
		{
			T[i*3+j] = gsl_matrix_get(iM,i,j);
		}
	}

	gsl_permutation_free (p);
    gsl_matrix_free (M);
	gsl_matrix_free (iM);
}		   