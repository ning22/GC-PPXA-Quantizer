/*	
    Copyright 2012 Anna Jezierska
	This software implements the vector quantization algorithm
	described in
	
	[1] "A Spatial Regularization Approach for Vector Quantization",
	 C. Chaux, A. Jezierska, J.-C. Pesquet, and H. Talbot,
	 Journal of Mathematical Imaging and Vision, vol. 41, pp. 23-38, 2011
	
	[2] "Image quantization under spatial smoothness constraints" ,
	 A. Jezierska,C. Chaux,  J.-C. Pesquet, and H. Talbot, 
	 International Conference on Image Processing (ICIP), Honk Kong, 26-29 September 2010. 
	
	Functions computing data linear transform using Principal Component Analysis. 
		
	This software was developed by Anna Jezierska version of 12.2011
	(anna.jezierska@univ-paris-est.fr).
	
	If you use this software for research purposes, you should cite
	the aforementioned paper in any resulting publication.
	
*/

#ifndef PCA
#define PCA


#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h> 
#include <gsl/gsl_blas.h>  

/*
	Function finds a linear operator T using Principal Component Analysis 
	
	PARAMETERS
    Input: data, size
    Output: T
	
	data - input data (3 channels)
	size - data's size
	T - linear operator of size 3x3
*/

void getPCA   (
			  double **data,
              int size,
              double * T
              );
			  
/*
	Function converts data = {ch1, ch2, ch3} using linear operator T
	
	PARAMETERS
    Input: data, size, T
    Output: data
	
	data - input data (3 channels)
	size - data's size
	T - linear operator of size 3x3
*/			  

void D2T (
			  double ** data,
              int size,
              double * T
			);	

			
/*
	Function finds inverse matrix T^-1 of size (n x n)
	
	PARAMETERS
    Input: T,n
    Output: T
	
*/	

void invT (
           double *T,
		   int n
		   );		  

#endif //PCA