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
	 
	Functions updating codebook using PPXA+ algorithm, described in details in:
	
	[3] "Proximal splitting methods in signal processing"
	P. L. Combettes and J.-C. Pesquet.  
	In H. H. Bauschke, R. Burachik, P. L. Combettes, V. Elser, D. R. Luke, and H. Wolkowicz,
	editors, Fixed-Point Algorithms for Inverse Problems in Science and Engineering. Springer-Verlag, New York, 2010.
	
	This software was developed by Anna Jezierska version of 12.2011
	(anna.jezierska@univ-paris-est.fr).
	
*/




#ifndef PPXA
#define PPXA

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h> 
#include <gsl/gsl_blas.h>  

#define ABS(x) ((x)>=0?(x):-(x))
#define MAX(a,b) (b>a)?a:b

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
            double *q,
            double * data, 
            int *id,
            int nCodevectors,
            int sData,
            double de,
            int typ
            );


#endif 