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
	 
	Functions computing data fidelity term for:
	- quadratic norm
	- mean absolute value criterion
	
	This software was developed by Anna Jezierska version of 12.2011
	(anna.jezierska@univ-paris-est.fr).
	
	If you use this software for research purposes, you should cite
	the aforementioned paper in any resulting publication.
	
*/

#ifndef DATA_FIDELITY_TERM_H
#define DATA_FIDELITY_TERM_H

#include <math.h>

#define ABS(x) ((x)>=0?(x):-(x))
#define SQ(x)	((x) * (x))

/*
	Quadratic term
	
	PARAMETERS
	Input: x,y, nChannels
	Output: res
*/

void D_SQ (double* x, double* y, int nChannels, double* res);

/*
	Mean absolute value criterion
	
	PARAMETERS
	Input: x,y, nChannels
	Output: res
*/

void D_L1 (double*, double*, int, double*);

#endif 
