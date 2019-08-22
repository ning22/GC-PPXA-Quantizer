/*
    Copyright 2012 Anna Jezierska
    This software implements the vector quantization algorithm
    described in:

    [1] "A Spatial Regularization Approach for Vector Quantization",
	C. Chaux, A. Jezierska, J.-C. Pesquet, and H. Talbot,
	Journal of Mathematical Imaging and Vision, vol. 41, pp. 23-38, 2011
	
	[2] "Image quantization under spatial smoothness constraints" ,
	A. Jezierska,C. Chaux,  J.-C. Pesquet, and H. Talbot, 
	International Conference on Image Processing (ICIP), Honk Kong, 26-29 September 2010. 
	
    This software was developed by Anna Jezierska version of 12.2011
	(anna.jezierska@univ-paris-est.fr).
	
	If you use this software for research purposes, you should cite
	the aforementioned paper in any resulting publication.	
*/

#ifndef ENERGY
#define ENERGY

#include <stdint.h>
#include <sys/types.h>

#include "data_fidelity_term.h"
#include "regularization_term.h"

/*
	Function returns the energy value E:
	Minimizaed energy E = fid_func + weight * reg_func
	
	PARAMETERS
	Input: data, q, id, sx, sy, nChannels, fid_func, reg_func, weight
	Output: energy
	
	data - pointer to the data
	q - pointer to codebook
	sx,sy,nChannels - width, height and channel's number of image, respectively
	fid_func, reg_func, weight - define energy function
	
*/

double energy ( 
				double **data, 
				double ** q, 
				int *id, 
				int sx, 
				int sy,
				int nChannels,
				void (*fid_func)(double*, double*, int, double*),
			    void (*reg_func)(double, double,double*),
				double weight
				);

#endif