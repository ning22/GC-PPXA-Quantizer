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
	 
	Functions computing cost between 2 nodes for anisotropic TV:
	- mean absolute value criterion
	- quadratic norm
	- Potts cost
	
	This software was developed by Anna Jezierskaversion of 12.2011
	(anna.jezierska@univ-paris-est.fr).
	
	If you use this software for research purposes, you should cite
	the aforementioned paper in any resulting publication.
	
*/

#ifndef REGULARIZATION_H
#define REGULARIZATION_H

#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#define ABS(x) ((x)>=0?(x):-(x))
#define SQ(x)	((x) * (x))
#define DELTA(x,y)	((x)==(y)?(0):1)

/*
	Mean absolute value criterion
	
	PARAMETERS
	Input: x,y
	Output: res
*/

void psi_L1 (double x, double y, double* res);

/*
	Quadratic term
	
	PARAMETERS
	Input: x,y
	Output: res
*/

void psi_SQ (double x, double y, double* res);

/*
	Potts cost
	
	PARAMETERS
	Input: x,y
	Output: res
*/

void psi_Potts (double x, double y, double* res);


#endif //REGULARIZATION_H
