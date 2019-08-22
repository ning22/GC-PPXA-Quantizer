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


#ifndef VQUAN
#define VQUAN

#include <iostream>
#include <string.h>
#include "median_cut_vquan.h"
#include "energy.h"
#include "ppxa.h"

#define DBL_MAX 1.7976931348623158e+308 /* max value */
#define SQ(x)	((x) * (x))

/*
	Two step alternating minimization algorithm for image quantization (spatial regularization approach).

	PARAMETERS

	Inputs:
	data - pointer to image's data
	sx,sy,nChannels - width, height and channel's number of image, respectively
	Q - number of quantization levels
	fid_func,reg_func, weight,typ - parameters corresponding to minimized energy E
	quantizer - algorithm calculating new partisions 
	centroid - function calculating center of mass (for given data fidelity)
	niter - max number of iterations

	Output:
	data
*/

void spataialRegAlgForVectorQuant (
									double ** data,
									int sx,
									int sy,
									int nChannels,
									int Q,
									void (*fid_func)(double*, double*, int, double*),
									void (*reg_func)(double, double,double*),
									double (*quantizer) (double **, double **, int *, int , int , int, int,
									                  void (*fid_func)(double*, double*, int, double*), 
													  void (*reg_func)(double, double,double*), 
					                                   double , double ),
									void (*centroid)(double **,double **,int *,int ,int ,int ),					   
									double weight,
									int niter,
									int typ,
									int prior
									);


#endif