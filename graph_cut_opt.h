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

    Functions updating label image using:
	- alpha expansion move (for non-convex regularization term)
    - gradient descent Murota's algorithm (for convex regularization term) 
	
	This software was developed by Anna Jezierska version of 12.2011
	(anna.jezierska@univ-paris-est.fr).
	
	If you use this software for research purposes, you should cite
	the aforementioned paper in any resulting publication.
	
*/

#ifndef ALPHA_EXPANSION
#define ALPHA_EXPANSION

#include <iostream>
#include <stdlib.h>
#include "energy.h"

#include "graph.h"

//#include "BKgraph.h"
//#include "BKgraph.hpp"
//#include "BKmaxflow.hpp"

/*
	Alpha-expansion algorithm for multilabel problems.
	
	PARAMETERS
	Input: data,q, id, sx,sy, nChannels, nCodevectors, fid_func, reg_func, weight, e_in
	Output: id
	
	data - pointer to the data
    q - pointer to codebook
    id - pointer to label image
	sx,sy,nChannels - width, height and channel's number of image, respectively
	nCodevectors - codevector's number
	fid_func, reg_func, weight - energy parameters
	e_in - initial energy
*/

double a_exp_multilabel( 
						double **data, 
						double **q, 
						int *id, 
						int sx, 
						int sy, 
						int nChannels,
						int nCodevectors,
						void (*fid_func)(double*, double*, int, double*), 
						void (*reg_func)(double, double,double*), 
						double weight,
						double e_in
						);

/*
	Gradient-descent Murota's algorithm for multilabel problems.
	
	PARAMETERS
	Input: data,q, id, sx,sy, nChannels, nCodevectors, fid_func, reg_func, weight, e_in
	Output: id
	
	data - pointer to the data
    q - pointer to codebook
    id - pointer to label image
	sx,sy,nChannels - width, height and channel's number of image, respectively
	nCodevectors - codevector's number
	fid_func, reg_func, weight - energy parameters
	e_in - initial energy
	
*/						
						
double grad_descent_multilabel( 
							double **data, 
							double **q, 
							int *id, 
							int sx, 
							int sy, 
							int nChannels,
							int nCodevectors,
							void (*fid_func)(double*, double*, int, double*), 
							void (*reg_func)(double, double,double*), 
							double weight,
							double e_in
						);						
						

#endif