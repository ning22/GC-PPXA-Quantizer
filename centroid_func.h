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
	 
	Functions computing centroid of mass for:
	- quadratic norm
	- mean absolute value criterion
	
	This software was developed by Anna Jezierska version of 12.2011
	(anna.jezierska@univ-paris-est.fr).
	
	If you use this software for research purposes, you should cite
	the aforementioned paper in any resulting publication.
	
*/


#ifndef CENTROID_FUNC
#define CENTROID_FUNC


#include <string.h>

/*
   Function computs centroid of regions for quadratic norm
   
   PARAMETERS
   Input: data,q,id,sData,nCodevectors,nChannels
   Output: q
   
   data - pointer to the data
   q - pointer to codebook
   id - pointer to label image
   sData - data's size
   nChannels - number of image channels
   nCodevectors - number of codevectors 
   
*/

void centroid_quadratic(
					double **data,
					double **q,
					int *id,
					int sData,
                    int nCodevectors,
					int nChannels
					);


/*
   Function computs centroid of regions for mean absolute value criterion
   
   PARAMETERS
   Input: data,q,id,sData,nCodevectors,nChannels
   Output: q
   
   data - pointer to the data
   q - pointer to codebook
   id - pointer to label image
   sData - data's size
   nChannels - number of image channels
   nCodevectors - codebook's number of codevectors 
   
*/
					
void centroid_absolute(
					double **data,
					double **q,
					int *id,
					int sData,
                    int nCodevectors,
					int nChannels
					);



#endif