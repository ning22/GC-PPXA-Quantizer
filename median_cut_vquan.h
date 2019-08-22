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
	 
	Median cut - standard color quantization algorithm, described e.g. in:

	[1] P. Heckbert. Color image quantization for frame buffer display.
     SIGGRAPH Comput. Graph., 16(3):297{307, Jul. 1982. 
	
	This software was developed by Anna Jezierska version of 12.2011
	(anna.jezierska@univ-paris-est.fr).
		
*/


#ifndef MEDIAN_CUT_VQUAN
#define MEDIAN_CUT_VQUAN

#include <iostream>
#include <math.h>
#include <list>

/*
    Function implements median cut quantization algorithm
	
	PARAMETERS
	Input: sData, nCodevectors, data, q
	Output: q
	
	data - pointer to the data
    q - pointer to codebook
	sData - data's size
	nCodevectors - number of codevectors 

*/

void mcGetCodevectors (
                        int sData,
                        int nCodevectors,
						double **data,
						double **q
                        );

#endif 