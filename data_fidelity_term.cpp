// Copyright 2012 Anna Jezierska

#include "data_fidelity_term.h"

/*
	Quadratic term
	
	PARAMETERS
	Input: x,y, nChannels
	Output: res
	
	x,y - data vectors
	nChannels - vector's size
	res - resulting scalar
	
*/

void D_SQ (double*x, double*y, int nChannels, double *res)
{
	double r = 0;
	for (int ch = 0; ch < nChannels; ch++)
	{
		r+= SQ(x[ch] - y[ch]);
	}
	*res = r/ nChannels;
}

/*
	Mean absolute value criterion
	
	PARAMETERS
	Input: x,y, nChannels
	Output: res
	
	x,y - data vectors
	nChannels - vector's size
	res - resulting scalar
*/

void D_L1 (double*x, double*y, int nChannels, double *res)
{
	double r = 0;
	for (int ch = 0; ch < nChannels; ch++)
	{
		r+= ABS(x[ch] - y[ch]);
	}
	*res = r/ nChannels;
}