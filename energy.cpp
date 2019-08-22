// Copyright 2012 Anna Jezierska
#include "energy.h"

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

double energy ( double **data, double ** q, int *id, int sx, int sy, int nChannels, void (*fid_func)(double*, double*, int, double*),
			    void (*reg_func)(double, double, double*), double weight)
{
   
	int sData = sx*sy;

     double energy = 0;
   
    /*prior term*/
   
     for (int j=0; j<sy; ++j)
     {
        int k=j;
        for (int i=0; i<(sx-1); ++i)
        {
            double x = (int) id[k];
            double y = (int) id[k+sx];
			double res = 0; 
			reg_func (x,y,&res);
            energy += res;
            k+=sx;
        }
     }
	 
	 
     for (int j=0; j<sx; ++j)
     {
        int k=j*sx;
        for (int i=0; i<(sy-1); ++i)
        {
			double x = (int) id[k];
            double y = (int) id[k+1];
			double res = 0; 
			reg_func (x,y,&res);
            energy += res;
            ++k;
        }
    }
	energy = energy *weight;
	 
	/*data fidelity term*/
    for (int i = 0; i<sData; ++i)
    {
        int cv = id[i];
        double *x = new double[nChannels];
		double *y = new double[nChannels];
		for (int ch = 0; ch < nChannels; ch ++)
		{
			x[ch] = data[ch][i];
			y[ch] = q[ch][cv];
		}
		double res = 0; 
        fid_func( x, y, nChannels,&res);
        energy += res;
		delete[] x;
		delete[] y;
    }
   
    return energy;
}



