#include "vquan.h"

/*
    Copyright 2012 Anna Jezierska
	Function returns codevectors sorted in
	incresing oreder with respect to channel 1.
	PARAMETERS:
	Input: q, Q, nChannels
	Output: q
	
	q - pointer to codebook
	Q - number of codevectors
	nChannels - number of data's channels

*/

void orderCodevectors (double **q, int Q, int nChannels)
{
    int ind = 0;
    double min = DBL_MAX;
	double *y =  new double[Q];
	
	double **qtemp = new double*[nChannels];
	for (int i =0; i < nChannels; i++)
	{
		qtemp[i] = new double [Q];
	}
	
    for (int i=0; i<Q; i++)
    {
		y[i]= q[0][i];		
    }
    for (int j=0; j<Q; j++)
    {
        for (int i=0; i< Q; i++)
        {
            if (y[i] < min)
            {
                min = y[i];
                ind = i;
            }
        }
        y[ind] = DBL_MAX;
        min = DBL_MAX;
		for (int i =0; i < nChannels; i++)
		{
			qtemp[i][j] = q[i][ind];
		}
		
    }
    for (int i=0; i<Q; i++)
    {
		for (int j =0; j <nChannels; j++ )
		{
			q[j][i] = qtemp[j][i];
		}
    }
    delete[] y;
	delete [] qtemp;
}
/*
	Function returns one when the elements of array q are not sorted in increasing order.
	
	PARAMETERS:
	Input: q, nCodevectors
	Output: 0 - decreasing order, 1 - increasing order	
	
	q - pointer to codebook
	nCodevectors - number of codevectors
*/
int isIncOrder (double *q, int nCodevectors)
{ 
    
    int result = 0;
    for (int i=1; i<nCodevectors; i++){ if (q[i]<q[i-1]) {result = 1; break;}}
    return result;
	
}

/*
	Function initialize variable id, based on codevectors q and the input data
	
	PARAMETERS:
	Input: q, data, sData, nChannels, nCodevectors
	Output: id	
*/

void  initID ( int* id, double **q, double **data, int sData, int nChannels, int nCodevectors)
{
	int res = 0;
	for (int i = 0; i < sData; i++)
	{
		double distance = DBL_MAX;
		for (int j = 0; j < nCodevectors; j++)
		{
			double d = 0;
			for (int ch = 0; ch <nChannels; ch ++ )
			{
				d+= SQ (data[ch][i] - q[ch][j]);
			}
			if (d<distance)
			{
				distance = d;
				res=j;
			}
		}
		id[i] = res;
	}
}

/*
	PARAMETERS:
	Input: sData, nCodevectors, q, id
	Output: data
*/

void setResult ( int sData, int nChannels, double **data, double ** q, int *id)
{
	for (int i =0; i < sData; i++)
	{
		int cv = id[i];
		for (int ch = 0 ; ch <  nChannels; ch++ )
		{
			data[ch][i] = q[ch][cv];
		}
	}
}
/*
	Codevector's initialization
	1 -channel images: uniform quantization
	3 -channel images: median cut quantization result
	
	PARAMETERS:
	Input: sData,  nChannels, Q, data, q
	Output: q
*/
void getCodevectors (int sData,int nChannels, int Q, double **data, double **q)
{
	if (nChannels == 3)
	{
		mcGetCodevectors (sData,Q,data,q); 
	}
	else if (nChannels == 1)
	{
	    double * p = data[0];
		double mind = 1000;
		double maxd = -1000;
	    // find min(sData) and max(data[0])
		for (int i = 0; i < sData; i++)
		{
		  if (*(p)> maxd) maxd = *(p);
		  if (*(p)< mind) mind = *(p);
		  p++;
		}
		double range = maxd - mind;
		double step = range/Q;
		double f = mind+step/2;
		for (int i = 0; i < Q; i++)
		{
			q[0][i] = f+i*step; 
		}
	}
}

/*
	Two step alternating minimization algorithm for image quantization (spatial regularization approach).
	
	PARAMETERS
	Input: data, sx, sy, nChannels, Q, fid_func, reg_func, quantizer, centroid, weight, niter, typ, prior
	Output: data
*/

void spataialRegAlgForVectorQuant (
									double ** data,
									int sx,
									int sy,
									int nChannels,
									int Q,
									void (*fid_func)(double*, double*,int, double*),
									void (*reg_func)(double, double,  double*),
									double (*quantizer) (double **, double **, int *, int , int , int, int,
													void (*fid_func)(double*, double*, int, double*), 
													void (*reg_func)(double, double,double*), 
													double , double ),
									void (*centroid)(double **,double **,int *,int ,int ,int ),				
									double weight,
									int niter,
									int typ,
									int prior
									)
{
    double **q = new double*[nChannels];
	for (int i =0; i < nChannels; i++)
	{
		q[i] = new double [Q];
	}

	int sData = sx*sy;
	int *id = new int[sData];
	
    //INITIALIZATION
 
	getCodevectors (sData,nChannels,Q,data,q); 
    orderCodevectors (q,Q,nChannels);  
    initID (id,q,data,sx*sy,nChannels,Q); 
	
	// MAIN LOOP
	
	for (int i =0; i < niter; i++)
	{
		
	    // Calculate energy
		double e = energy ( data, q, id, sx, sy, nChannels, fid_func, reg_func, weight);
		// Step 1: update id
		double e1  = quantizer (data, q, id, sx, sy,nChannels,Q,fid_func,reg_func, weight, e); 
		std::cout <<"Iteration "<< i << " energy "  <<   e1  <<   " improvement "  <<  e - e1 << std::endl;
		// Step 2: updating of the codebook q
        centroid(data,q,id,sData,Q,nChannels);
		int ord = isIncOrder(q[0], Q);
		if (ord != 0 && prior !=3)
		{
			codebook_ppxa_update(q[0],data[0],id,Q,sData,0,typ);                     
		}
		if (e1 == e)
		{
			break;
		}
		else if (i == niter-1)
		{
			std:: cout << "Max number of iterations reached" << std::endl;
		}
	}
	
	setResult ( sData, nChannels, data, q, id);
	
	// CLEAN MEMORY
	delete [] id;
	delete [] q;
}