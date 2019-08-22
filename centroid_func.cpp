// Copyright 2012 Anna Jezierska

#include "centroid_func.h"
#include <iostream>

// QUADRATIC NORM

/*
	Function compute the central of mass for given codevector q and lebeling id (quadratic norm)
	
	PARAMETERS
	Input: data, q, id, sData, nCodevectors
	Output: q
*/

void mean (double *data, double *q, int *id, int sData, int nCodevectors)
{
	double *t = new double [nCodevectors];
    int *st = new int [nCodevectors];
	memset (t,0,nCodevectors*sizeof(double));
	memset (st,0,nCodevectors*sizeof(int));
	
	int * idp = id;
	double * datap = data;
	for (int i =0; i < sData; i++ )
	{
		int index = *(idp);
		t[index] += *(datap);
		st[index]++;
		idp++;
		datap++;
	}
	double *tp = t;
	int *stp = st;
	for (int i = 0; i <nCodevectors; i++ )
	{
		if (*(stp) != 0)
		{
			q[i] = *(tp)/(*(stp));
		}
		tp++;
		stp++;
	}
	stp = st+1;
	for (int i = 1; i <nCodevectors-1; i++ )
	{
		if (*(stp) == 0)
		{
			q[i] = (q[i-1]+ q[i+1])/2;
		}
		stp++;
	}
	
	delete[] t;
	delete[] st;
}

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
   nCodevectors - codebook's number of codevectors 
   
*/

void centroid_quadratic ( 
					double **data,
					double **q,
					int *id,
					int sData,
                    int nCodevectors,
					int nChannels
                    )
{
    for (int i=0; i<nChannels; ++i)
    {
        mean(data[i], q[i], id, sData, nCodevectors);
    }

}


/*
CENTROID ABSOLUTE VALUE CRITERION
*/

/*
	Algorithm from N. Wirth's book, implementation by N. Devillard.
	This code in public domain.
*/

typedef double elem_type ;
#define ELEM_SWAP(a,b) { register elem_type t=(a);(a)=(b);(b)=t; }

/*
	Function : kth_smallest()
	In : array of elements, # of elements in the array, rank k
	Out : one element
	Job : find the kth smallest element in the array
	Notice : use the median() macro defined below to get the median.
	Reference:
	Author: Wirth, Niklaus
	Title: Algorithms + data structures = programs
	Publisher: Englewood Cliffs: Prentice-Hall, 1976
	Physical description: 366 p.
	Series: Prentice-Hall Series in Automatic Computation
*/

elem_type kth_smallest(elem_type a[], int n, int k)
{
	register int i,j,l,m ;
	register elem_type x ;
	l=0 ; m=n-1 ;
	while (l<m) 
	{
		x=a[k] ;
		i=l ;
		j=m ;
		do 
		{
			while (a[i]<x) i++ ;
			while (x<a[j]) j-- ;
			if (i<=j) 
			{
				ELEM_SWAP(a[i],a[j]) ;
				i++ ; j-- ;
			}
		} while (i<=j) ;
		if (j<k) l=i ;
		if (k<i) m=j ;
	}
	return a[k] ;
}

#define median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))

/*
	Function compute the central of mass for given codevector q and lebeling id (mean absolute value criterion)
	
	PARAMETERS
	Input: data, q, id, sData, nCodevectors
	Output: q
*/

void arrmedian (double *data, double *q, int *id, int sData, int nCodevectors)
{
    // find the number of members of each region
	int *st = new int [nCodevectors];
	memset (st,0,nCodevectors*sizeof(int));
	int * idp = id;
	for (int i =0; i < sData; i++ )
	{
		int index = *(idp);
		st[index]++;
		idp++;
	}
	
	// create the arrays
	double *t [nCodevectors];  
	for (int i=0 ; i<nCodevectors; ++i)    t[i] = new double [st[i]];
	
	// find members of each region
	double * datap = data;
	idp = id;
	int *ct = new int [nCodevectors];
	memset (ct,0,nCodevectors*sizeof(int));
	for ( int i = 0; i< sData ; ++i)
    {
	    int index1 = *(idp);
		int index2 = ct[index1];
		t[index1] [index2] = *(datap);
		ct[index1] = index2+1;
		idp++;
		datap++;
    }
	delete[] ct;
	for (int i=0 ; i<nCodevectors; ++i)  
	{
		q[i] = median(t[i], st[i]);
		delete[] t[i];
	}
	int *stp = st+1;
	for (int i = 1; i <nCodevectors-1; i++ )
	{
		if (*(stp) == 0)
		{
			q[i] = (q[i-1]+ q[i+1])/2;
		}
		stp++;
	}
	delete[] st;
}

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

void centroid_absolute ( 
					double **data,
					double **q,
					int *id,
					int sData,
                    int nCodevectors,
					int nChannels
                    )
{
    for (int i=0; i<nChannels; ++i)
    {
        arrmedian(data[i], q[i], id, sData, nCodevectors);
    }
}

#undef ELEM_SWAP


