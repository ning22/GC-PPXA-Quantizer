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
	
	This software was developed by Anna Jezierska version of 12.2011
	(anna.jezierska@univ-paris-est.fr).
	
	If you use this software for research purposes, you should cite
	the aforementioned paper in any resulting publication.
	
*/

#include <iostream>
#include <ctime>
#include <sstream>

#include "pca.h"
#include "vquan.h"
#include "data_fidelity_term.h"
#include "regularization_term.h"
#include "graph_cut_opt.h"
#include "centroid_func.h"


using namespace std;



static const char *usage =                                    \
    "usage: %s <size x> <size y> <nChannels> <ch1> <ch2> <ch3> <Q> <data fidelity> <reg term> <weight> <niter>\n"
	"      size x, size y - width and high of image, positive integer \n" 
	"      nChannels - number of channels, either 1 or 3\n" 
	"      ch1, ch2, ch3 - image data , .txt format\n"
	"      Q - number of codevectors, positive integer \n" 
	"      Minimizaed energy E = data fidelity (x,y) + weight reg term(x) \n"
	"      data fidelity - 1 -> L1, 2-> (L2)^2 \n" 
	"      reg term (anisotropic TV) - 1 -> L1, 2-> (L2)^2 , 3 -> Potts\n"
	"      weight - positive double \n"
	"      niter - max number of iterations ,positive integer";
	
static const char *pos_err = 
    "      Image size, number of codevectors and number of iterations is suppposed to be possitive integer\n" 
	"      size x, size y - width and high of image, positive integer \n" 
	"      Q - number of codevectors, positive integer \n" 
	"      niter - max number of iterations ,positive integer \n" ;



static const char *chan_err = 
    "      This software vesrion supports only 1-channel images and 3-channel images\n" ;
	

/*
Function read double value from txt file "filename" into varible x
*/	
	
void read_signal(double *x, int n, char *filename)
{
    int i;
    FILE *fp=fopen(filename, "rt");
	if (fp==NULL)
    {
       cout<<"Error opening file"<< filename << endl;
	   exit(1);
    }
	else
	{
		for(i=0; i<n; ++i){
			float t;
			fscanf(fp,"%f\n", &t);
			x[i] = (double)t;
		}
		fclose(fp);
	}
}	

/*
Function save variable x under txt file "filename"
*/	
void save_signal(double* x,int n,char* filename) 
{
    int i;
    FILE *fp=fopen(filename,"wt");
    for (i=0;i<n;i++) {
        fprintf(fp,"%e\n",x[i]);
    }
    fclose(fp);
}

	
	
int main(int argc, char *argv[])
{
	if (argc < 11) 
	{
        fprintf(stderr, usage, argv[0]);
        exit(1);
    }
	int sx =  atoi(argv[1]);
	int sy =  atoi(argv[2]);
	int nChannels =  atoi(argv[3]);
	int size =  sx*sy;
	
	int Q = atoi(argv[7]);
	int ifid = atoi(argv[8]); // function choice
	int ireg = atoi(argv[9]); // function choice
	double lambda = atof(argv[10]);
	int niter = atoi(argv[11]);
	
	
	if ( nChannels!=1 && nChannels!=3 ) 
	{
		fprintf(stderr, chan_err, argv[0]);
        exit(1);
	}
	
	
	if (Q<0 || sx < 1 || sy < 1 ) 
	{
		fprintf(stderr, pos_err, argv[0]);
        exit(1);
	}
	
	void (*fidelity)(double *, double*, int, double*);
	void (*regularization)(double, double, double*);
	double (*quantizer) (double **, double **, int *, int , int , int, int,
					   void (*fid_func)(double*, double*, int, double*), void (*reg_func)(double, double,double*), double , double );
	void (*centroid) (double **,double **,int *,int ,int ,int );				   
   
	switch(ifid)
	{
		case 1:
		fidelity = &(D_L1);
		centroid = &(centroid_absolute);
		break;
		case 2:
		fidelity = &(D_SQ);
		centroid = &(centroid_quadratic);
		break;
		default:
		cout << "UNKNOWN DATA FIDELITY TERM" << endl;
		exit(1);
		break;
	}
	
	switch(ireg)
	{
		case 1:
		regularization = &(psi_L1);
		quantizer =  &(grad_descent_multilabel);
		break;
		case 2:
		regularization = &(psi_SQ);
		quantizer =  &(grad_descent_multilabel);
		break;
		case 3:
		regularization = &(psi_Potts);
		quantizer =  &(a_exp_multilabel);
		break;
		default:
		cout << "UNKNOWN REGULARIZATION TERM" << endl;
		exit(1);
		break;
	}
	
	double **data = new double*[nChannels];
	char** ch_name = new char *[nChannels];
	
	for (int i =0; i < nChannels; i++)
	{
	    ch_name[i] = argv[i+4];
		data[i] = new double [size];
		read_signal(data[i],size, ch_name[i]);
	}
	
	double * T = new double [9]; // linear operator (PCA)
	
	////////////////////////////////////////
	
	if (nChannels == 3)
	{
		getPCA ( data,size, T);
		D2T (data,size,T);
	}
	
	
	spataialRegAlgForVectorQuant (data,sx,sy, nChannels, Q, fidelity, regularization, quantizer, centroid, lambda, niter, ifid, ireg);
	
	if (nChannels == 3)
	{
		invT(T,nChannels);
		D2T (data,size,T);
	}

	
	///////////////////////////////////////
	
	cout << "SAVE DATA" << endl;
	for (int i =0; i <nChannels; i++)
	{
	    ostringstream st;
        st <<  "out_ch" << i+1 << ".txt" ;
        string outstring = st.str();
		save_signal(data[i],size,(char*)outstring.c_str()); 
	}
	delete [] data;
	delete [] ch_name;
    return 0;
}
