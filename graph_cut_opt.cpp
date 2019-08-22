#include "graph_cut_opt.h"

/*
ALPHA EXPANSION
*/


void scramble_label_table(int * label_table, int label_table_size)
{
	for (int i =0; i < label_table_size; i++)
	{
		int r1 = rand()%label_table_size;
        int r2 = rand()%label_table_size;
		int temp = label_table[r1];
		label_table[r1] = label_table[r2];
        label_table[r2] = temp;
	}
}



/*
Performs alpha-expansion for  regular grid graph  (binary case, label specified by alpha)
*/

void a_exp_binary (
                        double *qalpha,
						int  ialpha,
						double **data, 
						double **q, 
						int *id, 
						int sx, 
						int sy, 
						int nChannels,
						void (*fid_func)(double*, double*, int, double*), 
						void (*reg_func)(double, double,double*), 
						double weight
					)
{
	 int sData = sx*sy; 
	 
	 typedef Graph<double,double,double> GraphAlpha;
     GraphAlpha *g  = new GraphAlpha(sData,2*sData); 
	 g -> add_node(sData); 
	//add t-link for pixel
	for (int i=0; i<sData; ++i)
	{
		int eid = id[i];
		double *d = new double[nChannels];
		double *a = new double[nChannels];
        for (int ch = 0; ch < nChannels; ch ++)
		{
			d[ch] = data[ch][i];
			a[ch] = q[ch][eid];
		}		
		double weight_alpha=0;
		double weight_alpha_neg =0;
		if (eid!=ialpha)
		{
			double E0 = 0;
			double E1 = 0;
			fid_func (d,a,nChannels,&E0);
			fid_func (d,qalpha,nChannels,&E1);
			weight_alpha_neg = E0 - E1;
		}
		if (weight_alpha_neg<0)
		{
			weight_alpha = -weight_alpha_neg;
			weight_alpha_neg = 0;
		}
		g -> add_tweights(i, weight_alpha, weight_alpha_neg);
		delete[] d;
		delete[] a;
	}
	/*add n-link edges -> vertical*/

	for (int j=0; j<sy; ++j)
	{
		int node2=j;
		for (int i=0; i<(sx-1); ++i)
		{
			int node1 = node2;
			node2 = node1+sx;
			double A = 0, B = 0, C = 0;
			reg_func (id[node1], id[node2], &A);
			reg_func (id[node1], ialpha, &B);
			reg_func (ialpha, id[node2], &C);
			double w = weight*(B+C-A);
			if (w !=0.0) g ->add_edge(node1,node2,w,0);
			g -> add_tweights(node2,0, weight*C);
			w = weight*(C-A);
			if (w!=0.0)
			{
				if (weight>0)
					g -> add_tweights(node1,w,0);
				else
					g -> add_tweights(node1,0,-w);
			}

		}
	}
	
	 /*add n-link  -> horizontal*/

	for (int j=0; j<sx; ++j)
	{
		int node2=j*sx;
		for (int i=0; i<(sy-1); ++i)
		{
			int node1 = node2;
			node2 = node1+1;
			double A = 0, B = 0, C = 0;
			reg_func (id[node1], id[node2], &A);
			reg_func (id[node1], ialpha, &B);
			reg_func (ialpha, id[node2], &C);
			double w = weight*(B+C-A);
			if (w != 0) g ->add_edge(node1,node2,w,0);
			g -> add_tweights(node2,0, weight*C);
			w = weight*(C-A);
			if (w!=0.0)
			{
				if (weight>0.0)
					g -> add_tweights(node1,w,0);
				else
					g -> add_tweights(node1,0,-w);
			}
		}
	}
    
	g->maxflow();
	
	for (int i=0;i <sData;i++)
	{
		 if (g->what_segment(i) == GraphAlpha::SINK)
		 {
			 id[i] = ialpha;
		 }
	}
	  
	g->~Graph();
}



/*
Performs alpha-expansion for  regular grid graph (multilabel case)
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
						)
{
   int* label_table = new int [nCodevectors]; // init label table
   for (int i =0; i < nCodevectors; i++) label_table[i] = i;
   scramble_label_table(label_table,nCodevectors); // random label order
   double old_energy = e_in;
   int change = 1;
   while (change)
   {
		for (int i =0; i <nCodevectors; i++)
		{
			int index = label_table[i];
			double *alpha = new double[nChannels];
			for (int ch = 0; ch < nChannels; ch ++)
			{
				alpha[ch] = q[ch][index];
			}
			a_exp_binary (alpha, index, data, q, id, sx, sy, nChannels, fid_func, reg_func, weight);
			delete[] alpha;
		}
		double e = energy ( data, q, id, sx, sy, nChannels, fid_func, reg_func, weight);
		
		if (e == old_energy) change = 0;
		else old_energy = e;
		
	}
	delete [] label_table;
	return old_energy;
}		

/*
MUROTA'S GRADIENT DESCENT ALGORITHM
*/	

/*
Performs Murota's gradient descent for  regular grid graph  (binary case, for one given step)
*/

void grad_descent_binary (
						int  istep,
						double **data, 
						double **q, 
						int *id, 
						int sx, 
						int sy, 
						int nChannels,
						void (*fid_func)(double*, double*, int, double*), 
						void (*reg_func)(double, double,double*), 
						double weight,
						int nCodevectors
					)			
{
	int change = 0;
	
    int sData = sx*sy; 

	double infval = nCodevectors*100;
	
	typedef Graph<double,double,double> GraphGD;
    GraphGD *g  = new GraphGD(sData,2*sData); 
	g -> add_node(sData); 
	//add t-link for pixel
	for (int i=0; i<sData; ++i)
    {
		int eid = id[i];
		
		double *d = new double[nChannels];
		double *a0 = new double[nChannels];
        for (int ch = 0; ch < nChannels; ch ++)
		{
			d[ch] = data[ch][i];
			a0[ch] = q[ch][eid];
		}
		
		int id_proposed = eid+istep;
		if (id_proposed < 0 ) id_proposed = 0;
		if (id_proposed > nCodevectors-1 ) id_proposed = nCodevectors-1;
		double *a1 = new double[nChannels];
		
        if (id_proposed < 0 || id_proposed > nCodevectors-1) 
		{ 
			for (int ch = 0; ch < nChannels; ch ++)
			{
				a1[ch] = infval;
			}
		} 
        else
        {	
			for (int ch = 0; ch < nChannels; ch ++)
			{
				a1[ch] = q[ch][id_proposed];
			}
		}
		
		double E0 = 0;
		double E1 = 0;
		fid_func (d,a0,nChannels,&E0);
		fid_func (d,a1,nChannels,&E1);
		
		double weight_alpha_neg = E0 - E1;
		double weight_alpha=0;
		
		if (weight_alpha_neg<0)
		{
			weight_alpha = -weight_alpha_neg;
			weight_alpha_neg = 0;
		}
		g -> add_tweights(i, weight_alpha, weight_alpha_neg);
		delete[] d;
		delete[] a0;
		delete[] a1;
		
    }
	
	//add n-link edges -> vertical
	for (int j=0; j<sy; ++j)
	{
		int node2=j;
		for (int i=0; i<(sx-1); ++i)
		{
			int node1 = node2;
			node2 = node1+sx;
			double A = 0, B = 0, C = 0;
			int id1 = id[node1];
			int id2 = id[node2];
			int id1_proposed = id1 + istep;
			int id2_proposed = id2 + istep;
			
			reg_func (id1, id2, &A);
			reg_func (id1, id2_proposed, &B);
			reg_func (id1_proposed, id2, &C);

			double w = weight*(B+C-A-A);
			if (w !=0.0) g ->add_edge(node1,node2,w,0);
			w = weight*(C-A);
			if (w!=0.0)
			{
				if (weight>0)
				{
					g -> add_tweights(node1,w,0);
					g -> add_tweights(node2,0, w);
				}
				else
				{
					g -> add_tweights(node1,0,-w);
					g -> add_tweights(node2,-w, 0);
				}
			}

		}
	}
	
	
	//add n-link  -> horizontal

	for (int j=0; j<sx; ++j)
	{
		int node2=j*sx;
		for (int i=0; i<(sy-1); ++i)
		{
			int node1 = node2;
			node2 = node1+1;
			double A = 0, B = 0, C = 0;
			int id1 = id[node1];
			int id2 = id[node2];
			int id1_proposed = id1 + istep;
			int id2_proposed = id2 + istep;
			
			reg_func (id1, id2, &A);
			reg_func (id1, id2_proposed, &B);
			reg_func (id1_proposed, id2, &C);
			
			double w = weight*(B+C-A-A);
			if (w != 0.0) g ->add_edge(node1,node2,w,0);
			w = weight*(C-A);
			if (weight>0.0)
			{
				g -> add_tweights(node1,w,0);
				g -> add_tweights(node2,0, w);
			}
			else
			{
				g -> add_tweights(node1,0,-w);
				g -> add_tweights(node2,-w, 0);
			}
		}
	}

	g->maxflow();
	
	for (int i=0;i <sData;i++)
	{
		 if (g->what_segment(i) == GraphGD::SINK)
		 {
			 id[i] = id[i] + istep ;
			 if (id[i] < 0) id[i] = 0;
			 else if (id[i] > nCodevectors-1) id[i] = nCodevectors-1;
		 }
	}
	g->~Graph();
}

/*
Performs Murota's gradient descent for regular grid graph  (full loop up to the convergence)
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
						)
{
	//init step
	int steps[] = {-1,1};
	int change = 1;
	int iter = 0;
	double old_energy = energy ( data, q, id, sx, sy, nChannels, fid_func, reg_func, weight);
	double e = old_energy;
	while(change)
	{		
		for (int i =0; i <2; i++)
		{
			iter++;
			int istep = steps[iter%2];
			grad_descent_binary (istep,data, q, id, sx, sy, nChannels, fid_func, reg_func, weight, nCodevectors);
			double etemp = energy ( data, q, id, sx, sy, nChannels, fid_func, reg_func, weight);	
			if (e>etemp) e= etemp;
		}		
		if (e == old_energy) change = 0;
		else old_energy = e;
	}
	//old_energy = energy ( data, q, id, sx, sy, fid_func, reg_func, weight);
	return old_energy;
	
}						

