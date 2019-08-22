//Copyright 2012 Anna Jezierska
#include "median_cut.h"

using namespace std;

void mcGetCodevectors (
                        int sData,
                        int nCodevectors,
						double **data,
						double **q
                        )
{
    int i = 0; // counter
	double *ch1 = data[0];
    double *ch2 = data[1];
    double *ch3 = data[2];
    double *Q_val_1 = q[0];
    double *Q_val_2 = q[1];
    double *Q_val_3 = q[2];
	
	
	Point* points = new Point[sData];
    for ( i =0; i<sData; i++)
    {
        points[i].x[0] = ch1[i];
        points[i].x[1] = ch2[i];
        points[i].x[2] = ch3[i];

    }
    list<Point> palette = medianCut(points, sData, nCodevectors);
    list<Point>::iterator iter;
    i=0;
    for (iter = palette.begin() ; iter != palette.end(); iter++)
    {
        Q_val_1[i] = (double)iter->x[0];
        Q_val_2[i] = (double)iter->x[1];
        Q_val_3[i] = (double)iter->x[2];
        i++;
    }
	delete points;

}
/*
void init_q_im_media_cut (
                            int ps,
                            int Q,
                            double *R,
                            double *G,
                            double *B,
                            double *Q_val_1,
                            double *Q_val_2,
                            double *Q_val_3,
                            double *q_1,
                            double *q_2,
                            double *q_3
                            )
{
    int32_t i =0, j=0, k=0;
    int32_t  distance = 0, d=0;
    for (i=0; i<ps; i++)
    {
        distance = 3*pow(256,2);
        for (k=0; k<Q; k++)
        {
            d = pow( ( R[i]- Q_val_1[k]),2)+
                pow( ( G[i]- Q_val_2[k]),2)+
                pow( ( B[i]- Q_val_3[k]),2);
            if (d<distance)
            {
                distance = d;
                j=k;
            }

        }
        q_1[i]=Q_val_1[j];
        q_2[i]=Q_val_2[j];
        q_3[i]=Q_val_3[j];
    }
}
*/
