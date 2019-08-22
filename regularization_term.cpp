//Copyright 2012 Anna Jezierska
#include "regularization_term.h"

void psi_L1 (double x, double y , double *res)
{
	*res = ABS(x-y);
}

void psi_SQ (double x, double y, double *res)
{
	*res = SQ(x-y);
}

void psi_Potts (double x, double y, double *res)
{
	*res = DELTA(x,y);
}

