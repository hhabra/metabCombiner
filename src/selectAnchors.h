#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

void updateLabels(SEXP labels, int f, double windx, double windy, 
				   double* rtx, double* rty, double* rts);

int findMaxQ(SEXP labels, double* Qx, double* Qy, double* rtx, int length);

SEXP selectAnchorsByID(SEXP labels, SEXP ids, SEXP rtx, SEXP rty, 
					   SEXP windx, SEXP windy);
					   
SEXP selectIterativeAnchors(SEXP labels, SEXP rtx, SEXP rty, SEXP Qx, SEXP Qy, 
							SEXP windx, SEXP windy);


