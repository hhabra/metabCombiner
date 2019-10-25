#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <float.h>

/*
 * File Description:
 * -----------------
 * This file contains a function for iteratively assigning anchor labels ("A") 
 * to abundant paired features in the two datasets. The anchor selection process 
 * continues until no features are labeled "P" (or potential anchors).
*/

/*
 * Function: findMaxQ
 * -------------------
 * Helper function for selectIterativeAnchors(). Searches for the most abundant
 * feature pair still labeled "P", as determined by Qx. Index of the most abundant 
 * complimentary feature (as determined by Qy) is chosen as an anchor match.
 *
*/
int findMaxQ(SEXP labels, double* Qx, double* Qy, double* rtx, int length)
{
	const char *label;
	double maxQ = -DBL_MAX, iQ; 
	int index = -1;
	
	for (int f = 0; f < length; f++){
		label = CHAR(STRING_ELT(labels, f));
		
		if(strcmp(label,"P") == 0){
			iQ = Qx[f];
			
			//finding largest abundance, maxQ
			if(iQ > maxQ){
				maxQ = iQ;
				index = f;
			} 
			
			else if(iQ == maxQ && rtx[f] == rtx[index] && Qy[f] > Qy[index]){
				index = f;
			}	
		}
	}
	
	return index;
}

/*
 * Function: updateLabels
 * -------------------------
 *
 * Helper function for selectIterativeAnchors. All features within windX units
 * in rtx or windY units in rty of the matched identity are eliminated as 
 * potential anchors by being labeled "N".
 * 

void updateLabels(SEXP labels, int f, double windX, double windY,
				   double* rtx, double* rty, double* rts)
{
	const char *label = CHAR(STRING_ELT(labels, f)); 
			
	if(strcmp(label,"P") == 0){
		if(fabs(rtx[f] - rts[0]) < windX || fabs(rty[f] - rts[1]) < windY){
			SET_STRING_ELT(labels, f, mkChar("N"));
		}
	}
			
	return;		
}
*/

/*
 * Function: selectIterativeAnchors
 * ---------------------------------
 * Iterative algorithm for selecting paired features whose retention times will
 * serve as ordered pairs (anchors) for a non-linear retention time model.
 *
 * Returns: updated anchor labels ("A" for anchor, "I" for identity, "N" for non-anchor
 * 
 * PARAMETERS:
 *  
 * labels: vector of characters labeled "P" (for potential anchors), "I" 
 * for Identities, "A" for anchors or "N" for non-anchors
 *
 * rtx: vector of numeric retention times corresponding to x-dataset features
 *
 * rty: vector of numeric retention times corresponding to y-dataset features
 *
 * Qx: vector of numeric abundance Q values corresponding to x-dataset features
 *
 * Qy: vector of numeric abundance Q values corresponding to y-dataset features  
 *
 * windX: numeric RT window value. Features within windX units of anchor in rtx 
 * are labeled "N" (for non-anchor).
 *
 * windY: numeric RT window value. Features within windY units of anchor in rty
 * are labeled "N" (for non-anchor)
 *
*/

SEXP selectIterativeAnchors(SEXP labels, SEXP rtx, SEXP rty, SEXP Qx, SEXP Qy, 
							SEXP windX, SEXP windY)
{						
	//duplicate copy of labels
	SEXP labels_copy = PROTECT(duplicate(labels));
	
	//obtaining C types
	double* rt_x = REAL(rtx);
	double* rt_y = REAL(rty);
	double* Q_x = REAL(Qx);
	double* Q_y = REAL(Qy);
	
	double wind_X = REAL(windX)[0];
	double wind_Y = REAL(windY)[0];
	
	//loop variables
	int length = LENGTH(labels);
	int nextAnchor;
	double rt_a[2];
		
	while(1){
		nextAnchor = findMaxQ(labels_copy, Q_x, Q_y, rt_x, length);
		
		if(nextAnchor < 0)
			break;
		
		SET_STRING_ELT(labels_copy, nextAnchor, mkChar("A"));
	
		rt_a[0] = rt_x[nextAnchor];
		rt_a[1] = rt_y[nextAnchor];
		
		//update labels: exclude features within RT window of labeled anchors
		for(int feat = 0; feat < length; feat++){
			updateLabels(labels_copy, feat, wind_X, wind_Y, rt_x, rt_y, rt_a);			
		}		
	}
		
	UNPROTECT(1);
	
	return labels_copy;
}




