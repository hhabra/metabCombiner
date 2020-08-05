#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

/*
 * File Description:
 * -----------------
 * This file contains two methods for the assignment of anchor labels to ordered
 * pair features whose retention times will be used to plot a nonlinear splines model.
 * The main method selectIterativeAnchors() uses an iterative approach based on features
 * with high abundances in both input datasets. An optional method, selectAnchorsByID(),
 * can be applied first for shared identified compounds. 
 *
*/



/*
 * Function: updateLabels
 * -------------------------
 *
 * Helper function. All features within windx units in rtx or windy units in rty of the
 * matched identity are eliminated as potential anchors by being labeled "N".
 *  
*/
void updateLabels(SEXP labels, int f, double windx, double windy, 
				   double* rtx, double* rty, double* rts)
{
	const char *label = CHAR(STRING_ELT(labels, f)); 
			
	if(strcmp(label,"P") == 0)
	{
		if(	fabs(rtx[f] - rts[0]) < windx || 
			fabs(rty[f] - rts[1]) < windy)
		{
			SET_STRING_ELT(labels, f, mkChar("N"));
		}
	}
			
	return;		
}

/*
 * Function: findMaxQ
 * -------------------
 * Helper function. Searches for the most abundant feature pair still labeled "P", as 
 * determined by Qx. Index of the most abundant  complimentary feature (as determined by 
 * Qy) is chosen as an anchor match.
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
			
			//determine most abundant counterpart of feature with maxQ
			else if(iQ == maxQ && rtx[f] == rtx[index] && Qy[f] > Qy[index]){
				index = f;
			}	
		}
	}
	
	return index;
}


/*
 * Function: selectAnchorsByID
 * -----------------------------
 * 
 * For each previously assigned identity pair, removes feature pairs labeled "P" 
 * (potential anchors) whose rtx or rty is within windx or windy units. 
 * 
 * Returns: updated labels ("I" for identity, "P" for potential, "N" for non-anchor)
 * 
 * PARAMETERS:
 * 
 * labels: vector of characters labeled "P" for potential anchors, "I" 
 * for Identities, "A" for anchors or "N" for non-anchors
 *
 * ids: vector of integers indicating row position of identities (labeled "I")
 *
 * rtx: vector of numeric retention times corresponding to x-dataset features
 *
 * rty: vector of numeric retention times corresponding to y-dataset features 
 *
 * windx: numeric retention time window value. Features within windx units in rtx 
 * are labeled "N" (for non-anchor).
 *
 * windy: numeric retention time window value. Features within windy units of anchor
 * in rty are labeled "N" (for non-anchor).
 *
*/ 
SEXP selectAnchorsByID(SEXP labels, SEXP ids, SEXP rtx, SEXP rty, 
					   SEXP windx, SEXP windy)
{

	//duplicate copy of labels
	SEXP labels_copy = PROTECT(duplicate(labels));

	//obtaining C types
	int* ids_c = INTEGER(ids);
	double* rtx_c = REAL(rtx);
	double* rty_c = REAL(rty);
	double windx_c = REAL(windx)[0];
	double windy_c = REAL(windy)[0];
	
	int length = LENGTH(labels);
	int numIds = LENGTH(ids);
	double rts[2];
	
	for(int id = 0; id < numIds; id++){
		int index = ids_c[id]-1;
	    rts[0] = rtx_c[index];
		rts[1] = rty_c[index];
		
		for(int feat = 0; feat < length; feat++){
			updateLabels(labels_copy, feat, windx_c, windy_c, rtx_c, rty_c, rts);
		}
	}
	
	UNPROTECT(1);
	
	return labels_copy;
}

/*
 * Function: selectIterativeAnchors
 * ---------------------------------
 * Iterative algorithm for selecting paired features whose retention times will
 * serve as ordered pairs (anchors) for a non-linear retention time model.
 *
 * Returns: updated anchor labels ("A" for anchor, "I" for identity, "N" for non-anchor)
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
 * windx: numeric RT window value. Features within windx units of anchor in rtx 
 * are labeled "N" (for non-anchor).
 *
 * windy: numeric RT window value. Features within windy units of anchor in rty
 * are labeled "N" (for non-anchor)
 *
*/
SEXP selectIterativeAnchors(SEXP labels, SEXP rtx, SEXP rty, SEXP Qx, SEXP Qy, 
							SEXP windx, SEXP windy)
{						
	//duplicate copy of labels
	SEXP labels_copy = PROTECT(duplicate(labels));
	
	//obtaining C types
	double* rtx_c = REAL(rtx);
	double* rty_c = REAL(rty);
	double* Qx_c = REAL(Qx);
	double* Qy_c = REAL(Qy);
	
	double windx_c = REAL(windx)[0];
	double windy_c = REAL(windy)[0];
	
	//loop variables
	int length = LENGTH(labels);
	int nextAnchor;
	double rts[2];
		
	while(1){
		nextAnchor = findMaxQ(labels_copy, Qx_c, Qy_c, rtx_c, length);
		
		if(nextAnchor < 0)
			break;
		
		SET_STRING_ELT(labels_copy, nextAnchor, mkChar("A"));
	
		rts[0] = rtx_c[nextAnchor];
		rts[1] = rty_c[nextAnchor];
		
		//update labels: exclude features within RT window of labeled anchors
		for(int feat = 0; feat < length; feat++){
			updateLabels(labels_copy, feat, windx_c, windy_c, rtx_c, rty_c, rts);			
		}		
	}
		
	UNPROTECT(1);
	
	return labels_copy;
}





