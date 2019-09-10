#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

/*
* File Description:
* -----------------
* This file contains an optionally applied function for assigning non-anchor
* labels "N" to feature pairs with retention time values falling within some 
* window of features previously labeled "I" (Identities).
*
* Function is called from the identityAnchorSelection function in R.
*
*/


/*
 * Function: updateLabels
 * -------------------------
 *
 * Helper function for selectAnchorsByID. All features within windX units
 * in rtx or windY units in rty of the matched identity are eliminated as 
 * potential anchors by being labeled "N".
 * 
 * 
 *
*/
void updateLabels(SEXP labels, int f, double windrX, double windrY, 
				   double* rtx, double* rty, double* rts)
{
	const char *label = CHAR(STRING_ELT(labels, f)); 
			
	if(strcmp(label,"P") == 0)
	{
		if(	fabs(rtx[f] - rts[0]) < windrX || 
			fabs(rty[f] - rts[1]) < windrY)
		{
			SET_STRING_ELT(labels, f, mkChar("N"));
		}
	}
			
	return;		
}


/*
 * Function: selectAnchorsByID
 * -----------------------------
 * 
 * For each previously assigned identity pair, removes feature pairs labeled "P" 
 * (potential anchors) whose rtx or rty is within windX or windY units. 
 * 
 * Returns: updated anchor labels.
 * 
 * PARAMETERS:
 * 
 * labels: vector of characters labeled "P" (for potential anchors), "I" 
 * for Identities, "A" for anchors or "N" for non-anchors
 *
 * Ids: vector of integers indicating position of identities (labeled "I")
 *
 * rtx: vector of numeric retention times corresponding to x-dataset features
 *
 * rty: vector of numeric retention times corresponding to y-dataset features 
 *
 * windX: numeric retention time window value. Features within windX units in rtx 
 * are labeled "N" (for non-anchor).
 *
 * windY: numeric retention time window value. Features within windY units of anchor
 * in rty are labeled "N" (for non-anchor).
 *
*/ 
SEXP selectAnchorsByID(SEXP labels, SEXP Ids, SEXP rtx, SEXP rty, 
					   SEXP windX, SEXP windY)
{

	//duplicate copy of labels
	SEXP labels_copy = PROTECT(duplicate(labels));

	//obtaining C types
	int* rIds = INTEGER(Ids);
	double* RTx = REAL(rtx);
	double* RTy = REAL(rty);
	double windrX = REAL(windX)[0];
	double windrY = REAL(windY)[0];
	
	int length = LENGTH(labels);
	int numIds = LENGTH(Ids);
	double rts[2];
	
		
	for(int Id = 0; Id < numIds; Id++){
		int index = rIds[Id];
	    rts[0] = RTx[index];
		rts[1] = RTy[index];
		
		for(int feat = 0; feat < length; feat++){
			updateLabels(labels_copy, feat, windrX, windrY, RTx, RTy, rts);
		}
	}
	
	UNPROTECT(1);
	
	return labels_copy;
}