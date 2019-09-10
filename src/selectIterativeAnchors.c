#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

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
* feature still labeled "P", as determined by rQx. Index of the most abundant 
* complimentary feature (as determined by rQy) is chosen as an anchor match; 
* others are deemed mismatches and labeled "N".
*/
int findMaxQ(SEXP labels, double* rQx, double* rQy, double* RTx, int length)
{
	const char *label;
	double maxQ = 0, iQ = 0; 
	int index = 0;
	
	for (int f = 0; f < length; f++){
		label = CHAR(STRING_ELT(labels, f));
		
		if(strcmp(label,"P") == 0){
			iQ = rQx[f];
			
			//finding largest abundance, maxQ
			if(iQ > maxQ){
				maxQ = iQ;
				index = f;
			} 
			
			//eliminating potential mismatches for feature with Q = maxQ
			else if(iQ == maxQ && RTx[f] == RTx[index] && rQy[f] > rQy[index]){
				//SET_STRING_ELT(labels, index, mkChar("N"));
				index = f;
			}
			
			else if(iQ == maxQ && RTx[f] == RTx[index] && rQy[f] < rQy[index]){
				//SET_STRING_ELT(labels, f, mkChar("N"));
			}
		}
	}
	
	return(index);
}

/*
* Function: stopCondition
* -------------------
* Helper function for selectIterativeAnchors(). In main function loop,
* determines if there are any features with "P" anchor labels.
* If so, stop condition not met and main function loop continues. If not, 
* stop condition for main loop has been reached.
*
*/
int stopCondition(SEXP labels, int length)
{
	const char *label;
	
	for(int feat = 0; feat < length; feat++){
		label = CHAR(STRING_ELT(labels, feat));
		
		if(strcmp(label, "P") == 0)
			return 0;
	}
	
	return 1;
}


/*
* Function: updateLabels
* -------------------------
*
* Helper function for selectIterativeAnchors. All features within windX units
* in RTx or windY units in RTy of the matched identity are eliminated as 
* potential anchors by being labeled "N".
* 
*/
void updateLabels(SEXP labels, int f, double windrX, double windrY,
				   double* RTx, double* RTy, double* rts)
{
	const char *label = CHAR(STRING_ELT(labels, f)); 
			
	if(strcmp(label,"P") == 0){
		if(fabs(RTx[f] - rts[0]) < windrX || fabs(RTy[f] - rts[1]) < windrY){
			SET_STRING_ELT(labels, f, mkChar("N"));
		}
	}
			
	return;		
}


/*
* Function: selectIterativeAnchors
* ---------------------------------
* Iterative algorithm for selecting paired features whose retention times will
* serve as ordered pairs (anchors) for a non-linear retention time model.
*
* returns: updated anchor labels
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
	double* RTx = REAL(rtx);
	double* RTy = REAL(rty);
	double* rQx = REAL(Qx);
	double* rQy = REAL(Qy);
	
	double windrX = REAL(windX)[0];
	double windrY = REAL(windY)[0];
	int length = LENGTH(labels);
	
	//loop variables
	int nextAnchor;
	double rts[2];
		
	while(!stopCondition(labels_copy, length)){
		nextAnchor = findMaxQ(labels_copy, rQx, rQy, RTx, length);
				
		SET_STRING_ELT(labels_copy, nextAnchor, mkChar("A"));
	
		rts[0] = RTx[nextAnchor];
		rts[1] = RTy[nextAnchor];
		
		//update labels: exclude features within RT window of labeled anchors
		for(int feat = 0; feat < length; feat++){
			updateLabels(labels_copy, feat, windrX, windrY, RTx, RTy, rts);			
		}		
	}
		
	UNPROTECT(1);
	
	return labels_copy;
}




