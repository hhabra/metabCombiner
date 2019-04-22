#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

/*
* File Description:
* ------------------
*
* Contains function to bin metabolomics features by m/z values, 
* to be used within mzGroup() function. Each m/z group must contain features 
* from both datasets.
*/

/*
* Function: binByMZ
* ------------------- 
*
* binByMZ performs continuous binning of features by m/z (pre-sorted). 
* Features with m/z values within 'gap' Daltons are considered for a potential 
* mzGroup; increment bin count if group contains features from both datasets.
* 
*/
SEXP binByMZ(SEXP mz, SEXP length, SEXP datasets, SEXP gap){
	//obtaining C types
	double *rmz = REAL(mz);
	int rlength = INTEGER(length)[0];	
	double rgap = REAL(gap)[0];
		
	//groups will contain the final returned vector; initialize to 0 vector
	SEXP groups = PROTECT(allocVector(INTSXP, rlength));
	int *rgroups = INTEGER(groups);
	memset(rgroups, 0, rlength * sizeof(int));    

	//initiating loop variables
	double diff;
	int cond1 = 0, cond2 = 0, bin = 1, start, count;
	const char *d1, *d2;
	
	//bin by MZ and assign new groups if:
	//condition 1: consecutive mz gap less than gap
	//condition 2: at least one feature from both datasets
	for(int i = 0; i < rlength-1; i++){
		diff = rmz[i+1] - rmz[i];
		
		//initiate a potential group
		if (!(cond1) && diff < rgap){
			cond1 = 1;
			start = i;
			count = 0;
		}
		
		//elongate group and check if cond2 met
		if(cond1 && diff < rgap){
			count++;
			
			if(!cond2){
				d1 = CHAR(STRING_ELT(datasets, i));
				d2 = CHAR(STRING_ELT(datasets, i+1));
				
				if(abs(strcmp(d1,d2)))    //compare strings
					cond2 = 1;
			}
		}
	
		//reset and update group count, if condition 2 met
		if(diff > rgap || i == (rlength - 1)){
			cond1 = 0;
			
			if(cond2){
				for(int j = start; j <= start+count; j++){
					rgroups[j] = bin;
				}
				
				bin++;
			}
			
			cond2 = 0;
		}
	}
	
	UNPROTECT(1);
	
	return(groups);
}