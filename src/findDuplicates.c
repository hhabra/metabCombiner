#define EPS 1e-6
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

/*
* File Description:
* ------------------
*
* Contains function for finding and labeling duplicate metabolomics features, to
* be removed in the R function adjustData().
*
*/

/*
 * Function: comparePair
 * ----------------------
 * Helper function for determineDuplicates. Used to compare pairs of features deemed to 
 * be potential duplicates and assigns labels if found to meet tolerances. 
 *
 * The label will be assigned to the feature with greater percent missing values (first)  
 * smaller counts (median/ mean abundance). If there is a tie in both, the second feature
 * is named a duplicate.
 *
 *
 *
*/
void comparePair(int i, int j, double *rt, double tolRT, double *missing, 
				 double *counts, int *duplicates)
{	
	//check if a feature is a pre-designated duplicate
	if(duplicates[i] == 1 || duplicates[j] == 1)
		return;
	
	//check if feature is within R/T radius
	if(fabs(rt[i] - rt[j]) > tolRT + EPS)
	    return;
	
	//duplicate condition detected if above statement is false
	if(missing[i] > missing[j])
	    duplicates[i] = 1;
	
	else if(missing[i] < missing[j])
	    duplicates[j] = 1;
	    
    else if (counts[i] > counts[j])
        duplicates[j] = 1;
        
    else if (counts[i] < counts[j])
        duplicates[i] = 1;
	 
	else
		duplicates[j] = 1;
		
	return;
}

/*
 * Function: findDuplicates
 * ------------------------------
 * Takes in a <m/z, rt> feature list sorted in m/z order and determines features that 
 * are deemed "duplicates" as defined by very small differences in m/z & rt (tolerances 
 * determined by tolMZ & tolRT. 
 *
 * Returns a vector of characters labeled 'o' (okay) or 'd' (duplicate)
*/
SEXP findDuplicates(SEXP mz, SEXP rt, SEXP tolMZ, SEXP tolRT, 
					SEXP missing, SEXP counts)
{
	//obtaining C types
    double *mz_r = REAL(mz);
    double *rt_r = REAL(rt);
    
    double tol_MZ = REAL(tolMZ)[0];
    double tol_RT = REAL(tolRT)[0];
    double *missing_r = REAL(missing);
    double *counts_r = REAL(counts);
    	
	int length = LENGTH(mz);
	
    SEXP duplicates = PROTECT(allocVector(INTSXP, length));
    int *duplicates_r = INTEGER(duplicates);
    memset(duplicates_r, 0, length * sizeof(int));
    
    for(int i = 0; i < length-1; i++){
        if(duplicates_r[i] == 1)
            continue;
		
        int k = 1;
        
        while(mz_r[i + k] - mz_r[i] < (tol_MZ + EPS)){
            comparePair(i, i+k, rt_r, tol_RT, missing_r, counts_r, duplicates_r);
			k++;
			if((i + k) >= length)      //handling boundary condition
                break; 
        }    
    }
         
    UNPROTECT(1);
    
    return duplicates;
}
