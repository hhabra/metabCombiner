#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

/*
* File Description:
* ------------------
*
* Contains function for finding and labeling duplicate metabolomics features, to be removed in
* the R function adjustData().
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
*/
void comparePair(int i, int j, double *RT, double rtolRT, double *rmissing, 
				 double *rcounts, int *rduplicates)
{
	double eps = 1e-6;
	if(rduplicates[i] == 1 || rduplicates[j] == 1)
		return;
		
	if(fabs(RT[i] - RT[j]) > rtolRT + eps)
	    return;
	
	//duplicate condition detected if above statement is false
	if(rmissing[i] > rmissing[j])
	    rduplicates[i] = 1;
	
	else if(rmissing[i] < rmissing[j])
	    rduplicates[j] = 1;
	    
    else if (rcounts[i] > rcounts[j])
        rduplicates[j] = 1;
        
    else if (rcounts[i] < rcounts[j])
        rduplicates[i] = 1;
	 
	else
		rduplicates[j] = 1;
		
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
SEXP findDuplicates(SEXP mz, SEXP rt, SEXP length, SEXP tolMZ, SEXP tolRT, 
					SEXP missing, SEXP counts)
{
	//obtaining C types
    double *rmz = REAL(mz);
    double *RT = REAL(rt);
    int rlength = INTEGER(length)[0];
    double rtolMZ = REAL(tolMZ)[0];
    double rtolRT = REAL(tolRT)[0];
    double *rmissing = REAL(missing);
    double *rcounts = REAL(counts);
    
    double eps = 1e-6;
	
    SEXP duplicates = PROTECT(allocVector(INTSXP, rlength));
    int *rduplicates = INTEGER(duplicates);
    memset(rduplicates, 0, rlength * sizeof(int));
    
    for(int i = 0; i < rlength-2; i++){
        if(rduplicates[i] == 1)
            continue;
        int k = 1;
        
        while(rmz[i + k] - rmz[i] < rtolMZ + eps){
            comparePair(i, i+k, RT, rtolRT, rmissing, rcounts, rduplicates);
            
            k = k+1;
            if((i + k) > rlength - 1)      //handling boundary condition
                break; 
        }    
    }
         
    UNPROTECT(1);
    
    return duplicates;
}
