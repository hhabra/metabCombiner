#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

/*
* File Description:
* ------------------
* Contains function to bin metabolomics features from 2 data sets by m/z values,
* to be used within mzGroup() function. Each m/z group must contain features
* from both datasets.
*/

/*
* Function: binByMZ
* -------------------
*
* binByMZ performs continuous binning of features by m/z (pre-sorted).
* Features with m/z values within 'gap' Daltons are considered for a potential
* mzGroup. Increment group count if mzGroup contains features from both datasets.
*
* PARAMETERS:
*
* mz: A pooled and sorted list of m/z values from two separate feature tables
*
* datasets: either 'x' or 'y'
*
* gap: m/z gap for binning
*/
SEXP binByMZ(SEXP mz, SEXP datasets, SEXP gap)
{
  	//obtaining C types
  	double *mz_c = REAL(mz);
  	double gap_c = REAL(gap)[0];
  	int length = LENGTH(mz);

  	//groups will contain the final returned vector; initialize to 0 vector
  	SEXP groups = PROTECT(allocVector(INTSXP, length));
  	int *groups_c = INTEGER(groups);
  	memset(groups_c, 0, length * sizeof(int));

  	//initiating loop variables
  	double diff;
  	int cond1 = 0, cond2 = 0, bin = 1, start, count;
  	const char *set1, *set2;

  	//bin by m/z and assign new groups if:
  	//condition 1: difference between consecutive m/z is less than gap
  	//condition 2: at least one feature found from both datasets
  	for(int i = 0; i < length-1; i++){
    		diff = mz_c[i+1] - mz_c[i];

    		//initiate a potential group
    		if (!(cond1) && diff < gap_c){
      			cond1 = 1;
      			start = i;
      			count = 0;
    		}

    		//elongate group and check if cond2 met
    		if(cond1 && diff < gap_c){
    			 count++;
    			 if(!cond2){
      				set1 = CHAR(STRING_ELT(datasets, i));
      				set2 = CHAR(STRING_ELT(datasets, i+1));

      				if(abs(strcmp(set1,set2)))    //compare strings
      					  cond2 = 1;
    			}
  		}

  		//reset and update group count, if condition 2 met
  		if(i == (length - 2) || diff > gap_c){
    			 cond1 = 0;
    			 if(cond2){
    				    for(int j = start; j <= start+count; j++){
    				        groups_c[j] = bin;
    			      }
    				    bin++;
    			  }
    			 cond2 = 0;
    		}
  	}
  	UNPROTECT(1);

	  return(groups);
}
