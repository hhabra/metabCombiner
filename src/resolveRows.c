#define EPS 1e-5
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include "labelRows.h"

/*
 * File Description:
 * ------------------
 * This file contains a function for automatically resolving alignment conflicts to
 * obtaining the most likely set of 1-1 matching feature pairs
 *
*/


/*
 * Function: sharedFeatCheck
 * ------------------------------
 * For a pair of rows within a subgroup, checks if they share an X or Y feature in common
 * or if a row is previously assigned as removable.
 *
*/
int sharedFeatCheck(SEXP labels, int ri, int rj, double* mzx, double *mzy, double *rtx, 
                	double *rty)
{	
	if(strcmp("REMOVE", CHAR(STRING_ELT(labels,rj))) == 0)
		return 0;
	
	if(fabs(mzx[ri] - mzx[rj]) < EPS && fabs(rtx[ri] - rtx[rj]) < EPS)
		return 0;

	if(fabs(mzy[ri] - mzy[rj]) < EPS && fabs(rty[ri] - rty[rj]) < EPS)
		return 0;

	return 1;

}

/*
 * Function: retOrderCheck
 * ------------------------------
 * For a pair of rows within a subgroup, checks if the retention time order of the
 * respective X & Y features is consistent
 *
*/
int retOrderCheck(int ri, int rj, double *rtx, double *rty)
{
	//RT order rearrangement detected
	if((rtx[ri] > rtx[rj] && rty[ri] < rty[rj]) || 
		(rtx[ri] < rtx[rj] && rty[ri] > rty[rj]))
		return 0;
	
	else
		return 1;
}

/*
 * Function: form_rowset
 * -------------------------------
 * This takes an initially constructed integer pointer (rowset) containing the head row
 * of a collection of 1-1 aligned feature pairs. Rows within this collection must pass
 * the shared feature and retention order checks to be included; passing these checks
 * updates the rowset and rowsetnum pointers.
 *
*/
void form_rowset(SEXP labels, int* rowset, int* rowsetnum, int end, double* mzx, 
				double *mzy, double *rtx, double *rty, int retOrder)
{
	int setSize = 1;  //total number in rowset
	int first = rowset[0];  //initial row number in rowset
	
	int setIndex, newIndex, shareCheck = 1, rtOrderCheck = 1; 
	
	for(newIndex = first+1; newIndex <= end; newIndex++){
		if(rowsetnum[newIndex] > 0)
			continue;
	    for(int i = 0; i < setSize; i++){
	    	setIndex = rowset[i];   //indices of rows in rowset
	    	
	    	shareCheck = sharedFeatCheck(labels, setIndex, newIndex, mzx, mzy, rtx, rty);
	    	
	    	if(retOrder && shareCheck)
	    		rtOrderCheck = retOrderCheck(setIndex, newIndex, rtx, rty);
	    	
	    	if(shareCheck == 0 || rtOrderCheck == 0)
	    		break;
	    }
	    	    
	    if(shareCheck && rtOrderCheck){
	    	rowset[setSize++] = newIndex;
	    	rowsetnum[newIndex] = rowsetnum[setIndex];
	    }	    
	}
}

/*
 * Function: score_rowset
 * ----------------------------
 *
 * Calculates the score of a collection of rows ("rowset") by the sum of individual 
 * constituent rows, then assigns to the 'resolveScore' array that is presented as  
 * an output column in the resulting table.
 *
 *
*/
double score_rowset(int* rowset, double* score, int nrow, double *resolveScore)
{
	int index = 0;
	double setScore = 0;
	
	while(index < nrow && rowset[index] != -1)
		setScore += score[rowset[index++]];
	
	index = 0;
	while(index < nrow && rowset[index] != -1)
		resolveScore[rowset[index++]] = setScore;	
	
	return setScore;
}


/*
 * Function: update_rowset_labels
 * ---------------------------------
 *
 * For the highest scoring rowset, updates the labels corresponding to the constituent
 * rows to "RESOLVED" from "CONFLICT"
 *
*/
void update_rowset_labels(SEXP labels, int* rowset, int nrow)
{
	int index = 0;
	
	while(index < nrow && rowset[index] != -1){
		if(strcmp("CONFLICT", CHAR(STRING_ELT(labels,rowset[index]))) == 0)
			SET_STRING_ELT(labels, rowset[index], mkChar("RESOLVED"));
		index++;
	}
}

/*
 * Function: resolve
 * -----------------------
 *
 * This function resolves a single subgroup by forming the highest-scoring combination
 * of 1-1 paired alignment rows. The sets of non-conflicting rows are assigned 
 * (form_rowset), then scored (score_rowset). The best-scoring set of rows are then
 * re-annotated (update_rowset_labels) before the process concludes.
 *
*/
void resolve(SEXP labels, int* rowsetnum, int start, int end, double* score, 
            double* mzx, double *mzy, double *rtx, double *rty, int retOrderOpt,
            double* resolveScore)
{
	int nrow = end - start + 1; 
	int **rowsets = (int**) calloc(nrow, sizeof(int*));   //tracks 1-1 row groupings
	int nsets = 0, bestset = 0;
	
	//forming the possible 1-1 feature pair sets
	for(int row = start; row <= end; row++){
		if(strcmp("REMOVE", CHAR(STRING_ELT(labels,row))) == 0 || rowsetnum[row] > 0)
			continue;		
		else{
			rowsetnum[row] = ++nsets;
			rowsets[nsets-1] = calloc(nrow, sizeof(int*));
			memset(rowsets[nsets-1], -1, sizeof(int) * nrow);
			rowsets[nsets-1][0] = row;
			
			form_rowset(labels, rowsets[nsets-1], rowsetnum, end, mzx, mzy, rtx, 
			            rty, retOrderOpt);
			
		}
	}
		
	//scoring the sets 
	double setscore = 0, bestscore = 0;
	for(int set = 0; set < nsets; set++){
		setscore = score_rowset(rowsets[set], score, nrow, resolveScore);
		if(bestscore < setscore){
			bestscore = setscore;
			bestset = set;
		} 
	}
	
	update_rowset_labels(labels, rowsets[bestset], nrow);
	
	//freeing data
	for(int n = 0; n < nrow; n++){
		free(rowsets[n]);
	}
	
	free(rowsets);
}  


/*
 * Function: resolveRows
 * -----------------------------
 *
 * This function attempts to find, within each previously computed subgroup of rows 
 * within the combined table report, the set of 1-1 paired alignments between features
 * in datasets X and Y that maximizes the sum of scores (resolveScore). Each potential
 * row is inserted into at most one collection of rows ("rowset"), subject to two 
 * conditions, the latter of which is optional: 1) the rows must share no features in
 * common, i.e. no feature is repeated more than once; 2) the retention order must be
 * consistent between the respective X and Y features. The resulting rowsets are scored
 * and rows belonging to the best-scoring rowset will be annotated as "RESOLVED". 
 * 
 * subgroup: subdivision of feature groupings based on conflict assignments
 *
 * mzx: m/z values from dataset X.
 *
 * mzy: m/z values from dataset Y.
 *
 * rtx: retention time values from dataset X.
 *
 * rty: retention time values from dataset Y.
 *
 * score: Calculated similarity between features from datasets X & Y.
 *
 * retOrder: integer option to expect consistent retention order when forming row sets
 *
 * resolveScore: a double pointer storing the rowset scores computed by the sum of 
 *                individual scores of the constituent rows
 *
 *
 *
*/
SEXP resolveRows(SEXP labels, SEXP subgroup, SEXP mzx, SEXP mzy, SEXP rtx, SEXP rty, 
                SEXP score, SEXP retOrder, SEXP resolveScore)
{
    SEXP labels_c = PROTECT(duplicate(labels));
    int* subgroup_c = INTEGER(subgroup);
    double* mzx_c = REAL(mzx);
	double* mzy_c = REAL(mzy);
	double* rtx_c = REAL(rtx);
	double* rty_c = REAL(rty);
    double* score_c = REAL(score);
    double* resolveScore_c = REAL(resolveScore);
    int retOrder_c = INTEGER(retOrder)[0];
    
    //loop variables
	int cursor = 0;
	int* start = calloc(1,sizeof(int));
	int* end = calloc(1,sizeof(int));
	int* rowsetnum = calloc(LENGTH(subgroup),sizeof(int));
		
	while(cursor < LENGTH(subgroup)){
		cursor = detectGroup(cursor, subgroup_c, start, end, LENGTH(subgroup));
				
		resolve(labels_c, rowsetnum, *start, *end, score_c, mzx_c, mzy_c, rtx_c, rty_c, 
		        retOrder_c, resolveScore_c);		
	}

	UNPROTECT(1);
	return labels_c;
}





