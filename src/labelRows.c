#define EPS 1e-6
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <stdbool.h>

/*
 * File Description:
 * ------------------
 * This file contains a function for determining whether alignment rows in the
 * metabCombiner report table are removable or in conflict with other rows.
 *
*/


/*
 * findGroup:
 * -------------
 *
 * Finds the starting and ending index of a group of alignments
 *
 * Returns: cursor value of next group or final index
 *
 * PARAMETERS:
 *
 * cursor: initial look-up index of group
 *
 * start: pointer to starting group index
 *
 * end: pointer to ending group index
 *
*/
int detectGroup(int cursor, int* group, int* start, int *end, int length)
{
	*start = cursor;    //start of the group

	int groupVal = group[cursor];

	while(cursor < length && group[cursor] == groupVal)
		cursor++;

	*end = cursor - 1;  //end of the group

	return cursor;
}

/*
 * balancedGroups:
 * ----------------
 *
 * Special handler for groups deemed "balanced", defined as a group with an equal
 * number of features from both datasets and no conflicting top-matches between
 * any pairs of features. If a m/z group meets this definition, all but the
 * top-matches
 *
 * Returns: Updated index of last row of group, pointed to by "end"
 *
 * PARAMETERS:
 *
 * start: pointer to starting group index
 *
 * end: pointer to ending group index
 *
*/
int balancedGroups(SEXP labels, int *start, int *end, int* rankX, int* rankY)
{
	int nrows = *end - *start + 1;
	int root = 1;

	//check if nrows is a square number; if not, group is unbalanced
	while(root * root < nrows)
		root++;
	if(root * root != nrows)
		return *end;

	int topMatches = 0;

	//check if there are appropriate number of top-matches
	for(int j = *start; j <= *end; j++){
		if(rankX[j] == 1 && rankY[j] == 1)    //a top-match
			topMatches++;

		else if (rankX[j] == 1 || rankX[j] == 1) //a conflicting top-match
			return *end;
	}

	if(topMatches != root)
		return *end;

	//balanced group detected; label non-identities and non-topMatches as removables
	int updateEnd = 1;

	for(int k = *end; k >= *start; k--){
		//signal to stop updating end index
		if(updateEnd && rankX[k] == 1 && rankY[k] == 1)
			updateEnd = 0;

		else if(rankX[k] != 1 || rankY[k] != 1){
			if(strcmp("", CHAR(STRING_ELT(labels,k))) == 0)
				SET_STRING_ELT(labels, k, mkChar("REMOVE"));

			if(updateEnd)
				(*end)--;

		}
	}

	return *end;
}

/*
 * filterScoreRankErr:
 * ---------------------
 *
 * Finds rows with score < minScore, rankX > maxRankX, rankY > maxRankY, or rterr >
 * maxRTerr. Updates and returns 'end' index pointer to lowest non-removed row (all other
 * rows need not be considered).
 *
*/
int filterScoreRankErr(SEXP labels, int *start, int *end, double* score, int* rankX,
						int* rankY, double minScore, int maxRankX, int maxRankY,
						double maxRTerr, double* rterr)
{

	int updateEnd = 1;

	for(int i = *end; i >= *start; i--){
		if(strcmp("", CHAR(STRING_ELT(labels,i))) == 0){
			if(rankX[i] > maxRankX ||  rankY[i] > maxRankY ||
				score[i] < minScore || rterr[i] > maxRTerr)
				SET_STRING_ELT(labels, i, mkChar("REMOVE"));

			else
				updateEnd = 0;

			if(updateEnd)
				(*end)--;
		}
	}

	return *end;
}


/*
 * Function: detect_con_score
 * ------------------------
 *
 * Score-based conflict detection. A pair of alignments are deemed in conflict if they
 * are within a limited score difference of the head alignment score. Those that fail to
 * meet this threshold are deemed removable.
 *
*/
void detect_con_score(SEXP labels, int* sub, int* alt, int ri, int rj, double* delta,
                     double* mz, double* rt, double* score, int* head, int type, int* max)
{
	double scoreGap = delta[0];

	if(head[ri] == 0)
		head[ri] = ri;

	double scoreHead = score[head[ri]];

	if(fabs(score[rj] - scoreHead) > scoreGap){
		if(strcmp("", CHAR(STRING_ELT(labels,rj))) == 0)
		    SET_STRING_ELT(labels, rj, mkChar("REMOVE"));

		return;
	}

	if(sub[ri] == 0){
	    if(strcmp("", CHAR(STRING_ELT(labels,ri))) == 0)
			SET_STRING_ELT(labels, ri, mkChar("CONFLICT"));

        sub[ri] = ++*max;
	}

	if(sub[rj] == 0){
		if(strcmp("", CHAR(STRING_ELT(labels,rj))) == 0)
		    SET_STRING_ELT(labels, rj, mkChar("CONFLICT"));
		sub[rj] = sub[ri];
		head[rj] = head[ri];
	}

	else if(sub[rj] != sub[ri])
		alt[rj] = sub[ri];
}



/*
 * Function: detect_con_mzrt
 * ------------------------
 *
 * mzrt-based conflict detection. A pair of alignment rows are deemed in conflict if the
 * lower-scoring alignment matches with one dataset feature while the other is within
 * some deviation in m/z & rt (as determined by conflict argument) of the higher-scoring
 * alignment. If they fall outside of this threshold, then the lower-scoring row is
 * deemed removable.
 *
*/
void detect_con_mzrt(SEXP labels, int* sub, int* alt, int ri, int rj, double* delta,
                     double* mz, double* rt, double* score, int* head, int type, int* max)
{
	double mzTol = (type == 1) ? delta[0] : delta[2];
	double rtTol = (type == 1) ? delta[1] : delta[3];

	if(fabs(mz[rj] - mz[ri]) > mzTol || fabs(rt[rj] - rt[ri]) > rtTol){
		if(strcmp("", CHAR(STRING_ELT(labels,rj))) == 0)
		    SET_STRING_ELT(labels, rj, mkChar("REMOVE"));

		return;
	}

	//unassigned subgroup
	if(sub[ri] == 0){
		if(strcmp("", CHAR(STRING_ELT(labels,ri))) == 0)
			SET_STRING_ELT(labels, ri, mkChar("CONFLICT"));

        sub[ri] = ++*max;
	}


	if(sub[rj] == 0){
		if(strcmp("", CHAR(STRING_ELT(labels,rj))) == 0)
		    SET_STRING_ELT(labels, rj, mkChar("CONFLICT"));
		sub[rj] = sub[ri];
	}

	else if(sub[rj] != sub[ri])
		alt[rj] = sub[ri];
}

/*
 * Function: findCons
 * ------------------------
 * This function loops through the features of a group that remain after score and rank
 * filtering, finding pairs of alignments which may be in conflict over a single feature.
 *
*/
void findCons(SEXP labels, int* sub, int* alt, int* max, int *start,
              int *end, double* delta, double *mzx, double *mzy,
              double *rtx, double *rty, double* score, int* head,
			  void (*detect)(SEXP, int*, int*, int, int, double*, double*, double*,
			                 double*, int*, int, int*))
{
	for(int ri = *start; ri <= *end; ri++){
		if(strcmp("REMOVE", CHAR(STRING_ELT(labels,ri))) == 0)
			continue;

		for(int rj = ri+1; rj <= *end; rj++){
			if(strcmp("REMOVE", CHAR(STRING_ELT(labels,rj))) == 0)
				continue;

			if(sub[rj] > 0 && alt[rj] > 0)
				continue;

			//determine if a pair of alignments are conflicting
			if(fabs(mzx[rj] - mzx[ri]) < EPS && fabs(rtx[rj] - rtx[ri]) < EPS)
				detect(labels, sub, alt, ri, rj, delta, mzy, rty, score, head, 2, max);

			if(fabs(mzy[rj] - mzy[ri]) < EPS && fabs(rty[rj] - rty[ri]) < EPS)
				detect(labels, sub, alt, ri, rj, delta, mzx, rtx, score, head, 1, max);

		}
	}
}


/*
 * Function: labelRows
 * --------------------------
 *
 * Called from the labelRows() function in R to find removable & conflicting rows from
 * metabCombiner report. Rows in the labels column shall be identified as either:
 *
 * 1) "IDENTITY" - Features with matching pre-determined identities. These rows not
 *    labeled as removable, even if score or rank restrictions are violated.
 *
 * 2) "CONFLICT" - Rows that have a conflicting match that is within score or m/z and rt
 *	  tolerances defined by "conflict", as well as meeting rank and score
 *	  restrictions defined by maxRankX, maxRankY, and minScore.
 *
 * 3) "REMOVE" - Rows that fail to meet at least one of the criteria (score <
 *    minScore, rankX > maxRankX, rankY > maxRankY, ...).
 *
 * 4) "" (empty) - feature pair alignments with no competitive matches
 *
 * Returns: Vector of labels.
 *
 * PARAMETERS:
 *
 * labels: report table row labels
 *
 * subgroup: subdivision of feature groupings based on conflict assignments
 *
 * alt: alternative subgrouping of features that conflict with multiple subgroups
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
 * rankX: Score ranking for X dataset features
 *
 * rankY: Score ranking for Y dataset features
 *
 * group: Integer feature group values (as determined by m/z)
 *
 * balanced: Boolean option to process "balanced" groups (see: balancedGroups()).
 *
 * delta: either a numeric constant score difference or a length 4 vector specifying m/z
 * and rt tolerances for determining if pairs of rows have a conflicting alignment.
 *
 * minScore: Minimum allowable feature similarity scores.
 *
 * maxRankX: Maximum allowable integer feature score rank for X dataset features
 *
 * maxRankY: Maximum allowable integer feature score rank for Y dataset features
 *
 * method: integer conflict detection method (1 = score, 2 = mzrt)
 *
 * maxRTerr: Maximum allowable retention time fitting error
 *
 * rterr: retention time fitting errors
 *
*/
SEXP labelRows(SEXP labels, SEXP subgroup, SEXP alt, SEXP mzx, SEXP mzy, SEXP rtx,
               SEXP rty, SEXP score, SEXP rankX, SEXP rankY, SEXP group, SEXP balanced,
			   SEXP delta, SEXP minScore, SEXP maxRankX, SEXP maxRankY, SEXP method,
			   SEXP maxRTerr, SEXP rterr)
{
	SEXP labels_c = PROTECT(duplicate(labels));

	//R variables in C
	int* subgroup_c = INTEGER(subgroup);
	int* alt_c = INTEGER(alt);

	double* mzx_c = REAL(mzx);
	double* mzy_c = REAL(mzy);
	double* rtx_c = REAL(rtx);
	double* rty_c = REAL(rty);
	double* score_c = REAL(score);
	int* group_c = INTEGER(group);
	int* rankX_c = INTEGER(rankX);
	int* rankY_c = INTEGER(rankY);
	bool balanced_c = LOGICAL(balanced);
	double* delta_c = REAL(delta);
	double minScore_c = REAL(minScore)[0];
	int maxRankX_c = INTEGER(maxRankX)[0];
	int maxRankY_c = INTEGER(maxRankY)[0];

	double maxRTerr_c = REAL(maxRTerr)[0];
	double* rterr_c = REAL(rterr);

	//choosing conflict detection method
	int method_c = INTEGER(method)[0];
	void (*detect_fun)(SEXP, int*, int*, int, int, double*, double*, double*, double*,
	                    int*, int, int*);

	int* head;   //head index of subgroup; special variable for score-conflict method

	if (method_c == 2){
	    detect_fun = detect_con_mzrt;
		head = calloc(1, sizeof(int));
	}

	else{
        detect_fun = detect_con_score;
        head = calloc(LENGTH(group), sizeof(int));
	}

	//loop variables
	int cursor = 0;
	int* start = calloc(1,sizeof(int));
	int* end = calloc(1,sizeof(int));
	int* maxSub = calloc(1, sizeof(int));    //current subgroup number

	//Main loop
	while(cursor < LENGTH(group)){
		cursor = detectGroup(cursor, group_c, start, end, LENGTH(group));

		if(balanced_c){
			*end = balancedGroups(labels_c, start, end, rankX_c, rankY_c);
		}

		*end = filterScoreRankErr(labels_c, start, end, score_c, rankX_c,
						        rankY_c, minScore_c, maxRankX_c, maxRankY_c,
						        maxRTerr_c, rterr_c);

		findCons(labels_c, subgroup_c, alt_c, maxSub, start, end, delta_c, mzx_c,
		        mzy_c, rtx_c, rty_c, score_c, head, detect_fun);
	}

	UNPROTECT(1);

	return labels_c;
}
