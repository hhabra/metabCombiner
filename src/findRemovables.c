#define EPS 1e-6
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <stdbool.h>

/*
 * File Description:
 * ------------------
 *
 *
 *
 *
*/



/*
 * findGroup:
 * -------------
 *
 * Finds the next group of paired features, starting from "cursor" index.
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
 * labels:
 *
 * start: pointer to starting group index
 *
 * end: pointer to ending group index
 *
 * rankX: vector of rankings for features in dataset X
 *
 * rankY: vector of rankings for features in dataset Y
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
 * filterScoreAndRank:
 * ---------------------
 *
 * Finds rows with score less than minScore and rankX & rankY in excess of
 * maxRankX & maxRankY, respectively. Updates end index pointer to the last row.
 *
 * PARAMETERS:
 *
 * labels:
 *
 * start:
 *
 * end:
 *
 * score:
 *
 * rankX:
 *
 * rankY:
 *
 * minScore:
 *
 * maxRankX:
 *
 * maxRankY:
 *
*/
int filterScoreAndRank(SEXP labels, int *start, int *end, double* score, int* rankX,
						int* rankY, double minScore, int maxRankX, int maxRankY)
{

	int updateEnd = 1;

	for(int i = *end; i >= *start; i--){
		if(strcmp("", CHAR(STRING_ELT(labels,i))) == 0){
			if(rankX[i] > maxRankX ||  rankY[i] > maxRankY || score[i] < minScore)
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
 * Function: detectConflict
 * ------------------------
 *
 * labels:
 *
 * i
 *
 *
 *
 *
*/
void detectConflict(SEXP labels, int row1, int row2, double mzTol, double rtTol,
				  double* mz, double *rt)
{
	if(fabs(mz[row2] - mz[row1]) > mzTol || fabs(rt[row2] - rt[row1]) > rtTol){
		SET_STRING_ELT(labels, row2, mkChar("REMOVE"));
		return;
	}

	if(strcmp("", CHAR(STRING_ELT(labels,row1))) == 0)
		SET_STRING_ELT(labels, row1, mkChar("CONFLICT"));

	if(strcmp("", CHAR(STRING_ELT(labels,row2))) == 0)
		SET_STRING_ELT(labels, row2, mkChar("CONFLICT"));

}
/*
 * Function: conflictCheck
 * ------------------------
 *
 * PARAMETERS:
 *
 * labels:
 *
 * start:
 *
 * end:
 *
 * mzx:
 *
 * mzy:
 *
 * rtx:
 *
 * rty:
 *
 * conflict:
*/
void conflictCheck(SEXP labels, int *start, int *end, double *mzx, double *mzy,
				 double *rtx, double *rty, double *conflict)
{
	double mzTolX = conflict[0];
	double rtTolX = conflict[1];
	double mzTolY = conflict[2];
	double rtTolY = conflict[3];

	for(int i = *start; i <= *end; i++){
		if(strcmp("REMOVE", CHAR(STRING_ELT(labels,i))) == 0)
			continue;

		for(int j = i+1; j <= *end; j++){
			if(strcmp("REMOVE", CHAR(STRING_ELT(labels,j))) == 0)
				continue;

			//determine if an conflicting alignment is found
			if(fabs(mzx[j] - mzx[i]) < EPS && fabs(rtx[j] - rtx[i]) < EPS){
				detectConflict(labels, i, j, mzTolY, rtTolY, mzy, rty);
			}

			else if(fabs(mzy[j] - mzy[i]) < EPS && fabs(rty[j] - rty[i]) < EPS){
				detectConflict(labels, i, j, mzTolX, rtTolX, mzx, rtx);
			}
		}
	}
}

/*
 * Function: findRemovables
 * --------------------------
 *
 * Called from the reduceTable() function in R to find removable rows from
 * Combiner report. Rows in the labels column are identified as either:
 *
 * 1) "IDENTITY" - Features with matching pre-determined identities. These rows not
 *    labeled as removable, even if score or rank restrictions are violated.
 *
 * 2) "CONFLICT" - Rows that have a conflicting match that is within m/z and rt
 *	  tolerances defined by "conflict", as well as meeting rank and score
 *	  restrictions defined by rankX, rankY, and minScore.
 *
 * 3) "REMOVE" - Rows that fail to meet at least one of the criteria (score <
 *    minScore, rankX > maxRankX, rankY > maxRankY, non-conflicting non-top match).
 *
 * 4) "" (empty) - Top-ranked feature alignments with no contested matches
 *
 * Returns: Vector of labels.
 *
 * PARAMETERS:
 *
 * labels: Combiner report table row labels (contain pre-assigned identity tags).
 *
 * mzx: m/z values from dataset X.
 *
 * mzy: m/z values from dataset Y.
 *
 * rtx: retention time values from dataset X.
 *
 * rty: retention time values from dataset Y.
 *
 * score: Calculated similarity scores between features from datasets X & Y.
 *
 * rankX: Score ranking for X dataset features
 *
 * rankY: Score ranking for Y dataset features
 *
 * group: Integer feature group values (as determined by m/z)
 *
 * balanced: Boolean option to process "balanced" groups (see: balancedGroups()).
 *
 * conflict: Length 4 vector specifying m/z and rt tolerances for determining
 * if pairs of rows have a conflicting alignment.
 *
 * minScore: Minimum allowable feature similarity scores.
 *
 * maxRankX: Maximum allowable integer feature score rank for X dataset features
 *
 * maxRankY: Maximum allowable integer feature score rank for Y dataset features
 *
*/
SEXP findRemovables(SEXP labels, SEXP mzx, SEXP mzy, SEXP rtx, SEXP rty,
					SEXP score, SEXP rankX, SEXP rankY, SEXP group, SEXP balanced,
					SEXP conflict, SEXP minScore, SEXP maxRankX, SEXP maxRankY)
{
	SEXP labels_copy = PROTECT(duplicate(labels));

	double* mz_x = REAL(mzx);
	double* mz_y = REAL(mzy);
	double* rt_x = REAL(rtx);
	double* rt_y = REAL(rty);
	double* score_r = REAL(score);
	int* group_r = INTEGER(group);
	int* rank_x = INTEGER(rankX);
	int* rank_y = INTEGER(rankY);
	bool balanced_r = LOGICAL(balanced);

	double* conflict_r = REAL(conflict);
	double min_Score = REAL(minScore)[0];
	int maxRank_X = INTEGER(maxRankX)[0];
	int maxRank_Y = INTEGER(maxRankY)[0];

	//loop variables
	int cursor = 0;
	int* start = malloc(sizeof(int));
	int* end = malloc(sizeof(int));

	while(cursor < LENGTH(group)){
		cursor = detectGroup(cursor, group_r, start, end, LENGTH(group));

		if(balanced_r){
			*end = balancedGroups(labels_copy, start, end, rank_x, rank_y);
		}

		*end = filterScoreAndRank(labels_copy, start, end, score_r, rank_x,
						          rank_y, min_Score, maxRank_X, maxRank_Y);


		conflictCheck(labels_copy, start, end, mz_x, mz_y, rt_x, rt_y, conflict_r);
	}

	UNPROTECT(1);

	return labels_copy;
}

