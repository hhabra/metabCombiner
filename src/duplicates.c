#define EPS 1e-6
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

/*
 * File Description:
 * ------------------
 * Contains functions for finding, grouping and labeling duplicate metabolomics
 * features (i.e. with similar m/z and RT values).
 *
 *
 * PARAMETERS:
 *
 * mz: numeric m/z value of metabolomics features
 *
 * rt: numeric retention time value of metabolomics features
 *
 * tolMZ: numeric m/z tolerance for duplicate features
 *
 * tolRT: numeric retention time tolerance for duplicate features
 *
 * missing: proportion of feature's values which are missing across samples
 *
 * counts: feature's median or mean abundance value
 *
*/


/*
 * Function: comparePair
 * ----------------------
 * Helper function for findDuplicates. Used to compare pairs of features
 * deemed to be potential duplicates and assigns labels if found to meet
 * tolerances. A label will be assigned to the feature with greater percent
 * missingness or otherwise smaller counts (median/ mean abundance). If there
 * is a tie in both, the second feature is named a duplicate.
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
 * Function: groupPair
 * ---------------------
 * Helper function for groupDuplicates, used to determine if a pair of features
 * are duplicates (based on m/z and RT distances). Since the data is pre-sorted
 * by proportion missingness and abundance, the second of the pair (j) can be
 * assigned as the removable copy.
 *
*/
void groupPair(int i, int j, double *mz, double *rt, double tolMZ, double tolRT,
               int *dupgroup, int* dupremove)
{
    if(fabs(mz[i] - mz[j]) > (tolMZ + EPS))
        return;
    if(fabs(rt[i] - rt[j]) > (tolRT + EPS))
        return;
    dupgroup[j] = dupgroup[i];
    dupremove[j] = 1;
    return;
}


/*
 * Function: binDuplicates
 * ---------------------
 * This function bins features by m/z for eventual grouping as potential
 * duplicate features (as defined by m/z and RT similarity).
 *
 */
SEXP binDuplicates(SEXP mz, SEXP tolMZ)
{
    double *mz_c = REAL(mz);
    double tolMZ_c = REAL(tolMZ)[0];
    int length = LENGTH(mz);
    SEXP dupbin = PROTECT(allocVector(INTSXP, length));
    int *dupbin_c = INTEGER(dupbin);
    memset(dupbin_c, 0, length * sizeof(int));
    dupbin_c[0] = 1;
    for(int i = 1; i < length; i++){
        dupbin_c[i] = dupbin_c[i-1];
        if(mz_c[i] - mz_c[i-1] > (tolMZ_c + EPS))
            dupbin_c[i]++;
    }
    UNPROTECT(1);
    return dupbin;
}


/*
 * Function: groupDuplicates
 * ---------------------------
 * Features are grouped in a pairwise manner with a "master copy" if their m/z
 * and RT values are within specified tolerances. Data is previously arranged
 * and sorted by m/z bins, proportion missingness and (descending) abundance.
 * The second of each pair are assigned to the same group as the master copy and
 * marked as removable.
 *
 */
SEXP groupDuplicates(SEXP mz, SEXP rt, SEXP tolMZ, SEXP tolRT, SEXP dupbin)
{
    double *mz_c = REAL(mz);
    double *rt_c = REAL(rt);
    double tolMZ_c = REAL(tolMZ)[0];
    double tolRT_c = REAL(tolRT)[0];
    int *dupbin_c = INTEGER(dupbin);
    int length = LENGTH(mz);
    SEXP dupgroup = PROTECT(allocVector(INTSXP, length));
    SEXP dupremove = PROTECT(allocVector(INTSXP, length));
    int *dupgroup_c = INTEGER(dupgroup);
    int *dupremove_c = INTEGER(dupremove);
    memset(dupgroup_c, 0, length * sizeof(int));
    memset(dupremove_c, 0, length * sizeof(int));
    int nextgroup = 0, k = 1;
    for(int i = 0; i <= length-1; i++){
        if(dupremove_c[i] == 1)
            continue;
        k = 1;
        if(dupgroup_c[i] == 0)
            dupgroup_c[i] = ++nextgroup;
        while((i + k) < length){
            if(dupbin_c[i] != dupbin_c[i+k])
                break;
            groupPair(i, i+k, mz_c, rt_c, tolMZ_c, tolRT_c, dupgroup_c,
                      dupremove_c);
            k++;
        }
    }
    SEXP dupdata = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dupdata, 0, dupgroup);
    SET_VECTOR_ELT(dupdata, 1, dupremove);
    UNPROTECT(3);
    return dupdata;
}


/*
 * Function: findDuplicates
 * ------------------------------
 * Takes in a <m/z, rt> feature list sorted in m/z order and determines features
 * that are deemed "duplicates" as defined by very small differences in m/z & rt
 * (tolerances determined by tolMZ & tolRT).
 *
 */
SEXP findDuplicates(SEXP mz, SEXP rt, SEXP tolMZ, SEXP tolRT,
                    SEXP missing, SEXP counts)
{
    double *mz_c = REAL(mz);
    double *rt_c = REAL(rt);
    double tolMZ_c = REAL(tolMZ)[0];
    double tolRT_c = REAL(tolRT)[0];
    double *missing_c = REAL(missing);
    double *counts_c = REAL(counts);
    int length = LENGTH(mz);
    SEXP duplicates = PROTECT(allocVector(INTSXP, length));
    int *duplicates_c = INTEGER(duplicates);
    memset(duplicates_c, 0, length * sizeof(int));
    for(int i = 0; i < length-1; i++){
        if(duplicates_c[i] == 1)
            continue;
        int k = 1;
        while((i + k) < length){
            if(mz_c[i + k] - mz_c[i] >= (tolMZ_c + EPS))
                break;
            comparePair(i, i+k, rt_c, tolRT_c, missing_c, counts_c, duplicates_c);
            k++;
        }
    }
    UNPROTECT(1);
    return duplicates;
}


/*
 * Function: merge_duplicate_values
 * --------------------------------
 * This function handles merging of feature abundances in the mergeDuplicates
 * subroutine. For each group of duplicate features, the maximum sample value
 * is determined as the respective entry in the merged row.
 *
 */
SEXP merge_duplicate_values(SEXP values, SEXP dupgroup, SEXP ngroup, SEXP n)
{
    double* values_c = REAL(values);
    int* dupgroup_c = INTEGER(dupgroup);
    int ngroup_c = INTEGER(ngroup)[0];
    int n_c = INTEGER(n)[0];
    SEXP merged = PROTECT(allocVector(REALSXP, ngroup_c * n_c));
    double* merged_c = REAL(merged);
    int index = -1, group = -1;

    for(int i = 0; i < LENGTH(values); i++){
        if(dupgroup_c[i] ==  group)
            merged_c[index] = fmax(merged_c[index], values_c[i]);
        else
            merged_c[++index] = values_c[i];
        group = dupgroup_c[i];
        if(index >= ngroup_c * n_c)
            break;
    }
    UNPROTECT(1);
    return merged;
}


