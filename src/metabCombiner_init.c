#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP binByMZ(SEXP, SEXP, SEXP);
extern SEXP binDuplicates(SEXP, SEXP);
extern SEXP findDuplicates(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP groupDuplicates(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP labelRows(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                      SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP merge_duplicate_values(SEXP, SEXP, SEXP, SEXP);

extern SEXP selectAnchorsByID(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP selectIterativeAnchors(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP write2file(SEXP, SEXP, SEXP, SEXP);
extern SEXP resolveRows(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
    {"binByMZ",                (DL_FUNC) &binByMZ,                 3},
    {"binDuplicates",         (DL_FUNC) &binDuplicates,            2},
    {"findDuplicates",         (DL_FUNC) &findDuplicates,          6},
    {"groupDuplicates",         (DL_FUNC) &groupDuplicates,        5},
    {"labelRows",              (DL_FUNC) &labelRows,              19},
    {"merge_duplicate_values",  (DL_FUNC) &merge_duplicate_values, 4},
    {"selectAnchorsByID",      (DL_FUNC) &selectAnchorsByID,       6},
    {"selectIterativeAnchors", (DL_FUNC) &selectIterativeAnchors,  7},
    {"write2file",             (DL_FUNC) &write2file,    4},
    {"resolveRows",             (DL_FUNC) &resolveRows,  9},
    {NULL, NULL, 0}
};

void R_init_metabCombiner(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
