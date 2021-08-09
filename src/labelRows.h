#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

int detectGroup(int cursor, int* group, int* start, int *end, int length);

int balancedGroups(SEXP labels, int *start, int *end, int* rankX, int* rankY);

int filterScoreRankErr(SEXP labels, int *start, int *end, double* score, int* rankX,
						int* rankY, double minScore, int maxRankX, int maxRankY,
						double maxRTerr, double* rterr);
						
						
void detect_con_score(SEXP labels, int* sub, int* alt, int ri, int rj, double* delta,
                double* mz, double* rt, double* score, int* head, int type, int* max);

void detect_con_mzrt(SEXP labels, int* sub, int* alt, int ri, int rj, double* delta,
                	double* mz, double* rt, double* score, int* head, int type, int* max);
                	
void findCons(SEXP labels, int* sub, int* alt, int* max, int *start,
              int *end, double* delta, double *mzx, double *mzy,
              double *rtx, double *rty, double* score, int* head,
			  void (*detect)(SEXP, int*, int*, int, int, double*, double*, double*,
			                 double*, int*, int, int*));
			                 
int sharedFeatCheck(SEXP labels, int ri, int rj, double* mzx, double *mzy, double *rtx, 
                	double *rty);
                	
int retOrderCheck(int ri, int rj, double *rtx, double *rty);

void form_rowset(SEXP labels, int* rowset, int* rowsetnum, int end, double* mzx, 
				double *mzy, double *rtx, double *rty, int retOrder);
				
double score_rowset(int* rowset, double* score, int max, double *solveScore);

void update_rowset_labels(SEXP labels, int* rowset, int nrow);

void resolve(SEXP labels, int* rowsetnum, int start, int end, double* score, 
               double* mzx, double *mzy, double *rtx, double *rty, int retOrderOpt,
               double* solveScore);
               
          	


