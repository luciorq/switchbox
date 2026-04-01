/* Registration of C routines */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void CalculateSignedScoreCore(int *situation, int *rowLenPtr, double *data1, int *colLen1Ptr, double *data2, int *colLen2Ptr,
                              double *output1, double *output3, double *output4, double *output5, double *output6, double *output7, double *output8);

void CalculateSignedScoreCoreTieHandler(int *situation, int *rowLenPtr, double *data1, int *colLen1Ptr, double *data2, int *colLen2Ptr,
                                        double *output1, double *output2, double *output3, double *output4);

void CalculateSignedScoreRestrictedPairsCore(int *situation, int *rowLenPtr, double *data1, int *colLen1Ptr, double *data2,
                                             int *colLen2Ptr, int *edges1, int *edges2, int *nopairsPtr,
                                             double *output1, double *output2, double *output3);

void CalculateSignedScoreRestrictedPairsCoreTieHandler(int *situation, int *rowLenPtr, double *data1, int *colLen1Ptr, double *data2,
                                                       int *colLen2Ptr, int *edges1, int *edges2, int *nopairsPtr,
                                                       double *output1, double *output2, double *output3, double *output4);


#if _MSC_VER >= 1000
__declspec(dllexport)
#endif

static const R_CMethodDef cMethods[] = {
    {"CalculateSignedScoreCore", (DL_FUNC)&CalculateSignedScoreCore, 13},
    {"CalculateSignedScoreCoreTieHandler", (DL_FUNC)&CalculateSignedScoreCoreTieHandler, 10},
    {"CalculateSignedScoreRestrictedPairsCore", (DL_FUNC)&CalculateSignedScoreRestrictedPairsCore, 12},
    {"CalculateSignedScoreRestrictedPairsCoreTieHandler", (DL_FUNC)&CalculateSignedScoreRestrictedPairsCoreTieHandler, 13},
    {NULL, NULL, 0},
};

void R_init_switchbox(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
