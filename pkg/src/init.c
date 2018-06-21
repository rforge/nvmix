/* Register routines with R ***************************************************/

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "evalf.h"
#include "evalfbonly.h"

static const R_CallMethodDef callMethods[] = {
	{"evalf_",      (DL_FUNC) &evalf_, 9},
	{"evalfbonly_", (DL_FUNC) &evalfbonly_, 8},
	{NULL, NULL, 0}
};

void R_init_student(DllInfo *dll)
{
	R_useDynamicSymbols(dll, FALSE);
	R_registerRoutines(dll, NULL, callMethods, NULL, NULL); /* s. WRE (2015, Section 5.4) */
}
