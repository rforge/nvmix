/* Register routines with R ***************************************************/

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "eval_int_t.h"
#include "eval_int_normal.h"
#include "eval_int_mix.h"

static const R_CallMethodDef callMethods[] = {
	{"eval_int_t_", (DL_FUNC) &eval_int_t_, 9},
    {"eval_int_normal_", (DL_FUNC) &eval_int_normal_, 8},
    {"eval_int_mix_", (DL_FUNC) &eval_int_normal_, 8},
	{NULL, NULL, 0}
};

void R_init_nvmix(DllInfo *dll)
{
	R_useDynamicSymbols(dll, FALSE);
	R_registerRoutines(dll, NULL, callMethods, NULL, NULL); /* s. WRE (2015, Section 5.4) */
}
