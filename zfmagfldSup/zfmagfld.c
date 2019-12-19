#include <stdlib.h>
#include <registryFunction.h>
#include <aSubRecord.h>
#include <epicsExport.h>

#include "zfmagfld.h"

static long matrix_multiply(aSubRecord *prec)
{
    return matrix_multiply_impl(prec);
}

epicsRegisterFunction(matrix_multiply);
