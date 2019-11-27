#include <string.h>
#include <stdlib.h>
#include <registryFunction.h>
#include <aSubRecord.h>
#include <menuFtype.h>
#include <errlog.h>
#include <epicsString.h>
#include <epicsExport.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

#include "zfmagfld.h"

long matrix_multiply_impl(aSubRecord *prec) 
{
    /*
     * Validate input data types
     */

    if (prec->fta != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect input argument type A", prec->name);
		return 1;
	}

    if (prec->ftb != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect input argument type B", prec->name);
		return 1;
	}

    if (prec->ftc != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect input argument type C", prec->name);
		return 1;
	}

    if (prec->ftd != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect input argument type D", prec->name);
		return 1;
	}

    if (prec->fte != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect input argument type E", prec->name);
		return 1;
	}

    if (prec->ftf != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect input argument type F", prec->name);
		return 1;
	}

    if (prec->ftg != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect input argument type G", prec->name);
		return 1;
	}

    if (prec->fth != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect input argument type H", prec->name);
		return 1;
	}

    if (prec->fti != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect input argument type I", prec->name);
		return 1;
	}

    if (prec->ftj != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect input argument type J", prec->name);
		return 1;
	}

    if (prec->ftk != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect input argument type K", prec->name);
		return 1;
	}

    if (prec->ftl != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect input argument type L", prec->name);
		return 1;
	}

    gsl_vector* data_vector = gsl_vector_alloc(3);

    gsl_vector_set(data_vector, 0, *(epicsFloat64*)prec->a);
    gsl_vector_set(data_vector, 1, *(epicsFloat64*)prec->b);
    gsl_vector_set(data_vector, 2, *(epicsFloat64*)prec->c);

    gsl_matrix *sensor_matrix = gsl_matrix_alloc(3, 3);

    gsl_matrix_set(sensor_matrix, 0, 0, *(epicsFloat64*)prec->d);
    gsl_matrix_set(sensor_matrix, 0, 1, *(epicsFloat64*)prec->e);
    gsl_matrix_set(sensor_matrix, 0, 2, *(epicsFloat64*)prec->f);
    gsl_matrix_set(sensor_matrix, 1, 0, *(epicsFloat64*)prec->g);
    gsl_matrix_set(sensor_matrix, 1, 1, *(epicsFloat64*)prec->h);
    gsl_matrix_set(sensor_matrix, 1, 2, *(epicsFloat64*)prec->i);
    gsl_matrix_set(sensor_matrix, 2, 0, *(epicsFloat64*)prec->j);
    gsl_matrix_set(sensor_matrix, 2, 1, *(epicsFloat64*)prec->k);
    gsl_matrix_set(sensor_matrix, 2, 2, *(epicsFloat64*)prec->l);

    gsl_vector* output_vector = gsl_vector_calloc(3); 

    gsl_blas_dgemv(CblasTrans, 1.0, sensor_matrix, data_vector, 0.0, output_vector);
    
    gsl_vector_free(data_vector);
    gsl_matrix_free(sensor_matrix);

    *(epicsFloat64*)prec->vala = gsl_vector_get(output_vector, 0);
    *(epicsFloat64*)prec->valb = gsl_vector_get(output_vector, 1);
    *(epicsFloat64*)prec->valc = gsl_vector_get(output_vector, 2);

    gsl_vector_free(output_vector);
    return 0; /* process output links */
}
