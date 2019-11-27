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

static gsl_matrix* populate_sensor_matrix(const double elem11, const double elem12, const double elem13, const double elem21, const double elem22, const double elem23, const double elem31, const double elem32, const double elem33)
{
    gsl_matrix *sensor_matrix = gsl_matrix_alloc(3, 3);

    gsl_matrix_set(sensor_matrix, 0, 0, elem11);
    gsl_matrix_set(sensor_matrix, 0, 1, elem12);
    gsl_matrix_set(sensor_matrix, 0, 2, elem13);
    gsl_matrix_set(sensor_matrix, 1, 0, elem21);
    gsl_matrix_set(sensor_matrix, 1, 1, elem22);
    gsl_matrix_set(sensor_matrix, 1, 2, elem23);
    gsl_matrix_set(sensor_matrix, 2, 0, elem31);
    gsl_matrix_set(sensor_matrix, 2, 1, elem32);
    gsl_matrix_set(sensor_matrix, 2, 2, elem33);

    return sensor_matrix;
}

static gsl_vector* populate_data_vector(const double readingX, const double readingY, const double readingZ)
{
    gsl_vector *data_vector = gsl_vector_alloc(3);

    gsl_vector_set(data_vector, 0, readingX);
    gsl_vector_set(data_vector, 1, readingY);
    gsl_vector_set(data_vector, 2, readingZ);

    return data_vector;
}

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

    //gsl_matrix data_vector = populate_data_vector((const double)*prec->fta, (const double)*prec->ftb, (const double)*prec->ftc);

    //double 1.0;

    //gsl_matrix sensor_matrix = populate_sensor_matrix((const double)*prec->ftd, (const double)*prec->fte, (const double)*prec->ftf, (const double)*prec->ftg, (const double)*prec->fth, (const double)*prec->fti, (const double)*prec->ftj, (const double)*prec->ftk, (const double)*prec->ftl);
    gsl_vector *data_vector = gsl_vector_alloc(3);

    gsl_vector_set(data_vector, 0, *(double*)prec->fta);
    gsl_vector_set(data_vector, 1, *(double*)prec->ftb);
    gsl_vector_set(data_vector, 2, *(double*)prec->ftc);

    gsl_matrix *sensor_matrix = gsl_matrix_alloc(3, 3);

    gsl_matrix_set(sensor_matrix, 0, 0, *(double*)prec->ftd);
    gsl_matrix_set(sensor_matrix, 0, 1, *(double*)prec->fte);
    gsl_matrix_set(sensor_matrix, 0, 2, *(double*)prec->ftf);
    gsl_matrix_set(sensor_matrix, 1, 0, *(double*)prec->ftg);
    gsl_matrix_set(sensor_matrix, 1, 1, *(double*)prec->fth);
    gsl_matrix_set(sensor_matrix, 1, 2, *(double*)prec->fti);
    gsl_matrix_set(sensor_matrix, 2, 0, *(double*)prec->ftj);
    gsl_matrix_set(sensor_matrix, 2, 1, *(double*)prec->ftk);
    gsl_matrix_set(sensor_matrix, 2, 2, *(double*)prec->ftl);

    gsl_vector* empty_vector = gsl_vector_calloc(3); 

    gsl_blas_dgemv(CblasTrans, 1.0, sensor_matrix, data_vector, 0.0, empty_vector);

    gsl_vector_free(empty_vector);
    gsl_vector_free(data_vector);
    gsl_matrix_free(sensor_matrix);


    *(epicsInt32*)prec->vala = *(epicsInt32*)prec->a;
    return 0; /* process output links */
}