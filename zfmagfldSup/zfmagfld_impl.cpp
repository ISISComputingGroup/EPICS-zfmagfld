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

/**
 * Deallocates the arrays used in to do the matrix multiplication
 */
void free_array_data(gsl_vector* data_vector, gsl_matrix* sensor_matrix, gsl_vector* field_vector){
    gsl_vector_free(data_vector);
    gsl_matrix_free(sensor_matrix);
    gsl_vector_free(field_vector);
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

    if (prec->ftva != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect output argument type A", prec->name);
		return 1;
	}
    
    if (prec->ftvb != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect output argument type B", prec->name);
		return 1;
	}

    if (prec->ftvc != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect output argument type C", prec->name);
		return 1;
	}

    // Allocate vector and assign values from inputs
    gsl_vector* data_vector = gsl_vector_alloc(3);

    gsl_vector_set(data_vector, 0, *(epicsFloat64*)prec->a);
    gsl_vector_set(data_vector, 1, *(epicsFloat64*)prec->b);
    gsl_vector_set(data_vector, 2, *(epicsFloat64*)prec->c);

    // Allocate matrix and assign values from inputs
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

    gsl_vector* field_vector = gsl_vector_calloc(3); 

    /*
     * dgemv is a BLAS standard function to multiply a matrix and a vector, implemented in gsl.
     * We need to multiply the static field matrix by the measured field with the offset subtracted:
     * field_vector = sensor_matrix * data_vector
     * 
     * CblasNoTrans means 'do not transpose the static field matrix'. It is provided by gsl.
     * 
     * 1.0 and 0.0 are scaling factors, used here to keep or remove terms in the equation dgemv solves.
     * 
     */
    gsl_blas_dgemv(CblasNoTrans, 1.0, sensor_matrix, data_vector, 0.0, field_vector);
    
    double fieldStrength = 0.0;
    double* field_strength = &fieldStrength;

    /* 
     * Here we are calculating the magnitude of the field strength vector:
     *    field strength^2 = field_x^2 + field_y^2 + field_z^2
     * 
     * The BLAS function ddot is used to calculate the dot product 
     * of the field vector with itself.
     */
    gsl_blas_ddot(field_vector, field_vector, field_strength);

    *(epicsFloat64*)prec->vala = gsl_vector_get(field_vector, 0); // X component
    *(epicsFloat64*)prec->valb = gsl_vector_get(field_vector, 1); // Y component
    *(epicsFloat64*)prec->valc = gsl_vector_get(field_vector, 2); // Z component
    *(epicsFloat64*)prec->vald = sqrt(*field_strength); // Magnitude of field vector

    free_array_data(data_vector, sensor_matrix, field_vector);
    return 0; /* process output links */
}
