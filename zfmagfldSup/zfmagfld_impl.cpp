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
#include <gsl/gsl_errno.h>

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

#include "zfmagfld.h"

/**
 * Multiplies the scaled, offset magnetometer data with the fixed sensor matrix.
 * Inputs:
 * (A, B, C) are the components of the measured, scaled field. Double precision
 * D through L are the values of the 3x3 fixed sensor matrix, in row-major order. Double precision
 * 
 * Outputs:
 * (VALA, VALB, VALC) are the components of the corrected field. Double precision
 * VALD is the magnitute of the corrected field strength. Double precision
 */
long matrix_multiply_impl(aSubRecord *prec) 
{
    // Disable default error handler for GSL as it terminates the process. Check return values instead.
    gsl_set_error_handler_off();

    // Validate input and output data types
    if (prec->fta != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type input A", prec->name);
		return 1;
	}

    if (prec->ftb != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type input B", prec->name);
		return 1;
	}

    if (prec->ftc != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type input C", prec->name);
		return 1;
	}

    if (prec->ftd != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type input D", prec->name);
		return 1;
	}

    if (prec->fte != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type input E", prec->name);
		return 1;
	}

    if (prec->ftf != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type input F", prec->name);
		return 1;
	}

    if (prec->ftg != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type input G", prec->name);
		return 1;
	}

    if (prec->fth != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type input H", prec->name);
		return 1;
	}

    if (prec->fti != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type input I", prec->name);
		return 1;
	}

    if (prec->ftj != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type input J", prec->name);
		return 1;
	}

    if (prec->ftk != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type input K", prec->name);
		return 1;
	}

    if (prec->ftl != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type input L", prec->name);
		return 1;
	}

    if (prec->ftva != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type output VALA", prec->name);
		return 1;
	}
    
    if (prec->ftvb != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type output VALB", prec->name);
		return 1;
	}

    if (prec->ftvc != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type output VALC", prec->name);
		return 1;
	}

    if (prec->ftvc != menuFtypeDOUBLE)
	{
        errlogSevPrintf(errlogMajor, "%s incorrect argument type output VALD", prec->name);
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

    // Allocate 3-element field vector, initialise to 0
    gsl_vector* field_vector = gsl_vector_calloc(3);

    // Allocate 1-element field strength variable, initialise to 0
    double fieldStrength = 0.0;
    double* field_strength = &fieldStrength;

    int status = 0;

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
    status = gsl_blas_dgemv(CblasNoTrans, 1.0, sensor_matrix, data_vector, 0.0, field_vector);

    if (status) {
        // Matrix multiplication failed
        errlogSevPrintf(errlogMajor, "GSL matrix multiplication (DGEMV) error %d, %s", status, gsl_strerror(status));
    } 

    if (!status) {
        /* 
         * DGEMV worked, calculate the magnitude of the field strength vector:
         *    field strength^2 = field_x^2 + field_y^2 + field_z^2
         * 
         * The BLAS function ddot is used to calculate the dot product 
         * of the field vector with itself.
         */
        status = gsl_blas_ddot(field_vector, field_vector, field_strength);

        if (status){
            // Dot product failed
            errlogSevPrintf(errlogMajor, "GSL dot product (DDOT) error %d, %s", status, gsl_strerror(status));
        }
    }

    *(epicsFloat64*)prec->vala = gsl_vector_get(field_vector, 0); // X component
    *(epicsFloat64*)prec->valb = gsl_vector_get(field_vector, 1); // Y component
    *(epicsFloat64*)prec->valc = gsl_vector_get(field_vector, 2); // Z component
    *(epicsFloat64*)prec->vald = sqrt(*field_strength); // Magnitude of field vector

    // Deallocate arrays
    gsl_vector_free(data_vector);
    gsl_matrix_free(sensor_matrix);
    gsl_vector_free(field_vector);

    return status; /* process output links or return error */
}
