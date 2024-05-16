/**
 *
 * @file sumsq_update.h
 *
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CPU auxiliary sumsq_update routine
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2020-03-03
 *
 */
#ifndef _sumsq_update_h_
#define _sumsq_update_h_

/**
 *******************************************************************************
 *
 * @brief Update the couple (scale, sumsq) with one element when
 * computing the Froebnius norm.
 *
 * The frobenius norm is equal to scale * sqrt( sumsq ), this method allows to
 * avoid overflow in the sum square computation.
 *
 *******************************************************************************
 *
 * @param[inout] scale
 *           On entry, the former scale
 *           On exit, the update scale to take into account the value
 *
 * @param[inout] sumsq
 *           On entry, the former sumsq
 *           On exit, the update sumsq to take into account the value
 *
 * @param[in] value
 *          The value to integrate into the couple (scale, sumsq)
 *
 *******************************************************************************/
static inline void
#if defined(PRECISION_d) || defined(PRECISION_z)
sumsq_update( int nb, double *scale, double *sumsq, const double *value )
{
    double absval = fabs(*value);
    double ratio;
    if ( absval != 0. ){
        if ( (*scale) < absval ) {
            ratio = (*scale) / absval;
            *sumsq = (double)nb + (*sumsq) * ratio * ratio;
            *scale = absval;
        } else {
            ratio = absval / (*scale);
            *sumsq = (*sumsq) + (double)nb * ratio * ratio;
        }
    }
}
#elif defined(PRECISION_s) || defined(PRECISION_c)
sumsq_update( int nb, float *scale, float *sumsq, const float *value )
{
    float absval = fabs(*value);
    float ratio;
    if ( absval != 0. ){
        if ( (*scale) < absval ) {
            ratio = (*scale) / absval;
            *sumsq = (float)nb + (*sumsq) * ratio * ratio;
            *scale = absval;
        } else {
            ratio = absval / (*scale);
            *sumsq = (*sumsq) + (float)nb * ratio * ratio;
        }
    }
}
#endif

/**
 *******************************************************************************
 *
 * @brief Update the couple (scale, sumsq) by adding another couple when
 * computing the Froebenius norm.
 *
 * The frobenius norm is equal to scale * sqrt( sumsq ), this method allows to
 * avoid overflow in the sum square computation.
 *
 *******************************************************************************
 *
 * @param[inout] scale
 *           On entry, the former scale
 *           On exit, the update scale to take into account the value
 *
 * @param[inout] sumsq
 *           On entry, the former sumsq
 *           On exit, the update sumsq to take into account the value
 *
 * @param[in] value
 *          The value to integrate into the couple (scale, sumsq)
 *
 *******************************************************************************/
static inline void
#if defined(PRECISION_d) || defined(PRECISION_z)
sumsq_update_2( const double *scalein, const double *sumsqin, double *scaleout, double *sumsqout )
{
    if (*scaleout >= 0.) {
        if ( (*scaleout) < (*scalein) ) {
            *sumsqout = *sumsqin  + (*sumsqout * (( *scaleout / *scalein ) * ( *scaleout / *scalein )));
            *scaleout = *scalein;
        } else {
            if ( (*scaleout) > 0. ){
                *sumsqout = *sumsqout + (*sumsqin * (( *scalein / *scaleout ) * ( *scalein / *scaleout )));
            }
        }
    }
}
#elif defined(PRECISION_s) || defined(PRECISION_c)
sumsq_update_2( const float *scalein, const float *sumsqin, float *scaleout, float *sumsqout )
{
    if (*scaleout >= 0.) {
        if ( (*scaleout) < (*scalein) ) {
            *sumsqout = *sumsqin  + (*sumsqout * (( *scaleout / *scalein ) * ( *scaleout / *scalein )));
            *scaleout = *scalein;
        } else {
            if ( (*scaleout) > 0. ){
                *sumsqout = *sumsqout + (*sumsqin * (( *scalein / *scaleout ) * ( *scalein / *scaleout )));
            }
        }
    }
}
#endif

#endif /* _sumsq_update_h_ */
