/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file auxcompute_z.h
 *
 * This file contains the declarations of computational auxiliary functions.
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/
#ifndef __AUXCOMPUTE_Z__
#define __AUXCOMPUTE_Z__

void HICMA_znormest( int M, int N, double _Complex *A, double *e, double _Complex *work);
#endif
