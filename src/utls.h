/*
 * =====================================================================================
 *
 *       Filename:  utls.h
 *
 *    Description:  pins utilities
 *
 *        Version:  1.0
 *        Created:  19/09/19 16:23:58
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#ifndef __UTLS_H
#define __UTLS_H

#include <math.h>
static float norm_cdf(int x, float p, int n);

static float norm_cdf(int x, float p, int n) {
    float mean = n * p;
    float sd = sqrt(n * p * (1 - p));
    return 0.5 * (1 + erf((x - mean)/(sd * sqrt(2))));
}

#endif
