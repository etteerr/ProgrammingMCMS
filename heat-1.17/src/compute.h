#ifndef COMPUTE_H
#define COMPUTE_H

#ifndef HM_VERBOSE
#define HM_VERBOSE 1
#endif

#include "input.h"
#include "output.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void do_compute(const struct parameters *p,
                struct results *r);

/**
 * hm_init_map
 *  returns 0 on succes
 */
int hm_init_map(struct parameters *p);

#endif
