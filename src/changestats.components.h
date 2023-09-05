/*
######################################################################
#
# changestats.components.h
#
# copyright (c) 2022, Carter T. Butts <buttsc@uci.edu>
# Last Modified 8/02/23
# Licensed under the GNU General Public License v3 or later
#
# Portions based on the statnet ergm.userterms package (Hunter et al.,
# 2013, JSS).
#
# Part of the R/ergm.components package
#
# This file contains headers related to the changescores for
# component and connectivity-related ergm terms.
#
######################################################################
*/

#ifndef CHANGESTATS_H
#define CHANGESTATS_H

/*IMPORTED HEADERS AND MACROS------------------------------------------------*/

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include <R.h>
/*#include <R_ext/Lapack.h>*/

#define DYADMATCH(i,j,k,l) ((MIN(i,j)==MIN(k,l))&&(MAX(i,j)==MAX(k,l)))
#define DDYADMATCH(i,j,k,l) (i==k)&&(j==l)


/*MEMORY STRUCTURES----------------------------------------------------------*/

/*Element for stacks 'n queues*/
typedef struct elementtype{
   Vertex v;
   struct elementtype *next,*prev;
} element;


/*UTILITY FUNCTION HEADERS---------------------------------------------------*/

void lBCDFS(Vertex v, Vertex root, Network *nwp, Vertex tail, Vertex head, int thpresent, int dyadic, Vertex *iter, char *visited, Vertex *miniter, Vertex *viter, double *bc);

double localBridgeCnt(Network *nwp, Vertex tail, Vertex head, int thpresent, int dyadic);


/*CHANGESTAT FUNCTION HEADERS------------------------------------------------*/

CHANGESTAT_FN(d_dimers_udir);

CHANGESTAT_FN(d_dimers_dir);

CHANGESTAT_FN(d_compsizesum);

CHANGESTAT_FN(d_components);

CHANGESTAT_FN(d_bridges);

#endif
