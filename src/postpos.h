/*------------------------------------------------------------------------------
* postpos.h : post-processing positioning header
*
*    Copyright (C) 2014 by Geospatial Information Authority of Japan,
*    All rights reserved.
*    
*    Released under the BSD, and GPL Licenses.
*
*
*  Original software: RTKLIB ver.2.4.2 p4
*
*    Copyright (C) 2007-2013 by T.Takasu, All rights reserved.
*
*
* references :
*
* history : 2015/01/08  1.0  new
*-----------------------------------------------------------------------------*/

#ifndef POSTPOS_H
#define POSTPOS_H
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#ifdef WIN32
#include <winsock2.h>
#include <windows.h>
#else
#include <pthread.h>
#endif
#ifdef __cplusplus
extern "C" {
#endif

#include "rtklib.h"

#define MAXPRCDAYS  100          /* max days of continuous processing */

extern pcvs_t pcvss;        /* receiver antenna parameters */
extern pcvs_t pcvsr;        /* satellite antenna parameters */
extern obs_t obss;          /* observation data */
extern nav_t navs;          /* navigation data */
extern  sbs_t sbss;          /* sbas messages */
extern lex_t lexs;          /* lex messages */
extern int nepoch;            /* number of observation epochs */
extern int iobsu;            /* current rover observation data index */
extern int iobsr;            /* current reference observation data index */
extern int isbs;            /* current sbas message index */
extern int ilex;            /* current lex message index */
extern int aborts;            /* abort status */
extern char proc_rov [64];   /* rover for current processing */
extern char proc_base[64];   /* base station for current processing */
extern char rtcm_file[1024]; /* rtcm data file */
extern char rtcm_path[1024]; /* rtcm data path */
extern rtcm_t rtcm;             /* rtcm control struct */
extern FILE *fp_rtcm;      /* rtcm data file pointer */
extern double diff_se;
extern int posepoch;            /* Poss. # of obs epochs */
extern gtime_t obsts;					/* obs start time */
extern gtime_t obste;					/* obs end time */

//static gtime_t tstart_={0};         // time start for progress-bar
//static gtime_t tend_  ={0};         // time end for progress-bar

//extern int openses(const prcopt_t *popt, const solopt_t *sopt,
//                   const filopt_t *fopt, nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr);
//extern void closeses(nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr);

extern int nextobsf(const obs_t *obs, int *i, int rcv);
extern int nextobsb(const obs_t *obs, int *i, int rcv);
extern int inputobs(obsd_t *obs, int solq, const prcopt_t *popt);
extern int readobsnav(gtime_t ts, gtime_t te, double ti, char **infile,
                      const int *index, int n, const prcopt_t *prcopt,
                      obs_t *obs, nav_t *nav, sta_t *sta);
extern void freeobsnav(obs_t *obs, nav_t *nav);
extern int checkbrk(const char *format, ...);
extern FILE *openfile(const char *outfile);


#ifdef __cplusplus
}
#endif
#endif /* RTKLIB_H */
