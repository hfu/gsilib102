/*------------------------------------------------------------------------------
* solution.h : solution header
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

#ifndef SOLUTION_H
#define SOLUTION_H

#ifdef __cplusplus
extern "C" {
#endif

/* constants and macros ------------------------------------------------------*/
#define MAXFIELD   64           /* max number of fields in a record */
#define MAXNMEA    256          /* max length of nmea sentence */

#define KNOT2M     0.514444444  /* m/knot */

#define SOL1 "backward"
#define SOL0 "forward"
#define SOL2 "combined"

extern const char rcsid[];

extern const int solq_nmea[];

/* solution option to field separator ----------------------------------------*/
extern const char *opt2sep(const solopt_t *opt);

#ifdef __cplusplus
}
#endif
#endif /* SOLUTION_H */
