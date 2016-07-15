/*------------------------------------------------------------------------------
* solution_cmn.c : solution common functions
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
*     [1] National Marine Electronic Association and International Marine
*         Electronics Association, NMEA 0183 version 4.10, August 1, 2012
*
* history : 2015/01/08  1.0  new
*-----------------------------------------------------------------------------*/

#include <ctype.h>
#include "rtklib.h"
#include "solution.h"

extern const char rcsid_sol[]="$Id: solution.c,v 1.1 2008/07/17 21:48:06 ttaka Exp $";

extern const int solq_nmea[]={  /* nmea quality flags to rtklib sol quality */
	/* nmea 0183 v.2.3 quality flags: */
    /*  0=invalid, 1=gps fix (sps), 2=dgps fix, 3=pps fix, 4=rtk, 5=float rtk */
    /*  6=estimated (dead reckoning), 7=manual input, 8=simulation */
    
    SOLQ_NONE ,SOLQ_SINGLE, SOLQ_DGPS, SOLQ_PPP , SOLQ_FIX,
    SOLQ_FLOAT,SOLQ_DR    , SOLQ_NONE, SOLQ_NONE, SOLQ_NONE
};

/* solution option to field separator ----------------------------------------*/
extern const char *opt2sep(const solopt_t *opt)
{
	if (!*opt->sep) return " ";
	else if (!strcmp(opt->sep,"\\t")) return "\t";
	return opt->sep;
}

