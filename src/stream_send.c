/*------------------------------------------------------------------------------
* stream_send.c : stream(send) functions
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

#include "rtklib.h"

/* send nmea request -----------------------------------------------------------
* send nmea gpgga message to stream
* args   : stream_t *stream I   stream
*          double *pos      I   position {x,y,z} (ecef) (m)
* return : none
*-----------------------------------------------------------------------------*/
extern void strsendnmea(stream_t *stream, const double *pos)
{
    sol_t sol={{0}};
    unsigned char buff[1024];
    int i,n;
    
    tracet(3,"strsendnmea: pos=%.3f %.3f %.3f\n",pos[0],pos[1],pos[2]);
    
    sol.stat=SOLQ_SINGLE;
    sol.time=utc2gpst(timeget());
    for (i=0;i<3;i++) sol.rr[i]=pos[i];
	n=outnmea_gga(buff,&sol);
    strwrite(stream,buff,n);
}
/* send receiver command -------------------------------------------------------
* send receiver commands to stream
* args   : stream_t *stream I   stream
*          char   *cmd      I   receiver command strings
* return : none
*-----------------------------------------------------------------------------*/
extern void strsendcmd(stream_t *str, const char *cmd)
{
    unsigned char buff[1024];
    const char *p=cmd,*q;
    char msg[1024],cmdend[]="\r\n";
    int n,m,ms;
    
    tracet(3,"strsendcmd: cmd=%s\n",cmd);
    
    for (;;) {
        for (q=p;;q++) if (*q=='\r'||*q=='\n'||*q=='\0') break;
        n=(int)(q-p); strncpy(msg,p,n); msg[n]='\0';
        
        if (!*msg||*msg=='#') { /* null or comment */
            ;
        }
        else if (*msg=='!') { /* binary escape */
            
            if (!strncmp(msg+1,"WAIT",4)) { /* wait */
                if (sscanf(msg+5,"%d",&ms)<1) ms=100;
                if (ms>3000) ms=3000; /* max 3 s */
                sleepms(ms);
            }
            else if (!strncmp(msg+1,"UBX",3)) { /* ublox */
                if ((m=gen_ubx(msg+4,buff))>0) strwrite(str,buff,m);
            }
            else if (!strncmp(msg+1,"STQ",3)) { /* skytraq */
                if ((m=gen_stq(msg+4,buff))>0) strwrite(str,buff,m);
            }
            else if (!strncmp(msg+1,"NVS",3)) { /* nvs */
                if ((m=gen_nvs(msg+4,buff))>0) strwrite(str,buff,m);
            }
            else if (!strncmp(msg+1,"LEXR",3)) { /* lex receiver */
                if ((m=gen_lexr(msg+5,buff))>0) strwrite(str,buff,m);
            }
        }
        else {
            strwrite(str,(unsigned char *)msg,n);
            strwrite(str,(unsigned char *)cmdend,2);
        }
        if (*q=='\0') break; else p=q+1;
    }
}
