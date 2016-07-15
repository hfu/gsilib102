/*------------------------------------------------------------------------------
* solution_out.c : functions to output solution
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

/* sqrt of covariance --------------------------------------------------------*/
static double sqvar(double covar)
{
	return covar<0.0?-sqrt(-covar):sqrt(covar);
}
/* solution to covariance ----------------------------------------------------*/
extern void soltocov(const sol_t *sol, double *P)
{
    P[0]     =sol->qr[0]; /* xx or ee */
    P[4]     =sol->qr[1]; /* yy or nn */
    P[8]     =sol->qr[2]; /* zz or uu */
    P[1]=P[3]=sol->qr[3]; /* xy or en */
    P[5]=P[7]=sol->qr[4]; /* yz or nu */
    P[2]=P[6]=sol->qr[5]; /* zx or ue */
}
/* output solution as the form of x/y/z-ecef ---------------------------------*/
static int outecef(unsigned char *buff, const char *s, const sol_t *sol,
				   const solopt_t *opt)
{
	const char *sep=opt2sep(opt);
	char *p=(char *)buff;

	trace(3,"outecef:\n");

	p+=sprintf(p,"%s%s%14.4f%s%14.4f%s%14.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s%6.1f",
			   s,sep,sol->rr[0],sep,sol->rr[1],sep,sol->rr[2],sep,sol->stat,sep,
			   sol->ns,sep,SQRT(sol->qr[0]),sep,SQRT(sol->qr[1]),sep,SQRT(sol->qr[2]),
			   sep,sqvar(sol->qr[3]),sep,sqvar(sol->qr[4]),sep,sqvar(sol->qr[5]),
			   sep,sol->age,sep,sol->ratio);
	return p-(char *)buff;
}
/* output solution as the form of x/y/z-ecef ---------------------------------*/
static int outrva(unsigned char *buff, const char *s, const sol_t *sol,
				   const solopt_t *opt)
{
	const char *sep=opt2sep(opt);
	char *p=(char *)buff;

	trace(3,"outrva:\n");

	p+=sprintf(p,"%s%s%14.4f%s%14.4f%s%14.4f%s%9.4f%s%9.4f%s%9.4f%s%8.4f%s%8.4f%s%8.4f%s%3d%s%3d%s"
	"%11.4f%s%11.4f%s%11.4f%s%11.4f%s%11.4f%s%11.4f%s"
	"%11.4f%s%11.4f%s%11.4f%s%11.4f%s%11.4f%s%11.4f%s"
	"%11.4f%s%11.4f%s%11.4f%s%11.4f%s%11.4f%s%11.4f%s"
	"%6.2f%s%6.1f",
			   s,sep,sol->rr[0],sep,sol->rr[1],sep,sol->rr[2],sep,sol->rr[3],sep,sol->rr[4],sep,sol->rr[5],sep,sol->rr[6],sep,sol->rr[7],sep,sol->rr[8],sep,
			   sol->stat,sep,sol->ns,sep,
			   SQRT(sol->qr[ 0]),sep,SQRT(sol->qr[ 1]),sep,SQRT(sol->qr[ 2]),sep,sqvar(sol->qr[ 3]),sep,sqvar(sol->qr[ 4]),sep,sqvar(sol->qr[ 5]),sep,
			   SQRT(sol->qr[ 6]),sep,SQRT(sol->qr[ 7]),sep,SQRT(sol->qr[ 8]),sep,sqvar(sol->qr[ 9]),sep,sqvar(sol->qr[10]),sep,sqvar(sol->qr[11]),sep,
			   SQRT(sol->qr[12]),sep,SQRT(sol->qr[13]),sep,SQRT(sol->qr[14]),sep,sqvar(sol->qr[15]),sep,sqvar(sol->qr[16]),sep,sqvar(sol->qr[17]),sep,
			   sol->age,sep,sol->ratio);
	return p-(char *)buff;
}
/* output solution as the form of lat/lon/height -----------------------------*/
static int outpos(unsigned char *buff, const char *s, const sol_t *sol,
				  const solopt_t *opt)
{
	double pos[3],dms1[3],dms2[3],P[9],Q[9];
	const char *sep=opt2sep(opt);
	char *p=(char *)buff;

	trace(3,"outpos  :\n");

	ecef2pos(sol->rr,pos);
	soltocov(sol,P);
	covenu(pos,P,Q);
	if (opt->height==1) { /* geodetic height */
		pos[2]-=geoidh(pos);
	}
	if (opt->degf) {
		deg2dms(pos[0]*R2D,dms1);
		deg2dms(pos[1]*R2D,dms2);
		p+=sprintf(p,"%s%s%4.0f%s%02.0f%s%08.5f%s%4.0f%s%02.0f%s%08.5f",s,sep,
				   dms1[0],sep,dms1[1],sep,dms1[2],sep,dms2[0],sep,dms2[1],sep,
				   dms2[2]);
	}
	else p+=sprintf(p,"%s%s%14.9f%s%14.9f",s,sep,pos[0]*R2D,sep,pos[1]*R2D);
	p+=sprintf(p,"%s%10.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s%6.1f",
			   sep,pos[2],sep,sol->stat,sep,sol->ns,sep,SQRT(Q[4]),sep,
			   SQRT(Q[0]),sep,SQRT(Q[8]),sep,sqvar(Q[1]),sep,sqvar(Q[2]),
			   sep,sqvar(Q[5]),sep,sol->age,sep,sol->ratio);
	return p-(char *)buff;
}
/* output solution as the form of e/n/u-baseline -----------------------------*/
static int outenu(unsigned char *buff, const char *s, const sol_t *sol,
				  const double *rb, const solopt_t *opt)
{
	double pos[3],rr[3],enu[3],P[9],Q[9];
	int i;
	const char *sep=opt2sep(opt);
	char *p=(char *)buff;

	trace(3,"outenu  :\n");

	for (i=0;i<3;i++) rr[i]=sol->rr[i]-rb[i];
	ecef2pos(rb,pos);
	soltocov(sol,P);
	covenu(pos,P,Q);
	ecef2enu(pos,rr,enu);
//	p+=sprintf(p,"%s%s%14.4f%s%14.4f%s%14.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s%6.1f\n",
	p+=sprintf(p,"%s%s%14.4f%s%14.4f%s%14.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s%6.1f",
			   s,sep,enu[0],sep,enu[1],sep,enu[2],sep,sol->stat,sep,sol->ns,sep,
			   SQRT(Q[0]),sep,SQRT(Q[4]),sep,SQRT(Q[8]),sep,sqvar(Q[1]),
			   sep,sqvar(Q[5]),sep,sqvar(Q[2]),sep,sol->age,sep,sol->ratio);
	return p-(char *)buff;
}
/* output solution in the form of nmea RMC sentence --------------------------*/
extern int outnmea_rmc(unsigned char *buff, const sol_t *sol)
{
	static double dirp=0.0;
	gtime_t time;
	double ep[6],pos[3],enuv[3],dms1[3],dms2[3],vel,dir,amag=0.0;
	char *p=(char *)buff,*q,sum,*emag="E";

	trace(3,"outnmea_rmc:\n");

	if (sol->stat<=SOLQ_NONE) {
		p+=sprintf(p,"$GPRMC,,,,,,,,,,,,");
		for (q=(char *)buff+1,sum=0;*q;q++) sum^=*q;
		p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
		return p-(char *)buff;
	}
	time=gpst2utc(sol->time);
	if (time.sec>=0.995) {time.time++; time.sec=0.0;}
	time2epoch(time,ep);
	ecef2pos(sol->rr,pos);
	ecef2enu(pos,sol->rr+3,enuv);
	vel=norm(enuv,3);
	if (vel>=1.0) {
		dir=atan2(enuv[0],enuv[1])*R2D;
		if (dir<0.0) dir+=360.0;
		dirp=dir;
	}
	else dir=dirp;
	deg2dms(fabs(pos[0])*R2D,dms1);
	deg2dms(fabs(pos[1])*R2D,dms2);
	p+=sprintf(p,"$GPRMC,%02.0f%02.0f%05.2f,A,%02.0f%010.7f,%s,%03.0f%010.7f,%s,%4.2f,%4.2f,%02.0f%02.0f%02d,%.1f,%s,%s",
			   ep[3],ep[4],ep[5],dms1[0],dms1[1]+dms1[2]/60.0,pos[0]>=0?"N":"S",
			   dms2[0],dms2[1]+dms2[2]/60.0,pos[1]>=0?"E":"W",vel/KNOT2M,dir,
			   ep[2],ep[1],(int)ep[0]%100,amag,emag,
			   sol->stat==SOLQ_DGPS||sol->stat==SOLQ_FLOAT||sol->stat==SOLQ_FIX?"D":"A");
	for (q=(char *)buff+1,sum=0;*q;q++) sum^=*q; /* check-sum */
	p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
	return p-(char *)buff;
}
/* output solution in the form of nmea GGA sentence --------------------------*/
extern int outnmea_gga(unsigned char *buff, const sol_t *sol)
{
	gtime_t time;
	double h,ep[6],pos[3],dms1[3],dms2[3],dop=1.0;
	int solq;
	char *p=(char *)buff,*q,sum;

	trace(3,"outnmea_gga:\n");

	if (sol->stat<=SOLQ_NONE) {
		p+=sprintf(p,"$GPGGA,,,,,,,,,,,,,,");
		for (q=(char *)buff+1,sum=0;*q;q++) sum^=*q;
		p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
		return p-(char *)buff;
	}
	for (solq=0;solq<8;solq++) if (solq_nmea[solq]==sol->stat) break;
	if (solq>=8) solq=0;
	time=gpst2utc(sol->time);
	if (time.sec>=0.995) {time.time++; time.sec=0.0;}
	time2epoch(time,ep);
	ecef2pos(sol->rr,pos);
	h=geoidh(pos);
	deg2dms(fabs(pos[0])*R2D,dms1);
	deg2dms(fabs(pos[1])*R2D,dms2);
	p+=sprintf(p,"$GPGGA,%02.0f%02.0f%05.2f,%02.0f%010.7f,%s,%03.0f%010.7f,%s,%d,%02d,%.1f,%.3f,M,%.3f,M,%.1f,",
			   ep[3],ep[4],ep[5],dms1[0],dms1[1]+dms1[2]/60.0,pos[0]>=0?"N":"S",
			   dms2[0],dms2[1]+dms2[2]/60.0,pos[1]>=0?"E":"W",solq,
			   sol->ns,dop,pos[2]-h,h,sol->age);
	for (q=(char *)buff+1,sum=0;*q;q++) sum^=*q; /* check-sum */
	p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
	return p-(char *)buff;
}
/* output solution in the form of nmea GSA sentence --------------------------*/
extern int outnmea_gsa(unsigned char *buff, const sol_t *sol,
					   const ssat_t *ssat)
{
	double azel[MAXSAT*2],dop[4];
	int i,sat,sys,nsat,prn[MAXSAT];
	char *p=(char *)buff,*q,sum;

	trace(3,"outnmea_gsa:\n");

	if (sol->stat<=SOLQ_NONE) {
		p+=sprintf(p,"$GPGSA,A,1,,,,,,,,,,,,,,,");
		for (q=(char *)buff+1,sum=0;*q;q++) sum^=*q;
		p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
		return p-(char *)buff;
	}
	/* GPGSA: gps/sbas */
	for (sat=1,nsat=0;sat<=MAXSAT&&nsat<12;sat++) {
		if (!ssat[sat-1].vs||ssat[sat-1].azel[1]<=0.0) continue;
		sys=satsys(sat,prn+nsat);
		if (sys!=SYS_GPS&&sys!=SYS_SBS) continue;
		if (sys==SYS_SBS) prn[nsat]+=33-MINPRNSBS;
		for (i=0;i<2;i++) azel[i+nsat*2]=ssat[sat-1].azel[i];
		nsat++;
	}
	if (nsat>0) {
		p+=sprintf(p,"$GPGSA,A,%d",sol->stat<=0?1:3);
		for (i=0;i<12;i++) {
			if (i<nsat) p+=sprintf(p,",%02d",prn[i]);
			else        p+=sprintf(p,",");
		}
		dops(nsat,azel,0.0,dop);
		p+=sprintf(p,",%3.1f,%3.1f,%3.1f,1",dop[1],dop[2],dop[3]);
		for (q=(char *)buff+1,sum=0;*q;q++) sum^=*q; /* check-sum */
		p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
	}
	/* GLGSA: glonass */
	for (sat=1,nsat=0;sat<=MAXSAT&&nsat<12;sat++) {
		if (!ssat[sat-1].vs||ssat[sat-1].azel[1]<=0.0) continue;
		if (satsys(sat,prn+nsat)!=SYS_GLO) continue;
		for (i=0;i<2;i++) azel[i+nsat*2]=ssat[sat-1].azel[i];
		nsat++;
	}
	if (nsat>0) {
		p+=sprintf(p,"$GLGSA,A,%d",sol->stat<=0?1:3);
		for (i=0;i<12;i++) {
			if (i<nsat) p+=sprintf(p,",%02d",prn[i]+64);
			else        p+=sprintf(p,",");
		}
		dops(nsat,azel,0.0,dop);
		p+=sprintf(p,",%3.1f,%3.1f,%3.1f,2",dop[1],dop[2],dop[3]);
		for (q=(char *)buff+1,sum=0;*q;q++) sum^=*q; /* check-sum */
		p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
	}
	/* GAGSA: galileo */
	for (sat=1,nsat=0;sat<=MAXSAT&&nsat<12;sat++) {
		if (!ssat[sat-1].vs||ssat[sat-1].azel[1]<=0.0) continue;
		if (satsys(sat,prn+nsat)!=SYS_GAL) continue;
		for (i=0;i<2;i++) azel[i+nsat*2]=ssat[sat-1].azel[i];
		nsat++;
	}
	if (nsat>0) {
		p+=sprintf(p,"$GAGSA,A,%d",sol->stat<=0?1:3);
		for (i=0;i<12;i++) {
			if (i<nsat) p+=sprintf(p,",%02d",prn[i]+64);
			else        p+=sprintf(p,",");
		}
		dops(nsat,azel,0.0,dop);
		p+=sprintf(p,",%3.1f,%3.1f,%3.1f,3",dop[1],dop[2],dop[3]);
		for (q=(char *)buff+1,sum=0;*q;q++) sum^=*q; /* check-sum */
		p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
	}
	return p-(char *)buff;
}
/* output solution in the form of nmea GSV sentence --------------------------*/
extern int outnmea_gsv(unsigned char *buff, const sol_t *sol,
					   const ssat_t *ssat)
{
	double az,el,snr;
	int i,j,k,n,sat,prn,sys,nmsg,sats[MAXSAT];
	char *p=(char *)buff,*q,*s,sum;

	trace(3,"outnmea_gsv:\n");

	if (sol->stat<=SOLQ_NONE) {
		p+=sprintf(p,"$GPGSV,1,1,0,,,,,,,,,,,,,,,,");
		for (q=(char *)buff+1,sum=0;*q;q++) sum^=*q;
        p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
        return p-(char *)buff;
    }
    /* GPGSV: gps/sbas */
    for (sat=1,n=0;sat<MAXSAT&&n<12;sat++) {
        sys=satsys(sat,&prn);
        if (sys!=SYS_GPS&&sys!=SYS_SBS) continue;
        if (ssat[sat-1].vs&&ssat[sat-1].azel[1]>0.0) sats[n++]=sat;
    }
    nmsg=n<=0?0:(n-1)/4+1;
    
    for (i=k=0;i<nmsg;i++) {
        s=p;
        p+=sprintf(p,"$GPGSV,%d,%d,%02d",nmsg,i+1,n);
        
        for (j=0;j<4;j++,k++) {
            if (k<n) {
                if (satsys(sats[k],&prn)==SYS_SBS) prn+=33-MINPRNSBS;
                az =ssat[sats[k]-1].azel[0]*R2D; if (az<0.0) az+=360.0;
                el =ssat[sats[k]-1].azel[1]*R2D;
                snr=ssat[sats[k]-1].snr[0]*0.25;
                p+=sprintf(p,",%02d,%02.0f,%03.0f,%02.0f",prn,el,az,snr);
            }
            else p+=sprintf(p,",,,,");
        }
        p+=sprintf(p,",1"); /* L1C/A */
        for (q=s+1,sum=0;*q;q++) sum^=*q; /* check-sum */
        p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
    }
    /* GLGSV: glonass */
    for (sat=1,n=0;sat<MAXSAT&&n<12;sat++) {
        if (satsys(sat,&prn)!=SYS_GLO) continue;
        if (ssat[sat-1].vs&&ssat[sat-1].azel[1]>0.0) sats[n++]=sat;
    }
    nmsg=n<=0?0:(n-1)/4+1;
    
    for (i=k=0;i<nmsg;i++) {
        s=p;
        p+=sprintf(p,"$GLGSV,%d,%d,%02d",nmsg,i+1,n);
        
        for (j=0;j<4;j++,k++) {
            if (k<n) {
                satsys(sats[k],&prn); prn+=64; /* 65-99 */
                az =ssat[sats[k]-1].azel[0]*R2D; if (az<0.0) az+=360.0;
                el =ssat[sats[k]-1].azel[1]*R2D;
                snr=ssat[sats[k]-1].snr[0]*0.25;
                p+=sprintf(p,",%02d,%02.0f,%03.0f,%02.0f",prn,el,az,snr);
            }
            else p+=sprintf(p,",,,,");
        }
        p+=sprintf(p,",1"); /* L1C/A */
        for (q=s+1,sum=0;*q;q++) sum^=*q; /* check-sum */
        p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
    }
    /* GAGSV: galileo */
    for (sat=1,n=0;sat<MAXSAT&&n<12;sat++) {
        if (satsys(sat,&prn)!=SYS_GAL) continue;
        if (ssat[sat-1].vs&&ssat[sat-1].azel[1]>0.0) sats[n++]=sat;
    }
    nmsg=n<=0?0:(n-1)/4+1;
    
    for (i=k=0;i<nmsg;i++) {
        s=p;
        p+=sprintf(p,"$GAGSV,%d,%d,%02d",nmsg,i+1,n);
        
        for (j=0;j<4;j++,k++) {
            if (k<n) {
                satsys(sats[k],&prn); /* 1-36 */
                az =ssat[sats[k]-1].azel[0]*R2D; if (az<0.0) az+=360.0;
                el =ssat[sats[k]-1].azel[1]*R2D;
                snr=ssat[sats[k]-1].snr[0]*0.25;
                p+=sprintf(p,",%02d,%02.0f,%03.0f,%02.0f",prn,el,az,snr);
            }
            else p+=sprintf(p,",,,,");
        }
        p+=sprintf(p,",7"); /* L1BC */
        for (q=s+1,sum=0;*q;q++) sum^=*q; /* check-sum */
        p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
    }
    return p-(char *)buff;
}
/* output processing options ---------------------------------------------------
* output processing options to buffer
* args   : unsigned char *buff IO output buffer
*          prcopt_t *opt    I   processign options
* return : number of output bytes
*-----------------------------------------------------------------------------*/
extern int outprcopts(unsigned char *buff, const prcopt_t *opt)
{
    const int sys[]={SYS_GPS,SYS_GLO,SYS_GAL,SYS_QZS,SYS_SBS,0};
    const char *s1[]={"single","dgps","kinematic","static","moving-base","fixed",
                 "ppp-kinematic","ppp-static","ppp-fixed","mb-static",""};
    const char *s2[]={"L1","L1+L2","L1+L2+L5","L1+L2+L5+L6","L1+L2+L5+L6+L7",
                      "L1+L2+L5+L6+L7+L8",""};
    const char *s3[]={"forward","backward","combined"};
    const char *s4[]={"off","broadcast","sbas","iono-free","estimation",
                      "ionex tec","qzs","lex","vtec_sf","vtec_ef","gtec",""};
    const char *s5[]={"off","saastamoinen","sbas","est ztd","est ztd+grad",""};
    const char *s6[]={"broadcast","precise","broadcast+sbas","broadcast+ssr apc",
                      "broadcast+ssr com","qzss lex",""};
    const char *s7[]={"gps","glonass","galileo","qzss","sbas",""};
    const char *s8[]={"off","continuous","instantaneous","fix and hold",""};
    const char *s9[]={"off","on","auto calib","Use IFB table",""};
	const char *s10[]={"L2P","L2C",""};
	const char *s11[]={"off","table","estimation",""};
	const char *s12[]={"user settings","table",""};
	const char *s13[]={"off","table","estimation","estimation (only P)","estimation (only L)",""};
	const char *s14[]={"off","table","estimation",""};
	const char *s15[]={"in all","exc. glonass",""};
    int i;
    char *p=(char *)buff;
    
    trace(3,"outprcopts:\n");
    
    p+=sprintf(p,"%s pos mode  : %s\n",COMMENTH,s1[opt->mode]);
    
    if (PMODE_DGPS<=opt->mode&&opt->mode<=PMODE_FIXED) {
    //   p+=sprintf(p,"%s freqs     : %s\n",COMMENTH,s2[opt->nfreq-1]);
		sprintf(p, "%s", freqstrs[opt->oprfrq[0]]);
        for(i=1;i<opt->nfreq;i++) sprintf( p, "+%s", freqstrs[opt->oprfrq[i]]);
		sprintf(p, "\n" );
    }
    if (opt->mode>PMODE_SINGLE) {
        p+=sprintf(p,"%s solution  : %s\n",COMMENTH,s3[opt->soltype]);
    }
    p+=sprintf(p,"%s elev mask : %.1f deg\n",COMMENTH,opt->elmin*R2D);
    if (opt->mode>PMODE_SINGLE) {
        p+=sprintf(p,"%s dynamics  : %s\n",COMMENTH,opt->dynamics?"on":"off");
        p+=sprintf(p,"%s tidecorr  : %s\n",COMMENTH,opt->tidecorr?"on":"off");
    }
    if (opt->mode<=PMODE_FIXED) {
        p+=sprintf(p,"%s ionos opt : %s\n",COMMENTH,s4[opt->ionoopt]);
    }
    p+=sprintf(p,"%s tropo opt : %s\n",COMMENTH,s5[opt->tropopt]);
    p+=sprintf(p,"%s ephemeris : %s\n",COMMENTH,s6[opt->sateph]);
    if (opt->navsys!=SYS_GPS) {
        p+=sprintf(p,"%s navi sys  :",COMMENTH);
        for (i=0;sys[i];i++) {
            if (opt->navsys&sys[i]) p+=sprintf(p," %s",s7[i]);
        }
        p+=sprintf(p,"\n");
    }
    if (PMODE_KINEMA<=opt->mode&&opt->mode<=PMODE_FIXED) {
        p+=sprintf(p,"%s amb res   : %s\n",COMMENTH,s8[opt->modear]);
        if (opt->navsys&SYS_GLO) {
            p+=sprintf(p,"%s amb glo   : %s\n",COMMENTH,s9[opt->glomodear]);
        }
        if (opt->thresar[0]>0.0) {
            p+=sprintf(p,"%s val thres : %.1f\n",COMMENTH,opt->thresar[0]);
        }
    }
    if (opt->mode==PMODE_MOVEB&&opt->baseline[0]>0.0) {
        p+=sprintf(p,"%s baseline  : %.4f %.4f m\n",COMMENTH,
                   opt->baseline[0],opt->baseline[1]);
    }
    for (i=0;i<2;i++) {
        if (opt->mode==PMODE_SINGLE||(i>=1&&opt->mode>PMODE_FIXED)) continue;
        p+=sprintf(p,"%s antenna%d  : %-21s (%7.4f %7.4f %7.4f)\n",COMMENTH,
                   i+1,opt->anttype[i],opt->antdel[i][0],opt->antdel[i][1],
                   opt->antdel[i][2]);
    }

	p+=sprintf(p,"%s L2 code   : %s\n",COMMENTH,s10[opt->l2cprior]);
	p+=sprintf(p,"%s pha shift : %s\n",COMMENTH,s11[opt->phasshft]);
	p+=sprintf(p,"%s err model : %s\n",COMMENTH,s12[opt->errmodel]);
    p+=sprintf(p,"%s isb       : %s\n",COMMENTH,s13[opt->isb]);
    p+=sprintf(p,"%s L2P-C bias: %s\n",COMMENTH,s14[opt->gl2bias]);

    return p-(char *)buff;
}

/* output solution header ------------------------------------------------------
* output solution header to buffer
* args   : unsigned char *buff IO output buffer
*          prcopt_t *popt   I   processing options
*          solopt_t *opt    I   solution options
* return : number of output bytes
*-----------------------------------------------------------------------------*/
extern int outsolheads(unsigned char *buff, const prcopt_t *popt, const solopt_t *opt)
{
    const char *s1[]={"WGS84","Tokyo"},*s2[]={"ellipsoidal","geodetic"};
    const char *s3[]={"GPST","UTC ","JST "},*sep=opt2sep(opt);
    char *p=(char *)buff;
	int timeu=opt->timeu<0?0:(opt->timeu>20?20:opt->timeu);

    trace(3,"outsolheads:\n");
    
    if (opt->posf==SOLF_NMEA) return 0;
    
    if (opt->outhead) {
        p+=sprintf(p,"%s (",COMMENTH);
        if      (opt->posf==SOLF_XYZ) p+=sprintf(p,"x/y/z-ecef=WGS84");
        else if (opt->posf==SOLF_RVA) p+=sprintf(p,"x/y/z-ecef=WGS84");
        else if (opt->posf==SOLF_ENU) p+=sprintf(p,"e/n/u-baseline=WGS84");
        else p+=sprintf(p,"lat/lon/height=%s/%s",s1[opt->datum],s2[opt->height]);
        p+=sprintf(p,",Q=1:fix,2:float,3:sbas,4:dgps,5:single,6:ppp,ns=# of satellites)\n");
    }
	p+=sprintf(p,"%s  %-*s%s",COMMENTH,(opt->timef?16:8)+timeu+1,s3[opt->times],sep);

	if (opt->posf==SOLF_LLH) { /* lat/lon/hgt */
        if (opt->degf) {
            p+=sprintf(p,"%16s%s%16s%s%10s%s%3s%s%3s%s%8s%s%8s%s%8s%s%8s%s%8s%s%8s%s%6s%s%6s",
                       "latitude(d'\")",sep,"longitude(d'\")",sep,"height(m)",sep,
                       "Q",sep,"ns",sep,"sdn(m)",sep,"sde(m)",sep,"sdu(m)",sep,
                       "sdne(m)",sep,"sdeu(m)",sep,"sdue(m)",sep,"age(s)",sep,"ratio");
        }
        else {
            p+=sprintf(p,"%14s%s%14s%s%10s%s%3s%s%3s%s%8s%s%8s%s%8s%s%8s%s%8s%s%8s%s%6s%s%6s",
                       "latitude(deg)",sep,"longitude(deg)",sep,"height(m)",sep,
                       "Q",sep,"ns",sep,"sdn(m)",sep,"sde(m)",sep,"sdu(m)",sep,
                       "sdne(m)",sep,"sdeu(m)",sep,"sdun(m)",sep,"age(s)",sep,"ratio");
        }
    }
	else if (opt->posf==SOLF_XYZ) { /* x/y/z-ecef */
        p+=sprintf(p,"%14s%s%14s%s%14s%s%3s%s%3s%s%8s%s%8s%s%8s%s%8s%s%8s%s%8s%s%6s%s%6s",
                   "x-ecef(m)",sep,"y-ecef(m)",sep,"z-ecef(m)",sep,"Q",sep,"ns",sep,
                   "sdx(m)",sep,"sdy(m)",sep,"sdz(m)",sep,"sdxy(m)",sep,
                   "sdyz(m)",sep,"sdzx(m)",sep,"age(s)",sep,"ratio");
    }
	else if (opt->posf==SOLF_RVA) { /* r/v/a */
		p+=sprintf(p,"%14s%s%14s%s%14s%s%9s%s%9s%s%9s%s%8s%s%8s%s%8s%s%3s%s%3s%s%11s%s%11s%s%11s%s%11s%s%11s%s%11s%s"
		"%11s%s%11s%s%11s%s%11s%s%11s%s%11s%s%11s%s%11s%s%11s%s%11s%s%11s%s%11s%s%6s%s%6s",
				   "x-ecef(m)",sep,"y-ecef(m)",sep,"z-ecef(m)",sep,
				   "vx(m/s)",sep,"vy(m/s)",sep,"vz(m/s)",sep,
				   "ax(m/s2)",sep,"ay(m/s2)",sep,"az(m/s2)",sep,"Q",sep,"ns",sep,
				   "sdx(m)",sep,"sdy(m)",sep,"sdz(m)",sep,"sdxy(m)",sep,"sdyz(m)",sep,"sdzx(m)",sep,
				   "sdvx(m/s)",sep,"sdvy(m/s)",sep,"sdvz(m/s)",sep,"sdvxy(m/s)",sep,"sdvyz(m/s)",sep,"sdvzx(m/s)",sep,
				   "sdax(m/s2)",sep,"sday(m/s2)",sep,"sdaz(m/s2)",sep,"sdaxy(m/s2)",sep,"sdayz(m/s2)",sep,"sdazx(m/s2)",sep,
				   "age(s)",sep,"ratio");
	}
    else if (opt->posf==SOLF_ENU) { /* e/n/u-baseline */
        p+=sprintf(p,"%14s%s%14s%s%14s%s%3s%s%3s%s%8s%s%8s%s%8s%s%8s%s%8s%s%8s%s%6s%s%6s",
                   "e-baseline(m)",sep,"n-baseline(m)",sep,"u-baseline(m)",sep,
                   "Q",sep,"ns",sep,"sde(m)",sep,"sdn(m)",sep,"sdu(m)",sep,
                   "sden(m)",sep,"sdnu(m)",sep,"sdue(m)",sep,"age(s)",sep,"ratio");
    }
    if(opt->tropout) {
		if((popt->mode==PMODE_DGPS) || (popt->mode==PMODE_KINEMA) || (popt->mode==PMODE_STATIC) || (popt->mode==PMODE_MOVEB) || (popt->mode==PMODE_FIXED)) {
			if(popt->tropopt==TROPOPT_EST) {
				p+=sprintf(p,"%s%14s%s%14s",sep,"ztd_r(m)",sep,"ztd_b(m)");
			}
			else if(popt->tropopt==TROPOPT_ESTG) {
				p+=sprintf(p,"%s%14s%s%14s%s%14s%s%14s%s%14s%s%14s"
					,sep,"ztd_r(m)",sep,"tgn_r(m)",sep,"tge_r(r)",sep,"ztd_b(m)",sep,"tgn_b",sep,"tge_b");
			}
		}
        else if((popt->mode==PMODE_PPP_KINEMA) || (popt->mode==PMODE_PPP_STATIC) || (popt->mode==PMODE_PPP_FIXED)) {
            if(popt->tropopt==TROPOPT_EST) {
                p+=sprintf(p,"%s%14s",sep,"ztd");
            }
            else if(popt->tropopt==TROPOPT_ESTG) {
                p+=sprintf(p,"%s%14s%s%14s%s%14s",sep,"ztd(rover)",sep,"tgrad",sep,"tgrad");
            }
        }
    }
    p+=sprintf(p,"\n");

    return p-(char *)buff;
}
/* output solution body --------------------------------------------------------
* output solution body to buffer
* args   : unsigned char *buff IO output buffer
*          sol_t  *sol      I   solution
*          double *rb       I   base station position {x,y,z} (ecef) (m)
*          prcopt_t *popt   I   processing options
*          solopt_t *opt    I   solution options
* return : number of output bytes
*-----------------------------------------------------------------------------*/
extern int outsols(unsigned char *buff, const sol_t *sol, const double *rb,
                   const prcopt_t *popt, const solopt_t *opt)
{
    gtime_t time,ts={0};
    double gpst;
    int week,timeu;
    const char *sep=opt2sep(opt);
//    char s[64];
	char s[64]={0};
    unsigned char *p=buff;
    int nt=0;
    int i;
    
    trace(3,"outsols :\n");
    
    if (opt->posf==SOLF_NMEA) {
        if (opt->nmeaintv[0]<0.0) return 0;
        if (!screent(sol->time,ts,ts,opt->nmeaintv[0])) return 0;
    }
    if (sol->stat<=SOLQ_NONE||(opt->posf==SOLF_ENU&&norm(rb,3)<=0.0)) {
        return 0;
    }
	timeu=opt->timeu<0?0:(opt->timeu>20?20:opt->timeu);
    
    time=sol->time;
    if (opt->times>=TIMES_UTC) time=gpst2utc(time);
	if (opt->times==TIMES_JST) time=timeadd(time,9*3600.0);
    
	if (opt->timef) time2str(time,s,timeu);
    else {
        gpst=time2gpst(time,&week);
        if (86400*7-gpst<0.5/pow(10.0,timeu)) {
            week++;
            gpst=0.0;
        }
        sprintf(s,"%4d%s%*.*f",week,sep,6+(timeu<=0?0:timeu+1),timeu,gpst);
    }

	if(opt->tropout) {
		if((popt->mode==PMODE_DGPS) || (popt->mode==PMODE_KINEMA) || (popt->mode==PMODE_STATIC) || (popt->mode==PMODE_MOVEB) || (popt->mode==PMODE_FIXED)) {
            nt = rtknt(popt);
        }
        else if((popt->mode==PMODE_PPP_KINEMA) || (popt->mode==PMODE_PPP_STATIC) || (popt->mode==PMODE_PPP_FIXED)) {
            nt = pppnt(popt);
        }
    }
	switch (opt->posf) {
        case SOLF_LLH:
			p+=outpos (p,s,sol,opt);
			for(i=0;i<nt;i++) p+=sprintf((char*)p,"%s% 14.4f",sep,sol->trop[i]);
            p+=sprintf((char*)p,"\n");
            break;
        case SOLF_XYZ:
            p+=outecef(p,s,sol,opt);
            for(i=0;i<nt;i++) p+=sprintf((char*)p,"%s% 14.4f",sep,sol->trop[i]);
            p+=sprintf((char*)p,"\n");
            break;
        case SOLF_ENU:
            p+=outenu(p,s,sol,rb,opt);
            for(i=0;i<nt;i++) p+=sprintf((char*)p,"%s% 14.4f",sep,sol->trop[i]);
            p+=sprintf((char*)p,"\n");
            break;
        case SOLF_NMEA:
            p+=outnmea_rmc(p,sol);
            p+=outnmea_gga(p,sol);
            break;
		case SOLF_RVA:
			p+=outrva(p,s,sol,opt);
            for(i=0;i<nt;i++) p+=sprintf((char*)p,"%s% 14.4f",sep,sol->trop[i]);
            p+=sprintf((char*)p,"\n");
            break;
	}
	return p-buff;
}
/* output solution extended ----------------------------------------------------
* output solution exteneded infomation
* args   : unsigned char *buff IO output buffer
*          sol_t  *sol      I   solution
*          ssat_t *ssat     I   satellite status
*          solopt_t *opt    I   solution options
* return : number of output bytes
* notes  : only support nmea
*-----------------------------------------------------------------------------*/
extern int outsolexs(unsigned char *buff, const sol_t *sol, const ssat_t *ssat,
                     const solopt_t *opt)
{
    gtime_t ts={0};
    unsigned char *p=buff;
    
    trace(3,"outsolexs:\n");
    
    if (opt->posf==SOLF_NMEA) {
        if (opt->nmeaintv[1]<0.0) return 0;
        if (!screent(sol->time,ts,ts,opt->nmeaintv[1])) return 0;
    }
    if (opt->posf==SOLF_NMEA) {
        p+=outnmea_gsa(p,sol,ssat);
        p+=outnmea_gsv(p,sol,ssat);
    }
    return p-buff;
}
/* output estimated isb  ----------------------------------------------------
* output estimated isb to file
* args   : FILE   *fp       I   output file pointer
*          sol_t  *sol      I   solution
*          ssat_t *ssat     I   satellite status
*          solopt_t *opt    I   solution options
* return : output size (bytes)
* notes  : 
*-----------------------------------------------------------------------------*/
extern void outisbs(unsigned char *buff, const sol_t *sol, const ssat_t *ssat,
                     const solopt_t *opt)
{
    const char *s1[]={"WGS84","Tokyo"},*s2[]={"ellipsoidal","geodetic"};
    const char *s3[]={"GPST","UTC ","JST "},*sep=opt2sep(opt);
    char *p=(char *)buff;

    p+=sprintf(p,"Inter System Bias Table                                       2013/12/10\n");
    p+=sprintf(p,"------------------------------------------------------------------------\n");
    p+=sprintf(p,"RECEIVER TYPE        S F O BIAS(ns)\n");
    p+=sprintf(p,"******************** * * * **********************\n");
}


/* output processing option ----------------------------------------------------
* output processing option to file
* args   : FILE   *fp       I   output file pointer
*          prcopt_t *opt    I   processing options
* return : none
*-----------------------------------------------------------------------------*/
extern void outprcopt(FILE *fp, const prcopt_t *opt)
{
    unsigned char buff[MAXSOLMSG+1];
    int n;
    
    trace(3,"outprcopt:\n");
    
    if ((n=outprcopts(buff,opt))>0) {
        fwrite(buff,n,1,fp);
    }
}
/* output solution header ------------------------------------------------------
* output solution heade to file
* args   : FILE   *fp       I   output file pointer
*          prcopt_t *popt   I   processing options
*          solopt_t *opt    I   solution options
* return : none
*-----------------------------------------------------------------------------*/
extern void outsolhead(FILE *fp, const prcopt_t *popt, const solopt_t *opt)
{
    unsigned char buff[MAXSOLMSG+1];
    int n;
    
    trace(3,"outsolhead:\n");
    
    if ((n=outsolheads(buff,popt,opt))>0) {
        fwrite(buff,n,1,fp);
    }
}
/* output solution body --------------------------------------------------------
* output solution body to file
* args   : FILE   *fp       I   output file pointer
*          sol_t  *sol      I   solution
*          double *rb       I   base station position {x,y,z} (ecef) (m)
*          prcopt_t *popt   I   processing options
*          solopt_t *opt    I   solution options
* return : none
*-----------------------------------------------------------------------------*/
extern void outsol(FILE *fp, const sol_t *sol, const double *rb,
                   const prcopt_t *popt, const solopt_t *opt)
{
//    unsigned char buff[MAXSOLMSG+1];
    unsigned char buff[MAXSOLMSG+1] = {0};
    int n;
    
    trace(3,"outsol  :\n");

	if ((n=outsols(buff,sol,rb,popt,opt))>0) {
	fprintf(stdout, "");
		fwrite(buff,n,1,fp);

	fflush(fp);
    }
}
/* output solution extended ----------------------------------------------------
* output solution exteneded infomation to file
* args   : FILE   *fp       I   output file pointer
*          sol_t  *sol      I   solution
*          ssat_t *ssat     I   satellite status
*          solopt_t *opt    I   solution options
* return : output size (bytes)
* notes  : only support nmea
*-----------------------------------------------------------------------------*/
extern void outsolex(FILE *fp, const sol_t *sol, const ssat_t *ssat,
                     const solopt_t *opt)
{
    unsigned char buff[MAXSOLMSG+1];
    int n;
    
    trace(3,"outsolex:\n");
    
    if ((n=outsolexs(buff,sol,ssat,opt))>0) {
        fwrite(buff,n,1,fp);
    }
}
/* output estimated isb  ----------------------------------------------------
* output estimated isb to file
* args   : FILE   *fp       I   output file pointer
*          ssat_t *ssat     I   satellite status
*          solopt_t *opt    I   solution options
*          sol_t  *sol      I   solution
*          mbs_t  *mbs      I   solution(MBS)
* return : output size (bytes)
* notes  : 
*-----------------------------------------------------------------------------*/
//extern void outisb(FILE *fp, const sol_t *sol, const ssat_t *ssat,
//                     const solopt_t *opt)
extern void outisbtable(const prcopt_t *popt, const solopt_t *sopt, const sol_t *sol, mbs_t *mbs)
{
	int i,f,fr;
	FILE *fp;
	char mode[2] = "";
	double isb;
	int ir;
	int sys;
	long sz;
	int nf=pppnf(popt);
	char recname[256]={0};
	char recname_base[256]={0};
	char sysch;

	if(    (popt->mode!=PMODE_SINGLE)
		&& (popt->mode!=PMODE_PPP_KINEMA) && (popt->mode!=PMODE_PPP_STATIC) && (popt->mode!=PMODE_PPP_FIXED)
		&& (popt->mode!=PMODE_KINEMA) && (popt->mode!=PMODE_STATIC)  && (popt->mode!=PMODE_FIXED)
		&& (popt->mode!=PMODE_MULTI)) {
		return;
	}

	if (sopt->isbout==ISBOUTOPT_OFF) return;
	mode[0]=sopt->isbout==ISBOUTOPT_NEW?'w':'a';

	if (!*sopt->isbfile) return;



	if((fp=fopen(sopt->isbfile,mode))==NULL){
		showmsg("isb output file open error: %s",sopt->isbfile);
		trace(1,"isb output file open error: %s\n",sopt->isbfile);
		return;
	}

	trace(3,"outisb:\n");
	fseek(fp, 0, SEEK_END);
	if ((ftell(fp)==0) || (sopt->isbout==ISBOUTOPT_NEW)) {
		fprintf(fp,"Inter System Bias Table                                                 \n");
		fprintf(fp,"----------------------------------------------------------------------\n");
		fprintf(fp,"RECEIVER TYPE        RECEIVER TYPE (BASE) S F O               BIAS(ns)\n");
		fprintf(fp,"******************** ******************** * * * **********************\n");
	}

	if((popt->mode==PMODE_KINEMA) || (popt->mode!=PMODE_STATIC) || (popt->mode!=PMODE_FIXED)) {
		if(strlen(popt->rectype[0])>0) {
			strcpy(recname, popt->rectype[0]);
		}
		else {
			strcpy(recname, stas[0].rectype);
		}
		if(strlen(popt->rectype[1])>0) {
			strcpy(recname_base, popt->rectype[1]);
		}
		else {
			strcpy(recname_base, stas[1].rectype);
		}

		for(i=ISYSGPS+1;i<=NSYS;i++) {
			sysch=syschar(i);
			if((popt->isb==ISBOPT_EST) || (popt->isb==ISBOPT_EST_L) || (popt->isb==ISBOPT_EST_0M)) {
				if(popt->ionoopt!=IONOOPT_IFLC) {
					for(f=0;f<nf;f++) {

						fr=popt->oprfrq[f];
						if (!(sysno(i)&sol->sys_prd[fr])) continue;
						fprintf(fp,"%-20s %-20s %c %d L % 22.9f\n", recname, recname_base, sysch, freqs[fr], (sol->isb[i-ISYSGPS][f][0]-sol->isb[0][f][0])/CLIGHT*1.0e9);
						recname[0] = '\0';
					}
				}
				else {
					if(i==ISYSGAL) f=2;
					else           f=1;
					if (!(sysno(i)&sol->sys_prd[0])) continue;
					fprintf(fp,"%-20s %-20s %c %d L % 22.9f\n", recname, recname_base , sysch, freqs[0], sol->isb[i-ISYSGPS][0][0]/CLIGHT*1.0e9);
					fprintf(fp,"%-20s %-20s %c %d L % 22.9f\n", ""     , ""           , sysch, freqs[f], sol->isb[i-ISYSGPS][1][0]/CLIGHT*1.0e9);
					recname[0] = '\0';
				}

			}
			if((popt->isb==ISBOPT_EST) || (popt->isb==ISBOPT_EST_P) || (popt->isb==ISBOPT_EST_0M)) {
				if(popt->ionoopt!=IONOOPT_IFLC) {
					for(f=0;f<nf;f++) {

						fr=popt->oprfrq[f];
						if (!(sysno(i)&sol->sys_prd[fr])) continue;
						fprintf(fp,"%-20s %-20s %c %d P % 22.9f\n", recname, recname_base, sysch, freqs[fr], sol->isb[i-ISYSGPS][f][1]/CLIGHT*1.0e9);
						recname[0] = '\0';
					}
				}
				else {
					if(i==ISYSGAL) f=2;
					if (!(sysno(i)&sol->sys_prd[0])) continue;
					fprintf(fp,"%-20s %-20s %c %d P % 22.9f\n", recname, recname_base, sysch, freqs[0], sol->isb[i-ISYSGPS][0][1]/CLIGHT*1.0e9);
					fprintf(fp,"%-20s %-20s %c %d P % 22.9f\n", ""     , ""          , sysch, freqs[f], sol->isb[i-ISYSGPS][1][1]/CLIGHT*1.0e9);
					recname[0] = '\0';
				}
			}
		}
	}
	else if(popt->mode!=PMODE_MULTI) {
		if(strlen(popt->rectype[0])>0) {
			strcpy(recname, popt->rectype[0]);
		}
		else {
			strcpy(recname, stas[0].rectype);
		}

		for(i=ISYSGPS;i<=NSYS;i++) {
			sysch=syschar(i);
			if((popt->mode!=PMODE_SINGLE) && (popt->mode!=PMODE_DGPS)) {
				if((popt->isb==ISBOPT_EST) || (popt->isb==ISBOPT_EST_L)) {
					if(popt->ionoopt!=IONOOPT_IFLC) {
						for(f=0;f<nf;f++) {
							fr=popt->oprfrq[f];
							if (!(sysno(i)&sol->sys_prd[fr])) continue;
							fprintf(fp,"%-20s %-20s %c %d L % 22.9f\n", recname, "", sysch, freqs[fr], sol->isb[i-ISYSGPS][f][0]/CLIGHT*1.0e9);
							recname[0] = '\0';
						}
					}
					else {
						if(i==ISYSGAL) f=2;
						else           f=1;
						if (!(sysno(i)&sol->sys_prd[0])) continue;
						fprintf(fp,"%-20s %-20s %c %d L % 22.9f\n", recname, "", sysch, freqs[0], sol->isb[i-ISYSGPS][0][0]/CLIGHT*1.0e9);
						fprintf(fp,"%-20s %-20s %c %d L % 22.9f\n", ""     , "", sysch, freqs[f], sol->isb[i-ISYSGPS][1][0]/CLIGHT*1.0e9);
						recname[0] = '\0';
					}

				}
			}
			if((popt->isb==ISBOPT_EST) || (popt->isb==ISBOPT_EST_P) || (popt->isb==ISBOPT_EST_0M)) {
				if(popt->ionoopt!=IONOOPT_IFLC) {
					for(f=0;f<nf;f++) {

						fr=popt->oprfrq[f];
						if (!(sysno(i)&sol->sys_prd[fr])) continue;
						fprintf(fp,"%-20s %-20s %c %d P % 22.9f\n", recname, "",sysch, freqs[fr], sol->isb[i-ISYSGPS][f][1]/CLIGHT*1.0e9);
						recname[0] = '\0';
					}
				}
				else {
					if(i==ISYSGAL) f=2;
					if (!(sysno(i)&sol->sys_prd[0])) continue;
					fprintf(fp,"%-20s %-20s %c %d P % 22.9f\n", recname, "", sysch, freqs[0], sol->isb[i-ISYSGPS][0][1]/CLIGHT*1.0e9);
					fprintf(fp,"%-20s %-20s %c %d P % 22.9f\n", ""     , "", sysch, freqs[f], sol->isb[i-ISYSGPS][1][1]/CLIGHT*1.0e9);
					recname[0] = '\0';
				}
			}
		}
	}
	else {
		for(ir=0;ir<mbs->stas.nsta;ir++) {
			for (sys=SYS_GPS,i=1;sys<=popt->navsys;sys<<=1) {
				if(!(popt->navsys & sys)) continue;
				if((popt->isb==ISBOPT_EST) || (popt->isb==ISBOPT_EST_L)) {
					isb = mbs->param.par[mbs->param.iisbL + (i - 1) * mbs->stas.nsta + ir];
					fprintf(fp,"%-20s %c   L % 22.9f\n", mbs->stas.sta[ir].rectype, syschar(sysind(sys)), isb/CLIGHT*1.0e9);
				}
				if((popt->isb==ISBOPT_EST) || (popt->isb==ISBOPT_EST_P) || (popt->isb==ISBOPT_EST_0M)) {
					isb = mbs->param.par[mbs->param.iisbP + (i - 1) * mbs->stas.nsta + ir];
					fprintf(fp,"%-20s %c   P % 22.9f\n", mbs->stas.sta[ir].rectype, syschar(sysind(sys)), isb/CLIGHT*1.0e9);
				}
				++i;
			}
		}
	}

	fclose(fp);
}

/* output estimated gl2  ----------------------------------------------------
* output estimated gl2 to file
* args   : FILE   *fp       I   output file pointer
*          ssat_t *ssat     I   satellite status
*          solopt_t *opt    I   solution options
*          sol_t  *sol      I   solution
*          mbs_t  *mbs      I   solution(MBS)
* return : output size (bytes)
* notes  :
*-----------------------------------------------------------------------------*/
extern void outgl2table(const prcopt_t *popt, const solopt_t *sopt, const sol_t *sol, mbs_t *mbs)
{
	int i;
	FILE *fp;
	char mode[2] = "";
	double l2b;
	int ir;
	int sys;
	long sz;

	if(popt->gl2bias!=GL2OPT_EST) return;

	if(sopt->gl2out==ISBOUTOPT_OFF) return;

	mode[0]=sopt->gl2out==GL2OUTOPT_NEW?'w':'a';

	if (!*sopt->gl2file) return;


	if((fp=fopen(sopt->gl2file,mode))==NULL){
		showmsg("gl2 output file open error: %s",sopt->gl2file);
		trace(1,"gl2 output file open error: %s\n",sopt->gl2file);
		return;
	}

	trace(3,"outl2pc:\n");

	fprintf(fp,"\n");
	fprintf(fp,"DIFFERENTIAL (P2-C2) CODE BIASES FOR SATELLITES AND RECEIVERS:\n");
	fprintf(fp,"\n");
	fprintf(fp,"PRN / STATION NAME        VALUE (NS)  RMS (NS) \n");
	fprintf(fp,"***   ****************    *****.***   *****.***\n");

	if(popt->mode!=PMODE_MULTI) {
		fprintf(fp,"G%2s   %-16s    %9.3f   %9s\n", "", stas[0].name, sol->gl2[0]/CLIGHT*1.0e9,"");
		if((popt->mode>=PMODE_DGPS) && (popt->mode<=PMODE_FIXED)) {
			fprintf(fp,"G%2s   %-16s    %9.3f   %9s\n", "", stas[1].name, sol->gl2[1]/CLIGHT*1.0e9,"");
		}
	}
	else {
		for(ir=0;ir<mbs->stas.nsta;ir++) {
			l2b = mbs->param.par[mbs->param.il2b + ir];
			fprintf(fp,"G%2s   %-16s    %9.3f   %9s\n", "", mbs->stas.sta[ir].name2, l2b/CLIGHT*1.0e9,"");
		}
    }

	fclose(fp);
}

