/*------------------------------------------------------------------------------
* solution_read.c : functions to input solution
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

/* separate fields -----------------------------------------------------------*/
static int tonum(char *buff, const char *sep, double *v)
{
    int n,len=(int)strlen(sep);
    char *p,*q;
    
    for (p=buff,n=0;n<MAXFIELD;p=q+len) {
        if ((q=strstr(p,sep))) *q='\0'; 
        if (*p) v[n++]=atof(p);
        if (!q) break;
    }
    return n;
}
/* convert ddmm.mm in nmea format to deg -------------------------------------*/
static double dmm2deg(double dmm)
{
    return floor(dmm/100.0)+fmod(dmm,100.0)/60.0;
}
/* convert time in nmea format to time ---------------------------------------*/
static void septime(double t, double *t1, double *t2, double *t3)
{
    *t1=floor(t/10000.0);
    t-=*t1*10000.0;
    *t2=floor(t/100.0);
    *t3=t-*t2*100.0;
}
/* covariance to solution ----------------------------------------------------*/
//static void covtosol(const double *P, sol_t *sol)
static void covtosol(const double *P, float *qr)
{
	qr[0]=(float)P[0]; /* xx or ee */
	qr[1]=(float)P[4]; /* yy or nn */
	qr[2]=(float)P[8]; /* zz or uu */
	qr[3]=(float)P[1]; /* xy or en */
	qr[4]=(float)P[5]; /* yz or nu */
	qr[5]=(float)P[2]; /* zx or ue */
}
/* decode nmea gprmc: recommended minumum data for gps -----------------------*/
static int decode_nmearmc(char **val, int n, sol_t *sol)
{
    double tod=0.0,lat=0.0,lon=0.0,vel=0.0,dir=0.0,date=0.0,ang=0.0,ep[6];
    double pos[3]={0};
    char act=' ',ns='N',ew='E',mew='E',mode='A';
    int i;
    
    trace(4,"decode_nmearmc: n=%d\n",n);
    
    for (i=0;i<n;i++) {
        switch (i) {
            case  0: tod =atof(val[i]); break; /* time in utc (hhmmss) */
            case  1: act =*val[i];      break; /* A=active,V=void */
            case  2: lat =atof(val[i]); break; /* latitude (ddmm.mmm) */
            case  3: ns  =*val[i];      break; /* N=north,S=south */
            case  4: lon =atof(val[i]); break; /* longitude (dddmm.mmm) */
            case  5: ew  =*val[i];      break; /* E=east,W=west */
            case  6: vel =atof(val[i]); break; /* speed (knots) */
            case  7: dir =atof(val[i]); break; /* track angle (deg) */
            case  8: date=atof(val[i]); break; /* date (ddmmyy) */
            case  9: ang =atof(val[i]); break; /* magnetic variation */
            case 10: mew =*val[i];      break; /* E=east,W=west */
            case 11: mode=*val[i];      break; /* mode indicator (>nmea 2) */
                                      /* A=autonomous,D=differential */
                                      /* E=estimated,N=not valid,S=simulator */
        }
    }
    if ((act!='A'&&act!='V')||(ns!='N'&&ns!='S')||(ew!='E'&&ew!='W')) {
        trace(2,"invalid nmea gprmc format\n");
        return 0;
    }
    pos[0]=(ns=='S'?-1.0:1.0)*dmm2deg(lat)*D2R;
    pos[1]=(ew=='W'?-1.0:1.0)*dmm2deg(lon)*D2R;
    septime(date,ep+2,ep+1,ep);
    septime(tod,ep+3,ep+4,ep+5);
    ep[0]+=ep[0]<80.0?2000.0:1900.0;
    sol->time=utc2gpst(epoch2time(ep));
    pos2ecef(pos,sol->rr);
    sol->stat=mode=='D'?SOLQ_DGPS:SOLQ_SINGLE;
    sol->ns=0;
    
    sol->type=0; /* postion type = xyz */
    
    trace(5,"decode_nmearmc: %s rr=%.3f %.3f %.3f stat=%d ns=%d vel=%.2f dir=%.0f ang=%.0f mew=%c mode=%c\n",
          time_str(sol->time,0),sol->rr[0],sol->rr[1],sol->rr[2],sol->stat,sol->ns,
          vel,dir,ang,mew,mode);
    
    return 1;
}
/* decode nmea gpgga: fix information ----------------------------------------*/
static int decode_nmeagga(char **val, int n, sol_t *sol)
{
    gtime_t time;
    double tod=0.0,lat=0.0,lon=0.0,hdop=0.0,alt=0.0,msl=0.0,ep[6],tt;
    double pos[3]={0};
    char ns='N',ew='E',ua=' ',um=' ';
    int i,solq=0,nrcv=0;
    
    trace(4,"decode_nmeagga: n=%d\n",n);
    
    for (i=0;i<n;i++) {
        switch (i) {
            case  0: tod =atof(val[i]); break; /* time in utc (hhmmss) */
            case  1: lat =atof(val[i]); break; /* latitude (ddmm.mmm) */
            case  2: ns  =*val[i];      break; /* N=north,S=south */
            case  3: lon =atof(val[i]); break; /* longitude (dddmm.mmm) */
            case  4: ew  =*val[i];      break; /* E=east,W=west */
            case  5: solq=atoi(val[i]); break; /* fix quality */
            case  6: nrcv=atoi(val[i]); break; /* # of satellite tracked */
            case  7: hdop=atof(val[i]); break; /* hdop */
            case  8: alt =atof(val[i]); break; /* altitude in msl */
            case  9: ua  =*val[i];      break; /* unit (M) */
            case 10: msl =atof(val[i]); break; /* height of geoid */
            case 11: um  =*val[i];      break; /* unit (M) */
        }
    }
    if ((ns!='N'&&ns!='S')||(ew!='E'&&ew!='W')) {
        trace(2,"invalid nmea gpgga format\n");
        return 0;
    }
    if (sol->time.time==0.0) {
        trace(2,"no date info for nmea gpgga\n");
        return 0;
    }
    pos[0]=(ns=='N'?1.0:-1.0)*dmm2deg(lat)*D2R;
    pos[1]=(ew=='E'?1.0:-1.0)*dmm2deg(lon)*D2R;
    pos[2]=alt+msl;
    
    time2epoch(sol->time,ep);
    septime(tod,ep+3,ep+4,ep+5);
    time=utc2gpst(epoch2time(ep));
    tt=timediff(time,sol->time);
    if      (tt<-43200.0) sol->time=timeadd(time, 86400.0);
    else if (tt> 43200.0) sol->time=timeadd(time,-86400.0);
    else sol->time=time;
    pos2ecef(pos,sol->rr);
    sol->stat=0<=solq&&solq<=8?solq_nmea[solq]:SOLQ_NONE;
    sol->ns=nrcv;
    
    sol->type=0; /* postion type = xyz */
    
    trace(5,"decode_nmeagga: %s rr=%.3f %.3f %.3f stat=%d ns=%d hdop=%.1f ua=%c um=%c\n",
          time_str(sol->time,0),sol->rr[0],sol->rr[1],sol->rr[2],sol->stat,sol->ns,
          hdop,ua,um);
    
    return 1;
}
/* decode nmea ---------------------------------------------------------------*/
static int decode_nmea(char *buff, sol_t *sol)
{
    char *p,*q,*val[MAXFIELD];
    int n=0;
    
    trace(4,"decode_nmea: buff=%s\n",buff);
    
    /* parse fields */
    for (p=buff;*p&&n<MAXFIELD;p=q+1) {
        if ((q=strchr(p,','))||(q=strchr(p,'*'))) {
            val[n++]=p; *q='\0';
        }
        else break;
    }
    /* decode nmea sentence */
    if (!strcmp(val[0],"$GPRMC")) {
        return decode_nmearmc(val+1,n-1,sol);
    }
    else if (!strcmp(val[0],"$GPGGA")) {
        return decode_nmeagga(val+1,n-1,sol);
    }
    return 0;
}
/* decode solution time ------------------------------------------------------*/
static char *decode_soltime(char *buff, const solopt_t *opt, gtime_t *time)
{
    double v[MAXFIELD];
    char *p,*q,s[64]=" ";
    int n,len;
    
    trace(4,"decode_soltime:\n");
    
    if (!strcmp(opt->sep,"\\t")) strcpy(s,"\t");
    else if (*opt->sep) strcpy(s,opt->sep);
    len=(int)strlen(s);
    
    /* yyyy/mm/dd hh:mm:ss or yyyy mm dd hh:mm:ss */
    if (sscanf(buff,"%lf/%lf/%lf %lf:%lf:%lf",v,v+1,v+2,v+3,v+4,v+5)>=6) {
        if (v[0]<100.0) {
            v[0]+=v[0]<80.0?2000.0:1900.0;
        }
        *time=epoch2time(v);
        if (opt->times==TIMES_UTC) {
            *time=utc2gpst(*time);
        }
        else if (opt->times==TIMES_JST) {
            *time=utc2gpst(timeadd(*time,-9*3600.0));
        }
        if (!(p=strchr(buff,':'))||!(p=strchr(p+1,':'))) return NULL;
        for (p++;isdigit((int)*p)||*p=='.';) p++;
        return p+len;
    }
    if (opt->posf==SOLF_GSIF) {
        if (sscanf(buff,"%lf %lf %lf %lf:%lf:%lf",v,v+1,v+2,v+3,v+4,v+5)<6) {
            return NULL;
        }
        *time=timeadd(epoch2time(v),-12.0*3600.0);
        if (!(p=strchr(buff,':'))||!(p=strchr(p+1,':'))) return NULL;
        for (p++;isdigit((int)*p)||*p=='.';) p++;
        return p+len;
    }
    /* wwww ssss */
    for (p=buff,n=0;n<2;p=q+len) {
        if ((q=strstr(p,s))) *q='\0'; 
        if (*p) v[n++]=atof(p);
        if (!q) break;
    }
    if (n>=2&&0.0<=v[0]&&v[0]<=3000.0&&0.0<=v[1]&&v[1]<604800.0) {
        *time=gpst2time((int)v[0],v[1]);
        return p;
    }
    return NULL;
}
/* decode x/y/z-ecef ---------------------------------------------------------*/
static int decode_solxyz(char *buff, const solopt_t *opt, sol_t *sol)
{
    double val[MAXFIELD],P[9]={0};
    int i=0,j,n;
    const char *sep=opt2sep(opt);
    
    trace(4,"decode_solxyz:\n");
    
    if ((n=tonum(buff,sep,val))<3) return 0;
    
    for (j=0;j<3;j++) {
        sol->rr[j]=val[i++]; /* xyz */
    }
    if (i<n) sol->stat=(unsigned char)val[i++];
    if (i<n) sol->ns  =(unsigned char)val[i++];
	if (i+3<n) {
		P[0]=val[i]*val[i]; i++; /* sdx */
        P[4]=val[i]*val[i]; i++; /* sdy */
		P[8]=val[i]*val[i]; i++; /* sdz */
        if (i+3<n) {
            P[1]=P[3]=NSQR(val[i]); i++; /* sdxy */
            P[5]=P[7]=NSQR(val[i]); i++; /* sdyz */
            P[2]=P[6]=NSQR(val[i]); i++; /* sdzx */
        }
		covtosol(P,sol->qr);
    }
	if (i<n) sol->age  =(float)val[i++];
    if (i<n) sol->ratio=(float)val[i];
    
	sol->type=0; /* postion type = xyz */
	sol->rva=0; /* postion type != rva */
    
    if (MAXSOLQ<sol->stat) sol->stat=SOLQ_NONE;
    return 1;
}
/* decode x/y/z-ecef ---------------------------------------------------------*/
static int decode_solrva(char *buff, const solopt_t *opt, sol_t *sol)
{
    double val[MAXFIELD],P[9]={0};
    int i=0,j,n;
    const char *sep=opt2sep(opt);

	trace(4,"decode_solrva:\n");

    if ((n=tonum(buff,sep,val))<3) return 0;

	for (j=0;j<9;j++) {
		sol->rr[j]=val[i++]; /* rva */
	}
	if (i<n) sol->stat=(unsigned char)val[i++];
	if (i<n) sol->ns  =(unsigned char)val[i++];
	if (i+3<n) {
		P[0]=val[i]*val[i]; i++; /* sdx */
		P[4]=val[i]*val[i]; i++; /* sdy */
		P[8]=val[i]*val[i]; i++; /* sdz */
		if (i+3<n) {
			P[1]=P[3]=NSQR(val[i]); i++; /* sdxy */
			P[5]=P[7]=NSQR(val[i]); i++; /* sdyz */
			P[2]=P[6]=NSQR(val[i]); i++; /* sdzx */
		}
		covtosol(P,sol->qr);
	}
	if (i+6<n) {
		P[0]=val[i]*val[i]; i++; /* sdx */
		P[4]=val[i]*val[i]; i++; /* sdy */
		P[8]=val[i]*val[i]; i++; /* sdz */
		P[1]=P[3]=NSQR(val[i]); i++; /* sdxy */
		P[5]=P[7]=NSQR(val[i]); i++; /* sdyz */
		P[2]=P[6]=NSQR(val[i]); i++; /* sdzx */
		covtosol(P,sol->qr+6);
	}
	if (i+6<n) {
		P[0]=val[i]*val[i]; i++; /* sdx */
		P[4]=val[i]*val[i]; i++; /* sdy */
		P[8]=val[i]*val[i]; i++; /* sdz */
		P[1]=P[3]=NSQR(val[i]); i++; /* sdxy */
		P[5]=P[7]=NSQR(val[i]); i++; /* sdyz */
		P[2]=P[6]=NSQR(val[i]); i++; /* sdzx */
		covtosol(P,sol->qr+12);
	}
	if (i<n) sol->age  =(float)val[i++];
    if (i<n) sol->ratio=(float)val[i];
    
	sol->type=0; /* postion type = xyz */
	sol->rva=1; /* postion type = rva */
    
    if (MAXSOLQ<sol->stat) sol->stat=SOLQ_NONE;
    return 1;
}
/* decode lat/lon/height -----------------------------------------------------*/
static int decode_solllh(char *buff, const solopt_t *opt, sol_t *sol)
{
    double val[MAXFIELD],pos[3],Q[9]={0},P[9];
    int i=0,n;
    const char *sep=opt2sep(opt);
    
    trace(4,"decode_solllh:\n");
    
    n=tonum(buff,sep,val);
    
    if (!opt->degf) {
        if (n<3) return 0;
        pos[0]=val[i++]*D2R; /* lat/lon/hgt (ddd.ddd) */
        pos[1]=val[i++]*D2R;
        pos[2]=val[i++];
    }
    else {
        if (n<7) return 0;
        pos[0]=dms2deg(val  )*D2R; /* lat/lon/hgt (ddd mm ss) */
        pos[1]=dms2deg(val+3)*D2R;
        pos[2]=val[6];
        i+=7;
    }
    pos2ecef(pos,sol->rr);
    if (i<n) sol->stat=(unsigned char)val[i++];
    if (i<n) sol->ns  =(unsigned char)val[i++];
    if (i+3<n) {
        Q[4]=val[i]*val[i]; i++; /* sdn */
        Q[0]=val[i]*val[i]; i++; /* sde */
        Q[8]=val[i]*val[i]; i++; /* sdu */
        if (i+3<n) {
            Q[1]=Q[3]=NSQR(val[i]); i++; /* sdne */
            Q[2]=Q[6]=NSQR(val[i]); i++; /* sdeu */
			Q[5]=Q[7]=NSQR(val[i]); i++; /* sdun */
		}
		covecef(pos,Q,P);
		covtosol(P,sol->qr);
    }
    if (i<n) sol->age  =(float)val[i++];
    if (i<n) sol->ratio=(float)val[i];
    
    sol->type=0; /* postion type = xyz */
    sol->rva=0; /* postion type != rva */

    if (MAXSOLQ<sol->stat) sol->stat=SOLQ_NONE;
    return 1;
}
/* decode e/n/u-baseline -----------------------------------------------------*/
static int decode_solenu(char *buff, const solopt_t *opt, sol_t *sol)
{
    double val[MAXFIELD],Q[9]={0};
    int i=0,j,n;
    const char *sep=opt2sep(opt);
    
    trace(4,"decode_solenu:\n");
    
    if ((n=tonum(buff,sep,val))<3) return 0;
    
    for (j=0;j<3;j++) {
        sol->rr[j]=val[i++]; /* enu */
    }
    if (i<n) sol->stat=(unsigned char)val[i++];
    if (i<n) sol->ns  =(unsigned char)val[i++];
    if (i+3<n) {
        Q[0]=val[i]*val[i]; i++; /* sde */
        Q[4]=val[i]*val[i]; i++; /* sdn */
        Q[8]=val[i]*val[i]; i++; /* sdu */
        if (i+3<n) {
            Q[1]=Q[3]=NSQR(val[i]); i++; /* sden */
            Q[5]=Q[7]=NSQR(val[i]); i++; /* sdnu */
            Q[2]=Q[6]=NSQR(val[i]); i++; /* sdue */
        }
		covtosol(Q,sol->qr);
    }
    if (i<n) sol->age  =(float)val[i++];
    if (i<n) sol->ratio=(float)val[i];
    
    sol->type=1; /* postion type = enu */
    sol->rva=0; /* postion type != rva */

    if (MAXSOLQ<sol->stat) sol->stat=SOLQ_NONE;
    return 1;
}
/* decode gsi f solution -----------------------------------------------------*/
static int decode_solgsi(char *buff, const solopt_t *opt, sol_t *sol)
{
    double val[MAXFIELD];
    int i=0,j;
    
    trace(4,"decode_solgsi:\n");
    
    if (tonum(buff," ",val)<3) return 0;
    
    for (j=0;j<3;j++) {
        sol->rr[j]=val[i++]; /* xyz */
    }
    sol->stat=SOLQ_FIX;
    return 1;
}
/* decode solution position --------------------------------------------------*/
static int decode_solpos(char *buff, const solopt_t *opt, sol_t *sol)
{
    sol_t sol0={{0}};
    char *p=buff;
    
    trace(4,"decode_solpos: buff=%s\n",buff);
    
    *sol=sol0;
    
    /* decode solution time */
    if (!(p=decode_soltime(p,opt,&sol->time))) {
        return 0;
    }
    /* decode solution position */
    switch (opt->posf) {
		case SOLF_XYZ : return decode_solxyz(p,opt,sol);
        case SOLF_LLH : return decode_solllh(p,opt,sol);
        case SOLF_ENU : return decode_solenu(p,opt,sol);
		case SOLF_RVA : return decode_solrva(p,opt,sol);
		case SOLF_GSIF: return decode_solgsi(p,opt,sol);
	}
    return 0;
}
/* decode reference position -------------------------------------------------*/
static void decode_refpos(char *buff, const solopt_t *opt, double *rb)
{
    double val[MAXFIELD],pos[3];
    int i,n;
    const char *sep=opt2sep(opt);
    
    trace(3,"decode_refpos: buff=%s\n",buff);
    
    if ((n=tonum(buff,sep,val))<3) return;
    
	if ((opt->posf==SOLF_XYZ) || (opt->posf==SOLF_RVA)) { /* xyz */
        for (i=0;i<3;i++) rb[i]=val[i];
    }
    else if (opt->degf==0) { /* lat/lon/hgt (ddd.ddd) */
        pos[0]=val[0]*D2R;
        pos[1]=val[1]*D2R;
        pos[2]=val[2];
        pos2ecef(pos,rb);
    }
    else if (opt->degf==1&&n>=7) { /* lat/lon/hgt (ddd mm ss) */
        pos[0]=dms2deg(val  )*D2R;
        pos[1]=dms2deg(val+3)*D2R;
        pos[2]=val[6];
        pos2ecef(pos,rb);
    }
}
/* decode solution -----------------------------------------------------------*/
static int decode_sol(char *buff, const solopt_t *opt, sol_t *sol, double *rb)
{
    char *p;
    
    trace(4,"decode_sol: buff=%s\n",buff);
    
    if (!strncmp(buff,COMMENTH,1)) { /* reference position */
        if (!strstr(buff,"ref pos")&&!strstr(buff,"slave pos")) return 0;
        if (!(p=strchr(buff,':'))) return 0;
        decode_refpos(p+1,opt,rb);
        return 0;
    }
    if (!strncmp(buff,"$GP",3)) { /* decode nmea */
        if (!decode_nmea(buff,sol)) return 0;
        
        /* for time update only */
        if (opt->posf!=SOLF_NMEA&&!strncmp(buff,"$GPRMC",6)) return 2;
    }
    else { /* decode position record */
        if (!decode_solpos(buff,opt,sol)) return 0;
    }
    return 1;
}
/* decode solution options ---------------------------------------------------*/
static void decode_solopt(char *buff, solopt_t *opt, char *posmod, int *lastresult)
{
    char *p,c;
    int i,j;
	char sol[10]={0};
	char *solname[3]={SOL0,SOL1,SOL2};
	trace(4,"decode_solhead: buff=%s\n",buff);
    
    if (strncmp(buff,COMMENTH,1)&&strncmp(buff,"+",1)) return;
    
    if      (strstr(buff,"GPST")) opt->times=TIMES_GPST;
    else if (strstr(buff,"UTC" )) opt->times=TIMES_UTC;
    else if (strstr(buff,"JST" )) opt->times=TIMES_JST;
    
    if ((p=strstr(buff,"pos mode"))) {
        for(i=0;i<20;i++){
        strncpy(&c,buff+14+i,1);
            if((c<='z'&&c>='a')||c=='-'){
                strncpy(posmod+i,buff+14+i,1);
            }else{
				*(posmod+i)='\0';
				if(!strcmp(posmod,"single")){
                	*lastresult=0;
				}
				break;
			}
		}
	}

	if ((p=strstr(buff,"solution"))) {
		for(i=0;i<10;i++){
			strncpy(&c,buff+14+i,1);
			if(c<='z'&&c>='a'){
				strncpy(sol+i,buff+14+i,1);
			}else{
				*(sol+i)='\0';
				for(j=0;j<3;j++){
					if(!strcmp(solname[j],sol)){
						*lastresult=j;
                        i=10;
						break;
					}
				}
			}
		}
	}
    
	if ((p=strstr(buff,"vx(m/s)"))) {
		opt->posf=SOLF_RVA;
		opt->degf=0;
		strncpy(opt->sep,p+9,1);        // 要チェック
		opt->sep[1]='\0';
	}
	else if ((p=strstr(buff,"x-ecef(m)"))) {
		opt->posf=SOLF_XYZ;
		opt->degf=0;
		strncpy(opt->sep,p+9,1);
		opt->sep[1]='\0';
	}
    else if ((p=strstr(buff,"latitude(d'\")"))) {
        opt->posf=SOLF_LLH;
        opt->degf=1;
        strncpy(opt->sep,p+14,1);
        opt->sep[1]='\0';
    }
    else if ((p=strstr(buff,"latitude(deg)"))) {
        opt->posf=SOLF_LLH;
        opt->degf=0;
        strncpy(opt->sep,p+13,1);
        opt->sep[1]='\0';
    }
    else if ((p=strstr(buff,"e-baseline(m)"))) {
        opt->posf=SOLF_ENU;
        opt->degf=0;
        strncpy(opt->sep,p+13,1);
        opt->sep[1]='\0';
    }
    else if ((p=strstr(buff,"+SITE/INF"))) { /* gsi f2/f3 solution */
        opt->times=TIMES_GPST;
        opt->posf=SOLF_GSIF;
        opt->degf=0;
        strcpy(opt->sep," ");
    }
}

/* read solution option ------------------------------------------------------*/
static void readsolopt(FILE *fp, solopt_t *opt, char *posmod, int *lr)
{
    char buff[MAXSOLMSG+1];
    int i;
    
    trace(3,"readsolopt:\n");
    
    for (i=0;fgets(buff,sizeof(buff),fp)&&i<100;i++) { /* only 100 lines */
        
        /* decode solution options */
        decode_solopt(buff,opt,posmod,lr);
    }
}

/* input solution data from stream ---------------------------------------------
* input solution data from stream
* args   : unsigned char data I stream data
*          gtime_t ts       I  start time (ts.time==0: from start)
*          gtime_t te       I  end time   (te.time==0: to end)
*          double tint      I  time interval (0: all)
*          int    qflag     I  quality flag  (0: all)
*          solbuf_t *solbuf IO solution buffer
* return : status (1:solution received,0:no solution,-1:disconnect received)
*-----------------------------------------------------------------------------*/
extern int inputsol(unsigned char data, gtime_t ts, gtime_t te, double tint,
                    int qflag, const solopt_t *opt, solbuf_t *solbuf)
{
    sol_t sol={{0}};
    int stat;
    
    trace(4,"inputsol: data=0x%02x\n",data);
    
    sol.time=solbuf->time;
    
    if (data=='$'||(!isprint(data)&&data!='\r'&&data!='\n')) { /* sync header */
        solbuf->nb=0;
    }
    solbuf->buff[solbuf->nb++]=data;
    if (data!='\n'&&solbuf->nb<MAXSOLMSG) return 0; /* sync trailer */
    
    solbuf->buff[solbuf->nb]='\0';
    solbuf->nb=0;
    
    /* check disconnect message */
    if (!strcmp((char *)solbuf->buff,MSG_DISCONN)) {
        trace(3,"disconnect received\n");
        return -1;
    }
    /* decode solution */
    if ((stat=decode_sol((char *)solbuf->buff,opt,&sol,solbuf->rb))>0) {
        solbuf->time=sol.time; /* update current time */
    }
    if (stat!=1||!screent(sol.time,ts,te,tint)||(qflag&&sol.stat!=qflag)) {
        return 0;
    }
    /* add solution to solution buffer */
    return addsol(solbuf,&sol);
}
/* read solution data --------------------------------------------------------*/
static int readsoldata(FILE *fp, gtime_t ts, gtime_t te, double tint, int qflag,
                      const solopt_t *opt, solbuf_t *solbuf)
{
    int c;
    
    trace(3,"readsoldata:\n");
    
    while ((c=fgetc(fp))!=EOF) {
        
        /* input solution */
        inputsol((unsigned char)c,ts,te,tint,qflag,opt,solbuf);
    }
    return solbuf->n>0;
}
/* compare solution data -----------------------------------------------------*/
static int cmpsol(const void *p1, const void *p2)
{
    sol_t *q1=(sol_t *)p1,*q2=(sol_t *)p2;
    double tt=timediff(q1->time,q2->time);
    return tt<-0.0?-1:(tt>0.0?1:0);
}
/* sort solution data --------------------------------------------------------*/
static int sort_solbuf(solbuf_t *solbuf)
{
    sol_t *solbuf_data;
    
    trace(4,"sort_solbuf: n=%d\n",solbuf->n);
    
    if (solbuf->n<=0) return 0;
    
    if (!(solbuf_data=(sol_t *)realloc(solbuf->data,sizeof(sol_t)*solbuf->n))) {
        trace(1,"sort_solbuf: memory allocation error\n");
        free(solbuf->data); solbuf->data=NULL; solbuf->n=solbuf->nmax=0;
        return 0;
    }
    solbuf->data=solbuf_data;
    qsort(solbuf->data,solbuf->n,sizeof(sol_t),cmpsol);
    solbuf->nmax=solbuf->n;
    solbuf->start=0;
    solbuf->end=solbuf->n-1;
    return 1;
}
/* read solutions data from solution files -------------------------------------
* read solution data from soluiton files
* args   : char   *files[]  I  solution files
*          int    nfile     I  number of files
*         (gtime_t ts)      I  start time (ts.time==0: from start)
*         (gtime_t te)      I  end time   (te.time==0: to end)
*         (double tint)     I  time interval (0: all)
*         (int    qflag)    I  quality flag  (0: all)
*          solbuf_t *solbuf O  solution buffer
* return : status (1:ok,0:no data or error)
*-----------------------------------------------------------------------------*/
extern int readsolt(char *files[], int nfile, gtime_t ts, gtime_t te,
					double tint, int qflag, solbuf_t *solbuf, char *posmod, int *lr)
{
    FILE *fp;
    solopt_t opt=solopt_default;
    int i;
    
    trace(3,"readsolt: nfile=%d\n",nfile);
    
    initsolbuf(solbuf,0,0);
    
    for (i=0;i<nfile;i++) {
        if (!(fp=fopen(files[i],"rb"))) {
            trace(1,"readsolt: file open error %s\n",files[i]);
            continue;
        }
        /* read solution options in header */
        readsolopt(fp,&opt,posmod,lr);
        rewind(fp);
        
        /* read solution data */
        if (!readsoldata(fp,ts,te,tint,qflag,&opt,solbuf)) {
            trace(1,"readsolt: no solution in %s\n",files[i]);
        }
        fclose(fp);
    }
    return sort_solbuf(solbuf);
}

extern int readsol(char *files[], int nfile, solbuf_t *sol)
{
    gtime_t time={0};
    char *posmod={0};
    int *lr={0};
    trace(3,"readsol: nfile=%d\n",nfile);
    
    return readsolt(files,nfile,time,time,0.0,0,sol,posmod,lr);
}

/* add solution data to solution buffer ----------------------------------------
* add solution data to solution buffer
* args   : solbuf_t *solbuf IO solution buffer
*          sol_t  *sol      I  solution data
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int addsol(solbuf_t *solbuf, const sol_t *sol)
{
    sol_t *solbuf_data;
    
    trace(4,"addsol:\n");
    
    if (solbuf->cyclic) { /* ring buffer */
        if (solbuf->nmax<=1) return 0;
        solbuf->data[solbuf->end]=*sol;
        if (++solbuf->end>=solbuf->nmax) solbuf->end=0;
        if (solbuf->start==solbuf->end) {
            if (++solbuf->start>=solbuf->nmax) solbuf->start=0;
        }
        else solbuf->n++;
        
        return 1;
    }
    if (solbuf->n>=solbuf->nmax) {
        solbuf->nmax=solbuf->nmax==0?8192:solbuf->nmax*2;
        if (!(solbuf_data=(sol_t *)realloc(solbuf->data,sizeof(sol_t)*solbuf->nmax))) {
            trace(1,"addsol: memory allocation error\n");
            free(solbuf->data); solbuf->data=NULL; solbuf->n=solbuf->nmax=0;
            return 0;
        }
        solbuf->data=solbuf_data;
    }
    solbuf->data[solbuf->n++]=*sol;
    return 1;
}
/* get solution data from solution buffer --------------------------------------
* get solution data by index from solution buffer
* args   : solbuf_t *solbuf I  solution buffer
*          int    index     I  index of solution (0...)
* return : solution data pointer (NULL: no solution, out of range)
*-----------------------------------------------------------------------------*/
extern sol_t *getsol(solbuf_t *solbuf, int index)
{
    trace(4,"getsol: index=%d\n",index);
    
    if (index<0||solbuf->n<=index) return NULL;
    if ((index=solbuf->start+index)>=solbuf->nmax) {
        index-=solbuf->nmax;
    }
    return solbuf->data+index;
}
/* initialize solution buffer --------------------------------------------------
* initialize position solutions
* args   : solbuf_t *solbuf I  solution buffer
*          int    cyclic    I  solution data buffer type (0:linear,1:cyclic)
*          int    nmax      I  initial number of solution data
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern void initsolbuf(solbuf_t *solbuf, int cyclic, int nmax)
{
    gtime_t time0={0};
    
    trace(3,"initsolbuf: cyclic=%d nmax=%d\n",cyclic,nmax);
    
    solbuf->n=solbuf->nmax=solbuf->start=solbuf->end=0;
    solbuf->cyclic=cyclic;
    solbuf->time=time0;
    solbuf->data=NULL;
    if (cyclic) {
        if (nmax<=2) nmax=2;
        if (!(solbuf->data=malloc(sizeof(sol_t)*nmax))) {
            trace(1,"initsolbuf: memory allocation error\n");
            return;
        }
        solbuf->nmax=nmax;
    }
}
/* free solution ---------------------------------------------------------------
* free memory for solution buffer
* args   : solbuf_t *solbuf I  solution buffer
* return : none
*-----------------------------------------------------------------------------*/
extern void freesolbuf(solbuf_t *solbuf)
{
    trace(3,"freesolbuf: n=%d\n",solbuf->n);
    
    free(solbuf->data);
    solbuf->n=solbuf->nmax=solbuf->start=solbuf->end=0;
    solbuf->data=NULL;
}
extern void freesolstatbuf(solstatbuf_t *solstatbuf)
{
    trace(3,"freesolstatbuf: n=%d\n",solstatbuf->n);
    
    solstatbuf->n=solstatbuf->nmax=0;
    free(solstatbuf->data);
    solstatbuf->data=NULL;
}
/* compare solution status ---------------------------------------------------*/
static int cmpsolstat(const void *p1, const void *p2)
{
    solstat_t *q1=(solstat_t *)p1,*q2=(solstat_t *)p2;
    double tt=timediff(q1->time,q2->time);
    return tt<-0.0?-1:(tt>0.0?1:0);
}
/* sort solution data --------------------------------------------------------*/
static int sort_solstat(solstatbuf_t *statbuf)
{
    solstat_t *statbuf_data;
    
    trace(4,"sort_solstat: n=%d\n",statbuf->n);
    
    if (statbuf->n<=0) return 0;
    
    if (!(statbuf_data=realloc(statbuf->data,sizeof(solstat_t)*statbuf->n))) {
        trace(1,"sort_solstat: memory allocation error\n");
        free(statbuf->data); statbuf->data=NULL; statbuf->n=statbuf->nmax=0;
        return 0;
    }
    statbuf->data=statbuf_data;
    qsort(statbuf->data,statbuf->n,sizeof(solstat_t),cmpsolstat);
    statbuf->nmax=statbuf->n;
    return 1;
}
/* decode solution status ----------------------------------------------------*/
static int decode_solstat(char *buff, solstat_t *stat)
{
    static const solstat_t stat0={{0}};
    double tow,az,el,resp,resc;
    int n,week,sat,frq,code=0,vsat,snr,fix,slip,lock,outc,slipc,rejc;
    char id[32]="",*p;
    int prmc;

	trace(4,"decode_solstat: buff=%s\n",buff);

	if (strstr(buff,"$SAT")!=buff) return 0;

	for (p=buff,prmc=1;*p;p++) {
		if (*p==',') {
			*p=' ';
			++prmc;
		}
	}
	if(18==prmc) {
		n=sscanf(buff,"$SAT%d%lf%s%d%d%lf%lf%lf%lf%d%d%d%d%d%d%d%d",
			&week,&tow,id,&frq,&code,&az,&el,&resp,&resc,&vsat,&snr,&fix,&slip,
			&lock,&outc,&slipc,&rejc);
	}
	else {
		n=sscanf(buff,"$SAT%d%lf%s%d%lf%lf%lf%lf%d%d%d%d%d%d%d%d",
			&week,&tow,id,&frq,&az,&el,&resp,&resc,&vsat,&snr,&fix,&slip,
			&lock,&outc,&slipc,&rejc);
	}
    
    if (n<15) {
        trace(2,"invalid format of solution status: %s\n",buff);
        return 0;
    }
    if ((sat=satid2no(id))<=0) {
        trace(2,"invalid satellite in solution status: %s\n",id);
        return 0;
    }
    *stat=stat0;
    stat->time=gpst2time(week,tow);
    stat->sat  =(unsigned char)sat;
    stat->frq  =(unsigned char)frq;
	stat->code =(unsigned char)code;
    stat->az   =(float)(az*D2R);
    stat->el   =(float)(el*D2R);
    stat->resp =(float)resp;
    stat->resc =(float)resc;
    stat->flag =(unsigned char)((vsat<<5)+(slip<<3)+fix);
    stat->snr  =(unsigned char)(snr*4.0+0.5);
    stat->lock =(unsigned short)lock;
    stat->outc =(unsigned short)outc;
    stat->slipc=(unsigned short)slipc;
    stat->rejc =(unsigned short)rejc;
    return 1;
}
/* add solution status data --------------------------------------------------*/
static void addsolstat(solstatbuf_t *statbuf, const solstat_t *stat)
{
    solstat_t *statbuf_data;
    
    trace(4,"addsolstat:\n");
    
    if (statbuf->n>=statbuf->nmax) {
        statbuf->nmax=statbuf->nmax==0?8192:statbuf->nmax*2;
        if (!(statbuf_data=(solstat_t *)realloc(statbuf->data,sizeof(solstat_t)*
                                                statbuf->nmax))) {
            trace(1,"addsolstat: memory allocation error\n");
            free(statbuf->data); statbuf->data=NULL; statbuf->n=statbuf->nmax=0;
            return;
        }
        statbuf->data=statbuf_data;
    }
    statbuf->data[statbuf->n++]=*stat;
}

/* read solution status data -------------------------------------------------*/
static int readsolstatdata(FILE *fp, gtime_t ts, gtime_t te, double tint,
                           solstatbuf_t *statbuf)
{
    solstat_t stat={{0}};
    char buff[MAXSOLMSG+1];
    
    trace(3,"readsolstatdata:\n");
    
    while (fgets(buff,sizeof(buff),fp)) {
        
        /* decode solution status */
        if (!decode_solstat(buff,&stat)) continue;
        
        /* add solution to solution buffer */
        if (screent(stat.time,ts,te,tint)) {
            addsolstat(statbuf,&stat);
        }
    }
    return statbuf->n>0;
}
/* read solution status --------------------------------------------------------
* read solution status from solution status files
* args   : char   *files[]  I  solution status files
*          int    nfile     I  number of files
*         (gtime_t ts)      I  start time (ts.time==0: from start)
*         (gtime_t te)      I  end time   (te.time==0: to end)
*         (double tint)     I  time interval (0: all)
*          solstatbuf_t *statbuf O  solution status buffer
* return : status (1:ok,0:no data or error)
*-----------------------------------------------------------------------------*/
extern int readsolstatt(char *files[], int nfile, gtime_t ts, gtime_t te,
                        double tint, solstatbuf_t *statbuf)
{
    FILE *fp;
    char path[1024];
    int i;
    
    trace(3,"readsolstatt: nfile=%d\n",nfile);
    
    statbuf->n=statbuf->nmax=0;
    statbuf->data=NULL;
    
    for (i=0;i<nfile;i++) {
        sprintf(path,"%s.stat",files[i]);
        if (!(fp=fopen(path,"r"))) {
            trace(1,"readsolstatt: file open error %s\n",path);
            continue;
        }
        /* read solution status data */
        if (!readsolstatdata(fp,ts,te,tint,statbuf)) {
            trace(1,"readsolt: no solution in %s\n",path);
        }
        fclose(fp);
    }
    return sort_solstat(statbuf);
}
extern int readsolstat(char *files[], int nfile, solstatbuf_t *statbuf)
{
    gtime_t time={0};
    
    trace(3,"readsolstat: nfile=%d\n",nfile);
    
    return readsolstatt(files,nfile,time,time,0.0,statbuf);
}

