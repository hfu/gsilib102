/*------------------------------------------------------------------------------
* postcmn.c : post-processing positioning common functions
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

#include "postpos.h"

static const char rcsid[]="$Id:";

pcvs_t pcvss={0};        /* receiver antenna parameters */
pcvs_t pcvsr={0};        /* satellite antenna parameters */
obs_t obss={0};          /* observation data */
nav_t navs={0};          /* navigation data */
sbs_t sbss={0};          /* sbas messages */
lex_t lexs={0};          /* lex messages */
int nepoch=0;            /* number of observation epochs */
int iobsu =0;            /* current rover observation data index */
int iobsr =0;            /* current reference observation data index */
int isbs  =0;            /* current sbas message index */
int ilex  =0;            /* current lex message index */
int aborts=0;            /* abort status */
char proc_rov [64]="";   /* rover for current processing */
char proc_base[64]="";   /* base station for current processing */
char rtcm_file[1024]=""; /* rtcm data file */
char rtcm_path[1024]=""; /* rtcm data path */
rtcm_t rtcm;             /* rtcm control struct */
FILE *fp_rtcm=NULL;      /* rtcm data file pointer */
double diff_se=0;
int posepoch=0;            /* Poss. # of obs epochs */
gtime_t obsts;					/* obs start time */
gtime_t obste;					/* obs end time */
//int revs=0;            /* analysis direction (0:forward,1:backward) */

/* search next observation data index ----------------------------------------*/
extern int nextobsf(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;
    
    for (;*i<obs->n;(*i)++) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i+n<obs->n;n++) {
        tt=timediff(obs->data[*i+n].time,obs->data[*i].time);
        if (obs->data[*i+n].rcv!=rcv||tt>DTTOL) break;
    }
    return n;
}

extern int nextobsb(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;
    
    for (;*i>=0;(*i)--) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i-n>=0;n++) {
        tt=timediff(obs->data[*i-n].time,obs->data[*i].time);
        if (obs->data[*i-n].rcv!=rcv||tt<-DTTOL) break;
    }
    return n;
}

/* input obs data, navigation messages and sbas correction -------------------*/
extern int inputobs(obsd_t *obs, int solq, const prcopt_t *popt)
{
    gtime_t time={0};
    char path[1024];
    int i,nu,nr,n=0;
    
    trace(3,"infunc  : revs=%d iobsu=%d iobsr=%d isbs=%d\n",revs,iobsu,iobsr,isbs);
    
    if (0<=iobsu&&iobsu<obss.n) {
        settime((time=obss.data[iobsu].time));
        if (checkbrk("processing : %s Q=%d",time_str(time,0),solq)) {
            aborts=1; showmsg("aborted"); return -1;
        }
    }
    if (!revs) { /* input forward data */
        if ((nu=nextobsf(&obss,&iobsu,1))<=0) return -1;
        if (popt->intpref) {
            for (;(nr=nextobsf(&obss,&iobsr,2))>0;iobsr+=nr)
                if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)>-DTTOL) break;
        }
        else {
            for (i=iobsr;(nr=nextobsf(&obss,&i,2))>0;iobsr=i,i+=nr)
                if (timediff(obss.data[i].time,obss.data[iobsu].time)>DTTOL) break;
        }
        nr=nextobsf(&obss,&iobsr,2);
        for (i=0;i<nu&&n<MAXOBS;i++) obs[n++]=obss.data[iobsu+i];
        for (i=0;i<nr&&n<MAXOBS;i++) obs[n++]=obss.data[iobsr+i];
        iobsu+=nu;
        
        /* update sbas corrections */
        while (isbs<sbss.n) {
			time=gpst2time(sbss.msgs[isbs].week,(double)sbss.msgs[isbs].tow);
            
            if (getbitu(sbss.msgs[isbs].msg,8,6)!=9) { /* except for geo nav */
                sbsupdatecorr(sbss.msgs+isbs,&navs);
            }
            if (timediff(time,obs[0].time)>-1.0-DTTOL) break;
            isbs++;
        }
        /* update lex corrections */
        while (ilex<lexs.n) {
            if (lexupdatecorr(lexs.msgs+ilex,&navs,&time)) {
                if (timediff(time,obs[0].time)>-1.0-DTTOL) break;
            }
            ilex++;
        }
        /* update rtcm corrections */
        if (*rtcm_file) {
            
            /* open or swap rtcm file */
            reppath(rtcm_file,path,obs[0].time,"","");
            
            if (strcmp(path,rtcm_path)) {
                strcpy(rtcm_path,path);
                
                if (fp_rtcm) fclose(fp_rtcm);
                fp_rtcm=fopen(path,"rb");
                if (fp_rtcm) {
                    rtcm.time=obs[0].time;
                    input_rtcm3f(&rtcm,fp_rtcm);
                    trace(2,"rtcm file open: %s\n",path);
                }
            }
            if (fp_rtcm) {
                while (timediff(rtcm.time,obs[0].time)<0.0) {
                    if (input_rtcm3f(&rtcm,fp_rtcm)<-1) break;
                }
                for (i=0;i<MAXSAT;i++) navs.ssr[i]=rtcm.ssr[i];
            }
        }
    }
    else { /* input backward data */
        if ((nu=nextobsb(&obss,&iobsu,1))<=0) return -1;
        if (popt->intpref) {
            for (;(nr=nextobsb(&obss,&iobsr,2))>0;iobsr-=nr)
                if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)<DTTOL) break;
        }
        else {
            for (i=iobsr;(nr=nextobsb(&obss,&i,2))>0;iobsr=i,i-=nr)
                if (timediff(obss.data[i].time,obss.data[iobsu].time)<-DTTOL) break;
        }
        nr=nextobsb(&obss,&iobsr,2);
        for (i=0;i<nu&&n<MAXOBS;i++) obs[n++]=obss.data[iobsu-nu+1+i];
        for (i=0;i<nr&&n<MAXOBS;i++) obs[n++]=obss.data[iobsr-nr+1+i];
        iobsu-=nu;
        
        /* update sbas corrections */
        while (isbs>=0) {
            time=gpst2time(sbss.msgs[isbs].week,(double)sbss.msgs[isbs].tow);
            
            if (getbitu(sbss.msgs[isbs].msg,8,6)!=9) { /* except for geo nav */
                sbsupdatecorr(sbss.msgs+isbs,&navs);
            }
            if (timediff(time,obs[0].time)<1.0+DTTOL) break;
            isbs--;
        }
        /* update lex corrections */
        while (ilex>=0) {
            if (lexupdatecorr(lexs.msgs+ilex,&navs,&time)) {
                if (timediff(time,obs[0].time)<1.0+DTTOL) break;
            }
            ilex--;
        }
    }
    return n;
}


/* read obs and nav data -----------------------------------------------------*/
extern int readobsnav(gtime_t ts, gtime_t te, double ti, char **infile,
                      const int *index, int n, const prcopt_t *prcopt,
					  obs_t *obs, nav_t *nav, sta_t *sta)
{
    int i,j,ind=0,nobs=0,rcv=1;

	gtime_t obss;					/* obs start time */
	gtime_t obse;					/* obs end time */

	trace(3,"readobsnav: ts=%s n=%d\n",time_str(ts,0),n);

	obs->data=NULL; obs->n =obs->nmax =0;
	nav->eph =NULL; nav->n =nav->nmax =0;
	nav->geph=NULL; nav->ng=nav->ngmax=0;
	if(nav->seph){free(nav->seph);} nav->seph=NULL; nav->ns=nav->nsmax=0;
	nepoch=0;
    
    for (i=0;i<n;i++) {
        if (checkbrk("")) return 0;
        
        if (index[i]!=ind) {
            if (obs->n>nobs) rcv++;
            ind=index[i]; nobs=obs->n; 
        }
        /* read rinex obs and nav file */
        if (readrnxt(infile[i],rcv,ts,te,ti,prcopt->rnxopt[rcv<=1?0:1],obs,nav,
                     rcv<=2?sta+rcv-1:NULL)<0) {
            checkbrk("error : insufficient memory");
            trace(1,"insufficient memory\n");
            return 0;
        }
    }
    if (obs->n<=0) {
        checkbrk("error : no obs data");
        trace(1,"no obs data\n");
        return 0;
    }
    if (nav->n<=0&&nav->ng<=0&&nav->ns<=0) {
        checkbrk("error : no nav data");
        trace(1,"no nav data\n");
        return 0;
    }
    /* sort observation data */
    nepoch=sortobs(obs);
    
    /* delete duplicated ephemeris */
    uniqnav(nav);
    
    /* set time span for progress display */
//    if (ts.time==0||te.time==0) {
		for (i=0;   i<obs->n;i++) if (obs->data[i].rcv==1) break;
		for (j=obs->n-1;j>=0;j--) if (obs->data[j].rcv==1) break;
		if (i<j) {
			if     (ts.time==0)                       obss=obs->data[i].time;
			else if(obs->data[i].time.time < ts.time) obss=ts;
			else                                      obss=obs->data[i].time;

			if     (te.time==0)                       obse=obs->data[j].time;
			else if(te.time < obs->data[j].time.time) obse=te;
			else                                      obse=obs->data[j].time;

			settspan(obss, obse);

			diff_se = settspan_d(obse,obss);
			posepoch = (int)(diff_se/interval) + 1;
		}
//	}
	return 1;
}

/* free obs and nav data -----------------------------------------------------*/
extern void freeobsnav(obs_t *obs, nav_t *nav)
{
	int i;
    trace(3,"freeobsnav:\n");
	free(obs->data); obs->data=NULL; obs->n =obs->nmax =0;
    free(nav->eph ); nav->eph =NULL; nav->n =nav->nmax =0;
    free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
	free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
    if (nav->ifbs) {free(nav->ifbs ); nav->ifbs =NULL; nav->nifbs=nav->nifbsmax=0;}
	if (nav->errs) {free(nav->errs ); nav->errs =NULL; nav->nerrs=nav->nerrsmax=0;}
	if (nav->dcb) {free(nav->dcb ); nav->dcb =NULL; nav->nd=nav->ndmax=0;}
	if (nav->bipm) {free(nav->bipm); nav->bipm=NULL;nav->nbipm=nav->nbipmmax=0;}
	for (i=0;i<MAXSAT;i++) {
		nav->cbias[i][0]=0.0;
		nav->cbias[i][1]=0.0;
		nav->cbias[i][2]=0.0;
	}

}

/* show message and check break ----------------------------------------------*/
extern int checkbrk(const char *format, ...)
{
    va_list arg;
    char buff[1024],*p=buff;
    if (!*format) return showmsg("");
    va_start(arg,format);
    p+=vsprintf(p,format,arg);
    va_end(arg);
    if (*proc_rov&&*proc_base) sprintf(p," (%s-%s)",proc_rov,proc_base);
    else if (*proc_rov ) sprintf(p," (%s)",proc_rov );
    else if (*proc_base) sprintf(p," (%s)",proc_base);
    return showmsg(buff);
}

/* open output file for append -----------------------------------------------*/
extern FILE *openfile(const char *outfile)
{
    trace(3,"openfile: outfile=%s\n",outfile);
    
    return !*outfile?stdout:fopen(outfile,"a");
}