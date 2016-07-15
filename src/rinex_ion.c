/*------------------------------------------------------------------------------
* rinex_ion.c : rinex(ionospheric delay) functions
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

static FILE *fp_ion=NULL;        /* ion rinex file pointer */
static FILE *fp_ion_data=NULL;   /* ion rinex file pointer */
static int us[MAXSAT]={0};
static prcopt_t *popt=NULL;
static solopt_t *sopt=NULL;
static const char fniondata[]={"irnx.data"};
static gtime_t tstart;
static gtime_t tend;
static int ini=0;

static rnxopt_t ropt;

extern int openrnxobsh_ion(const prcopt_t *p_opt, const char *file)
{
	int i,j;
	char date[21];

    popt=p_opt;

	pctime_rnx(date);

	ropt.rnxver=3.00;
	ropt.navsys=popt->navsys;
	ropt.obstype;        /* observation type */
	ropt.freqtype;       /* frequency type */
	ropt.mask[6][64];   /* code mask {GPS,GLO,GAL,QZS,SBS,CMP} */
	ropt.staid [32];    /* station id for rinex file name */
	strcpy(ropt.prog,PROGRAM_NAME_VER);
	strcpy(ropt.runby,date);
	strcpy(ropt.marker,stas[0].name);
	strcpy(ropt.markerno,stas[0].marker);
	strcpy(ropt.name[0],stas[0].observer);
	strcpy(ropt.name[1],stas[0].agency);

	strcpy(ropt.rec[0],stas[0].recsno);
	strcpy(ropt.rec[1],stas[0].rectype);
	strcpy(ropt.rec[2],stas[0].recver);

	strcpy(ropt.ant[0],stas[0].antsno);
	strcpy(ropt.ant[1],stas[0].antdes);
	strcpy(ropt.ant[2],"");
	for(i=0;i<3;++i) ropt.apppos[i] = stas[0].pos[i];
	for(i=0;i<3;++i) ropt.antdel[i]= stas[0].del[i];

	for(i=0;i<MAXCOMMENT;++i) {
		ropt.comment[i][0]=NULL; /* comments */
	}
	ropt.outiono=0;        /* output iono correction */
	ropt.outtime=0;        /* output time system correction */
	ropt.outleaps=0;       /* output leap seconds */

	if(ropt.navsys&SYS_GPS) {
		ropt.nobs[0]=1;
		for(j=0;j<ropt.nobs[0];++j) {
			sprintf(ropt.tobs[0][j],"I%1d\0", j+1);
		}
	}
	if(ropt.navsys&SYS_GLO) {
		ropt.nobs[1]=1;
		for(j=0;j<ropt.nobs[1];++j) {
			sprintf(ropt.tobs[1][j],"I%1d\0", j+1);
		}
	}
	if(ropt.navsys&SYS_GAL) {
		ropt.nobs[2]=1;
		for(j=0;j<ropt.nobs[2];++j) {
			sprintf(ropt.tobs[2][j],"I%1d\0", j+1);
		}
	}
	if(ropt.navsys&SYS_QZS) {
		ropt.nobs[i]=1;
		for(j=0;j<ropt.nobs[3];++j) {
			sprintf(ropt.tobs[3][j],"I%1d\0", j+1);
		}
	}
	if(ropt.navsys&SYS_SBS) {
		ropt.nobs[4]=1;
		for(j=0;j<ropt.nobs[4];++j) {
			sprintf(ropt.tobs[4][j],"I%1d\0", j+1);
		}
	}
	if(ropt.navsys&SYS_CMP) {
		ropt.nobs[5]=1;
		for(j=0;j<ropt.nobs[5];++j) {
			sprintf(ropt.tobs[5][j],"I%1d\0", j+1);
		}
	}

	if (!(fp_ion=fopen(file,"w"))) {
		trace(1,"rtkopenirinex: file open error file=%s\n",file);
		return 0;
	}

	if (!(fp_ion_data=fopen(fniondata,"w+"))) {
		trace(1,"rtkopenirinex: file open error file=%s\n",fniondata);
		fclose(fp_ion);
		fp_ion=NULL;
		return 0;
	}

	ini=0;

	return 1;
}

extern int closernxobsh_ion()
{
	char buf[256]={0};

	if (fp_ion == NULL) return 0;

	outrnxobsh(fp_ion,&ropt,NULL);


	if (fp_ion_data) {
		rewind(fp_ion_data);
		while(NULL!=fgets(buf,sizeof(buf), fp_ion_data)) {
			fprintf(fp_ion, "%s", buf);
		}
		fclose(fp_ion_data);
	}

	if (fp_ion) fclose(fp_ion);
	popt=NULL;

	remove(fniondata);

	return 1;
}

extern int outrnxobsb_ion(const obsd_t *obs, const nav_t *nav, const double *ion, int n, int flag)
{
    const char *mask;
    double ep[6];
    char sats[MAXOBS][4]={""};
    int i,j,k,m,ns,sys,ind[MAXOBS],s[MAXOBS]={0};

	if (ropt.rnxver<=2.99) { /* ver.2 */
		return 0;
	}

	trace(3,"outrnxobsb: n=%d\n",n);

	time2epoch(obs[0].time,ep);

	for (i=ns=0;i<n&&ns<MAXOBS;i++) {
		sys=satsys(obs[i].sat,NULL);
		if (!(sys&ropt.navsys)||ropt.exsats[obs[i].sat-1]) continue;
		if (!sat2code(obs[i].sat,sats[ns])) continue;
		switch (sys) {
			case SYS_GPS: s[ns]=0; break;
			case SYS_GLO: s[ns]=1; break;
			case SYS_GAL: s[ns]=2; break;
			case SYS_QZS: s[ns]=3; break;
			case SYS_SBS: s[ns]=4; break;
			case SYS_CMP: s[ns]=5; break;
		}
		if (!ropt.nobs[s[ns]]) continue;
		ind[ns++]=i;
	}
	fprintf(fp_ion_data,"> %04.0f %2.0f %2.0f %2.0f %2.0f%11.7f  %d%3d%21s\n",
		ep[0],ep[1],ep[2],ep[3],ep[4],ep[5],flag,ns,"");
	for (i=0;i<ns;i++) {
		sys=satsys(obs[ind[i]].sat,NULL);
		fprintf(fp_ion_data,"%-3s",sats[i]);
		m=s[i];
		mask=ropt.mask[s[i]];
		outrnxobsf(fp_ion_data, ion[obs[ind[i]].sat-1]/nav->lam[obs[ind[i]].sat-1][0], -1);

		if (fprintf(fp_ion_data,"\n")==EOF) return 0;
	}

	if(ini==0) {
		ropt.tstart=obs->time;
		ini=1;
	}
	ropt.tend=obs->time;
	return 1;
}
