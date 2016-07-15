/*------------------------------------------------------------------------------
* crinex.c : clock rinex functions
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

static FILE *fp_crx=NULL;        /* clock rinex file pointer */
static FILE *fp_crx_data=NULL;   /* clock rinex file pointer */
static int us[MAXSAT]={0};
static prcopt_t *popt=NULL;
static solopt_t *sopt=NULL;
static mbs_t *mbs=NULL;
static sta_t *sta=NULL;

extern int rtkopencrinex(const prcopt_t *p_opt, const solopt_t *s_opt, const mbs_t *mbs0, const char *file, const sta_t* sta0)
{
	char s[64];
	int ntype;
	char date[21];
	long indHead;
	int i;
	double d;

	if(!(s_opt->recclout)) return 1;

	trace(3,"rtkopencrinex: file=%s\n",file);

	if (!(fp_crx=fopen(file,"w"))) {
		trace(1,"rtkopencrinex: file open error file=%s\n",file);
		return 0;
	}
	if (!(fp_crx_data=tmpfile())) {
		trace(1,"rtkopencrinex: file open error file=%s\n",file);
		fclose(fp_crx);
		fp_crx = NULL;
		return 0;
	}

	popt = p_opt;
	sopt = s_opt;
	mbs = mbs0;
	sta = sta0;

	indHead = ftell(fp_crx);
	time2str(utc2gpst(pctime()),s,0);
	pctime_rnx(date);

	fprintf(fp_crx,"     3.00           CLOCK DATA          Mixed               RINEX VERSION / TYPE\n");
	fprintf(fp_crx,"%-20sGSI                 %-20sPGM / RUN BY / DATE\n", PROGRAM_NAME_VER, date);
	fprintf(fp_crx,"   GPS                                                      TIME SYSTEM ID\n");

	if((popt->mode == PMODE_MULTI) && (sopt->recclout) && (sopt->satclout)) {
		fprintf(fp_crx,"     2    AR    AS                                          # / TYPES OF DATA\n");
	}
	else if((popt->mode == PMODE_MULTI) && (sopt->satclout)) {
		fprintf(fp_crx,"     1    AS                                                # / TYPES OF DATA\n");
	}
	else if(sopt->recclout) {
		fprintf(fp_crx,"     1    AR                                                # / TYPES OF DATA\n");
	}

	fprintf(fp_crx,"%-3s  %-55sANALYSIS CENTER\n",ANAL_CENTER_ABB, ANAL_CENTER_FULL);

	fflush(fp_crx);
	return 1;

}
extern void rtkclosecrinex(void)
{
	int i;
	int ns=0;

	char id[32];
	char buf[256]={0};

	if (!fp_crx) return;

	trace(3,"rtkclosecrinex:\n");

	if((popt->mode == PMODE_MULTI) && (mbs)) {
		fprintf(fp_crx,"%6d    %-50s# OF SOLN STA / TRF\n", mbs->stas.nsta, "");
		for(i=0;i<mbs->stas.nsta;i++) {
			fprintf(fp_crx,"%-4s %20s%11.0f %11.f %11.0fSOLN STA NAME / NUM\n"
				, mbs->stas.sta[i].name2,""
				, mbs->stas.sta[i].pos[0]*1000+0.5
				, mbs->stas.sta[i].pos[1]*1000+0.5
				, mbs->stas.sta[i].pos[2]*1000+0.5);
		}
	}
	else {
		fprintf(fp_crx,"%6d    %-50s# OF SOLN STA / TRF\n", 1, "");
		fprintf(fp_crx,"%-4s %20s%11.0f %11.0f %11.0fSOLN STA NAME / NUM\n"
			, sta[0].name,""
			, sta[0].pos[0]*1000+0.5
			, sta[0].pos[1]*1000+0.5
			, sta[0].pos[2]*1000+0.5);
	}

	if((popt) && (sopt)) {
		if((popt->mode == PMODE_MULTI) && (sopt->satclout)) {
			for(i=0;i<MAXSAT;++i) {
				if(us[i]==1) ++ns;
			}
			fprintf(fp_crx,"%6d                                                      # OF SOLN SATS\n", ns);
            ns=0;
			for(i=0;i<MAXSAT;++i) {
				if(us[i]==1) {
                    ++ns;
					satno2id(i+1, id);
					fprintf(fp_crx,"%-3s ", id);
					if(ns%15==0) fprintf(fp_crx,"PRN LIST\n");
				}
			}
			if(ns%15!=0) fprintf(fp_crx,"%-*sPRN LIST\n", 60-4*(ns%15), "");
		}
	}
	fprintf(fp_crx,"                                                            END OF HEADER\n");

	if (fp_crx_data) {
		rewind(fp_crx_data);
		while(NULL!=fgets(buf,sizeof(buf), fp_crx_data)) {
			fprintf(fp_crx, "%s", buf);
		}
		fclose(fp_crx_data);
	}

	if (fp_crx) fclose(fp_crx);
	if (fp_crx_data) fclose(fp_crx_data);
	fp_crx=NULL;
	popt=NULL;
	sopt=NULL;
	mbs=NULL;
	sta=NULL;
}

extern void outsolcrinexrec(gtime_t *t, char *id, double clc)
{
	double ep[6];

	if (!fp_crx) return;

	time2epoch(*t,ep);

	fprintf(fp_crx_data, "AR %-4s %04.0f %02.0f %02.0f %02.0f %02.0f %09.6f  1   %19.12f\n"
		,id,ep[0],ep[1],ep[2],ep[3],ep[4],ep[5], clc);
	fflush(fp_crx);
}

extern void outsolcrinexsat(gtime_t *t, int sat, double clc)
{
	double ep[6];
	char id[32];

	us[sat-1]=1;

	if (!fp_crx) return;

	time2epoch(*t,ep);
	satno2id(sat, id);

	fprintf(fp_crx_data, "AS %-3s  %04.0f %02.0f %02.0f %02.0f %02.0f %09.6f  1   %19.12f\n"
		,id,ep[0],ep[1],ep[2],ep[3],ep[4],ep[5], clc);

}

