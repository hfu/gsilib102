/*------------------------------------------------------------------------------
* qc.c :
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

#include <stdarg.h>
#include "rtklib.h"
#include "postpos.h"
#include "gsiqc.h"

#define MAXFILE     8                   /* max number of input files */

#define QCDOOPT  "0:off,1:on"
#define POINT_AVG     1    /* Points in MP moving avg */
#define PLOT_WIDTH    72    /* Width of ASCII summary plot */
#define GAP         10.0    /* Report data gap greater than(min) */
#define HIST_VER      19    /* Minimum signal to noise for L2 */

static double tick=0;
static double total_h=0;
static int outflag=0;
//static gtime_t obsts;					/* obs start time */
//static gtime_t obste;					/* obs end time */


static char *MP_Name[]={"12", "21", "15", "51", "25", "52", "16", "61", "26", "62", "56", "65", 
						"15b", "5b1", "25b", "5b2", "55b", "5b5", "65b", "5b6", "15a", "5a1", 
						"25a", "5a2", "55a", "5a5", "65a", "5a6", "5b5a", "5a5b"};

//typedef struct {
//	int num;
//	int del;
//	int abcomp;
//	int belcomp;
//	double rms;
//	double ele;
//	double constmp;
//	double val[POINT_AVG];
//	double msec;
//	unsigned char slip;
//	unsigned char msecslip;
//	unsigned char calcflag;
//} sv_mp;

typedef struct {
	double inition;
	double beforeion;
	double beforeiod;
	gtime_t time;
	unsigned char slip;
	unsigned char comp;
} sv_ion;

typedef struct {
	int above;
	int below;
} freq_lli;

typedef struct {
	int hor;						/* above 0°(include no obs) */
	double sum_ele_hor;				/* total elevation (above 0°(include no obs)) */
	double ele_hor;					/* mean elevation (above 0°(include no obs)) */
	int mask;						/* above mask(include no obs) */
	double sum_ele_mask;			/* total elevation (above mask(include no obs)) */
	double ele_mask;				/* mean elevation (above mask(include no obs)) */
	int reprt;						/* above mask(with obs) */
	int compl;						/* above mask &&  code && phase */
	int *count;
	int belowmask;					/* below mask(with obs) */
	
	sv_ion ion[NFREQ-1];
	int iodslipab;
	int iodslipbe;
	int iodmpslipab;
	int iodmpslipbe;
	freq_lli lli[NFREQ+NEXOBS];
//	sv_mp mp[(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1)];
} sv_sum;

typedef struct {
	int n[NSYS];
	int slps[NSYS];
	double rms[NSYS];
} hist_mp;

//typedef struct {
//	hist_mp mp[(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1)];
//	int S1n;
//	double totalS1;
//	double meanS1;
//	int S2n;
//	double totalS2;
//	double meanS2;
//	double *valueS1;
//	double *valueS2;
//	double sds1;
//	double sds2;
//	int ionslps;
//	int ionn;
//	double ionrms;
//} hist_data;

typedef struct {
	int n[NSYS+1];
	int above[NSYS+1];
	int below[NSYS+1];
	double rms[NSYS+1];
	double ele[NSYS+1];
} satvalue_mp;

//typedef struct {
//	int totalsat;
//	int totalobs;
//	int totalmask[NSYS+1];
//	int totalcomp[NSYS+1];		/* abode mask &&  code && phase */
//	int totaldel;				/* below mask && (!code | !phase) */
//	int totalbelow[NSYS+1];
//	int totalbelowcomp;
//	int totaliodslipab;
//	int totaliodslipbe;
//	int totaliodmpslipab;
//	int totaliodmpslipbe;
//	int totalslipaab[NSYS];
//	satvalue_mp mp[(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1)];
//
//	freq_lli lli[NFREQ];
//
//	char* mpsys[NSYS+1];
//	unsigned char obssys[NSYS+1];
//	int nsys;
//	int clkslip;
//	int msecmp;
//	int totalmp;
//} satvalue;

sv_sum sv[MAXSAT]={0};


static qc_param param={
	400,
	10,
	4,
	50,
	0,
	0,
	10,
	25,
	0.01,
};

opt_t qcopts[]={
	{"max_iod",             1, (void *)&param.iodsign, "cm/min", 1},
	{"data_gap",            1, (void *)&param.gap,        "min", 1},
	{"mp_threshold",        1, (void *)&param.mpsigma,      "m", 1},
	{"point_mp_moving_ave", 0, (void *)&param.ponitmp,       "", 1},
	{"min_s1",              1, (void *)&param.mins1,         "", 1},
	{"min_s2",              1, (void *)&param.mins2,         "", 1},
	{"ele_mask",            1, (void *)&param.elmin,      "deg", 1},
	{"ele_threshold",       1, (void *)&param.elcomp,     "deg", 1},
	{"tolerance_clck",      1, (void *)&param.tclock,      "ms", 1},
	{"outflag",             3, (void *)&param.outflag,"0:file&stdout,1:file,2:stdout", 1},
};

/* print help ----------------------------------------------------------------*/
static void printhelp(void)
{
    //for (i=0;i<sizeof(help)/sizeof(*help);i++) fprintf(stderr,"%s\n",help[i]);
    exit(0);
}


static int postpos_qc(gtime_t ts, gtime_t te, double ti, double tu,
                   const prcopt_t *popt, const solopt_t *sopt,
				   const filopt_t *fopt, char **infile, int n, char *outfile,
				   const char *rov, const char *base, const qc_param param);


/* open procssing session ----------------------------------------------------*/
extern int openses_qc(const prcopt_t *popt, const solopt_t *sopt,
                   const filopt_t *fopt, nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr)
{
    char *ext;
    int i,j;
    
    trace(3,"openses :\n");
    
//    for (i=ISYSGPS;i<=NSYS;i++) {
//        for (j=1;j<=MAXFREQ;j++) {
//            setcodepri( sysno(i), j, popt->codepri[i][j-1]);
//        }
//    }
//
//
//    /*L2Cのプライオリティ設定*/
//	//setrnxcodepri(1,popt->l2cprior ==0?"CPYWMNDSLX":"SLXCPYWMND");
//	if (popt->l2cprior ==1) {
//	  setcodepri(SYS_GPS, 2, "SLXCPYWMND");
//	}
//    /* read satellite antenna parameters */
//    if (*fopt->satantp&&!(readpcv(fopt->satantp,pcvs))) {
//		showmsg("error : no sat ant pcv in %s",fopt->satantp);
//		trace(1,"sat antenna pcv read error: %s\n",fopt->satantp);
//        return 0;
//    }
//    /* read receiver antenna parameters */
//    if (*fopt->rcvantp&&!(readpcv(fopt->rcvantp,pcvr))) {
//		showmsg("error : no rec ant pcv in %s",fopt->rcvantp);
//        trace(1,"rec antenna pcv read error: %s\n",fopt->rcvantp);
//        return 0;
//    }
//    /* read dcb parameters */
//    if (*fopt->dcb) {
//        readdcb(fopt->dcb,nav);
//    }
//
//    /* read isb parameters */
//    if (*fopt->isb) {
//        readisb(fopt->isb,nav);
//    }
//    
//	/*1/4波長シフトテーブルを読み込み*/
//	if(*popt->mopt.ifpcs){
//		if((readL2C(popt->mopt.ifpcs,nav))!=0){
//			showmsg("error : no 1/4cycle phase correction file");
//			trace(1,"no 1/4cycle phase correction file\n");
//			return 0;
//		}
//	}
//	
//	/*GLONASS IFBテーブルを使用？*/
//	if(popt->glomodear==GLO_ARMODE_IFB && popt->navsys & SYS_GLO){
//		/*IFBテーブルの読み込み*/
//		if((readifb(popt->mopt.ififb,nav))!=0){
//			showmsg("error : no GLONASS IFB table file");
//			trace(1,"no GLONASS IFB table file\n");
//			return 0;
//		}
//	}
//	/*観測誤差モデルを使用？*/
//	if(popt->errmodel==ERRMODEL_TABLE){
//		/*観測誤差モデルの読み込み*/
//		if((readerr(popt->mopt.iferr,nav))!=0){
//			showmsg("error : no error model file");
//        	trace(1,"no error model file\n");
//        	return 0;
//		}
//	}
//	
//	/*時系補正「Correction」モード */
//	if(popt->tsyscorr==TSYSCORR_CORR){
//		/*GLONASS時系変換パラメータの読み込み*/
//		if((readcirt(fopt->cirtfile,nav))!=0){
//			showmsg("error : no BIPM Circular T file");
//			trace(1,"no BIPM Circular T file\n");
//			return 0;
//		}
//	}
//
//    /* read ionosphere data file */
//    if (*fopt->iono&&(ext=strrchr(fopt->iono,'.'))) {
//        if (strlen(ext)==4&&(ext[3]=='i'||ext[3]=='I')) {
//            readtec(fopt->iono,nav,0);
//        }
//#ifdef EXTSTEC
//        else if (!strcmp(ext,".stec")||!strcmp(ext,".STEC")) {
//            stec_read(fopt->iono,nav);
//        }
//#endif
//    }
//    /* open geoid data */
//    if (sopt->geoid>0&&*fopt->geoid) {
//        if (!opengeoid(sopt->geoid,fopt->geoid)) {
//            showmsg("error : no geoid data %s",fopt->geoid);
//            trace(2,"no geoid data %s\n",fopt->geoid);
//        }
//    }
//    /* read erp data */
//    if (*fopt->eop) {
//        if (!readerp(fopt->eop,&nav->erp)) {
//            showmsg("error : no erp data %s",fopt->eop);
//            trace(2,"no erp data %s\n",fopt->eop);
//        }
//    }
    return 1;
}

/* close procssing session ---------------------------------------------------*/
extern void closeses_qc(nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr)
{
    trace(3,"closeses:\n");
    
    ///* free antenna parameters */
    //free(pcvs->pcv); pcvs->pcv=NULL; pcvs->n=pcvs->nmax=0;
    //free(pcvr->pcv); pcvr->pcv=NULL; pcvr->n=pcvr->nmax=0;
    //
    ///* close geoid data */
    //closegeoid();
    //
    ///* free erp data */
    //free(nav->erp.data); nav->erp.data=NULL; nav->erp.n=nav->erp.nmax=0;
    
    /* close solution statistics and debug trace */
//    rtkclosestat();
////	rtkclosecrinex();
//	closernxobsh_ion();
////    rtkclosesnx(NULL);
//    traceclose();
}


static void fprintf_qc(FILE *fp, char* out){
	if(outflag==0){
		fprintf(fp, out);
		fprintf(stdout, out);
	}
	else if(outflag==2){
		fprintf(fp, out);
	}
	else if(outflag==2){
		fprintf(stdout, out);
	}

}

static void calcion(obsd_t *obs, double iodsign, const nav_t *nav, double gap){
	int i;
	double lambd[NFREQ]={0},freq[NFREQ]={0},ion[NFREQ]={0};
	double dift=0,dt=0;
	gtime_t t={0};

	for(i=0;i<NFREQ;i++){
		lambd[i]=satwavelen(obs->sat,i,nav);
		freq[i]=(lambd[i]!=0)?CLIGHT/lambd[i]:0;
	}

#if NFREQ >= 2
	if(obs->L[0]>0 && obs->L[1]>0 && freq[1]*freq[1]-freq[0]*freq[0]!=0){
		ion[0]=freq[0]*freq[0]/(freq[1]*freq[1]-freq[0]*freq[0])*(obs->L[0]*lambd[0]-obs->L[1]*lambd[1]);
		sv[obs->sat-1].ion[0].comp=1;
	}
#if NFREQ >= 3
	if(obs->L[0]>0 && obs->L[2]>0 && freq[1]*freq[1]*(freq[2]*freq[2]-freq[0]*freq[0])!=0){
		ion[1]=freq[0]*freq[0]*freq[2]*freq[2]/(freq[1]*freq[1]*(freq[2]*freq[2]-freq[0]*freq[0]))*(obs->L[0]*lambd[0]-obs->L[2]*lambd[2]);
		sv[obs->sat-1].ion[1].comp=1;
	}
#if NFREQ >= 4
	if(obs->L[0]>0 && obs->L[3]>0 && freq[1]*freq[1]*(freq[3]*freq[3]-freq[0]*freq[0])!=0){
		ion[2]=freq[0]*freq[0]*freq[3]*freq[3]/(freq[1]*freq[1]*(freq[3]*freq[3]-freq[0]*freq[0]))*(obs->L[0]*lambd[0]-obs->L[3]*lambd[3]);
		sv[obs->sat-1].ion[2].comp=1;
	}
#if NFREQ >= 5
	if(obs->L[0]>0 && obs->L[4]>0 && freq[1]*freq[1]*(freq[4]*freq[4]-freq[0]*freq[0])!=0){
		ion[3]=freq[0]*freq[0]*freq[4]*freq[4]/(freq[1]*freq[1]*(freq[4]*freq[4]-freq[0]*freq[0]))*(obs->L[0]*lambd[0]-obs->L[4]*lambd[4]);
		sv[obs->sat-1].ion[3].comp=1;
	}
#if NFREQ >= 6
	if(obs->L[0]>0 && obs->L[5]>0 && freq[1]*freq[1]*(freq[5]*freq[5]-freq[0]*freq[0])!=0){
		ion[4]=freq[0]*freq[0]*freq[5]*freq[5]/(freq[1]*freq[1]*(freq[5]*freq[5]-freq[0]*freq[0]))*(obs->L[0]*lambd[0]-obs->L[5]*lambd[5]);
		sv[obs->sat-1].ion[4].comp=1;
	}
#endif
#endif
#endif
#endif
#endif
	for(i=0;i<NFREQ-1;i++){
		obs->ion[i]=0;
		obs->iod[i]=0;
		dt=timediff(t,sv[obs->sat-1].ion[i].time);
		dift=timediff(obs->time,sv[obs->sat-1].ion[i].time);

		if(dt!=0&&dift> gap * 60){
			sv[obs->sat-1].ion[i].inition=0;
			sv[obs->sat-1].ion[i].beforeion=0;
		}
		if(sv[obs->sat-1].ion[i].inition!=0 && ion[i]!=0){
			// ION
			obs->ion[i]=ion[i]-sv[obs->sat-1].ion[i].inition;
			trace(4,"calcion: sat=%2d ion0=%8.3f ion1=%8.3f\n",obs->sat,sv[obs->sat-1].ion[i].beforeion,obs->ion[i]);
			if(dift!=0){
				//IOD
				obs->iod[i]=(obs->ion[i]-sv[obs->sat-1].ion[i].beforeion) * 60.0 / dift;
			}
		}
		if(sv[obs->sat-1].ion[i].inition==0) sv[obs->sat-1].ion[i].inition=ion[i];
		// slip判定
		if(fabs(obs->iod[i])>=iodsign){
			trace(2,"calcion: slip detected rcv=%03d sat=%03d iod=%.3f (%8.3f->%8.3f)\n",
				  obs->rcv,obs->sat,fabs(obs->iod[i]),sv[obs->sat-1].ion[i].beforeion,obs->ion[i]);

			sv[obs->sat-1].ion[i].slip=1;
			sv[obs->sat-1].ion[i].inition+=(obs->ion[i]-sv[obs->sat-1].ion[i].beforeion);
			obs->ion[i]=sv[obs->sat-1].ion[i].beforeion;
			obs->iod[i]=sv[obs->sat-1].ion[i].beforeiod;
	/*		if(ion!=0 && (obs->LLI[0]==1 || obs->LLI[1]==1)){
				obs->ion = sv[obs->sat-1].beforeiod;
				obs->iod = sv[obs->sat-1].iod;
			}
	*/	}
		else if(obs->ion!=0) {
			sv[obs->sat-1].ion[i].beforeion=obs->ion[i];
			sv[obs->sat-1].ion[i].beforeiod=obs->iod[i];
		}
		if(ion[i]!=0) {
			sv[obs->sat-1].ion[i].time=obs->time;
		}
	}
}

extern void detionslip(ssat_t *ssat, obsd_t *obs, int n, const nav_t *nav, const prcopt_t *opt){

	int i,j;

	trace(3,"detionslip: n=%d\n",n);
	for(i=0;i<n&&i<MAXOBS;i++){
		calcion(&obs[i], opt->thresslip, nav, GAP);

		if(NFREQ>=2 && sv[obs[i].sat-1].ion[0].slip==1){
			for (j=0;j<opt->nfreq;j++) ssat[obs[i].sat-1].slip[j]|=1;
			sv[obs[i].sat-1].ion[0].slip=0;
		}
	}

}
//
//static void calcmp(obsd_t *obs, qc_param param){
//
//	int i,j,n;
//	double lambd[NFREQ]={0},freq[NFREQ]={0},c[(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1)]={0};
//	double total,before;
//	double msec=0;
//    int nfreq=NFREQ;
//
//	for(i=0;i<NFREQ;i++){
//		lambd[i]=satwavelen(obs->sat,i,&navs);
//		freq[i]=(lambd[i]!=0)?CLIGHT/lambd[i]:0;
//	}
//#if NFREQ >= 2
//	c[0]=(freq[0]*freq[0]-freq[1]*freq[1]!=0)?freq[0]*freq[0]/(freq[0]*freq[0]-freq[1]*freq[1]):0;
//	c[1]=(freq[0]*freq[0]-freq[1]*freq[1]!=0)?freq[1]*freq[1]/(freq[0]*freq[0]-freq[1]*freq[1]):0;
//#if NFREQ >= 3
//	c[2]=(freq[0]*freq[0]-freq[2]*freq[2]!=0)?freq[0]*freq[0]/(freq[0]*freq[0]-freq[2]*freq[2]):0;
//	c[3]=(freq[0]*freq[0]-freq[2]*freq[2]!=0)?freq[2]*freq[2]/(freq[0]*freq[0]-freq[2]*freq[2]):0;
//	c[4]=(freq[1]*freq[1]-freq[2]*freq[2]!=0)?freq[1]*freq[1]/(freq[1]*freq[1]-freq[2]*freq[2]):0;
//	c[5]=(freq[1]*freq[1]-freq[2]*freq[2]!=0)?freq[2]*freq[2]/(freq[1]*freq[1]-freq[2]*freq[2]):0;
//#if NFREQ >= 4
//	c[6]=(freq[0]*freq[0]-freq[3]*freq[3]!=0)?freq[0]*freq[0]/(freq[0]*freq[0]-freq[3]*freq[3]):0;
//	c[7]=(freq[0]*freq[0]-freq[3]*freq[3]!=0)?freq[3]*freq[3]/(freq[0]*freq[0]-freq[3]*freq[3]):0;
//	c[8]=(freq[1]*freq[1]-freq[3]*freq[3]!=0)?freq[1]*freq[1]/(freq[1]*freq[1]-freq[3]*freq[3]):0;
//	c[9]=(freq[1]*freq[1]-freq[3]*freq[3]!=0)?freq[3]*freq[3]/(freq[1]*freq[1]-freq[3]*freq[3]):0;
//	c[10]=(freq[2]*freq[2]-freq[3]*freq[3]!=0)?freq[2]*freq[2]/(freq[2]*freq[2]-freq[3]*freq[3]):0;
//	c[11]=(freq[2]*freq[2]-freq[3]*freq[3]!=0)?freq[3]*freq[3]/(freq[2]*freq[2]-freq[3]*freq[3]):0;
//#if NFREQ >= 5
//
//	c[12]=(freq[0]*freq[0]-freq[4]*freq[4]!=0)?freq[0]*freq[0]/(freq[0]*freq[0]-freq[4]*freq[4]):0;
//	c[13]=(freq[0]*freq[0]-freq[4]*freq[4]!=0)?freq[4]*freq[4]/(freq[0]*freq[0]-freq[4]*freq[4]):0;
//	c[14]=(freq[1]*freq[1]-freq[4]*freq[4]!=0)?freq[1]*freq[1]/(freq[1]*freq[1]-freq[4]*freq[4]):0;
//	c[15]=(freq[1]*freq[1]-freq[4]*freq[4]!=0)?freq[4]*freq[4]/(freq[1]*freq[1]-freq[4]*freq[4]):0;
//	c[16]=(freq[2]*freq[2]-freq[4]*freq[4]!=0)?freq[2]*freq[2]/(freq[2]*freq[2]-freq[4]*freq[4]):0;
//	c[17]=(freq[2]*freq[2]-freq[4]*freq[4]!=0)?freq[4]*freq[4]/(freq[2]*freq[2]-freq[4]*freq[4]):0;
//	c[18]=(freq[3]*freq[3]-freq[4]*freq[4]!=0)?freq[3]*freq[3]/(freq[3]*freq[3]-freq[4]*freq[4]):0;
//	c[19]=(freq[3]*freq[3]-freq[4]*freq[4]!=0)?freq[4]*freq[4]/(freq[3]*freq[3]-freq[4]*freq[4]):0;
//#if NFREQ >= 6
//	c[20]=(freq[0]*freq[0]-freq[5]*freq[5]!=0)?freq[0]*freq[0]/(freq[0]*freq[0]-freq[5]*freq[5]):0;
//	c[21]=(freq[0]*freq[0]-freq[5]*freq[5]!=0)?freq[5]*freq[5]/(freq[0]*freq[0]-freq[5]*freq[5]):0;
//	c[22]=(freq[1]*freq[1]-freq[5]*freq[5]!=0)?freq[1]*freq[1]/(freq[1]*freq[1]-freq[5]*freq[5]):0;
//	c[23]=(freq[1]*freq[1]-freq[5]*freq[5]!=0)?freq[5]*freq[5]/(freq[1]*freq[1]-freq[5]*freq[5]):0;
//	c[24]=(freq[2]*freq[2]-freq[5]*freq[5]!=0)?freq[2]*freq[2]/(freq[2]*freq[2]-freq[5]*freq[5]):0;
//	c[25]=(freq[2]*freq[2]-freq[5]*freq[5]!=0)?freq[5]*freq[5]/(freq[2]*freq[2]-freq[5]*freq[5]):0;
//	c[26]=(freq[3]*freq[3]-freq[5]*freq[5]!=0)?freq[3]*freq[3]/(freq[3]*freq[3]-freq[5]*freq[5]):0;
//	c[27]=(freq[3]*freq[3]-freq[5]*freq[5]!=0)?freq[5]*freq[5]/(freq[3]*freq[3]-freq[5]*freq[5]):0;
//	c[28]=(freq[4]*freq[4]-freq[5]*freq[5]!=0)?freq[4]*freq[4]/(freq[4]*freq[4]-freq[5]*freq[5]):0;
//	c[29]=(freq[4]*freq[4]-freq[5]*freq[5]!=0)?freq[5]*freq[5]/(freq[4]*freq[4]-freq[5]*freq[5]):0;
//#endif
//#endif
//#endif
//#endif
//#endif
//
//	for(i=0;i<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1);i++){
//		obs->MP[i]=0;
//	}
//
//	switch (NFREQ){
//		case 6:
//			if(obs->P[0]>0&&obs->L[0]>0&&obs->L[5]>0&&c[21]!=0) { obs->MP[20]=obs->P[0]-(2*c[21]+1)*obs->L[0]*lambd[0]+ 2*c[21]   *obs->L[5]*lambd[5];
//			sv[obs->sat-1].mp[20].calcflag=1;}
//			if(obs->P[5]>0&&obs->L[0]>0&&obs->L[5]>0&&c[20]!=0) { obs->MP[21]=obs->P[5]- 2*c[20]   *obs->L[0]*lambd[0]+(2*c[20]-1)*obs->L[5]*lambd[5];
//				sv[obs->sat-1].mp[21].calcflag=1;}
//			if(obs->P[1]>0&&obs->L[1]>0&&obs->L[5]>0&&c[23]!=0) { obs->MP[22]=obs->P[1]-(2*c[23]+1)*obs->L[1]*lambd[1]+ 2*c[23]   *obs->L[5]*lambd[5];
//				sv[obs->sat-1].mp[22].calcflag=1;}
//			if(obs->P[5]>0&&obs->L[1]>0&&obs->L[5]>0&&c[22]!=0) { obs->MP[23]=obs->P[5]- 2*c[22]   *obs->L[1]*lambd[1]+(2*c[22]-1)*obs->L[5]*lambd[5];
//				sv[obs->sat-1].mp[23].calcflag=1;}
//			if(obs->P[2]>0&&obs->L[2]>0&&obs->L[5]>0&&c[25]!=0) { obs->MP[24]=obs->P[2]-(2*c[25]+1)*obs->L[2]*lambd[2]+ 2*c[25]   *obs->L[5]*lambd[5];
//				sv[obs->sat-1].mp[24].calcflag=1;}
//			if(obs->P[5]>0&&obs->L[2]>0&&obs->L[5]>0&&c[24]!=0) { obs->MP[25]=obs->P[5]- 2*c[24]   *obs->L[2]*lambd[2]+(2*c[24]-1)*obs->L[5]*lambd[5];
//				sv[obs->sat-1].mp[25].calcflag=1;}
//			if(obs->P[3]>0&&obs->L[3]>0&&obs->L[5]>0&&c[27]!=0) { obs->MP[26]=obs->P[3]-(2*c[27]+1)*obs->L[3]*lambd[3]+ 2*c[27]   *obs->L[5]*lambd[5];
//				sv[obs->sat-1].mp[26].calcflag=1;}
//			if(obs->P[5]>0&&obs->L[3]>0&&obs->L[5]>0&&c[26]!=0) { obs->MP[27]=obs->P[5]- 2*c[26]   *obs->L[3]*lambd[3]+(2*c[26]-1)*obs->L[5]*lambd[5];
//				sv[obs->sat-1].mp[27].calcflag=1;}
//			if(obs->P[4]>0&&obs->L[4]>0&&obs->L[5]>0&&c[28]!=0) { obs->MP[28]=obs->P[4]-(2*c[29]+1)*obs->L[4]*lambd[4]+ 2*c[29]   *obs->L[5]*lambd[5];
//				sv[obs->sat-1].mp[28].calcflag=1;}
//			if(obs->P[5]>0&&obs->L[4]>0&&obs->L[5]>0&&c[29]!=0) { obs->MP[29]=obs->P[5]- 2*c[28]   *obs->L[4]*lambd[4]+(2*c[28]-1)*obs->L[5]*lambd[5];
//				sv[obs->sat-1].mp[29].calcflag=1;}
//		case 5:
//			if(obs->P[0]>0&&obs->L[0]>0&&obs->L[4]>0&&c[13]!=0) { obs->MP[12]=obs->P[0]-(2*c[13]+1)*obs->L[0]*lambd[0]+ 2*c[13]   *obs->L[4]*lambd[4];
//				sv[obs->sat-1].mp[12].calcflag=1;}
//			if(obs->P[4]>0&&obs->L[0]>0&&obs->L[4]>0&&c[12]!=0) { obs->MP[13]=obs->P[4]- 2*c[12]   *obs->L[0]*lambd[0]+(2*c[12]-1)*obs->L[4]*lambd[4];
//				sv[obs->sat-1].mp[13].calcflag=1;}
//			if(obs->P[1]>0&&obs->L[1]>0&&obs->L[4]>0&&c[15]!=0) { obs->MP[14]=obs->P[1]-(2*c[15]+1)*obs->L[1]*lambd[1]+ 2*c[15]   *obs->L[4]*lambd[4];
//				sv[obs->sat-1].mp[14].calcflag=1;}
//			if(obs->P[4]>0&&obs->L[1]>0&&obs->L[4]>0&&c[14]!=0) { obs->MP[15]=obs->P[4]- 2*c[14]   *obs->L[1]*lambd[1]+(2*c[14]-1)*obs->L[4]*lambd[4];
//				sv[obs->sat-1].mp[15].calcflag=1;}
//			if(obs->P[2]>0&&obs->L[2]>0&&obs->L[4]>0&&c[17]!=0) { obs->MP[16]=obs->P[2]-(2*c[17]+1)*obs->L[2]*lambd[2]+ 2*c[17]   *obs->L[4]*lambd[4];
//				sv[obs->sat-1].mp[16].calcflag=1;}
//			if(obs->P[4]>0&&obs->L[2]>0&&obs->L[4]>0&&c[16]!=0) { obs->MP[17]=obs->P[4]- 2*c[16]   *obs->L[2]*lambd[2]+(2*c[16]-1)*obs->L[4]*lambd[4];
//				sv[obs->sat-1].mp[17].calcflag=1;}
//			if(obs->P[3]>0&&obs->L[3]>0&&obs->L[4]>0&&c[19]!=0) { obs->MP[18]=obs->P[3]-(2*c[19]+1)*obs->L[3]*lambd[3]+ 2*c[19]   *obs->L[4]*lambd[4];
//				sv[obs->sat-1].mp[18].calcflag=1;}
//			if(obs->P[4]>0&&obs->L[3]>0&&obs->L[4]>0&&c[18]!=0) { obs->MP[19]=obs->P[4]- 2*c[18]   *obs->L[3]*lambd[3]+(2*c[18]-1)*obs->L[4]*lambd[4];
//				sv[obs->sat-1].mp[19].calcflag=1;}
//		case 4:
//			if(obs->P[0]>0&&obs->L[0]>0&&obs->L[3]>0&&c[7]!=0)  {  obs->MP[6]=obs->P[0]-(2*c[7]+1) *obs->L[0]*lambd[0]+ 2*c[7]    *obs->L[3]*lambd[3];
//				sv[obs->sat-1].mp[6].calcflag=1;}
//			if(obs->P[3]>0&&obs->L[0]>0&&obs->L[3]>0&&c[6]!=0)  {  obs->MP[7]=obs->P[3]- 2*c[6]    *obs->L[0]*lambd[0]+(2*c[6]-1) *obs->L[3]*lambd[3];
//				sv[obs->sat-1].mp[7].calcflag=1;}
//			if(obs->P[1]>0&&obs->L[1]>0&&obs->L[3]>0&&c[9]!=0)  {  obs->MP[8]=obs->P[1]-(2*c[9]+1) *obs->L[1]*lambd[1]+ 2*c[9]    *obs->L[3]*lambd[3];
//				sv[obs->sat-1].mp[8].calcflag=1;}
//			if(obs->P[3]>0&&obs->L[1]>0&&obs->L[3]>0&&c[8]!=0)  {  obs->MP[9]=obs->P[3]- 2*c[8]    *obs->L[1]*lambd[1]+(2*c[8]-1) *obs->L[3]*lambd[3];
//				sv[obs->sat-1].mp[9].calcflag=1;}
//			if(obs->P[2]>0&&obs->L[2]>0&&obs->L[3]>0&&c[10]!=0) { obs->MP[10]=obs->P[2]-(2*c[11]+1)*obs->L[2]*lambd[2]+ 2*c[11]   *obs->L[3]*lambd[3];
//				sv[obs->sat-1].mp[10].calcflag=1;}
//			if(obs->P[3]>0&&obs->L[2]>0&&obs->L[3]>0&&c[11]!=0) { obs->MP[11]=obs->P[3]- 2*c[10]   *obs->L[2]*lambd[2]+(2*c[10]-1)*obs->L[3]*lambd[3];
//				sv[obs->sat-1].mp[11].calcflag=1;}
//		case 3:
//			if(obs->P[0]>0&&obs->L[0]>0&&obs->L[2]>0&&c[3]!=0) { obs->MP[2]=obs->P[0]-(2*c[3]+1)*obs->L[0]*lambd[0]+ 2*c[3]   *obs->L[2]*lambd[2];
//				sv[obs->sat-1].mp[2].calcflag=1;}
//			if(obs->P[2]>0&&obs->L[0]>0&&obs->L[2]>0&&c[2]!=0) { obs->MP[3]=obs->P[2]- 2*c[2]   *obs->L[0]*lambd[0]+(2*c[2]-1)*obs->L[2]*lambd[2];
//				sv[obs->sat-1].mp[3].calcflag=1;}
//			if(obs->P[1]>0&&obs->L[1]>0&&obs->L[2]>0&&c[5]!=0) { obs->MP[4]=obs->P[1]-(2*c[5]+1)*obs->L[1]*lambd[1]+ 2*c[5]   *obs->L[2]*lambd[2];
//				sv[obs->sat-1].mp[4].calcflag=1;}
//			if(obs->P[2]>0&&obs->L[1]>0&&obs->L[2]>0&&c[4]!=0) { obs->MP[5]=obs->P[2]- 2*c[4]   *obs->L[1]*lambd[1]+(2*c[4]-1)*obs->L[2]*lambd[2];
//				sv[obs->sat-1].mp[5].calcflag=1;}
//		case 2:
//			if(obs->P[0]>0&&obs->L[0]>0&&obs->L[1]>0&&c[1]!=0) { obs->MP[0]=obs->P[0]-(2*c[1]+1)*obs->L[0]*lambd[0]+ 2*c[1]   *obs->L[1]*lambd[1];
//				sv[obs->sat-1].mp[0].calcflag=1;}
//			if(obs->P[1]>0&&obs->L[0]>0&&obs->L[1]>0&&c[0]!=0) { obs->MP[1]=obs->P[1]- 2*c[0]   *obs->L[0]*lambd[0]+(2*c[0]-1)*obs->L[1]*lambd[1];
//				sv[obs->sat-1].mp[1].calcflag=1;}
//			break;
//	}
//
//	for(i=0;i<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1);i++){
//		sv[obs->sat-1].mp[i].slip=0;
//		sv[obs->sat-1].mp[i].msecslip=0;
//		n=0;
//		total=0;
//
//		for(j=1;j<param.ponitmp;j++){
//			sv[obs->sat-1].mp[i].val[j-1]=sv[obs->sat-1].mp[i].val[j];
//			if(obs->MP[i]!=0 && sv[obs->sat-1].mp[i].val[j-1]!=0){
//				total+=sv[obs->sat-1].mp[i].val[j-1];
//				n++;
//			}
//		}
//		sv[obs->sat-1].mp[i].val[param.ponitmp-1]=obs->MP[i];
//		if(obs->MP[i]!=0){
//			before=n>0?total/n:0;
//			sv[obs->sat-1].mp[i].msec=obs->MP[i]/CLIGHT*1000;
//			//if(fabs(sv[obs->sat-1].mp[i].msec)>=1-param.tclock){
//			if(fabs(sv[obs->sat-1].mp[i].msec)>=1+param.tclock){
//				sv[obs->sat-1].mp[i].msecslip=1;
//			}
//			obs->MP[i]-=(total+obs->MP[i])/(n+1);
//			if(fabs(obs->MP[i])>=param.mpsigma){
//				if(sv[obs->sat-1].mp[i].constmp==0){
//					sv[obs->sat-1].mp[i].constmp=sv[obs->sat-1].mp[i].val[param.ponitmp-1]-before;
//					sv[obs->sat-1].mp[i].val[param.ponitmp-1]=before;
//					obs->MP[i]=before-(total+before)/(n+1);
//					sv[obs->sat-1].mp[i].slip=1;
//				}
//				else{
//					obs->MP[i]=sv[obs->sat-1].mp[i].val[param.ponitmp-1]-sv[obs->sat-1].mp[i].constmp;
//					obs->MP[i]-=(total+obs->MP[i])/(n+1);
//					if(fabs(obs->MP[i])>=param.mpsigma){
//						sv[obs->sat-1].mp[i].constmp=sv[obs->sat-1].mp[i].val[param.ponitmp-1]-before;
//						sv[obs->sat-1].mp[i].val[param.ponitmp-1]=before;
//						obs->MP[i]=before-(total+before)/(n+1);
//						sv[obs->sat-1].mp[i].slip=1;
//					}
//					else{
//						sv[obs->sat-1].mp[i].val[param.ponitmp-1]-=sv[obs->sat-1].mp[i].constmp;
//					}
//				}
//			}
//			else{
//				sv[obs->sat-1].mp[i].constmp=0;
//			}
//		}
//	}
//
//}

//
//static void initvec(hist_data *ele_hist,int obssn){
//
//	int i,j,k;
//	gtime_t ts={0};
//
//	for(i=0;i<MAXSAT;i++){
//		if(i<NSATGPS) sv[i].count = (int*)calloc(obstype[0] * 2,sizeof(int));
//		else if(i<NSATGPS+MAXPRNGLO) sv[i].count = (int*)calloc(obstype[1] * 2,sizeof(int));
//		else if(i<NSATGPS+MAXPRNGLO+NSATGAL) sv[i].count = (int*)calloc(obstype[2] * 2,sizeof(int));
//		else if(i<NSATGPS+MAXPRNGLO+NSATGAL+NSATQZS) sv[i].count = (int*)calloc(obstype[3] * 2,sizeof(int));
//		else if(i<NSATGPS+MAXPRNGLO+NSATGAL+NSATQZS+NSATCMP) sv[i].count = (int*)calloc(obstype[5] * 2,sizeof(int));
//		else sv[i].count = (int*)calloc(obstype[4] * 2,sizeof(int));
//		for(j=0;j<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1);j++){
//			sv[i].mp[j].abcomp=0;
//			sv[i].mp[j].belcomp=0;
//			sv[i].mp[j].calcflag=0;
//			sv[i].mp[j].constmp=0;
//			sv[i].mp[j].del=0;
//			sv[i].mp[j].ele=0;
//			sv[i].mp[j].num=0;
//			sv[i].mp[j].rms=0;
//			sv[i].mp[j].slip=0;
//			for(k=0;k<POINT_AVG;k++){
//				sv[i].mp[j].val[k]=0;
//			}
//		}
//	}
//	for(i=0;i<HIST_VER;i++){
//		ele_hist[i].valueS1 = (double *)calloc(obssn,sizeof(double));
//		ele_hist[i].valueS2 = (double *)calloc(obssn,sizeof(double));
//		for(j=0;j<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1);j++){
//			for(k=0;k<NSYS;k++){
//				ele_hist[i].mp[j].n[k]=0;
//				ele_hist[i].mp[j].rms[k]=0;
//				ele_hist[i].mp[j].slps[k]=0;
//			}
//		}
//	}
//}

static int signalkind(obsd_t obs){
	int kind=100;
	// A/S
	if(obs.LLIcode[0]==4||obs.LLIcode[1]==4||obs.LLIcode[2]==4){
		// L1
		if(obs.L[0]!=0){
			if(obs.P[0]!=0){
				// A/S on , L1 , C1
				if(obs.code[0]==CODE_L1C){
					if(obs.L[1]!=0 && obs.P[1]!=0){
						if(obs.code[1]>=CODE_L2C && obs.code[1]<=CODE_L2X){
							kind=ASCII_SYMBOL24_NO;
						}
						else{
							kind=ASCII_SYMBOL26_NO;
						}
					}
					else if(obs.L[1]==0 && obs.P[1]!=0 && obs.code[1]>=CODE_L2C && obs.code[1]<=CODE_L2X){
						kind=ASCII_SYMBOL17_NO;
					}
					else if(obs.L[2]!=0 && obs.P[2]!=0 && obs.code[2]>=CODE_L5I && obs.code[2]<=CODE_L5X){
						kind=ASCII_SYMBOL25_NO;
					}
					else if(obs.L[2]==0 && obs.P[2]!=0 && obs.code[2]>=CODE_L5I && obs.code[2]<=CODE_L5X){
						kind=ASCII_SYMBOL18_NO;
					}
					else{
						kind=ASCII_SYMBOL22_NO;
					}
				}
				// A/S on , L1 , P1
				else{
					if(obs.L[1]!=0 && obs.P[1]!=0){
						if(obs.code[1]>=CODE_L2C && obs.code[1]<=CODE_L2X){
							kind=ASCII_SYMBOL24_NO;
						}
						else{
							kind=ASCII_SYMBOL26_NO;
						}
					}
					else if(obs.L[1]==0 && obs.P[1]!=0 && obs.code[1]>=CODE_L2C && obs.code[1]<=CODE_L2X){
						kind=ASCII_SYMBOL17_NO;
					}
					else{
						kind=ASCII_SYMBOL23_NO;
					}
				}
			}
		}
		// no L1
		else{
			if(obs.code[0]==CODE_L1C && obs.P[0]!=0){
				kind=ASCII_SYMBOL21_NO;
			}
		}
	}
	else{
		// L1
		if(obs.L[0]!=0){
			if(obs.P[0]!=0){
				// no A/S , L1 , C1
				if(obs.code[0]==CODE_L1C){
					if(obs.L[1]!=0 && obs.P[1]!=0){
						if(obs.code[1]>=CODE_L2C && obs.code[1]<=CODE_L2X){
							kind=ASCII_SYMBOL24_NO;
						}
						else{
							kind=ASCII_SYMBOL19_NO;
						}
					}
					else if(obs.L[1]==0 && obs.P[1]!=0 && obs.code[1]>=CODE_L2C && obs.code[1]<=CODE_L2X){
						kind=ASCII_SYMBOL17_NO;
					}
					else if(obs.L[2]!=0 && obs.P[2]!=0 && obs.code[2]>=CODE_L5I && obs.code[2]<=CODE_L5X){
						kind=ASCII_SYMBOL25_NO;
					}
					else if(obs.L[2]==0 && obs.P[2]!=0 && obs.code[2]>=CODE_L5I && obs.code[2]<=CODE_L5X){
						kind=ASCII_SYMBOL18_NO;
					}
					else{
						kind=ASCII_SYMBOL15_NO;
					}
				}
				// no A/S , L1 , P1
				else{
					if(obs.L[1]!=0 && obs.P[1]!=0){
						if(obs.code[1]>=CODE_L2C && obs.code[1]<=CODE_L2X){
							kind=ASCII_SYMBOL24_NO;
						}
						else{
							kind=ASCII_SYMBOL20_NO;
						}
					}
					else if(obs.L[1]==0 && obs.P[1]!=0 && obs.code[1]>=CODE_L2C && obs.code[1]<=CODE_L2X){
						kind=ASCII_SYMBOL17_NO;
					}
					else{
						kind=ASCII_SYMBOL16_NO;
					}
				}
			}
		}
		// no L1
		else{
			if(obs.code[0]==CODE_L1C && obs.P[0]!=0){
				kind=ASCII_SYMBOL14_NO;
			}
		}
	}
	return kind;
}
//static int checksv(obsd_t obs, double elevation, int *repcp, double elmin, int abmask[][PLOT_WIDTH],
//						satvalue *satval, int qc[][PLOT_WIDTH],int qc_index, double *totals1,
//						double *totals2,int *s1n, int *s2n, double *sn[], int index,
//						hist_data *ele_hist,int dn[][PLOT_WIDTH], double iodsign, double sigma,
//						double elcomp, qc_param param, double *maxmsec, double *minmsec){
//
//	int h_index=0,k,kmax=0,sys_index=0,i,ionslip=0,clckflg;
//
//	// include code | phase
//	if(obs.L[0]>0 || obs.L[1]>0 || obs.P[0]>0 || obs.P[1]>0){
//		(*repcp)++;
//	}
//	// abode horizon
//	if(elevation > 0){
//		sv[index].hor++;
//		sv[index].sum_ele_hor += elevation*R2D;
//		sv[index].ele_hor = sv[index].sum_ele_hor / sv[index].hor;
//		// abode mask
//		if(elevation >= elmin){
//			sv[index].mask++;
//			sv[index].reprt++;
//			sv[index].sum_ele_mask += elevation*R2D;
//			sv[index].ele_mask = sv[index].sum_ele_mask / sv[index].mask;
//			abmask[index][qc_index]=1;
//			if(obs.L[0]>0 && obs.L[1]>0 && obs.P[0]>0 && obs.P[1]>0){
//				sv[index].compl++;
//				for(i=0;i<NFREQ;i++){
//					if(obs.LLI[i]==1){
//						if(elevation*R2D >= elcomp){
//							sv[index].lli[i].above++;
//						}
//						else{
//							sv[index].lli[i].below++;
//						}
//					}
//				}
//			}
//			else{
//				(satval->totaldel)++;
//			}
//			if(index<NSATGPS) kmax = obstype[0];
//			else if(index<NSATGPS+MAXPRNGLO) kmax = obstype[1];
//			else if(index<NSATGPS+MAXPRNGLO+NSATGAL) kmax = obstype[2];
//			else if(index<NSATGPS+MAXPRNGLO+NSATGAL+NSATQZS) kmax = obstype[3];
//			else if(index<NSATGPS+MAXPRNGLO+NSATGAL+NSATQZS+NSATCMP) kmax = obstype[5];
//			else kmax = obstype[4];
//			for(k=0;k<kmax;k++){
//				if(obs.value[k]!=0)sv[index].count[k]++;
//				else sv[index].count[k + kmax]++;
//			}
//			for(k=0;k<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1);k++){
//				if(sv[index].mp[k].calcflag==1){
//					sv[index].mp[k].num++;
//					sv[index].mp[k].rms+=obs.MP[k]*obs.MP[k];
//					sv[index].mp[k].ele+=elevation;
//					if(sv[index].mp[k].slip==1){
//						if(elevation*R2D >= elcomp){
//							sv[index].mp[k].abcomp++;
//						}
//						else{
//							sv[index].mp[k].belcomp++;
//						}
//					}
//				}
//				else sv[index].mp[k].del++;
//			}
//			// ION slip check
//			ionslip=0;
//			for(k=0;k<NFREQ-1;k++){
//				if(sv[index].ion[k].slip==1){
//					ionslip=1;
//					break;
//				}
//			}
//			if(ionslip==1){
//				sv[index].iodslipab++;
//				sv[index].iodmpslipab++;
//			}
//			else{
//				for(k=0;k<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1)/2;k++){
//					if(sv[index].mp[k*2].slip==1 && sv[index].mp[2*k+1].slip==1){
//						sv[index].iodmpslipab++;
//						break;
//					}
//				}
//			}
//		}
//	}
//	if(elevation < elmin){
//		sv[index].belowmask++;
//		if(elevation!=0){
//			if(obs.L[0]>0 && obs.L[1]>0 && obs.P[0]>0 && obs.P[1]>0){
//				(satval->totalbelowcomp)++;
//				qc[index][qc_index]=MIN(qc[index][qc_index], ASCII_SYMBOL12_NO);
//			}
//			else{
//				if(obs.L[0]>0 || obs.L[1]>0 || obs.P[0]>0 || obs.P[1]>0){
//					qc[index][qc_index]=MIN(qc[index][qc_index], ASCII_SYMBOL13_NO);
//				}
//				else if(obs.L[0]==0 && obs.L[1]==0 && obs.P[0]==0 && obs.P[1]==0){
//					qc[index][qc_index]=MIN(qc[index][qc_index], ASCII_SYMBOL29_NO);
//				}
//				dn[(index)*2+1][qc_index]++;
//			}
//		}
//		ionslip=0;
//		for(i=0;i<NFREQ-1;i++){
//			if(sv[index].ion[i].slip==1){
//				ionslip=1;
//				break;
//			}
//		}
//		if(ionslip==1){
//			sv[index].iodslipbe++;
//			sv[index].iodmpslipbe++;
//		}
//		else{
//			for(i=0;i<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1)/2;i++){
//				if(sv[index].mp[2*i].slip==1 && sv[index].mp[2*i+1].slip==1){
//					sv[index].iodmpslipbe++;
//					break;
//				}
//			}
//		}
//		if(nav_sat[index]==0){
//			if(obs.L[0]>0 && obs.L[1]>0 && obs.P[0]>0 && obs.P[1]>0){
//				for(i=0;i<NFREQ;i++){
//					if(obs.LLI[i]==1){
//						sv[index].lli[i].below++;
//					}
//				}
//			}
//			for(i=0;i<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1);i++){
//				if(sv[index].mp[i].calcflag==1){
//					sv[index].mp[i].num++;
//					sv[index].mp[i].rms+=obs.MP[i]*obs.MP[i];
//					if(sv[index].mp[i].slip==1){
//						sv[index].mp[i].belcomp++;
//					}
//				}
//			}
//		}
//	}
//	//if(elevation== 0) qc[index][qc_index]=-5;
//
//	if(index<NSATGPS) sys_index=NSYSGPS-1;
//	else if(index<NSATGPS+NSATGLO) sys_index=NSYSGPS+NSYSGLO-1;
//	else if(index<NSATGPS+NSATGLO+NSATGAL) sys_index=NSYSGPS+NSYSGLO+NSYSGAL-1;
//	else if(index<NSATGPS+NSATGLO+NSATGAL+NSATQZS) {
//		sys_index=NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS-1;
//		if(NSYSGPS+NSYSQZS==2) sys_index=NSYSGPS-1;
//	}
//	else if(index<NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP) sys_index=NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS+NSYSCMP-1;
//	else sys_index=NSYS;
//	if(sys_index<0 || sys_index>NSYS) sys_index=NSYS;
//
//	if(elevation <= 0){
//		h_index = 0;
//	}
//	else if(elevation == PI/2){
//		h_index=18;
//	}
//	else if(elevation < PI/2){
//		h_index = (int)((elevation*R2D)/5) + 1;
//	}
//	if(obs.SN[0]>param.mins1){
//		ele_hist[h_index].valueS1[ele_hist[h_index].S1n] = obs.SN[0];
//		ele_hist[h_index].S1n++;
//		ele_hist[h_index].totalS1 += obs.SN[0];
//		ele_hist[h_index].meanS1 = ele_hist[h_index].totalS1 / ele_hist[h_index].S1n;
//	}
//	if(obs.SN[1]>param.mins2){
//		ele_hist[h_index].valueS2[ele_hist[h_index].S2n] = obs.SN[1];
//		ele_hist[h_index].S2n++;
//		ele_hist[h_index].totalS2 += obs.SN[1];
//		ele_hist[h_index].meanS2 = ele_hist[h_index].totalS2 / ele_hist[h_index].S2n;
//	}
//	clckflg=0;
//	for(i=0;i<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1);i++){
//		if(sv[index].mp[i].calcflag==1){
//			ele_hist[h_index].mp[i].rms[sys_index]+=obs.MP[i]*obs.MP[i];
//			if(sv[index].mp[i].slip==1){
//				ele_hist[h_index].mp[i].slps[sys_index]++;
//				satval->totalmp++;
//				switch(i){
//					case 0:
//						qc[index][qc_index] = MIN(qc[index][qc_index],ASCII_SYMBOL5_NO);
//						break;
//					case 1:
//						qc[index][qc_index] = MIN(qc[index][qc_index],ASCII_SYMBOL6_NO);
//						break;
//					case 2:
//						qc[index][qc_index] = MIN(qc[index][qc_index],ASCII_SYMBOL7_NO);
//						break;
//					case 3:
//						qc[index][qc_index] = MIN(qc[index][qc_index],ASCII_SYMBOL8_NO);
//						break;
//
//				}
//			}
//			else{
//				ele_hist[h_index].mp[i].n[sys_index]++;
//			}
//			if(sv[index].mp[i].msecslip==1){
//				qc[index][qc_index] = MIN(qc[index][qc_index],ASCII_SYMBOL2_NO);
//				(*maxmsec)=MAX(*maxmsec, sv[index].mp[i].msec);
//				(*minmsec)=MIN(*minmsec, sv[index].mp[i].msec);
//			}
//			else{
//				clckflg=1;
//			}
//		}
//	}
//	for(i=0;i<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1)/2;i++){
//		if(sv[index].mp[2*i].slip==1 && sv[index].mp[2*i+1].slip==1){
//			qc[index][qc_index] = MIN(qc[index][qc_index],ASCII_SYMBOL4_NO);
//			break;
//		}
//	}
//	for(i=0;i<NFREQ-1;i++){
//		if(sv[index].ion[i].comp==1 && sv[index].ion[i].slip==0){
//			ele_hist[h_index].ionn++;
//			ele_hist[h_index].ionrms+=obs.ion[i]*obs.ion[i];
//		}
//		else if(sv[index].ion[i].slip==1){
//			ele_hist[h_index].ionslps++;
//			ele_hist[h_index].ionrms+=obs.ion[i]*obs.ion[i];
//			qc[index][qc_index] = MIN(qc[index][qc_index],ASCII_SYMBOL3_NO);
//		}
//	}
//	if(obs.LLI[0]==1 || obs.LLI[1]==1 || obs.LLI[2]==1){
//		qc[index][qc_index] = ASCII_SYMBOL10_NO;
//	}
//	if(qc[index][qc_index]>=0){
//		qc[index][qc_index] = MAX(qc[index][qc_index],signalkind(obs));
//	}
//	if(obs.SN[0]>param.mins1){
//		(*totals1)+=obs.SN[0];
//		sn[0][*s1n]=obs.SN[0];
//		(*s1n)++;
//	}
//	if(obs.SN[1]>param.mins2){
//		(*totals2)+=obs.SN[1];
//		sn[1][*s2n]=obs.SN[1];
//		(*s2n)++;
//	}
//
//	return clckflg;
//}


static void checksv_noobs(double elevation, int index, int dn[][PLOT_WIDTH], int abmask[][PLOT_WIDTH],
						  int qc[][PLOT_WIDTH], double elmin, int qc_index){

	if(elevation > 0){
		sv[index].hor++;
		sv[index].sum_ele_hor += elevation*R2D;
		sv[index].ele_hor = sv[index].sum_ele_hor / sv[index].hor;
		qc[index][qc_index] = ASCII_SYMBOL29_NO;
		if(elevation >= elmin){
			sv[index].mask++;
			sv[index].sum_ele_mask += elevation*R2D;
			sv[index].ele_mask = sv[index].sum_ele_mask / sv[index].mask;
			dn[(index)*2][qc_index]++;
			abmask[index][qc_index]=1;
			qc[index][qc_index]=MIN(qc[index][qc_index], ASCII_SYMBOL9_NO);
		}
		else{
			qc[index][qc_index]=MIN(qc[index][qc_index], ASCII_SYMBOL29_NO);
		}
	}
}

/* write header to output file -----------------------------------------------*/
static int outhead_qc(const char *outfile, char **infile, int n,
                   const prcopt_t *popt, const solopt_t *sopt)
{
    FILE *fp=stdout;
    
    trace(3,"outhead: outfile=%s n=%d\n",outfile,n);
    
    if (*outfile) {
        createdir(outfile);
        
        if (!(fp=fopen(outfile,"w"))) {
            showmsg("error : open output file %s",outfile);
            return 0;
        }
    }
    /* output header */
	fprintf(fp,"version : gsiqc\n\n");
    
    if (*outfile) fclose(fp);
    
    return 1;
}


static void outsvplot(FILE *fp, int qc[][PLOT_WIDTH], int dn[][PLOT_WIDTH], double elmin, int abmask[][PLOT_WIDTH],
					  double *ep1, double *ep2, int *clk){
	int tick_sum=0,i,j,k;
	char id[64];
	char out[1024];
   	const char *abbreviated_name[]={
		"","JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC",""
	};

	fprintf_qc(fp," SV+");
	for(i=0;i<PLOT_WIDTH;i++){
		if(tick_sum == (int)(total_h/PLOT_WIDTH * i / tick)){
			fprintf_qc(fp,"-");
		}
		else{
			fprintf_qc(fp,"|");
		}
		tick_sum = (int)(total_h/PLOT_WIDTH * i / tick);
	}
	fprintf_qc(fp,"+ SV\n");

	for(j=0;j<MAXSAT;j++){
		if(obs_sat[j]==0 && nav_sat[j]==0) continue;
		satno2id(j+1,id);
		sprintf(out,"%s|",id);
		fprintf_qc(fp,out);
		for(i=0;i<PLOT_WIDTH;i++){
			switch(qc[j][i]){
				case ASCII_SYMBOL1_NO:fprintf_qc(fp,ASCII_SYMBOL1);break;
				case ASCII_SYMBOL2_NO:fprintf_qc(fp,ASCII_SYMBOL2);break;
				case ASCII_SYMBOL3_NO:fprintf_qc(fp,ASCII_SYMBOL3);break;
				case ASCII_SYMBOL4_NO:fprintf_qc(fp,ASCII_SYMBOL4);break;
				case ASCII_SYMBOL5_NO:fprintf_qc(fp,ASCII_SYMBOL5);break;
				case ASCII_SYMBOL6_NO:fprintf_qc(fp,ASCII_SYMBOL6);break;
				case ASCII_SYMBOL7_NO:fprintf_qc(fp,ASCII_SYMBOL7);break;
				case ASCII_SYMBOL8_NO:fprintf_qc(fp,ASCII_SYMBOL8);break;
				case ASCII_SYMBOL9_NO:fprintf_qc(fp,ASCII_SYMBOL9);break;
				case ASCII_SYMBOL10_NO:fprintf_qc(fp,ASCII_SYMBOL10);break;
				case ASCII_SYMBOL11_NO:fprintf_qc(fp,ASCII_SYMBOL11);break;
				case ASCII_SYMBOL12_NO:fprintf_qc(fp,ASCII_SYMBOL12);break;
				case ASCII_SYMBOL13_NO:fprintf_qc(fp,ASCII_SYMBOL13);break;
				case ASCII_SYMBOL14_NO:fprintf_qc(fp,ASCII_SYMBOL14);break;
				case ASCII_SYMBOL15_NO:fprintf_qc(fp,ASCII_SYMBOL15);break;
				case ASCII_SYMBOL16_NO:fprintf_qc(fp,ASCII_SYMBOL16);break;
				case ASCII_SYMBOL17_NO:fprintf_qc(fp,ASCII_SYMBOL17);break;
				case ASCII_SYMBOL18_NO:fprintf_qc(fp,ASCII_SYMBOL18);break;
				case ASCII_SYMBOL19_NO:fprintf_qc(fp,ASCII_SYMBOL19);break;
				case ASCII_SYMBOL20_NO:fprintf_qc(fp,ASCII_SYMBOL20);break;
				case ASCII_SYMBOL21_NO:fprintf_qc(fp,ASCII_SYMBOL21);break;
				case ASCII_SYMBOL22_NO:fprintf_qc(fp,ASCII_SYMBOL22);break;
				case ASCII_SYMBOL23_NO:fprintf_qc(fp,ASCII_SYMBOL23);break;
				case ASCII_SYMBOL24_NO:fprintf_qc(fp,ASCII_SYMBOL24);break;
				case ASCII_SYMBOL25_NO:fprintf_qc(fp,ASCII_SYMBOL25);break;
				case ASCII_SYMBOL26_NO:fprintf_qc(fp,ASCII_SYMBOL26);break;
				case ASCII_SYMBOL27_NO:fprintf_qc(fp,ASCII_SYMBOL27);break;
				case ASCII_SYMBOL28_NO:fprintf_qc(fp,ASCII_SYMBOL28);break;
				case ASCII_SYMBOL29_NO:fprintf_qc(fp,ASCII_SYMBOL29);break;
				default:fprintf_qc(fp," ");break;
			}
		}
		sprintf(out,"|%s\n",id);
		fprintf_qc(fp,out);
	}
	fprintf_qc(fp,"-dn|");
	for(i=0;i<PLOT_WIDTH;i++){
		tick_sum=0;
		for(k=0;k<MAXSAT;k++){
			if(dn[k*2+1][i]>0){
				tick_sum--;
			}
			if(dn[MAXSAT*2][i]!=0 && dn[k*2][i]==dn[MAXSAT*2][i]){
				tick_sum++;
			}
		}
		if(tick_sum<0){
			fprintf_qc(fp,"+");
		}
		else if(tick_sum>0 && tick_sum<10){
			sprintf(out,"%d",tick_sum);
			fprintf_qc(fp,out);
		}
		else if(tick_sum>=10){
			sprintf(out,"%c",(char)(tick_sum+87));
			fprintf_qc(fp,out);
		}
		else{
			fprintf_qc(fp," ");
		}
	}
	fprintf_qc(fp,"|-dn\n");
	fprintf_qc(fp,"+dn|");
	for(i=0;i<PLOT_WIDTH;i++){
		tick_sum=0;
		for(k=0;k<MAXSAT;k++){
			if(dn[MAXSAT*2][i]!=0 && dn[k*2+1][i]==dn[MAXSAT*2][i]){
				tick_sum--;
			}
			if(dn[k*2][i]>0){
				tick_sum++;
			}
		}
		if(tick_sum<0){
			fprintf_qc(fp,"+");
		}
		else if(tick_sum>0 && tick_sum<10){
			sprintf(out,"%d",tick_sum);
			fprintf_qc(fp,out);
		}
		else if(tick_sum>=10){
			sprintf(out,"%c",(char)(tick_sum+87));
			fprintf_qc(fp,out);
		}
		else{
			fprintf_qc(fp," ");
		}
	}
	fprintf_qc(fp,"|+dn\n");
	sprintf(out,"+%2.0f|",elmin*R2D);
	fprintf_qc(fp,out);
	for(i=0;i<PLOT_WIDTH;i++){
		tick_sum=0;
		for(k=0;k<MAXSAT;k++){
			tick_sum += abmask[k][i];
		}
		if(tick_sum>0 && tick_sum<10){
			sprintf(out,"%d",tick_sum);
			fprintf_qc(fp,out);
		}
		else if(tick_sum>=10){
			sprintf(out,"%c",(char)(tick_sum+87));
			fprintf_qc(fp,out);
		}
		else{
			fprintf_qc(fp," ");
		}
	}
	sprintf(out,"|+%2.0f\n",elmin*R2D);
	fprintf_qc(fp,out);
	sprintf(out,"Clk|");
	fprintf_qc(fp,out);
	for(i=0;i<PLOT_WIDTH;i++){
		switch(clk[i]){
			case 2: fprintf_qc(fp,"-");break;
			case 1: fprintf_qc(fp,"+");break;
			default: fprintf_qc(fp," ");break;
		}
	}
	fprintf_qc(fp,"|Clk\n");
	fprintf_qc(fp,"   +");
	tick_sum=0;
	for(i=0;i<PLOT_WIDTH;i++){
		if(tick_sum == (int)(total_h/PLOT_WIDTH * i / tick)){
			fprintf_qc(fp,"-");
		}
		else{
			fprintf_qc(fp,"|");
		}
		tick_sum = (int)(total_h/PLOT_WIDTH * i / tick);
	}
	fprintf_qc(fp,"+   \n");

	time2epoch(obsts,ep1);
    time2epoch(obste,ep2);
	sprintf(out,"%02.0f:%02.0f:%06.3f                                                        %02.0f:%02.0f:%06.3f\n",
		ep1[3],ep1[4],ep1[5],ep2[3],ep2[4],ep2[5]);
	fprintf_qc(fp,out);
	sprintf(out,"%04.0f %s %02.0f                                                          %04.0f %s %02.0f\n",
		ep1[0],abbreviated_name[(int)ep1[1]],ep1[2],ep2[0],abbreviated_name[(int)ep2[1]],ep2[2]);
	fprintf_qc(fp,out);

}

static void outobsnav(FILE *fp,char *sysname,int start, int end){
	int i,count=0;
	char out[1024];
	sprintf(out,"%s SVs w/o OBS :",sysname);
	fprintf_qc(fp,out);
	for(i = start; i < end; i++){
		if(obs_sat[i]==0){
			count++;
			if(count==13 || count==25) fprintf_qc(fp,"\n                         ");
			sprintf(out," %2d ",i+1-start);
			fprintf_qc(fp,out);
		}
	}
	count=0;
	sprintf(out,"\n%s SVs w/o NAV :",sysname);
	fprintf_qc(fp,out);
	for(i = start; i < end; i++){
		if(obs_sat[i]==1 && nav_sat[i]==0){
			count++;
			if(count==13 || count==25) fprintf_qc(fp,"\n                         ");
			sprintf(out," %2d ",i+1-start);
			fprintf_qc(fp,out);
		}
	}
}
static void outsvnumber(FILE *fp,int totalsat){

	char out[1024];

	sprintf(out,"Total satellites w/ obs : %d\n",totalsat);
	fprintf_qc(fp,out);

	if(NSYSGPS){
		outobsnav(fp,"NAVSTAR GPS",0,NSATGPS);
	}
	if(NSYSGLO){
		fprintf_qc(fp,"\n");
		outobsnav(fp,"    GLONASS",NSATGPS,NSATGPS+NSATGLO);
	}
	if(NSYSGAL){
		fprintf_qc(fp,"\n");
		outobsnav(fp,"    GALILEO",NSATGPS+NSATGLO,NSATGPS+NSATGLO+NSATGAL);
	}
	if(NSYSQZS){
		fprintf_qc(fp,"\n");
		outobsnav(fp,"       QZSS",NSATGPS+NSATGLO+NSATGAL,NSATGPS+NSATGLO+NSATGAL+NSATQZS);
	}
	if(NSYSCMP){
		fprintf_qc(fp,"\n");
		outobsnav(fp,"     BeiDou",NSATGPS+NSATGLO+NSATGAL+NSATQZS,NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP);
	}
	fprintf_qc(fp,"\n");
}
//static void calctotalsatvalue(satvalue *satval){
//	int i=0,index=0,mpindex=0,j;
//	double mptotal[(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1)][NSYS+1]={0};
//	double mpele[(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1)][NSYS+1]={0};
//	char *sysname="";
//	unsigned char obssys;
//
//	satval->nsys=NSYS;
//	if(NSYSGPS+NSYSQZS==2)satval->nsys=NSYS-1;
//	for(i = 0; i < MAXSAT; i++){
//		if(sv[i].hor>0) {
//			satval->totalobs+=sv[i].hor;
//		}
//		index=NSYS;
//		mpindex=NSYS;
//		sysname="";
//		obssys=0;
//		if(i<NSATGPS) {
//			mpindex=index=NSYSGPS-1;
//			sysname="GPS";
//			obssys=obs_sys[0];
//			if(NSYSGPS+NSYSQZS==2) {
//				sysname="GPS + QZSS";
//			}
//		}
//		else if(i<NSATGPS+NSATGLO) {
//			mpindex=index=NSYSGPS+NSYSGLO-1;
//			sysname="GLONASS";
//			obssys=obs_sys[1];
//		}
//		else if(i<NSATGPS+NSATGLO+NSATGAL) {
//			mpindex=index=NSYSGPS+NSYSGLO+NSYSGAL-1;
//			sysname="Galileo";
//			obssys=obs_sys[2];
//		}
//		else if(i<NSATGPS+NSATGLO+NSATGAL+NSATQZS) {
//			mpindex=index=NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS-1;
//			sysname="QZSS";
//			obssys=obs_sys[3];
//			if(NSYSGPS+NSYSQZS==2) {
//				mpindex=NSYSGPS-1;
//				sysname="GPS + QZSS";
//				obssys=(obs_sys[3]||obs_sys[0]);
//			}
//		}
//		else if(i<NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP) {
//			mpindex=index=NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS+NSYSCMP-1;
//			sysname="BeiDou";
//			obssys=obs_sys[5];
//		}
//
//		if(index<0 || index>NSYS){
//			mpindex=index=NSYS;
//			sysname="";
//			obssys=0;
//		}
//
//		satval->totalmask[NSYS]+=sv[i].mask;
//		satval->totalcomp[NSYS]+=sv[i].compl;
//		if(mpindex!=NSYS){
//			satval->totalmask[mpindex]+=sv[i].mask;
//			satval->totalcomp[mpindex]+=sv[i].compl;
//		}
//
//		satval->totalbelow[index] += sv[i].belowmask;
//		satval->mpsys[mpindex] = sysname;
//		satval->obssys[mpindex]=obssys;
//
//		if(sv[i].iodslipab>0) satval->totaliodslipab+=sv[i].iodslipab;
//		if(sv[i].iodslipbe>0) satval->totaliodslipbe+=sv[i].iodslipbe;
//		if(sv[i].iodmpslipab>0) {
//			satval->totaliodmpslipab+=sv[i].iodmpslipab;
//			if(mpindex!=NSYS) satval->totalslipaab[mpindex]+=sv[i].iodmpslipab;
//		}
//		if(sv[i].iodmpslipbe>0) satval->totaliodmpslipbe+=sv[i].iodmpslipbe;
//		if(obs_sat[i]==1){
//			(satval->totalsat)++;
//		}
//		for(j=0;j<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1);j++){
//			satval->mp[j].n[mpindex]+=sv[i].mp[j].num;
//			mptotal[j][mpindex]+=sv[i].mp[j].rms;
//			mpele[j][mpindex]+=sv[i].mp[j].ele;
//			satval->mp[j].above[mpindex]+=sv[i].mp[j].abcomp;
//			satval->mp[j].below[mpindex]+=sv[i].mp[j].belcomp;
//		}
//		for(j=0;j<NFREQ;j++){
//			satval->lli[j].above+=sv[i].lli[j].above;
//			satval->lli[j].below+=sv[i].lli[j].below;
//		}
//	}
//	for(i=0;i<=NSYS;i++) {
//		for(j=0;j<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1);j++){
//			satval->mp[j].rms[i]=satval->mp[j].n[i]>0?sqrt(mptotal[j][i]/satval->mp[j].n[i]):0.0;
//			satval->mp[j].ele[i]=satval->mp[j].n[i]>0?sqrt(mpele[j][i]/satval->mp[j].n[i]):0.0;
//		}
//	}
//}

static void initobs(obsd_t *obs){
	int i,j;
	for(i=0;i<MAXSAT;i++){
		for (j=0;j<NFREQ+NEXOBS;j++) {
			obs[i].P[j]=obs[i].L[j]=obs[i].SN[j]=0.0; obs[i].D[j]=0.0f;
			obs[i].SNR[j]=obs[i].LLI[j]=obs[i].code[j]=0;
		}
		for (j=0;j<MAXOBSTYPE;j++) {
			obs[i].value[j]=0;
		}
		for (j=0;j<NFREQ-1;j++) {
			obs[i].ion[j]=0;
			obs[i].iod[j]=0;
			sv[i].ion[j].slip=0;
		}
	}
}
//static void outsnsum(FILE *fp,hist_data *ele_hist){
//	double maxmeanS1=0, maxmeanS2=0, totals1, totals2, sd, mean, maxmean[2]={0,0};
//	int tick=1, i, j, k,n;
//	char out[1024];
//
//	for(i=0;i<HIST_VER;i++){
//		totals1=0;
//		for(k=0;k<ele_hist[i].S1n;k++){
//			totals1+=pow(ele_hist[i].valueS1[k]-ele_hist[i].meanS1,2);
//		}
//		ele_hist[i].sds1=0;
//		if(ele_hist[i].S1n!=0) ele_hist[i].sds1=sqrt(totals1/ele_hist[i].S1n);
//		totals2=0;
//		for(k=0;k<ele_hist[i].S2n;k++){
//			totals2+=pow(ele_hist[i].valueS2[k]-ele_hist[i].meanS2,2);
//		}
//		ele_hist[i].sds2=0;
//		if(ele_hist[i].S2n!=0) ele_hist[i].sds2=sqrt(totals2/ele_hist[i].S2n);
//	}
//	for(i=0;i<HIST_VER;i++){
//		if(maxmean[0]<ele_hist[i].meanS1) maxmean[0]=ele_hist[i].meanS1;
//		if(maxmean[0]<ele_hist[i].sds1) maxmean[0]=ele_hist[i].sds1;
//		if(maxmean[1]<ele_hist[i].meanS2) maxmean[1]=ele_hist[i].meanS2;
//		if(maxmean[1]<ele_hist[i].sds2) maxmean[1]=ele_hist[i].sds2;
//	}
//	for(j=0;j<2;j++){
//		while(tick*40<maxmean[j]){
//			tick++;
//		}
//		sprintf(out,"\nS/N L%d summary (per elevation bin):\n",j+1);
//		fprintf_qc(fp,out);
//		sprintf(out,"elev (deg)  tot SN%d sig    mean         %2d|0      %2d|0      %2d|0      %2d|0\n",j+1,tick,tick*2,tick*3,tick*4);
//		fprintf_qc(fp,out);
//		for(i=18;i>=0;i--){
//			if(i!=0){
//				sprintf(out," %2d - %2d",(i-1)*5,i*5);
//				fprintf_qc(fp,out);
//			}
//			else{
//				fprintf_qc(fp,"    <  0");
//			}
//			switch(j){
//				case 0:	n=ele_hist[i].S1n; sd=ele_hist[i].sds1; mean=ele_hist[i].meanS1; break;
//				case 1:	n=ele_hist[i].S2n; sd=ele_hist[i].sds2; mean=ele_hist[i].meanS2; break;
//				default: n=0;sd=0;mean=0;
//			}
//			sprintf(out,"  %5d %7.3f %8.3f ", n, sd, mean);
//			fprintf_qc(fp,out);
//			for(k=0;k<40;k++){
//				if(k+1 <= (sd+tick*0.5)/tick && k+1 <= (mean+tick*0.5)/tick){
//					fprintf_qc(fp,"#");
//				}
//				else if(k+1 <= (sd+tick*0.5)/tick){
//					fprintf_qc(fp,"=");
//				}
//				else if(k+1 <= (mean+tick*0.5)/tick){
//					fprintf_qc(fp,"|");
//				}
//			}
//			fprintf_qc(fp,"\n");
//		}
//	}
//
//}

static void outsvcomp(FILE *fp, int *totalbelow, double elmin){
	int i, j, k, totalno, start, end, index=0;
	char nav_mark, id[64];
	char out[1024];

	for(j=0;j<NUMSYS;j++){
		if(j==0) {start=0;end=NSATGPS;
			if(NSYSGPS==0||obs_sys[0]==0)continue;}
		else if(j==1) {start=end;end+=NSATGLO;
			if(NSYSGLO==0||obs_sys[1]==0)continue;}
		else if(j==2) {start=end;end+=NSATGAL;
			if(NSYSGAL==0||obs_sys[2]==0)continue;}
		else if(j==3) {start=end;end+=NSATQZS;
			if(NSYSQZS==0||obs_sys[3]==0)continue;}
		else if(j==5) {start-=NSATCMP;end+=NSATCMP;
			if(NSYSCMP==0||obs_sys[5]==0)continue;}
		else{start=end+NSATCMP;end=MAXSAT;
			if(obs_sys[4]==0)continue;}
		fprintf_qc(fp,"\n SV  #+hor <ele> #+mask <ele> #reprt #compl");
		for(i=0;i<obstype[j];i++){
			sprintf(out," %5s ",type[j][i]);
			fprintf_qc(fp,out);
		}
		fprintf_qc(fp,"\n--- ------ ----- ------ ----- ------ ------");
		for(i=0;i<obstype[j];i++){
			 fprintf_qc(fp," ------");
		}
		fprintf_qc(fp,"\n");
		for(i=start;i<end;i++){
			if(obs_sat[i]||sv[i].hor>0){
				nav_mark=' ';
				if(nav_sat[i]==0) nav_mark='*';
				satno2id(i+1,id);
				sprintf(out,"%s%c%6d %5.2f %6d %5.2f %6d %6d ",
					id,nav_mark,sv[i].hor,sv[i].ele_hor,sv[i].mask, sv[i].ele_mask,sv[i].reprt,sv[i].compl);
				fprintf_qc(fp,out);
				for(k=0;k<obstype[j];k++){
					sprintf(out,"%6d ",sv[i].count[k]);
					fprintf_qc(fp,out);
				}
				fprintf_qc(fp,"\n");
			}
		}
	}
	fprintf_qc(fp,"   * = SV with no NAV info\n");

	for(j=0;j<NUMSYS;j++){
		if(j==0) {start=0;end=NSATGPS;
			if(NSYSGPS==0||obs_sys[0]==0)continue;
			fprintf_qc(fp,"\n-- GPS --\n");
			index=NSYSGPS-1;}
		else if(j==1) {start=end;end+=NSATGLO;
			if(NSYSGLO==0||obs_sys[1]==0)continue;
			fprintf_qc(fp,"\n-- GLONASS --\n");
			index=NSYSGPS+NSYSGLO-1;}
		else if(j==2) {start=end;end+=NSATGAL;
			if(NSYSGAL==0||obs_sys[2]==0)continue;
			fprintf_qc(fp,"\n-- Galileo --\n");
			index=NSYSGPS+NSYSGLO+NSYSGAL-1;}
		else if(j==3) {start=end;end+=NSATQZS;
			if(NSYSQZS==0||obs_sys[3]==0)continue;
			fprintf_qc(fp,"\n-- QZSS --\n");
			index=NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS-1;}
		else if(j==5) {start-=NSATCMP;end+=NSATCMP;
			if(NSYSCMP==0||obs_sys[5]==0)continue;
			fprintf_qc(fp,"\n-- BeiDou --\n");
			index=NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS+NSYSCMP-1;}
		else{start=end+NSATCMP;end=MAXSAT;
			if(obs_sys[4]==0)continue;
			fprintf_qc(fp,"\n-- SBAS --\n");
			index=NSYS;}
		if(index<0)index=NSYS;
		sprintf(out,"Obs below mask ( %5.2f deg)   : %6d\n", elmin*R2D, totalbelow[index]);
		fprintf_qc(fp,out);
		for(i=0;i<obstype[j];i++){
			totalno=0;
			for(k=start;k<end;k++){
				totalno=sv[k].count[i+obstype[j]];
			}
			sprintf(out,"Obs above mask w/ no %-8s : %6d\n",type[j][i],totalno);
			fprintf_qc(fp,out);
		}
	}
}
//static void outrmshist(FILE *fp, hist_data *ele_hist, int dtype,int index){
//	int tick1=1,tick2=1,k,i,n=0,slps=0;
//	double maxval=0,maxval1=0,maxval2=0,rms=0,slper=0,scale=0.5;
//	char out[1024];
//
//	for(i=0;i<HIST_VER;i++){
//		switch(dtype){
//			case 0:	rms=ele_hist[i].ionn+ele_hist[i].ionslps>0?sqrt(ele_hist[i].ionrms/(ele_hist[i].ionn+ele_hist[i].ionslps)):0;
//				slper=ele_hist[i].ionn>0?ele_hist[i].ionslps*100.0/ele_hist[i].ionn:0;
//				scale=0.5;
//				break;
//			default:rms=ele_hist[i].mp[dtype-1].n[index]+ele_hist[i].mp[dtype-1].slps[index]>0?
//						sqrt(ele_hist[i].mp[dtype-1].rms[index]/(ele_hist[i].mp[dtype-1].n[index]+ele_hist[i].mp[dtype-1].slps[index])):0;
//				slper=ele_hist[i].mp[dtype-1].n[index]>0?ele_hist[i].mp[dtype-1].slps[index]*100.0/ele_hist[i].mp[dtype-1].n[index]:0;
//				scale=0.05;
//				break;
//		}
//		if(maxval1<slper) maxval1=slper;
//		if(maxval2<rms) maxval2=rms;
//	}
//
//	while(tick1*0.5*40<maxval1){
//		tick1++;
//	}
//	while(tick2*scale*40<maxval2){
//		tick2++;
//	}
//	switch(dtype){
//		case 0:
//			sprintf(out,"elev (deg)  tot slps <ION rms, m>      %2d=%%%%      %2.0f|m      %2d=%%%%      %2.0f|m\n",
//				tick1*5,tick2*scale*20,tick1*5*3,tick2*2*scale*20);
//			fprintf_qc(fp,out);
//			break;
//		default:
//			sprintf(out,"elev (deg)  tot slps <MP%-3s rms, m>     %2d=%%%%      %2.0f|m      %2d=%%%%      %2.0f|m\n",
//				MP_Name[dtype-1],tick1*5,tick2*scale*20,tick1*5*3,tick2*2*scale*20);
//			fprintf_qc(fp,out);
//			break;
//	}
//	for(i=18;i>=0;i--){
//		if(i!=0){
//			sprintf(out," %2d - %2d ",(i-1)*5,i*5);
//			fprintf_qc(fp,out);
//		}
//		else{
//			fprintf_qc(fp,"    <  0*");
//		}
//		switch(dtype){
//			case 0:	n=ele_hist[i].ionn; slps=ele_hist[i].ionslps; rms=n+slps>0?sqrt(ele_hist[i].ionrms/(n+slps)):0; break;
//			default: n=ele_hist[i].mp[dtype-1].n[index]; slps=ele_hist[i].mp[dtype-1].slps[index]; rms=n+slps>0?sqrt(ele_hist[i].mp[dtype-1].rms[index]/(n+slps)):0; break;
//		}
//		slper=n>0?slps*100.0/n:0;
//		sprintf(out," %5d %4d %10.6f  ", n+slps, slps, rms);
//		fprintf_qc(fp,out);
//		for(k=0;k<40;k++){
//			if(k+1 <= (slper+tick1*0.5*0.5)/(tick1*0.5) && k+1 <= (rms+tick2*0.5*scale)/(tick2*scale)){
//				fprintf_qc(fp,"#");
//			}
//			else if(k+1 <= (slper+tick1*0.5*0.5)/(tick1*0.5)){
//				fprintf_qc(fp,"=");
//			}
//			else if(k+1 <= (rms+tick2*0.5*scale)/(tick2*scale)){
//				fprintf_qc(fp,"|");
//			}
//		}
//		fprintf_qc(fp,"\n");
//	}
//	fprintf_qc(fp,"     * or unknown elevation\n");
//}

//static void outmpsum(FILE *fp, int dtype, double elmin, double elcomp, satvalue satval){
//	int i,mpab=0,mpbel=0,obs=0;
//	char nav_mark;
//	char id[64];
//	double ele=0,rms=0;
//	char out[1024];
//
//	sprintf(out,"MP%-3s RMS summary (per SV):\n",MP_Name[dtype]);
//	fprintf_qc(fp,out);
//	fprintf_qc(fp,"                                      slips  L1 rx  L2 rx  slips  L1 rx  L2 rx\n");
//	sprintf(out," SV  obs>%2.0f  # del <elev> MP%-3s rms [m]  < %.0f   < %.0f   < %.0f   > %.0f   > %.0f   > %.0f\n",
//		elmin*R2D, MP_Name[dtype], elcomp, elcomp, elcomp, elcomp, elcomp, elcomp);
//	fprintf_qc(fp,out);
//	for(i=0;i<MAXSAT;i++){
//		if(obs_sat[i]){
//			nav_mark=' ';
//			if(nav_sat[i]==0) nav_mark='*';
//			satno2id(i+1,id);
//			sprintf(out,"%s%c%6d %6d %7.2f %9.6f %6d %6d %6d %6d %6d %6d \n",
//				id,nav_mark,sv[i].mp[dtype].num+sv[i].mp[dtype].del,sv[i].mp[dtype].del,
//				sv[i].mp[dtype].num!=0?sv[i].mp[dtype].ele/sv[i].mp[dtype].num*R2D:0,
//				sv[i].mp[dtype].num!=0?sqrt(sv[i].mp[dtype].rms/sv[i].mp[dtype].num):0,
//				sv[i].mp[dtype].belcomp,sv[i].lli[0].below,	sv[i].lli[1].below,
//				sv[i].mp[dtype].abcomp,sv[i].lli[0].above,sv[i].lli[1].above);
//			fprintf_qc(fp,out);
//		}
//	}
//	fprintf_qc(fp,"   * = SV with no NAV info\n");
//
//	for(i=0;i<satval.nsys;i++){
//		if(satval.obssys[i]==0) continue;
//		sprintf(out,"\n-- %s --\n",satval.mpsys[i]);
//		fprintf_qc(fp,out);
//
//		rms=satval.mp[dtype].rms[i]; ele=satval.mp[dtype].ele[i]*R2D; obs=satval.mp[dtype].n[i];
//		mpab=satval.mp[dtype].above[i]; mpbel=satval.mp[dtype].below[i];
//
//		sprintf(out,"mean MP%-3s rms        : %.6f m\n", MP_Name[dtype],rms);
//		fprintf_qc(fp,out);
//		sprintf(out,"total mean elevation  : %.2f degrees\n",ele);
//		fprintf_qc(fp,out);
//		sprintf(out,"# MP%-3s obs > %2.0f*     : %6d\n", MP_Name[dtype], elmin*R2D, obs);
//		fprintf_qc(fp,out);
//		sprintf(out,"# qc MP%-3s slips < %.0f : %6d\n", MP_Name[dtype], elcomp, mpbel);
//		fprintf_qc(fp,out);
//		sprintf(out,"# Rvr L1   slips < %.0f : %6d\n", elcomp, satval.lli[0].below);
//		fprintf_qc(fp,out);
//		sprintf(out,"# Rvr L2   slips < %.0f : %6d\n", elcomp, satval.lli[1].below);
//		fprintf_qc(fp,out);
//		sprintf(out,"# qc MP%-3s slips > %.0f : %6d\n", MP_Name[dtype], elcomp, mpab);
//		fprintf_qc(fp,out);
//		sprintf(out,"# Rvr L1   slips > %.0f : %6d\n", elcomp, satval.lli[0].above);
//		fprintf_qc(fp,out);
//		sprintf(out,"# Rvr L2   slips > %.0f : %6d\n", elcomp, satval.lli[1].above);
//		fprintf_qc(fp,out);
//		fprintf_qc(fp,"   * or SV with no NAV info\n");
//	}
//}

//static void outrmsmp(FILE *fp, satvalue satval, int ponitmp){
//	int i,j;
//	char out[1024];
//
//	for(i=0;i<satval.nsys;i++){
//		if(satval.obssys[i]==0) continue;
//		sprintf(out,"\n-- %s --\n",satval.mpsys[i]);
//		fprintf_qc(fp,out);
//		for(j=0;j<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1);j++){
//			sprintf(out,"Moving average MP%-3s    : %9.6f m\n",MP_Name[j], satval.mp[j].rms[i]);
//			fprintf_qc(fp,out);
//		}
//	}
//	sprintf(out,"\nPoints in MP moving avg : %d\n", ponitmp);
//	fprintf_qc(fp,out);
//
//}
//static void outslipnum(FILE *fp, double iodsign, double elmin, satvalue satval, double interval){
//	char out[1024];
//
//	sprintf(out,"IOD signifying a slip   : >%.1f cm/minute\n",iodsign);
//	fprintf_qc(fp,out);
//	sprintf(out,"IOD slips <  %4.1f deg*  : %6d\n",elmin*R2D,satval.totaliodslipbe);
//	fprintf_qc(fp,out);
//	sprintf(out,"IOD slips >  %4.1f deg   : %6d\n",elmin*R2D,satval.totaliodslipab);
//	fprintf_qc(fp,out);
//	sprintf(out,"IOD or MP slips <  %4.1f*: %6d\n",elmin*R2D,satval.totaliodmpslipbe);
//	fprintf_qc(fp,out);
//	sprintf(out,"IOD or MP slips >  %4.1f : %6d\n",elmin*R2D,satval.totaliodmpslipab);
//	fprintf_qc(fp,out);
//	fprintf_qc(fp," * or unknown elevation\n");
//}
//static void outsum(FILE *fp, double *ep1, double *ep2, satvalue satval){
//	int i,j;
//	char out[1024];
//
//	fprintf_qc(fp,"                  first epoch    last epoch    hrs   dt  #expt  #have   %% ");
//	for(i=0;i<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1);i++){
//		sprintf(out," mp%-3s",MP_Name[i]);
//		fprintf_qc(fp,out);
//	}
//	fprintf_qc(fp," o/slps\n");
//
//	for(i=0;i<satval.nsys;i++){
//		if(satval.obssys[i]==0) continue;
//		sprintf(out,"SUM %10s  %2.0f %2.0f %2.0f %02.0f:%02.0f %2.0f %2.0f %2.0f %02.0f:%02.0f %5.2f %3.0f %6d %6d %3.0f",
//			satval.mpsys[i],fmod(ep1[0],100),ep1[1],ep1[2],ep1[3],ep1[4],fmod(ep2[0],100),ep2[1],ep2[2],ep2[3],ep2[4],total_h,interval,
//			satval.totalmask[i], satval.totalcomp[i], satval.totalmask[i]>0?(double)satval.totalcomp[i]*100.0/(double)satval.totalmask[i]:0);
//			fprintf_qc(fp,out);
//		for(j=0;j<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1);j++){
//			sprintf(out," %5.2f", satval.mp[j].rms[i]);
//			fprintf_qc(fp,out);
//		}
//		sprintf(out," %6d\n",satval.totalslipaab[i]>0?satval.totalcomp[i]/satval.totalslipaab[i]:0);
//		fprintf_qc(fp,out);
//	}
//	fprintf_qc(fp,"\n\n");
//}
//
//static void outantenna(FILE *fp, char** infile, int filenum, double* ep1, double* ep2, double* rr, double interval){
//	int count=0,i;
//	double antpos[3],dms1[3],dms2[3],dist[3];
//	char out[1024];
//   	const char *abbreviated_name[]={
//		"","JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC",""
//	};
//
//	fprintf_qc(fp,"\n*********************\n");
//	fprintf_qc(fp,"QC of RINEX  file(s) : ");
//	if(infile[0]) {
//		sprintf(out,"%s",infile[0]);
//		fprintf_qc(fp,out);
//	}
//	fprintf_qc(fp,"\ninput RnxNAV file(s) : ");
//	for(i=1;i<filenum;i++){
//		if(infile[i]) {
//			if(count>0){
//				fprintf_qc(fp,"\n                       ");
//			}
//			sprintf(out,"%s",infile[i]);
//			fprintf_qc(fp,out);
//			count++;
//		}
//	}
//	fprintf_qc(fp,"\n*********************\n");
//
//	sprintf(out,"\n4-character ID          : %s",stas[0].name);
//	fprintf_qc(fp,out);
//	if(stas[0].marker[0]) {
//		sprintf(out," (# = %s)",stas[0].marker);
//		fprintf_qc(fp,out);
//	}
//	sprintf(out,"\nReceiver type           : %s",stas[0].rectype);
//	fprintf_qc(fp,out);
//	if(stas[0].recsno[0]) {
//		sprintf(out," (# = %s)",stas[0].recsno);
//		fprintf_qc(fp,out);
//	}
//	if(stas[0].recver[0]) {
//		sprintf(out," (fw = %s)",stas[0].recver);
//		fprintf_qc(fp,out);
//	}
//	sprintf(out,"\nAntenna type            : %s",stas[0].antdes);
//	fprintf_qc(fp,out);
//	if(stas[0].antsno[0]){
//		sprintf(out," (# = %s)",stas[0].antsno);
//		fprintf_qc(fp,out);
//	}
//	fprintf_qc(fp,"\n");
//
//	sprintf(out,"\nTime of start of window : %04.0f %s %02.0f  %02.0f:%02.0f:%06.3f\n",ep1[0],abbreviated_name[(int)ep1[1]],ep1[2],ep1[3],ep1[4],ep1[5]);
//	fprintf_qc(fp,out);
//	sprintf(out,"Time of  end  of window : %04.0f %s %02.0f  %02.0f:%02.0f:%06.3f\n",ep2[0],abbreviated_name[(int)ep2[1]],ep2[2],ep2[3],ep2[4],ep2[5]);
//	fprintf_qc(fp,out);
//	sprintf(out,"Time line window length : %.2f hour(s), ticked every %.1f hour(s)\n", total_h,tick);
//	fprintf_qc(fp,out);
//
//	sprintf(out,"  antenna WGS 84 (xyz)  : %.4f %.4f %.4f (m)\n",rr[0],rr[1],rr[2]);
//	fprintf_qc(fp,out);
//
//	ecef2pos(rr,antpos);
//	antpos[2]-=geoidh(antpos);
//    deg2dms(antpos[0]*R2D,dms1);
//    deg2dms(antpos[1]*R2D,dms2);
//
//	sprintf(out,"  antenna WGS 84 (geo)  : N %.0f deg %.0f' %.2f\"  E %.0f deg %.0f' %.2f\"\n",dms1[0],dms1[1],dms1[2],dms2[0],dms2[1],dms2[2]);
//	fprintf_qc(fp,out);
//	sprintf(out,"  antenna WGS 84 (geo)  :  %10.6f deg   %10.6f deg\n",antpos[0]*R2D,antpos[1]*R2D);
//	fprintf_qc(fp,out);
//	sprintf(out,"          WGS 84 height : %.4f m\n",antpos[2]);
//	fprintf_qc(fp,out);
//
//	dist[0] = rr[0] - stas[0].pos[0];
//	dist[1] = rr[1] - stas[0].pos[1];
//	dist[2] = rr[2] - stas[0].pos[2];
//	sprintf(out,"|qc - header| position  : %.4f m\n",norm(dist,3));
//	fprintf_qc(fp,out);
//	sprintf(out,"Observation interval    : %.4f seconds\n",interval);
//	fprintf_qc(fp,out);
//
//}
//static void outtotalobs(FILE *fp, int repeatnum, satvalue satval, double elmin, int duplnum){
//	char out[1024];
//
//	sprintf(out,"Poss. # of obs epochs   : %6d\n",posepoch);
//	fprintf_qc(fp,out);
//	sprintf(out,"Epochs w/ observations  : %6d\n",nepoch);
//	fprintf_qc(fp,out);
//	sprintf(out,"Epochs repeated         : %6d  (%.2f%%%%)\n",repeatnum,(nepoch==0)?0.0:(double)((double)repeatnum/nepoch*100));
//	fprintf_qc(fp,out);
//	sprintf(out,"Possible obs >   0.0 deg: %6d\n",satval.totalobs);
//	fprintf_qc(fp,out);
//	sprintf(out,"Possible obs >  %4.1f deg: %6d\n",elmin*R2D,satval.totalmask[NSYS]);
//	fprintf_qc(fp,out);
//	sprintf(out,"Complete obs >  %4.1f deg: %6d\n",elmin*R2D,satval.totalcomp[NSYS]);
//	fprintf_qc(fp,out);
//	sprintf(out," Deleted obs >  %4.1f deg: %6d\n",elmin*R2D,satval.totaldel);
//	fprintf_qc(fp,out);
//	sprintf(out,"  Masked obs <  %4.1f deg: %6d\n",elmin*R2D,satval.totalbelowcomp);
//	fprintf_qc(fp,out);
//	sprintf(out,"Obs w/ SV duplication   : %6d  (within non-repeated epochs)\n",duplnum);
//	fprintf_qc(fp,out);
//}
//static void outparameters(FILE *fp, qc_param param){
//	char out[1024];
//
//	sprintf(out,"Processing parameters are:\n");
//	fprintf_qc(fp,out);
//	sprintf(out,"Maximum ionospheric rate (L1)      : %.2f cm/min\n",param.iodsign);
//	fprintf_qc(fp,out);
//	sprintf(out,"Report data gap greater than       : %.2f min\n",param.gap);
//	fprintf_qc(fp,out);
//	sprintf(out,"Multipath slip sigma threshold     : %.2f m\n",param.mpsigma);
//	fprintf_qc(fp,out);
//	sprintf(out,"Points in MP moving averages       : %d\n",param.ponitmp);
//	fprintf_qc(fp,out);
//	sprintf(out,"Minimum signal to noise for L1     : %.2f\n",param.mins1);
//	fprintf_qc(fp,out);
//	sprintf(out,"Minimum signal to noise for L2     : %.2f\n",param.mins2);
//	fprintf_qc(fp,out);
//	sprintf(out,"Elevation mask (cutoff)            : %.2f degrees\n",param.elmin);
//	fprintf_qc(fp,out);
//	sprintf(out,"Elevation comparison threshold     : %.2f degrees\n",param.elcomp);
//	fprintf_qc(fp,out);
//	sprintf(out,"Width of ASCII summary plot        : %d\n",PLOT_WIDTH);
//	fprintf_qc(fp,out);
//	sprintf(out,"Tolerance for 1-ms clock slips     : %e ms\n\n\n\n",param.tclock);
//	fprintf_qc(fp,out);
//}

//
//static int singlepos(const prcopt_t *opt, sol_t *sol, ssat_t *ssat, const obsd_t *obs, int n, const nav_t *nav)
//{
////	prcopt_t *opt=&rtk->opt;
////    sol_t solb={{0}};
//    int nu;
////	char msg[128]="";
//
//	//trace(3,"rtkpos  : time=%s n=%d\n",time_str(obs[0].time,3),n);
//	//trace(4,"obs=\n"); traceobs(4,obs,n);
//
//	/* count rover/base station observations */
//    for (nu=0;nu   <n&&obs[nu   ].rcv==1;nu++) ;
//   
//    //time=sol->time; /* previous epoch */
//    
//    /* rover position by single point positioning */
//	if (!pntpos(obs,nu,nav,opt,sol,NULL,ssat,msg)) {
//        //errmsg(rtk,"rover point pos error (%s)\n",msg);
//        
//   //     if (!opt->dynamics) {
//			////outsolstat(rtk, obs);
//   //         return 0;
//   //     }
//    }
//
//    //if (time.time!=0) rtk->tt=timediff(rtk->sol.time,time);
//    
////    /* single point positioning */
////    if (opt->mode==PMODE_SINGLE) {
////        //outsolstat(rtk);
//////        outsolstat(rtk, obs);
////        pntposoutsolstat(rtk, obs, statlevel, fp_stat);
////		return 1;
////	}
//
//	return 1;
//}
//
///* process positioning -------------------------------------------------------*/
//static void procpos2(char *outfile[], const prcopt_t *popt, const solopt_t *sopt,
//                    int mode, char **infile,int filenum,const qc_param param)
//{
//    FILE *fp[OUTFILE_MAX];
//	gtime_t time={0},beforetime={0};
//    sol_t sol={{0}};
//    ssat_t ssat[MAXSAT]={{0}};
////    rtk_t rtk;
//    obsd_t obs[MAXSAT]={{0}};
//	double rb[3]={0},ep1[6],ep2[6],*sn[2],means1=0,means2=0,sds1=0,sds2=0;
//    char buff[MAXSOLMSG+1], id[64];
//    int i,nobs,n,solstatic,pri[]={0,1,2,3,4,5,1,6},j=0,k,rev_clk_slip,week=0;
//	int mpno[]={0,12,15,21,25,51,52};
//	int qc[MAXSAT][PLOT_WIDTH]={0},qc_index=0,dn[MAXSAT*2+1][PLOT_WIDTH]={0};
//	int abmask[MAXSAT][PLOT_WIDTH]={0},clk[PLOT_WIDTH]={0};
//	int tick_sum=0,flag,n_back,s1n=0,s2n=0,repcp=0,del;
//   	const char *abbreviated_name[]={
//		"","JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC",""
//	};
//	char out[1024];
//
//    double *rs,*dts,*var,*azel_,*resp,rr[3],pos[3],e[3],totals1=0,totals2=0,elevation,dift=0,spant=0;
//	double maxmsec, minmsec,sec;
//    int vsat[MAXSAT]={0},svh[MAXSAT];
//    char msg[128]="";
//	hist_data ele_hist[HIST_VER]={0};
//	satvalue satval={0};
//	double rev_clock=0, diff_rc=0, diff_rct=0;
//    int nu;
//
//    trace(3,"procpos : mode=%d\n",mode);
//
//    solstatic=sopt->solstatic&&
//			  (popt->mode==PMODE_STATIC||popt->mode==PMODE_PPP_STATIC);
//
//	for(i=0;i<OUTFILE_MAX;i++){
//		fp[i]  = openfile(outfile[i]);
//	}
//
////	rtkinit(&rtk,popt);
//    rtcm_path[0]='\0';
//	outflag=param.outflag;
//
//	sn[0] = (double *)calloc(obss.n+1,sizeof(double));
//	sn[1] = (double *)calloc(obss.n+1,sizeof(double));
//
//	initvec(ele_hist,obss.n+1);
//	total_h = (interval + diff_se)  / 3600;
//	if(total_h > 6){
//		tick = 3;
//	}
//	else{
//		tick = 1;
//	}
//
//	while ((nobs=inputobs(obs,sol.stat,popt))>=0) {
//
//		/* exclude satellites */
//        for (i=n=0;i<nobs;i++) {
//            if ((satsys(obs[i].sat,NULL)&popt->navsys)&&
//				popt->exsats[obs[i].sat-1]!=1) obs[n++]=obs[i];
//		}
//		n_back = n;
//		for(i=0;i<MAXSAT;i++){
//			flag=0;
//			for(k=0;k<n;k++){
//				if(i+1==obs[k].sat){
//					flag=1;
//					break;
//				}
//			}
//			if(flag==0){
//				obs[n_back++].sat=i+1;
//			}
//		}
//
//        if (n<=0) continue;
//
//
//        for (nu=0;nu<n&&obs[nu].rcv==1;nu++) ;
//        pntpos(obs,nu,&navs,popt,&sol,NULL,ssat,msg);
//
//		if(j!=0){
//			dift=timediff(obs->time,beforetime);
//		}
//		spant+=dift;
//		beforetime=obs->time;
//		j++;
//
//		qc_index = (int)(spant * PLOT_WIDTH / (3600.0 * total_h));
//		if(qc_index>=PLOT_WIDTH)qc_index=71;
//		dn[MAXSAT*2][qc_index]++;
//
//		if(dift!=0 && dift-interval>0){
//			clk[qc_index]=2;
//		}
//		else if(dift!=0 && dift-interval<0){
//			clk[qc_index]=MAX(clk[qc_index],1);
//		}
//
//		// obsなしSVのel計算
//	    rs=mat(6,n_back); dts=mat(2,n_back); var=mat(1,n_back); azel_=zeros(2,n_back); resp=mat(1,n_back);
//		satposs(obs[0].time, obs,n_back,&navs,popt->sateph,rs,dts,var,svh);
//		for(k=0;k<3;k++){
//			rr[k]=sol.rr[k];
//		}
//	    ecef2pos(rr,pos);
//		for(i=n;i<MAXSAT;i++){
//			azel_[(obs[i].sat-1)*2]=0;
//			azel_[(obs[i].sat-1)*2+1]=0;
//			if(geodist(rs+i*6,rr,e)<=0){
//				continue;
//			}
//			satazel(pos,e,azel_+(obs[i].sat-1)*2);
//		}
//
//		for(i=1;i<OUTFILE_MAX;i++){
//			if(fp[i]){
//				sprintf(buff,"%2d",n);
//				fwrite(buff,3,1,fp[i]);
//			}
//		}
//
//		//clock slip
//		rev_clk_slip=0;
//		maxmsec=-99999;
//		minmsec=99999;
//		for(i = 0; i < n;i++){
//
//			elevation = ssat[obs[i].sat-1].azel[1];
//			// MP計算
//			calcmp(&obs[i], param);
//			// ION計算
//			calcion(&obs[i], param.iodsign*0.01, &navs, param.gap);
//			for(k=1;k<OUTFILE_MAX;k++){
//				if(fp[k]){
//					satno2id(obs[i].sat,id);
//					sprintf(buff,"%s ",id);
//					fwrite(buff,4,1,fp[k]);
//				}
//			}
//			// 各観測データをチェック
//			rev_clk_slip+=checksv(obs[i], elevation, &repcp, popt->elmin, abmask, &satval, qc, qc_index,
//				&totals1, &totals2, &s1n, &s2n, sn, obs[i].sat-1, ele_hist, dn, param.iodsign*0.01, param.mpsigma,
//				param.elcomp, param, &maxmsec, &minmsec);
//		}
//		if(rev_clk_slip==0 && fabs(maxmsec-minmsec)<=param.tclock){
//			satval.clkslip++;
//			for(i = 0; i < n;i++){
//				qc[obs[i].sat-1][qc_index]=MIN(qc[obs[i].sat-1][qc_index], ASCII_SYMBOL1_NO);
//				for(j=0;j<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1);j++){
//					sv[obs[i].sat-1].mp[j].constmp=0;
//				}
//			}
//		}
//		else{
//			for(i = 0; i < n;i++){
//				for(j=0;j<(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1);j++){
//					if(sv[obs[i].sat-1].mp[j].calcflag==1 && sv[obs[i].sat-1].mp[j].msecslip==1){
//						satval.msecmp++;
//					}
//				}
//			}
//		}
//
//		for(i = n; i < MAXSAT;i++){
//			// MP計算(初期化)
//			calcmp(&obs[i], param);
//
//			checksv_noobs(azel_[(obs[i].sat-1)*2+1],obs[i].sat-1,dn,abmask,qc,popt->elmin,qc_index);
//		}
//
//		for(i=1;i<OUTFILE_MAX;i++){
//			if(fp[i]){
//				fprintf(fp[i],"\n");
//				for(k = 0; k < n;k++){
//					switch(i){
//						case OUTFILE_AZ: fprintf(fp[i],"%10.4f ",ssat[obs[k].sat-1].azel[0]*R2D);break;
//						case OUTFILE_EL: fprintf(fp[i],"%10.4f ",ssat[obs[k].sat-1].azel[1]*R2D);break;
//						case OUTFILE_ION:
//							if(NFREQ>=2) fprintf(fp[i],"%8.3f%c ",obs[k].ion[0],sv[obs[k].sat-1].ion[0].slip==1?'I':' ');break;
//						case OUTFILE_IOD:
//							if(NFREQ>=2) fprintf(fp[i],"%8.3f%c ",obs[k].iod[0],sv[obs[k].sat-1].ion[0].slip==1?'I':' ');break;
//						case OUTFILE_ION5:
//							if(NFREQ>=3) fprintf(fp[i],"%8.3f%c ",obs[k].ion[1],sv[obs[k].sat-1].ion[1].slip==1?'I':' ');break;
//						case OUTFILE_IOD5:
//							if(NFREQ>=3) fprintf(fp[i],"%8.3f%c ",obs[k].iod[1],sv[obs[k].sat-1].ion[1].slip==1?'I':' ');break;
//						case OUTFILE_MP12:
//							if((NFREQ+NEXOBS)*(NFREQ+NEXOBS-1)>=1) fprintf(fp[i],"%8.3f%c ",obs[k].MP[0],sv[obs[k].sat-1].mp[0].slip==1?'M':' ');
//							break;
//						case OUTFILE_MP21:
//							if((NFREQ+NEXOBS)*(NFREQ+NEXOBS-1)>=2) fprintf(fp[i],"%8.3f%c ",obs[k].MP[1],sv[obs[k].sat-1].mp[1].slip==1?'M':' ');
//							break;
//						case OUTFILE_MP15:
//							if((NFREQ+NEXOBS)*(NFREQ+NEXOBS-1)>=3) fprintf(fp[i],"%8.3f%c ",obs[k].MP[2],sv[obs[k].sat-1].mp[2].slip==1?'M':' ');
//							break;
//						case OUTFILE_MP51:
//							if((NFREQ+NEXOBS)*(NFREQ+NEXOBS-1)>=4) fprintf(fp[i],"%8.3f%c ",obs[k].MP[3],sv[obs[k].sat-1].mp[3].slip==1?'M':' ');
//							break;
//						case OUTFILE_MP25:
//							if((NFREQ+NEXOBS)*(NFREQ+NEXOBS-1)>=5) fprintf(fp[i],"%8.3f%c ",obs[k].MP[4],sv[obs[k].sat-1].mp[4].slip==1?'M':' ');
//							break;
//						case OUTFILE_MP52:
//							if((NFREQ+NEXOBS)*(NFREQ+NEXOBS-1)>=6) fprintf(fp[i],"%8.3f%c ",obs[k].MP[5],sv[obs[k].sat-1].mp[5].slip==1?'M':' ');
//							break;
//					}
//				}
//				fprintf(fp[i],"\n");
//			}
//		}
//
//		initobs(obs);
//    }
//
//	if(fp[OUTFILE_S]){
//		// SVPlotを出力
//		outsvplot(fp[OUTFILE_S],qc,dn,popt->elmin,abmask,ep1,ep2,clk);
//		// rinexfile,time,anntena情報等を出力
//		outantenna(fp[OUTFILE_S],infile,filenum,ep1,ep2,sol.rr,interval);
//
//		calctotalsatvalue(&satval);
//		// satellites
//		outsvnumber(fp[OUTFILE_S],satval.totalsat);
//
//		outtotalobs(fp[OUTFILE_S],repeatnum,satval,popt->elmin,duplnum);
//		// Moving average MP
//		outrmsmp(fp[OUTFILE_S],satval, param.ponitmp);
//
//		if(s1n!=0) means1=totals1/s1n;
//		if(s2n!=0) means2=totals2/s2n;
//		totals1=0;
//		totals2=0;
//		for(i=0;i<s1n;i++){
//			totals1+=pow(sn[0][i]-means1,2);
//		}
//		for(i=0;i<s2n;i++){
//			totals2+=pow(sn[1][i]-means2,2);
//		}
//		if(s1n!=0) sds1=sqrt(totals1/s1n);
//		if(s2n!=0) sds2=sqrt(totals2/s2n);
//		sprintf(out,"Mean S1 S2              : %.2f (sd=%.2f n=%d) %.2f (sd=%.2f n=%d)\n",means1,sds1,s1n,means2,sds2,s2n);
//		fprintf_qc(fp[OUTFILE_S],out);
//		sec=time2gpst(obsts,&week);
//		sprintf(out,"Freq no. and timecode   : %d %d ffffff\n",NFREQ, week*7+(int)(sec/86400));
//		fprintf_qc(fp[OUTFILE_S],out);
//		sprintf(out,"Report gap > than       : %.2f minute(s)\n",param.gap);
//		fprintf_qc(fp[OUTFILE_S],out);
//		sprintf(out,"epochs w/ msec clk slip : %d\n",satval.clkslip);
//		fprintf_qc(fp[OUTFILE_S],out);
//		sprintf(out,"other msec mp events    : %d (: %d)   {expect ~= 1:50}\n",satval.msecmp, satval.totalmp);
//		fprintf_qc(fp[OUTFILE_S],out);
//
//		outslipnum(fp[OUTFILE_S], param.iodsign, popt->elmin, satval, interval);
//
//		outsum(fp[OUTFILE_S], ep1, ep2, satval);
//
//		outparameters(fp[OUTFILE_S], param);
//
//		sprintf(out,"Observations start   : %04.0f %s %02.0f  %02.0f:%02.0f:%06.3f\n",ep1[0],abbreviated_name[(int)ep1[1]],ep1[2],ep1[3],ep1[4],ep1[5]);
//		fprintf_qc(fp[OUTFILE_S],out);
//		sprintf(out,"Observations  end    : %04.0f %s %02.0f  %02.0f:%02.0f:%06.3f\n",ep2[0],abbreviated_name[(int)ep2[1]],ep2[2],ep2[3],ep2[4],ep2[5]);
//		fprintf_qc(fp[OUTFILE_S],out);
//		sprintf(out,"Observation interval : %.4f second(s)\n",interval);
//		fprintf_qc(fp[OUTFILE_S],out);
//
//		outsvcomp(fp[OUTFILE_S], satval.totalbelow, popt->elmin);
//
//		del=satval.totaldel;
//		for(i=0;i<NSYS;i++){
//			del+=satval.totalbelow[i];
//		}
//		sprintf(out,"\nObs reported w/ code | phase : %6d\n",repcp);
//		fprintf_qc(fp[OUTFILE_S],out);
//		sprintf(out,"Obs deleted (any reason)     : %6d\n",del);
//		fprintf_qc(fp[OUTFILE_S],out);
//		sprintf(out,"Obs complete                 : %6d\n",satval.totalcomp[NSYS]);
//		fprintf_qc(fp[OUTFILE_S],out);
//
//		fprintf_qc(fp[OUTFILE_S],"\n");
//		outrmshist(fp[OUTFILE_S], ele_hist,0,0);
//		fprintf_qc(fp[OUTFILE_S],"\n\n");
//
//		for(j=1;j<=(NFREQ+NEXOBS)*(NFREQ+NEXOBS-1);j++){
//			fprintf_qc(fp[OUTFILE_S],"\n\n");
//			outmpsum(fp[OUTFILE_S],j-1,popt->elmin,param.elcomp,satval);
//			fprintf_qc(fp[OUTFILE_S],"\n");
//			for(i=0;i<satval.nsys;i++){
//				if(satval.obssys[i]==0) continue;
//				sprintf(out,"\n--- %s ---\n",satval.mpsys[i]);
//				fprintf_qc(fp[OUTFILE_S],out);
//				outrmshist(fp[OUTFILE_S], ele_hist,j,i);
//			}
//		}
//		fprintf_qc(fp[OUTFILE_S],"\n");
//		outsnsum(fp[OUTFILE_S], ele_hist);
//
//	}
//
//    //if (mode==0&&solstatic&&time.time!=0.0 && fp[OUTFILE_S]) {
//    //    sol.time=time;
//    //    outsol(fp[OUTFILE_S],&sol,rb,popt,sopt);
//    //}
//	for(i=0;i<OUTFILE_MAX;i++){
//		if(fp[i]) fclose(fp[i]);
//	}
//	for(i=0;i<HIST_VER;i++){
//		free(ele_hist[i].valueS1);
//		free(ele_hist[i].valueS2);
//	}
//	for(i=0;i<MAXSAT;i++){
//		free(sv[i].count);
//	}
//	free(sn[0]);
//	free(sn[1]);
//
////    rtkfree(&rtk);
//}
//
///* execute processing session for qc ---------------------------------------------*/
//static int execses_qc(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
//                   const solopt_t *sopt, const filopt_t *fopt, int flag,
//                   char **infile, const int *index, int n, char *outfile, const qc_param param)
//{
//    FILE *fp;
//    prcopt_t popt_=*popt;
////    char tracefile[1024],statfile[1024];
////	fcb_t *fcb;
//    const char *s1[]={"GPST","UTC","JST"};
//	char ofile[OUTFILE_MAX][1024],s[5],s2[32],*ofiles[OUTFILE_MAX];
//	double ep[6],t1;
//	int i,w1;
//	const char *dname[]={"","Azimuth","Elevation","ION","IOD","ION","IOD","MP12","MP21","MP15","MP51","MP25","MP52","S/N1","S/N2"};
//	const char *fname[]={"",".azi",".ele",".ion",".iod",".ion5",".iod5",".mp12",".mp21",".mp15",".mp51",".mp25",".mp52",".sn1","sn2"};
//
//    trace(3,"execses : n=%d outfile=%s\n",n,outfile);
//
//	for(i=0;i<MAXSAT;i++){
//		obs_sat[i]=0;
//		nav_sat[i]=0;
//	}
//	for(i=0;i<NUMSYS;i++){
//		obs_sys[i]=0;
//	}
//    /* read obs and nav data */
//    if (!readobsnav(ts,te,ti,infile,index,n,&popt_,&obss,&navs,stas)) return 0;
//
//    t1=time2gpst(obsts,&w1);
//    if (sopt->times>=1) obsts=gpst2utc(obsts);
//    if (sopt->times>=1) obste=gpst2utc(obste);
//    if (sopt->times==2) obsts=timeadd(obsts,9*3600.0);
//    if (sopt->times==2) obste=timeadd(obste,9*3600.0);
//    time2epoch(obsts,ep);
//    time2str(obsts,s2,1);
//	sprintf(s,".%02.0fS",fmod(ep[0],100));
//
//	for(i=0;i < OUTFILE_MAX;i++){
//		memset(ofile[i], 0, 1024);
//		strcpy(ofile[i],outfile);
//	}
//
//	for(i=0;i<OUTFILE_MAX;i++){
//		if(i==0) strcat(ofile[i],s);
//		else strcat(ofile[i],fname[i]);
//		ofiles[i]=ofile[i];
//	}
//
//	/* write header to output file */
//    if (flag&&!outhead_qc(ofiles[0],infile,n,&popt_,sopt)) {
//        freeobsnav(&obss,&navs);
//        return 0;
//    }
//	for(i=1;i<OUTFILE_MAX;i++){
//		if(fp = fopen(ofile[i],"w")){
//			fprintf(fp,"%s\nT_SAMP    %.4f\nSTART_TIME : %s %s (week%04d %8.1fs)\n",dname[i],interval,s2,s1[sopt->times],w1,t1);
//			fclose(fp);
//		}
//	}
//
//    iobsu=iobsr=isbs=ilex=revs=aborts=0;
//
//	procpos2(ofiles,&popt_,sopt,0,infile,n, param);
//
//    /* free obs and nav data */
//    freeobsnav(&obss,&navs);
//
//    return aborts?1:0;
//}
//
//static int postpos_qc(gtime_t ts, gtime_t te, double ti, double tu,
//                   const prcopt_t *popt, const solopt_t *sopt,
//				   const filopt_t *fopt, char **infile, int n, char *outfile,
//				   const char *rov, const char *base, const qc_param param)
//{
//    gtime_t tts,tte,ttte;
//    double tunit,tss;
//    int i,j,k,nf,stat=0,week,flag=1,index[MAXINFILE]={0};
//    char *ifile[MAXINFILE],ofile[1024],*ext;
//
//    trace(3,"postpos : ti=%.0f tu=%.0f n=%d outfile=%s\n",ti,tu,n,outfile);
//
////   /* 単一基線結果ファイル出力先をメイン画面から取得 */
////	strcpy( fsb, outfile );
//
//    /* open processing session */
//    if (!openses_qc(popt,sopt,fopt,&navs,&pcvss,&pcvsr)) return -1;
//
//    if (ts.time!=0&&te.time!=0&&tu>=0.0) {
//        if (timediff(te,ts)<0.0) {
//            showmsg("error : no period");
//            closeses_qc(&navs,&pcvss,&pcvsr);
//            return 0;
//        }
//        for (i=0;i<MAXINFILE;i++) {
//            if (!(ifile[i]=(char *)malloc(1024))) {
//                for (;i>=0;i--) free(ifile[i]);
//                closeses_qc(&navs,&pcvss,&pcvsr);
//                return -1;
//            }
//        }
//        if (tu==0.0||tu>86400.0*MAXPRCDAYS) tu=86400.0*MAXPRCDAYS;
//        settspan(ts,te);
//        tunit=tu<86400.0?tu:86400.0;
//        tss=tunit*(int)floor(time2gpst(ts,&week)/tunit);
//
//        for (i=0;;i++) { /* for each periods */
//            tts=gpst2time(week,tss+i*tu);
//            tte=timeadd(tts,tu-DTTOL);
//            if (timediff(tts,te)>0.0) break;
//            if (timediff(tts,ts)<0.0) tts=ts;
//            if (timediff(tte,te)>0.0) tte=te;
//
//            strcpy(proc_rov ,"");
//            strcpy(proc_base,"");
//            if (checkbrk("reading    : %s",time_str(tts,0))) {
//                stat=1;
//                break;
//            }
//            for (j=k=nf=0;j<n;j++) {
//
//                ext=strrchr(infile[j],'.');
//
//                if (ext&&(!strcmp(ext,".rtcm3")||!strcmp(ext,".RTCM3"))) {
//                    strcpy(ifile[nf++],infile[j]);
//                }
//                else {
//                    /* include next day precise ephemeris or rinex brdc nav */
//                    ttte=tte;
//                    if (ext&&(!strcmp(ext,".sp3")||!strcmp(ext,".SP3")||
//                              !strcmp(ext,".eph")||!strcmp(ext,".EPH"))) {
//                        ttte=timeadd(ttte,3600.0);
//                    }
//                    else if (strstr(infile[j],"brdc")) {
//                        ttte=timeadd(ttte,7200.0);
//                    }
//                    nf+=reppaths(infile[j],ifile+nf,MAXINFILE-nf,tts,ttte,"","");
//                }
//                while (k<nf) index[k++]=j;
//
//                if (nf>=MAXINFILE) {
//                    trace(2,"too many input files. trancated\n");
//                    break;
//                }
//            }
//            if (!reppath(outfile,ofile,tts,"","")&&i>0) flag=0;
//
//            /* execute processing session */
//            stat=execses_qc(tts,tte,ti,popt,sopt,fopt,flag,ifile,index,nf,ofile,param);
//
//            if (stat==1) break;
//        }
//        for (i=0;i<MAXINFILE;i++) free(ifile[i]);
//    }
//    else if (ts.time!=0) {
//        for (i=0;i<n&&i<MAXINFILE;i++) {
//            if (!(ifile[i]=(char *)malloc(1024))) {
//                for (;i>=0;i--) free(ifile[i]); return -1;
//            }
//            reppath(infile[i],ifile[i],ts,"","");
//            index[i]=i;
//        }
//        reppath(outfile,ofile,ts,"","");
//
//        /* execute processing session */
//        stat=execses_qc(ts,te,ti,popt,sopt,fopt,1,ifile,index,n,ofile,param);
//
//        for (i=0;i<n&&i<MAXINFILE;i++) free(ifile[i]);
//    }
//    else {
//        for (i=0;i<n;i++) index[i]=i;
//
//        /* execute processing session */
//        stat=execses_qc(ts,te,ti,popt,sopt,fopt,1,infile,index,n,outfile,param);
//    }
//    /* close processing session */
//    closeses_qc(&navs,&pcvss,&pcvsr);
//
//    return stat;
//}

