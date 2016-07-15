/*------------------------------------------------------------------------------
* rtkpos.c : precise positioning
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

static const char rcsid[]="$Id:$";

/* constants/macros ----------------------------------------------------------*/
#define VAR_POS     SQR(30.0) /* initial variance of receiver pos (m^2) */
#define VAR_VEL     SQR(10.0) /* initial variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(10.0) /* initial variance of receiver acc ((m/ss)^2) */
#define VAR_HWBIAS  SQR(1.0)  /* initial variance of h/w bias ((m/MHz)^2) */
#define VAR_GRA     SQR(0.001)/* initial variance of gradient (m^2) */
#define VAR_ISBL    SQR(0.01)/* initial variance of isb(phase) (m^2) */
#define VAR_ISBP    SQR(1.0)/* initial variance of isb(pseudorange) (m^2) */
#define VAR_GL2     SQR(1)/* initial variance of GPS L2P-L2C (m^2) */

#define INIT_ZWD    0.15     /* initial zwd (m) */

#define PRN_HWBIAS  1E-6     /* process noise of h/w bias (m/MHz/sqrt(s)) */
#define GAP_RESION  120      /* gap to reset ionosphere parameters (epochs) */
#define MAXACC      30.0     /* max accel for doppler slip detection (m/s^2) */

#define VAR_HOLDAMB 0.001    /* constraint to hold ambiguity (cycle^2) */

#define TTOL_MOVEB  (1.0+2*DTTOL)
                             /* time sync tolerance for moving-baseline (s) */

#define FIX_THRES_WLNL  0.999       	/* fix threshold (cycle): p0=0.9999 */
#define FIX_THRES       0.129       	/* fix threshold (cycle): p0=0.9999 */
#define MIN_FIX_CNT     10          /* min fix count to fix ambiguity */



/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nfreq)
#define NP(opt)     ((opt)->dynamics==0?3:9)
#define NI(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL(opt)     ((opt)->glomodear!=GLO_ARMODE_AUTO?0:NFREQGLO)

#define NS0(opt)     ((opt)->isb<ISBOPT_EST?0:(opt)->isb==ISBOPT_EST?(NSYS)*2:(NSYS))
#define NS(opt)     (NS0(opt)*NF(opt))


#define N2(opt)     ((opt)->gl2bias!=GL2OPT_EST?0:2)
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt)+NS(opt)+N2(opt))
#define NX(opt)     (NR(opt)+NB(opt))

/* state variable index */
#define II(s,opt)   (NP(opt)+(s)-1)                 /* ionos (s:satellite no) */
#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) /* tropos (r:0=rov,1:ref) */
#define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))   /* receiver h/w bias */

#define IS(sys,type,f,opt) (NP(opt)+NI(opt)+NT(opt)+NL(opt)+NS0(opt)*(f)+((opt)->isb==ISBOPT_EST?(NSYS)*(type):0)+(sys-NSYSGPS))   /* isb */

#define I2(r,opt) (NP(opt)+NI(opt)+NT(opt)+NL(opt)+NS(opt)+(r))   /* gl2bias (r:0=rov,1:ref) */
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1) /* phase bias (s:satno,f:freq) */

/* function prototypes for antest -----------------------------------------------*/
extern int sphharmonic(const double *azel, const prcopt_t *opt, double *sphCoef);

/* global variables ----------------------------------------------------------*/
static int statlevel=0;          /* rtk status output level (0:off) */
static FILE *fp_stat=NULL;       /* rtk status file pointer */
static FILE *fp_antest=NULL;     /* antenna file pointer */
static FILE *fp_phaseerr=NULL;   /* phase-error file pointer */
static int sph_start=0;          /* start index for sphharmonic */
static int base_satno[NFREQ];    /* no of base sat */
static FILE *fp_ion=NULL;        /* rtk rinex ion file pointer */

/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    rtk->x[i]=xi;
    for (j=0;j<rtk->nx;j++) {
        rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=i==j?var:0.0;
    }
}

/* initialize rtk control ------------------------------------------------------
* initialize rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          prcopt_t *opt    I   positioning options (see rtklib.h)
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkinit(rtk_t *rtk, const prcopt_t *opt)
{
	sol_t sol0={{0}};
	ssat_t ssat0={0};
	int i,j,k;
	passd_t pass0={{0}};
	ambc_t ambc0={0};
	/*rtk_t構造体に新たに追加した、gsi_t構造体の内容を初期化するコードを追加する。*/
	gsi_t gsi={0};
	int nx = opt->mode<=PMODE_FIXED?rtknx(opt):pppnx(opt);

	trace(3,"rtkinit :\n");

	rtk->sol=sol0;
	for (i=0;i<6;i++) rtk->rb[i]=0.0;
#ifndef RTKPOS4PCV
	rtk->nx=opt->mode<=PMODE_FIXED?rtknx(opt):pppnx(opt);
#else
	rtk->nx=nx;
	if (opt->antestmode == ANTEST_MODE_PCV) {
		double azel[2] = {0};
		rtk->nx += sphharmonic(azel, opt, NULL);
	}
#endif
	rtk->na=opt->mode<=PMODE_FIXED?rtknr(opt):0;
	rtk->tt=0.0;
	rtk->x=zeros(rtk->nx,1);
	rtk->P=zeros(rtk->nx,rtk->nx);
	rtk->xa=zeros(rtk->na,1);
	rtk->Pa=zeros(rtk->na,rtk->na);
	rtk->nfix=rtk->neb=0;
	for (i=0;i<MAXSAT;i++) {
		rtk->ssat[i]=ssat0;
		rtk->pass[i]=pass0;
		rtk->ambc[i]=ambc0;

		for(j=0;j<NFREQ;j++){
			for(k=0;k<2;k++){
				rtk->passtime[k][i][j].n=0;
				rtk->passtime[k][i][j].nmax=0;
				rtk->passtime[k][i][j].st=NULL;
				rtk->passtime[k][i][j].et=NULL;
			}
		}
	}
	rtk->ux=izeros(rtk->nx,1);
#ifndef RTKPOS4PCV
	rtk->gsi.tint=0;
	rtk->gsi.minsat=0;
	for (i=0;i<MAXSAT;i++) {
		rtk->gsi.ssat[0][i]=0;
		rtk->gsi.ssat[1][i]=0;
	}
	for (i=0;i<3;i++){
		rtk->gsi.rb[i]=0.0;
		rtk->gsi.rr[i]=0.0;
	}
	for (i=0;i<6;i++){
		rtk->gsi.qr[i]=0.0;
		rtk->gsi.rr[i]=0.0;
	}
	rtk->gsi.ndata=0;
	rtk->gsi.ave=0.0;
	rtk->gsi.var=0.0;
	rtk->gsi.ndop=0;
	rtk->gsi.rdop=0.0;
	rtk->gsi.nrej=0;
	rtk->gsi.ratio=0.0;
	rtk->gsi.sys=0;
#endif
	for (i=0;i<MAXERRMSG;i++) rtk->errbuf[i]=0;
	rtk->opt=*opt;
    /*rtk_t構造体に新たに追加した、gsi_t構造体の内容を初期化するコードを追加する。*/
	rtk->gsi=gsi;

	/*L2Cのプライオリティ設定*/
	//setrnxcodepri(1,opt->l2cprior ==0?"CPYWMNDSLX":"SLXCPYWMND");
	if (opt->l2cprior ==1) {
		setcodepri(SYS_GPS, 2, "SLXCPYWMND");
	}
	if (opt->codepri[1][0][0] != '0x00') {
		setcodepri(SYS_GLO, 1, opt->codepri[1][0]);
	}
	if (opt->codepri[1][1][0] != '0x00') {
		setcodepri(SYS_GLO, 2, opt->codepri[1][1]);
	}
#ifdef RTKPOS4PCV
    if (opt->antestmode == ANTEST_MODE_PCV) {
        /* 推定パラメータ＆共分散行列の初期化(球面調和関数の箇所のみ)
		* ※衛星バイアス部分は他で初期化される
        */
        sph_start = nx;
        for (i = sph_start; i < rtk->nx; i++) initx(rtk, 1e-20, rtk->opt.antestvar, i);
    }
#endif
}

/* free rtk control ------------------------------------------------------------
* free memory for rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkfree(rtk_t *rtk)
{
	int i,j,k;
	trace(3,"rtkfree :\n");

	rtk->nx=rtk->na=0;
    free(rtk->x ); rtk->x =NULL;
    free(rtk->P ); rtk->P =NULL;
    free(rtk->xa); rtk->xa=NULL;
	free(rtk->Pa); rtk->Pa=NULL;

	for (i=0;i<MAXSAT;i++) {
		for(j=0;j<NFREQ;j++){
			for(k=0;k<2;k++){
				rtk->passtime[k][i][j].n=0;
				rtk->passtime[k][i][j].nmax=0;
				if(rtk->passtime[k][i][j].st) {free(rtk->passtime[k][i][j].st);rtk->passtime[k][i][j].st=NULL;};
				if(rtk->passtime[k][i][j].et) {free(rtk->passtime[k][i][j].et);rtk->passtime[k][i][j].et=NULL;};
			}
		}
	}
	free(rtk->ux);

}

static void stat_stock(rtk_t *rtk, const int stat, const int nf, const int nu, const obsd_t *obs,
					   const int nv, const double *v, const int *vflg)
{
	/* No.0 変数宣言 */
	int i,n;
	double d;
	int ti;
	double azel[MAXOBS*2],dop[4];

	/* No.1 No.1xは、statがいずれの場合でも更新する。 */
	/* No.11 データ間隔の更新 */
	ti=ROUND(rtk->tt);
	if (ti!=0) {
		if (rtk->gsi.tint==0 || ti<rtk->gsi.tint) rtk->gsi.tint=ti;
	}

	trace(5,"tint: rtk->tt[%f]ti[%d]tint[%f]\n",rtk->tt,ti,rtk->gsi.tint);

	/* No.12 最小衛星数の更新 */
	if (rtk->gsi.minsat==0 || rtk->sol.ns<rtk->gsi.minsat) {
		rtk->gsi.minsat=rtk->sol.ns;
	}

	/* No.13 有効衛星の更新 */
//	for (i=0;i<MAXSAT;i++) {
//		for(j=0;j<NFREQ;j++){
//			if (rtk->ssat[i].vsat[j]) rtk->gsi.ssat[i]=1;
//		}
//	}

	/* No.15 RMSの計算 */
	for (i=0;i<nv;i++) {
		if (vflg[i]&16) continue;
		rtk->gsi.ndata++;
		d=v[i]-rtk->gsi.ave;
		rtk->gsi.ave+=d/rtk->gsi.ndata;
		if (rtk->gsi.ndata>1)
			rtk->gsi.var+=(d*d-rtk->gsi.var)/rtk->gsi.ndata;

		trace(5,"stat_stock RMS ndata[%d] v[%d]=[%f] ave[%f] d[%f] var[%f]\n",rtk->gsi.ndata,i,v[i],rtk->gsi.ave,d, rtk->gsi.var);

	}

	/* No.16 DOPの計算 */
	/* 新点のobsに対応したazelの中から、rtk->ssat[sat-1].vsat[0]=1のazelを抽出して、dops関数に渡す。 */
	for (i=n=0;i<nu;i++) {
		if (rtk->ssat[obs[i].sat-1].vsat[0]==0) continue;
		memcpy(azel+n*2,rtk->ssat[obs[i].sat-1].azel+1,sizeof(double)*2);
		trace(5,"dops: i[%d] azel=[%f]\n",i,azel+n*2);
		n++;
	}
	dops(n,azel,rtk->opt.elmin,dop);
	rtk->gsi.ndop++;
	rtk->gsi.rdop+=(dop[1]-rtk->gsi.rdop)/rtk->gsi.ndop;

	trace(5,"dops: ndops[%d]rdop[%f]dop[%f] n[%d]\n",rtk->gsi.ndop,rtk->gsi.rdop,dop[1],n);

	/* No.2 No.2xは、stat==SOLQ_FIXのとき更新する。 */
	if (stat==SOLQ_FIX) {
		/* No.21 基準点の保存 */
		memcpy(rtk->gsi.rb,rtk->rb,sizeof(double)*3);

		/* No.22 新点の保存 */
		memcpy(rtk->gsi.rr,rtk->sol.rr,sizeof(double)*3);

		/* No.23 RATIO */
		rtk->gsi.ratio=rtk->sol.ratio;

		/* No.24 共分散の保存 */
		memcpy(rtk->gsi.qr,rtk->sol.qr,sizeof(float)*6);
	trace(5,"qr: [%e][%e][%e][%e][%e][%e]\n",rtk->gsi.qr[0],rtk->gsi.qr[1],rtk->gsi.qr[2],rtk->gsi.qr[3],rtk->gsi.qr[4],rtk->gsi.qr[5]);

        rtk->gsi.nfixTotal++;
        for (i=0;i<3;i++) rtk->gsi.rr_fixAve[i] -= (rtk->gsi.rr_fixAve[i] - rtk->gsi.rr[i])/rtk->gsi.nfixTotal;
	}
#if 0
	/* パス情報の設定 */
	for(i=0;i<nu;i++){
		for(j=0;j<NFREQ;j++){
			/* パス情報確認 */
			if(!rtk->ssat[obs[i].sat-1].vsat[j]) continue;

			/* パス情報格納領域確認・確保 */
			if(rtk->passtime[obs[i].sat-1][j].n==0){

				/* 領域を確保(NUM_PASS_ADD分)、初期化 */
				if(addpasstime(&rtk->passtime[obs[i].sat-1][j],obs[i].time)!=0){
					trace(1,"stat_stock : malloc error");
					return;
				}
				continue;
			}

			/* 最後のパスデータの取得 */
			/* 最後のパスの終了時刻と、データ時刻の差を計算 */
			difftime = fabs(timediff( obs[i].time, rtk->passtime[obs[i].sat-1][j].et[rtk->passtime[obs[i].sat-1][j].n-1] ));

			if( difftime < MIN_ARC_GAP ){
				/* パス情報保存 */
				/* 最後のパスの終了時刻と、データ時刻の差がパス分割時間間隔(MIN_ARC_GAP)より小さい場合、
				最後のパスの終了時刻＝データ時刻とする */
				rtk->passtime[obs[i].sat-1][j].et[rtk->passtime[obs[i].sat-1][j].n-1] = obs[i].time;

			}else{
				/* 最後のパスの終了時刻と、データ時刻の差がパス分割時間間隔(MIN_ARC_GAP)以上の場合、
				次の配列に開始終了時刻を設定する。*/
				if(addpasstime(&rtk->passtime[obs[i].sat-1][j],obs[i].time)!=0){
					trace(1,"stat_stock : malloc error");
					return;
				}
				continue;
			}
		}
	}
#endif
}

/* open solution status file ---------------------------------------------------
* open solution status file and set output level
* args   : char     *file   I   rtk status file
*          int      level   I   rtk status level (0: off)
* return : status (1:ok,0:error)
* output : solution status file record format
*
*   $POS,week,tow,stat,posx,posy,posz,posxf,posyf,poszf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          posx/posy/posz    : position x/y/z ecef (m) float
*          posxf/posyf/poszf : position x/y/z ecef (m) fixed
*
*   $VELACC,week,tow,stat,vele,veln,velu,acce,accn,accu,velef,velnf,veluf,accef,accnf,accuf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          vele/veln/velu    : velocity e/n/u (m/s) float
*          acce/accn/accu    : acceleration e/n/u (m/s^2) float
*          velef/velnf/veluf : velocity e/n/u (m/s) fixed
*          accef/accnf/accuf : acceleration e/n/u (m/s^2) fixed
*
*   $CLK,week,tow,stat,clk1,clk2,clk3,clk4
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          clk1     : receiver clock bias GPS (m)
*          clk2     : receiver clock bias GLONASS (m)
*          clk3,clk4: reserved
*
*   $ION,week,tow,stat,sat,az,el,ion,ion-fixed
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          sat      : satellite id
*          az/el    : azimuth/elevation angle(deg)
*          ion      : vertical ionospheric delay L1 (m) float
*          ion-fixed: vertical ionospheric delay L1 (m) fixed
*
*   $TROP,week,tow,stat,rcv,ztd,ztdf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          rcv      : receiver (1:rover,2:base station)
*          ztd      : zenith total delay (m) float
*          ztdf     : zenith total delay (m) fixed
*
*   $HWBIAS,week,tow,stat,frq,bias,biasf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          frq      : frequency (1:L1,2:L2,...)
*          bias     : h/w bias coefficient (m/MHz) float
*          biasf    : h/w bias coefficient (m/MHz) fixed
*
*   $SAT,week,tow,sat,frq,az,el,resp,resc,vsat,snr,fix,slip,lock,outc,slipc,rejc
*          week/tow : gps week no/time of week (s)
*          sat/frq  : satellite id/frequency (1:L1,2:L2,...)
*          az/el    : azimuth/elevation angle (deg)
*          resp     : pseudorange residual (m)
*          resc     : carrier-phase residual (m)
*          vsat     : valid data flag (0:invalid,1:valid)
*          snr      : signal strength (dbHz)
*          fix      : ambiguity flag  (0:no data,1:float,2:fixed,3:hold,4:ppp)
*          slip     : cycle-slip flag (bit1:slip,bit2:parity unknown)
*          lock     : carrier-lock count
*          outc     : data outage count
*          slipc    : cycle-slip count
*          rejc     : data reject (outlier) count
*
*-----------------------------------------------------------------------------*/
extern int rtkopenstat(const char *file, int level)
{
    trace(3,"rtkopenstat: file=%s level=%d\n",file,level);
    
    if (level<=0) return 0;
	if (!(fp_stat=fopen(file,"w"))) {
        trace(1,"rtkopenstat: file open error file=%s\n",file);
        return 0;
    }
    statlevel=level;
    return 1;
}
/* close solution status file --------------------------------------------------
* close solution status file
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkclosestat(void)
{
    trace(3,"rtkclosestat:\n");
    
    if (fp_stat) fclose(fp_stat);
    fp_stat=NULL;
    statlevel=0;
}
/* output solution status ----------------------------------------------------*/
//static void outsolstat(rtk_t *rtk)
static void outsolstat(rtk_t *rtk, const obsd_t *obs)
{
	ssat_t *ssat;
    double tow,pos[3],vel[3],acc[3],vela[3]={0},acca[3]={0},xa[3];
	int i,j,f,k,week,est,nfreq,nf=NF(&rtk->opt),fr;
	char id[32];
    int type,sys;
    
    if (statlevel<=0||!fp_stat) return;
    
    trace(3,"outsolstat:\n");
    
    est=rtk->opt.mode>=PMODE_DGPS;
    nfreq=est?nf:1;
    tow=time2gpst(rtk->sol.time,&week);
    
    /* receiver position */
    if (est) {
        for (i=0;i<3;i++) xa[i]=i<rtk->na?rtk->xa[i]:0.0;
        fprintf(fp_stat,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
                rtk->sol.stat,rtk->x[0],rtk->x[1],rtk->x[2],xa[0],xa[1],xa[2]);
    }
    else {
        fprintf(fp_stat,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
                rtk->sol.stat,rtk->sol.rr[0],rtk->sol.rr[1],rtk->sol.rr[2],
                0.0,0.0,0.0);
    }
    /* receiver velocity and acceleration */
    if (est&&rtk->opt.dynamics) {
        ecef2pos(rtk->sol.rr,pos);
        ecef2enu(pos,rtk->x+3,vel);
		ecef2enu(pos,rtk->x+6,acc);
        if (rtk->na>=6) ecef2enu(pos,rtk->xa+3,vela);
        if (rtk->na>=9) ecef2enu(pos,rtk->xa+6,acca);
        fprintf(fp_stat,"$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
                week,tow,rtk->sol.stat,vel[0],vel[1],vel[2],acc[0],acc[1],acc[2],
                vela[0],vela[1],vela[2],acca[0],acca[1],acca[2]);
    }
	else {
		ecef2pos(rtk->sol.rr,pos);
        ecef2enu(pos,rtk->sol.rr+3,vel);
        fprintf(fp_stat,"$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
                week,tow,rtk->sol.stat,vel[0],vel[1],vel[2],
                0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    }
    /* receiver clocks */
    fprintf(fp_stat,"$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
            week,tow,rtk->sol.stat,1,rtk->sol.dtr[0]*1E9,
			(rtk->sol.dtr[0]+rtk->sol.dtr[1])*1E9,0.0,0.0);
    
    /* ionospheric parameters */
	if (est&&rtk->opt.ionoopt==IONOOPT_EST) {
		for (i=0;i<MAXSAT;i++) {
			ssat=rtk->ssat+i;
			if (!ssat->vs) continue;
			satno2id(i+1,id);
			j=II(i+1,&rtk->opt);
			xa[0]=j<rtk->na?rtk->xa[j]:0.0;
			fprintf(fp_stat,"$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.7f,%.7f\n",
					week,
					tow,
					rtk->sol.stat,
					id,
					ssat->azel[0]*R2D,
					ssat->azel[1]*R2D,
					rtk->x[j],
					xa[0]);
		}
    }
    /* tropospheric parameters */
    if (est&&(rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG)) {
        for (i=0;i<2;i++) {
            j=IT(i,&rtk->opt);
            xa[0]=j<rtk->na?rtk->xa[j]:0.0;
            fprintf(fp_stat,"$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n",week,tow,rtk->sol.stat,
                    i+1,rtk->x[j],xa[0]);
        }
    }
    /* receiver h/w bias */
    if (est&&rtk->opt.glomodear==GLO_ARMODE_AUTO) {
        for (i=0;i<nfreq;i++) {
            j=IL(i,&rtk->opt);
            xa[0]=j<rtk->na?rtk->xa[j]:0.0;
            fprintf(fp_stat,"$HWBIAS,%d,%.3f,%d,%d,%.4f,%.4f\n",week,tow,rtk->sol.stat,
					rtk->opt.oprfrq[i]+1,rtk->x[j],xa[0]);
		}
	}
	/* ISB */
	if (   (est&&rtk->opt.isb==ISBOPT_EST)   || (est&&rtk->opt.isb==ISBOPT_EST_P)
		|| (est&&rtk->opt.isb==ISBOPT_EST_L) || (est&&rtk->opt.isb==ISBOPT_EST_0M)) {

		for (sys=ISYSGPS;sys<=rtk->opt.navsys;sys++) {
			if(!(rtk->opt.navsys & sysno(sys))) continue;
			for (type=0;type<2;type++) {
				if (type==0 && rtk->opt.isb==ISBOPT_EST_P) continue;
				if (type==1 && rtk->opt.isb==ISBOPT_EST_L) continue;
				for (k=0;k<nfreq;k++) {
					f=rtk->opt.oprfrq[k];
					if(f<0) continue;
					if(f==2) fr=5;
					else fr=f+1;

					j=IS(sys,type,k,&rtk->opt);
					xa[0]=j<rtk->na?rtk->xa[j]:0.0;
					fprintf(fp_stat,"$ISB,%d,%.3f,%d,%d,%d,%d,%.4f,%.4f\n",week,tow,rtk->sol.stat,
							sys,fr,type,rtk->x[j],xa[0]);
				}
			}
		}
	}
	/* debug */
	for (sys=ISYSGPS;sys<=rtk->opt.navsys;sys++) {
		if(!(rtk->opt.navsys & sysno(sys))) continue;

		for (k=0;k<nfreq;k++) {
			f=rtk->opt.oprfrq[k];
			if(f<0) continue;
			if(f==2) fr=5;
			else fr=f+1;

			fprintf(fp_stat,"$resave,%d,%.3f,%d,%d,%d,%.15e,%.15e\n",week,tow,rtk->sol.stat,
					sys,fr,rtk->estisb.ave[sys-1][f],rtk->estisb.ave[sys-1][f]-rtk->estisb.ave[ISYSGPS-1][f]);
		}
	}

	/* receiver L2P-L2C */
	if (rtk->opt.gl2bias==GL2OPT_EST) {
		for (i=0;i<2;i++) {
			if(rtk->opt.mode==PMODE_DGPS) continue;
			j=I2(i,&rtk->opt);
			xa[0]=j<rtk->na?rtk->xa[j]:0.0;
			fprintf(fp_stat,"$RCVL2B,%d,%.3f,%d,%d,%.4f,%.4f\n",week,tow,rtk->sol.stat,
					i+1,rtk->x[j],xa[0]);
		}
	}
    if (rtk->sol.stat==SOLQ_NONE||statlevel<=1) return;
    
    /* residuals and status */
	for (i=0,k=0;i<MAXSAT;i++) {
		ssat=rtk->ssat+i;
		if (!ssat->vs) continue;
		while((obs[k].sat!=i+1)) {
			++k;
		}
		satno2id(i+1,id);
		for (j=0;j<nfreq;j++) {
			f=rtk->opt.oprfrq[j];
			fprintf(fp_stat,"$SAT,%d,%.3f,%s,%d,%d,%.1f,%.1f,%.4f,%.4f,%d,%.0f,%d,%d,%d,%d,%d,%d\n",
					week,tow,id,f+1,obs[k].code[f],ssat->azel[0]*R2D,ssat->azel[1]*R2D,
					ssat->resp [f],ssat->resc[f],  ssat->vsat[f],ssat->snr[f]*0.25,
					ssat->fix  [f],ssat->slip[f]&3,ssat->lock[j],ssat->outc[f],
					ssat->slipc[f],ssat->rejc[f]);
		}
		++k;
	}
}

/* open antenna result file ---------------------------------------------------
* open antenna result file
* args   : char     *file   I   antenna result file
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int rtkopenantest(const char *file)
{
    trace(3,"rtkopenantest: file=%s\n",file);

    if (!(fp_antest=fopen(file,"w"))) {
        trace(1,"rtkopenantest: file open error file=%s\n",file);
        return 0;
    }
    return 1;
}
/* close antenna result file --------------------------------------------------
* close antenna result file
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkcloseantest(void)
{
    trace(3,"rtkcloseantest:\n");

    if (fp_antest) fclose(fp_antest);
    fp_antest=NULL;
}
/* open phase-error file ---------------------------------------------------
* open phase-error file
* args   : char     *file   I   phase-error file
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int rtkopenphase(const char *file)
{
    trace(3,"rtkopenphase: file=%s\n",file);

    if (!(fp_phaseerr=fopen(file,"w"))) {
        trace(1,"rtkopenphase: file open error file=%s\n",file);
        return 0;
    }
    return 1;
}
/* close phase-error file --------------------------------------------------
* close phase-error file
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkclosephase(void)
{
    trace(3,"rtkclosephase:\n");

    if (fp_phaseerr) fclose(fp_phaseerr);
    fp_phaseerr=NULL;
}
/* output antenna result ----------------------------------------------------*/
static void outantest(rtk_t *rtk)
{
    int i;
    int use_freq = rtk->opt.antestfreq;
    char time_str[64];
    char satid[8];

    /* ファイルがオープンしていないなら何もしない */
    if (!fp_antest && !fp_phaseerr) return;

    /* 時刻を取得 */
    time2str(rtk->sol.time, time_str, 0);

    if (rtk->opt.antestmode == ANTEST_MODE_GLO) {
        for (i = 0; i < MAXSAT; i++) {
            if (!rtk->ssat[i].vsat[use_freq]) continue;

            /* [format]
             * time,satid,Q,bias
             */
            satno2id(i+1, satid);
            fprintf(fp_antest,"%s,%s,%d,%8.4f\n",
                    time_str, satid, rtk->sol.stat, rtk->x[IB(i,use_freq,&rtk->opt)]);
        }

    } else if (rtk->opt.antestmode == ANTEST_MODE_PCV
			   && revs == 1) {  /* backward only */
        /* [format]
         * time,Q,pcv...
         */
        fprintf(fp_antest,"%s,%d",
                time_str, rtk->sol.stat);
        for (i=sph_start; i<rtk->nx; i++) {
            fprintf(fp_antest,",%.8f",rtk->x[i]);
        }
        fprintf(fp_antest,"\n");

        /* [format]
         * time, satid, satno, base_satno, Q, az, el, phase-error
         */
        if (fp_phaseerr) {
            for (i=0;i<MAXSAT;i++) {
                if (!rtk->ssat[i].vsat[use_freq]) continue;

                satno2id(i+1, satid);
                fprintf(fp_phaseerr,"%s,%s,%d,%d,%d,%12.6f,%12.6f,%12.6f\n",
                        time_str, satid, i+1, base_satno[use_freq], rtk->sol.stat,
                        rtk->ssat[i].azel[0], rtk->ssat[i].azel[1], rtk->ssat[i].resc[use_freq]);
            }
        }
    }
}

/* save error message --------------------------------------------------------*/
static void errmsg(rtk_t *rtk, const char *format, ...)
{
    char buff[256],tstr[32];
    int n;
    va_list ap;
    time2str(rtk->sol.time,tstr,2);
    n=sprintf(buff,"%s: ",tstr+11);
    va_start(ap,format);
    n+=vsprintf(buff+n,format,ap);
    va_end(ap);
    n=n<MAXERRMSG-rtk->neb?n:MAXERRMSG-rtk->neb;
    memcpy(rtk->errbuf+rtk->neb,buff,n);
    rtk->neb+=n;
    trace(2,"%s",buff);
}
/* single-differenced observable ---------------------------------------------*/
static double sdobs(const obsd_t *obs, int i, int j, int f)
{
	double pi=f<NFREQ?obs[i].L[f]:obs[i].P[f-NFREQ];
	double pj=f<NFREQ?obs[j].L[f]:obs[j].P[f-NFREQ];
	return pi==0.0||pj==0.0?0.0:pi-pj;
}
/* single-differenced geometry-free linear combination of phase --------------*/
static double gfobs_L1L2(const obsd_t *obs, int i, int j, const double *lam)
{
	double pi=sdobs(obs,i,j,0)*lam[0],pj=sdobs(obs,i,j,1)*lam[1];
    return pi==0.0||pj==0.0?0.0:pi-pj;
}
static double gfobs_L1L5(const obsd_t *obs, int i, int j, const double *lam)
{
    double pi=sdobs(obs,i,j,0)*lam[0],pj=sdobs(obs,i,j,2)*lam[2];
    return pi==0.0||pj==0.0?0.0:pi-pj;
}
//static double varerr(int sys, double el, double bl, double dt, int f,
static double varerr(int sys, double el, double bl, double dt, int indf, int type,
    const prcopt_t *opt, const obsd_t *obs, const nav_t *nav, int r)
{
    double a, b, c=opt->err[3]*bl/1E4;
    double d=CLIGHT*opt->sclkstab*dt,fact=1.0;
	int i;

	if (type==1) fact=opt->eratio[indf];
	if (fact<=0.0)  fact=opt->eratio[0];

    /*衛星システムから分散荷重係数fの決定*/
	fact *= (sys == SYS_GLO)? EFACT_GLO: ((sys == SYS_SBS)? EFACT_SBS: EFACT_GPS);
    
	/* 観測誤差モデルのデフォルト値を設定 */
	a=opt->err[1];
	b=opt->err[2];

	/*観測誤差モデルテーブルを使用する？*/
    if(opt->errmodel == ERRMODEL_TABLE){
        /* 移動局(=0)/基準局(=1)について検索 */
		if(r==0 || r==1) {
			trace(4,"errtbl_before: r = %d,a = %f,b = %f\n",r,a,b);
			search_errtbl(opt->rectype[r], opt->anttype[r], sys, obs->code[0], 0, nav, &a, &b);
			trace(4,"errtbl_after: r = %d,a = %f,b = %f\n",r,a,b);
		}
	}

	/*疑似距離 かつDCBを補正できない？*/
	if ((type==1)&&(has_dcb(obs->sat, obs->code[0], nav, NULL)==0)) {
		/*設定値Code Error Ratioでfを荷重する */
		fact *= opt->dcbratio;
	}

    /*電離層フリーモード*/
	if(opt->ionoopt==IONOOPT_IFLC){
        /* 分散荷重係数 f=3.0*f */
        fact *= 3.0;
	}

#if TEST_CMPMODEL == 0		//	プロトタイプと一致
	return SQR(fact)*(SQR(a)+SQR(b/sin(el)));
#elif TEST_CMPMODEL == 1	//	RTKLIBと同一
	return SQR(fact)*(SQR(a)+SQR(b/sin(el))+SQR(a)+SQR(b/sin(el))+2.0*SQR(c))+SQR(d);
#else						//	正解式
    return SQR(fact)*(SQR(a)+SQR(b/sin(el))+SQR(c))+SQR(d)/2.0;
#endif
}
/* baseline length -----------------------------------------------------------*/
static double baseline(const double *ru, const double *rb, double *dr)
{
    int i;
    for (i=0;i<3;i++) dr[i]=ru[i]-rb[i];
    return norm(dr,3);
}
/* select common satellites between rover and reference station --------------*/
static int selsat(const obsd_t *obs, double *azel, int nu, int nr,
				  const prcopt_t *opt, int *sat, int *iu, int *ir)
{
	int i,j,k=0;

	trace(3,"selsat  : nu=%d nr=%d\n",nu,nr);

	for (i=0,j=nu;i<nu&&j<nu+nr;i++,j++) {
		if      (obs[i].sat<obs[j].sat) j--;
		else if (obs[i].sat>obs[j].sat) i--;
		else if (azel[1+j*2]>=opt->elmin) { /* elevation at base station */
			sat[k]=obs[i].sat; iu[k]=i; ir[k++]=j;
			trace(4,"(%2d) sat=%3d iu=%2d ir=%2d\n",k-1,obs[i].sat,i,j);
		}
    }
    return k;
}
/* temporal update of position/velocity/acceleration -------------------------*/
static void udpos(rtk_t *rtk, double tt)
{
    double *F,*FP,*xp,pos[3],Q[9]={0},Qv[9],var=0.0;
    int i,j;
    
    trace(3,"udpos   : tt=%.3f\n",tt);
    
    /* fixed mode */
    if (rtk->opt.mode==PMODE_FIXED) {
        for (i=0;i<3;i++) initx(rtk,rtk->opt.ru[i],1E-8,i);
        return;
    }
    /* initialize position for first epoch */
    if (norm(rtk->x,3)<=0.0) {
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        if (rtk->opt.dynamics) {
            for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
            for (i=6;i<9;i++) initx(rtk,1E-6,VAR_ACC,i);
        }
    }
    /* static mode */
    if (rtk->opt.mode==PMODE_STATIC) return;
    
    /* kinmatic mode without dynamics */
    if (!rtk->opt.dynamics) {
		for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        return;
    }
    /* check variance of estimated postion */
    for (i=0;i<3;i++) var+=rtk->P[i+i*rtk->nx]; var/=3.0;
    
    if (var>VAR_POS) {
        /* reset position with large variance */
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
        for (i=6;i<9;i++) initx(rtk,1E-6,VAR_ACC,i);
        trace(2,"reset rtk position due to large variance: var=%.3f\n",var);
        return;
    }
    /* state transition of position/velocity/acceleration */
    F=eye(rtk->nx); FP=mat(rtk->nx,rtk->nx); xp=mat(rtk->nx,1);
    
    for (i=0;i<6;i++) {
        F[i+(i+3)*rtk->nx]=tt;
    }
    /* x=F*x, P=F*P*F+Q */
    matmul("NN",rtk->nx,1,rtk->nx,1.0,F,rtk->x,0.0,xp);
    matcpy(rtk->x,xp,rtk->nx,1);
    matmul("NN",rtk->nx,rtk->nx,rtk->nx,1.0,F,rtk->P,0.0,FP);
    matmul("NT",rtk->nx,rtk->nx,rtk->nx,1.0,FP,F,0.0,rtk->P);
    
    /* process noise added to only acceleration */
    Q[0]=Q[4]=SQR(rtk->opt.prn[3]); Q[8]=SQR(rtk->opt.prn[4]);
    ecef2pos(rtk->x,pos);
    covecef(pos,Q,Qv);
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
        rtk->P[i+6+(j+6)*rtk->nx]+=Qv[i+j*3];
    }
    free(F); free(FP); free(xp);
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void udion(rtk_t *rtk, double tt, double bl, const int *sat, int ns)
{
    double el,fact;
    int i,j;
    
    trace(3,"udion   : tt=%.1f bl=%.0f ns=%d\n",tt,bl,ns);
    
    for (i=1;i<=MAXSAT;i++) {
		j=II(i,&rtk->opt);
        if (rtk->x[j]!=0.0&&
            rtk->ssat[i-1].outc[0]>GAP_RESION&&rtk->ssat[i-1].outc[1]>GAP_RESION)
            rtk->x[j]=0.0;
    }
    for (i=0;i<ns;i++) {
		j=II(sat[i],&rtk->opt);
        
        if (rtk->x[j]==0.0) {
            initx(rtk,1E-6,SQR(rtk->opt.std[1]*bl/1E4),j);
        }
        else {
            /* elevation dependent factor of process noise */
            el=rtk->ssat[sat[i]-1].azel[1];
            fact=cos(el);
            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[1]*bl/1E4*fact)*tt;
        }
    }
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop(rtk_t *rtk, double tt, double bl)
{
    int i,j,k;
    
    trace(3,"udtrop  : tt=%.1f\n",tt);
    
    for (i=0;i<2;i++) {
        j=IT(i,&rtk->opt);
        
        if (rtk->x[j]==0.0) {
            initx(rtk,INIT_ZWD,SQR(rtk->opt.std[2]),j); /* initial zwd */
            
            if (rtk->opt.tropopt>=TROPOPT_ESTG) {
                for (k=0;k<2;k++) initx(rtk,1E-6,VAR_GRA,++j);
            }
        }
        else {
            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[2])*tt;
            
            if (rtk->opt.tropopt>=TROPOPT_ESTG) {
                for (k=0;k<2;k++) {
                    rtk->P[++j*(1+rtk->nx)]+=SQR(rtk->opt.prn[2]*0.3)*fabs(rtk->tt);
                }
            }
        }
    }
}
/* temporal update of isb ------------------------------------*/
static void udisb(rtk_t *rtk, double tt)
{
	int i,j,k,type;
	int nfreq=NF(&rtk->opt);
    
	trace(3,"udisb: tt=%.1f\n",tt);

    for (i=ISYSGPS+1;i<=NSYS;i++) {
        
        for (type=0;type<2;type++) {
			if (type==0 && rtk->opt.isb==ISBOPT_EST_P ) continue;
			if (type==0 && rtk->opt.isb==ISBOPT_EST_0M) continue;
			if (type==1 && rtk->opt.isb==ISBOPT_EST_L ) continue;
			for (k=0;k<nfreq;k++) {
				j=IS(i,type,k,&rtk->opt);

				if (rtk->x[j]==0.0) {
					if (type==0){
						initx(rtk,1E-6,VAR_ISBL,j);
					}
					else if (type==1){
						initx(rtk,1E-6,VAR_ISBP,j);
					}
				}
				else {
					if (type==0){
						rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[5])*rtk->tt;
					}
					else if (type==1){
						rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[6])*rtk->tt;
					}
				}
			}
        }
    }
}
/* temporal update of receiver h/w biases ------------------------------------*/
static void udrcvbias(rtk_t *rtk, double tt)
{
    int i,j;
    
    trace(3,"udrcvbias: tt=%.1f\n",tt);
    
    for (i=0;i<NFREQGLO;i++) {
        j=IL(i,&rtk->opt);
        
        if (rtk->x[j]==0.0) {
            initx(rtk,1E-6,VAR_HWBIAS,j);
        }
        /* hold to fixed solution */
        else if (rtk->nfix>=rtk->opt.minfix&&rtk->sol.ratio>rtk->opt.thresar[0]) {
            initx(rtk,rtk->xa[j],rtk->Pa[j+j*rtk->na],j);
        }
        else {
            rtk->P[j+j*rtk->nx]+=SQR(PRN_HWBIAS)*tt;
        }
    }
}
/* temporal update of receiver GPS L2C L2P biases ------------------------------------*/
static void udgl2bias(rtk_t *rtk, double tt)
{
    int i,j;

	trace(3,"udgl2: tt=%.1f\n",tt);
	for (i=0;i<2;i++) {
		j=I2(i,&rtk->opt);
		if (rtk->x[j]==0.0) {
		//    initx(rtk,1E-6,VAR_HWBIAS,j);
			initx(rtk,1E-6,VAR_GL2,j);
		}
		else {
			rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[7])*tt;
		}
	}
}
/* detect cycle slip by LLI --------------------------------------------------*/
static void detslp_ll(rtk_t *rtk, const obsd_t *obs, int i, int rcv)
{
    unsigned char slip,LLI1,LLI2,LLI;
	int f,sat=obs[i].sat;
    
    trace(3,"detslp_ll: i=%d rcv=%d\n",i,rcv);

//	for (f=0;f<rtk->opt.nfreq;f++) {
	for (f=0;f<NFREQ;f++) {

		if (obs[i].L[f]==0.0) continue;

		/* restore previous LLI */
		LLI1=(rtk->ssat[sat-1].slip[f]>>6)&3;
		LLI2=(rtk->ssat[sat-1].slip[f]>>4)&3;
		LLI=rcv==1?LLI1:LLI2;

		/* detect slip by cycle slip flag */
		slip=(rtk->ssat[sat-1].slip[f]|obs[i].LLI[f])&3;

        if (obs[i].LLI[f]&1) {
            errmsg(rtk,"slip detected (sat=%2d rcv=%d LLI%d=%x)\n",
                   sat,rcv,f+1,obs[i].LLI[f]);
        }
        /* detect slip by parity unknown flag transition */
        if (((LLI&2)&&!(obs[i].LLI[f]&2))||(!(LLI&2)&&(obs[i].LLI[f]&2))) {
            errmsg(rtk,"slip detected (sat=%2d rcv=%d LLI%d=%x->%x)\n",
				   sat,rcv,f+1,LLI,obs[i].LLI[f]);
            slip|=1;
        }
        /* save current LLI and slip flag */
        if (rcv==1) rtk->ssat[sat-1].slip[f]=(obs[i].LLI[f]<<6)|(LLI2<<4)|slip;
        else        rtk->ssat[sat-1].slip[f]=(obs[i].LLI[f]<<4)|(LLI1<<6)|slip;
    }
}
/* detect cycle slip by L1-L2 geometry free phase jump -----------------------*/
static void detslp_gf_L1L2(rtk_t *rtk, const obsd_t *obs, int i, int j,
						   const nav_t *nav)
{
	int sat=obs[i].sat;
	int f;
	int nf;
	double g0,g1;

	trace(3,"detslp_gf_L1L2: i=%d j=%d\n",i,j);

	nf=0;
	for (f=0;f<rtk->opt.nfreq;f++) {
		if     (rtk->opt.oprfrq[f]==0) ++nf;
		else if(rtk->opt.oprfrq[f]==1) ++nf;
	}
	if (nf<2) return;

	if ((g1=gfobs_L1L2(obs,i,j,nav->lam[sat-1]))==0.0) return;

	g0=rtk->ssat[sat-1].gf;
	rtk->ssat[sat-1].gf=g1;
        
	if (g0!=0.0&&fabs(g1-g0)>rtk->opt.thresslip) {
        
		rtk->ssat[sat-1].slip[0]|=1;
		rtk->ssat[sat-1].slip[1]|=1;
        
        errmsg(rtk,"slip detected (sat=%2d GF_L1_L2=%.3f %.3f)\n",sat,g0,g1);
    }
}
/* detect cycle slip by L1-L5 geometry free phase jump -----------------------*/
static void detslp_gf_L1L5(rtk_t *rtk, const obsd_t *obs, int i, int j,
                           const nav_t *nav)
{
    int sat=obs[i].sat;
	int f;
	int nf;
    double g0,g1;
    
    trace(3,"detslp_gf_L1L5: i=%d j=%d\n",i,j);

	nf=0;
	for (f=0;f<rtk->opt.nfreq;f++) {
		if     (rtk->opt.oprfrq[f]==0) ++nf;
		else if(rtk->opt.oprfrq[f]==2) ++nf;
	}
	if (nf<2) return;

	if ((g1=gfobs_L1L5(obs,i,j,nav->lam[sat-1]))==0.0) return;

    g0=rtk->ssat[sat-1].gf2; rtk->ssat[sat-1].gf2=g1;
        
    if (g0!=0.0&&fabs(g1-g0)>rtk->opt.thresslip) {
        
        rtk->ssat[sat-1].slip[0]|=1;
        rtk->ssat[sat-1].slip[2]|=1;
        
        errmsg(rtk,"slip detected (sat=%2d GF_L1_L5=%.3f %.3f)\n",sat,g0,g1);
    }
}
/* detect cycle slip by doppler and phase difference -------------------------*/
static void detslp_dop(rtk_t *rtk, const obsd_t *obs, int i, int rcv,
                       const nav_t *nav)
{
    /* detection with doppler disabled because of clock-jump issue (v.2.3.0) */
#if 0
    int f,sat=obs[i].sat;
    double tt,dph,dpt,lam,thres;
    
    trace(3,"detslp_dop: i=%d rcv=%d\n",i,rcv);
    
    for (f=0;f<rtk->opt.nfreq;f++) {
        if (obs[i].L[f]==0.0||obs[i].D[f]==0.0||rtk->ph[rcv-1][sat-1][f]==0.0) {
            continue;
        }
		if (fabs(tt=timediff(obs[i].time,rtk->pt[rcv-1][sat-1][f]))<DTTOL) continue;
        if ((lam=nav->lam[sat-1][f])<=0.0) continue;
        
        /* cycle slip threshold (cycle) */
        thres=MAXACC*tt*tt/2.0/lam+rtk->opt.err[4]*fabs(tt)*4.0;
        
		/* phase difference and doppler x time (cycle) */
		dph=obs[i].L[f]-rtk->ph[rcv-1][sat-1][f];
        dpt=-obs[i].D[f]*tt;
        
        if (fabs(dph-dpt)<=thres) continue;
        
        rtk->slip[sat-1][f]|=1;
        
		errmsg(rtk,"slip detected (sat=%2d rcv=%d L%d=%.3f %.3f thres=%.3f)\n",
               sat,rcv,f+1,dph,dpt,thres);
    }
#endif
}
/* temporal update of phase biases -------------------------------------------*/
static void udbias(rtk_t *rtk, double tt, const obsd_t *obs, const int *sat,
                   const int *iu, const int *ir, int ns, const nav_t *nav)
{
    double cp,pr,cp1,cp2,pr1,pr2,*bias,offset,lami,lam1,lam2,C1,C2;
    int i,j,f,k,slip,reset,nf=NF(&rtk->opt);
    
    trace(3,"udbias  : tt=%.1f ns=%d\n",tt,ns);
    
    for (i=0;i<ns;i++) {
        
		/* detect cycle slip by LLI */
		for (f=0;f<NFREQ;f++) rtk->ssat[sat[i]-1].slip[f]&=0xFC;
		detslp_ll(rtk,obs,iu[i],1);
		detslp_ll(rtk,obs,ir[i],2);
        
        /* detect cycle slip by geometry-free phase jump */
		detslp_gf_L1L2(rtk,obs,iu[i],ir[i],nav);
        detslp_gf_L1L5(rtk,obs,iu[i],ir[i],nav);
        
        /* detect cycle slip by doppler and phase difference */
        detslp_dop(rtk,obs,iu[i],1,nav);
        detslp_dop(rtk,obs,ir[i],2,nav);
    }
	for (k=0;k<nf;k++) {
		f=rtk->opt.oprfrq[k];
		/* reset phase-bias if instantaneous AR or expire obs outage counter */
        for (i=1;i<=MAXSAT;i++) {
            
            reset=++rtk->ssat[i-1].outc[f]>(unsigned int)rtk->opt.maxout;
            
         // if (rtk->opt.modear==ARMODE_INST&&rtk->x[IB(i,f,&rtk->opt)]!=0.0) {
			if (rtk->opt.modear==ARMODE_INST&&rtk->x[IB(i,k,&rtk->opt)]!=0.0) {
		 //     initx(rtk,0.0,0.0,IB(i,f,&rtk->opt));
				initx(rtk,0.0,0.0,IB(i,k,&rtk->opt));
			}
		    else if (reset&&rtk->x[IB(i,k,&rtk->opt)]!=0.0) {
		 // else if (reset&&rtk->x[IB(i,f,&rtk->opt)]!=0.0) {
             // initx(rtk,0.0,0.0,IB(i,f,&rtk->opt));
                initx(rtk,0.0,0.0,IB(i,k,&rtk->opt));
                trace(3,"udbias : obs outage counter overflow (sat=%3d L%d n=%d)\n",
                      i,f+1,rtk->ssat[i-1].outc[f]);
            }
            if (rtk->opt.modear!=ARMODE_INST&&reset) {
                rtk->ssat[i-1].lock[f]=-rtk->opt.minlock;
            }
        }
        /* reset phase-bias if detecting cycle slip */
        for (i=0;i<ns;i++) {
            j=IB(sat[i],k,&rtk->opt);
            rtk->P[j+j*rtk->nx]+=rtk->opt.prn[0]*rtk->opt.prn[0]*tt;
            slip=rtk->ssat[sat[i]-1].slip[f];
            if (rtk->opt.ionoopt==IONOOPT_IFLC) slip|=rtk->ssat[sat[i]-1].slip[1];
            if (rtk->opt.modear==ARMODE_INST||!(slip&1)) continue;
            rtk->x[j]=0.0;
            rtk->ssat[sat[i]-1].lock[f]=-rtk->opt.minlock;
        }
        bias=zeros(ns,1);
        
        /* estimate approximate phase-bias by phase - code */
        for (i=j=0,offset=0.0;i<ns;i++) {
            
            if (rtk->opt.ionoopt!=IONOOPT_IFLC) {
                cp=sdobs(obs,iu[i],ir[i],f); /* cycle */
                pr=sdobs(obs,iu[i],ir[i],f+NFREQ);
                lami=nav->lam[sat[i]-1][f];
                if (cp==0.0||pr==0.0||lami<=0.0) continue;
                
                bias[i]=cp-pr/lami;
            }
            else {
                cp1=sdobs(obs,iu[i],ir[i],0);
                cp2=sdobs(obs,iu[i],ir[i],1);
                pr1=sdobs(obs,iu[i],ir[i],NFREQ);
                pr2=sdobs(obs,iu[i],ir[i],NFREQ+1);
                lam1=nav->lam[sat[i]-1][0];
                lam2=nav->lam[sat[i]-1][1];
                if (cp1==0.0||cp2==0.0||pr1==0.0||pr2==0.0||lam1<=0.0||lam2<=0.0) continue;
                
                C1= SQR(lam2)/(SQR(lam2)-SQR(lam1));
                C2=-SQR(lam1)/(SQR(lam2)-SQR(lam1));
                bias[i]=(C1*lam1*cp1+C2*lam2*cp2)-(C1*pr1+C2*pr2);
            }
            if (rtk->x[IB(sat[i],k,&rtk->opt)]!=0.0) {
                offset+=bias[i]-rtk->x[IB(sat[i],k,&rtk->opt)];
                j++;
            }
        }
        /* correct phase-bias offset to enssure phase-code coherency */
        if (j>0) {
            for (i=1;i<=MAXSAT;i++) {
				if (rtk->x[IB(i,k,&rtk->opt)]!=0.0) rtk->x[IB(i,k,&rtk->opt)]+=offset/j;
            }
        }
        /* set initial states of phase-bias */
        for (i=0;i<ns;i++) {
			if (bias[i]==0.0||rtk->x[IB(sat[i],k,&rtk->opt)]!=0.0) continue;
            initx(rtk,bias[i],SQR(rtk->opt.std[0]),IB(sat[i],k,&rtk->opt));
        }
        free(bias);
    }
}

/* temporal update of states --------------------------------------------------*/
static void udstate(rtk_t *rtk, const obsd_t *obs, const int *sat,
                    const int *iu, const int *ir, int ns, const nav_t *nav)
{
    double tt=fabs(rtk->tt),bl,dr[3];
    int i;
    
    trace(3,"udstate : ns=%d\n",ns);
    
    /* temporal update of position/velocity/acceleration */
    udpos(rtk,tt);
    
    /* temporal update of ionospheric parameters */
    if (rtk->opt.ionoopt>=IONOOPT_EST) {
        bl=baseline(rtk->x,rtk->rb,dr);
        udion(rtk,tt,bl,sat,ns);
    }
    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt>=TROPOPT_EST) {
        udtrop(rtk,tt,bl);
    }
    /* temporal update of isb */
	if (   (rtk->opt.isb==ISBOPT_EST)   || (rtk->opt.isb==ISBOPT_EST_P)
		|| (rtk->opt.isb==ISBOPT_EST_L) || (rtk->opt.isb==ISBOPT_EST_0M)) {
		udisb(rtk,tt);
    }

    /* temporal update of eceiver h/w bias */
    if (rtk->opt.glomodear==GLO_ARMODE_AUTO&&(rtk->opt.navsys&SYS_GLO)) {
        udrcvbias(rtk,tt);
    }
    /* temporal update of phase-bias */
    if (rtk->opt.mode>PMODE_DGPS) {
        udbias(rtk,tt,obs,sat,iu,ir,ns,nav);
	}
	/* temporal update of GPS L2P-L2C bias */
    if (rtk->opt.gl2bias == GL2OPT_EST) {
        for(i=0;i<NFREQ;++i) {
            if(rtk->opt.oprfrq[i]==ind_L2) {
                udgl2bias(rtk,tt);
                break;
            }
        }
	}
}

double chk_L2Csft(const obsd_t *obs, const prcopt_t *opt, const sta_t* stas){

    double y;
	if(!isL2C(obs->code[1])) {
	    y=0.0;
	}
	else if(opt->phasshft==0){
	    y=0.0;
	}
    else if(opt->phasshft==1){
        y=stas->phaseshift;
    }
    else if(opt->phasshft==2){
	    /* 観測データのシフト情報確認 */
	    if(obs->phasecorr!=0.0){
            y=obs->phasecorr;
		    /*波長補正済みの場合、通常-0.25*/
		}
        else {
	        /*波長補正されていない受信機へ、1/4を設定する。*/
			y=PHASE_CYCLE;
        }
    }
	return y;
}

void setL2Csft(const int phasshft, const char* rectype, const sft_t *sft, sta_t *sta){

    int i;
    sta->phaseshift = 0.0;
	if(phasshft==1){
	    /*波長補正されていない受信機へ、1/4を設定する。*/
        sta->phaseshift = PHASE_CYCLE;
        /*補正済み受信機か判定する。*/
        for(i=0;i<(sft->n);i++){
			if ( strcmp(rectype, sft->rectyp[i]) == 0 ) {
				/*	already shifted	*/
				sta->phaseshift -= sft->bias[i];
				break;
            }
	    }
    }
    return;
}

void setdcb(const nav_t *nav, sta_t *sta){

	int i,j;
	for(i=0;i<nav->nd;i++){
	    if(    (strcmp(sta->rectype,nav->dcb[i].sta_name)==0)
		    || (strcmp(sta->name   ,nav->dcb[i].sta_name)==0)) {
            for(j=0;j<NSYS;j++) {
                sta->dcb[j][0] = nav->dcb[i].cbias_s[j][0];
                sta->dcb[j][1] = nav->dcb[i].cbias_s[j][1];
				sta->dcb[j][2] = nav->dcb[i].cbias_s[j][2];
            }
		    break;
	    }
	}
    return;
}

/* undifferenced phase/code residual for satellite ---------------------------*/
static void zdres_sat(int base, double r, const obsd_t *obs, const nav_t *nav,
                      const double *azel, const double *dant,
                      const prcopt_t *opt, double *y)
{
	const double *lam=nav->lam[obs->sat-1];
	double f1,f2,C1,C2,dant_if,ydif;
    int i,nf=NF(opt),j;
    int sys;
    double ydifIdb[NFREQ][2] = {0};
	double ifb_if[2]={0};
    int f;
    int i0,i1;

	ydif=0.0;
	/*L2C位相波の場合、1/4波長補正チェック*/
	if (isL2C(obs->code[1])){
		ydif = chk_L2Csft(obs, opt, &stas[obs->rcv-1]);
	}

	chk_isb(satsys(obs->sat,NULL), opt, &stas[obs->rcv-1], ydifIdb);

	if (opt->ionoopt==IONOOPT_IFLC) { /* iono-free linear combination */
        i0=opt->oprfrq[0];
        i1=opt->oprfrq[1];
		if (lam[i0]==0.0||lam[i1]==0.0) return;

		if (testsnr(base,i0,azel[1],obs->SNR[i0]*0.25,&opt->snrmask)||
            testsnr(base,i1,azel[1],obs->SNR[i1]*0.25,&opt->snrmask)) return;
        
        f1=CLIGHT/lam[i0];
        f2=CLIGHT/lam[i1];
        C1= SQR(f1)/(SQR(f1)-SQR(f2));
        C2=-SQR(f2)/(SQR(f1)-SQR(f2));
        dant_if=C1*dant[i0]+C2*dant[i1];
		ifb_if[0]=C1*ydifIdb[i0][0]+C2*ydifIdb[i1][0];
        ifb_if[1]=C1*ydifIdb[i0][1]+C2*ydifIdb[i1][1];

		if (obs->L[i0]!=0.0&&obs->L[i1]!=0.0) {
			y[0]=(C1*obs->L[i0]*lam[i0]+C2*(obs->L[i1]-ydif)*lam[i1])-(r+dant_if+ifb_if[i0]);
        }
        if (obs->P[i0]!=0.0&&obs->P[i1]!=0.0) {
            y[1]=(C1*obs->P[i0]        +C2*obs->P[i1]               )-(r+dant_if+ifb_if[i1]);
        }
    }
	else {
		for (i=0;i<nf;i++) {

#ifdef RTKPOS4PCV
            /* check estimate PCV mode */
            if ((opt->antestmode != ANTEST_MODE_NONE) && (opt->antestfreq != i)) {
				/* アンテナ位相特性推定の場合、使用周波数と違うものは使わない */
                continue;
            }
#endif

			f=opt->oprfrq[i];
            if (lam[f]==0.0) continue;
            
            /* check snr mask */
            if (testsnr(base,f,azel[1],obs->SNR[f]*0.25,&opt->snrmask)) {
                continue;
            }
			/* residuals = observable - pseudorange */
			if (isL2C(obs->code[f])){
				if (obs->L[f]!=0.0) y[i     ]=(obs->L[f]+ydif)*lam[f]-(r+dant[f]+ydifIdb[f][0]);
				if (obs->P[f]!=0.0) y[i+nf  ]= obs->P[f]             -(r+dant[f]+ydifIdb[f][1]);
				if (obs->dPpc!=0.0) y[i+nf+1]= obs->dPpc;
			}
			else {
				if (obs->L[f]!=0.0) y[i     ]=(obs->L[f])*lam[f]-(r+dant[f]+ydifIdb[f][0]);
				if (obs->P[f]!=0.0) y[i+nf  ]= obs->P[f]        -(r+dant[f]+ydifIdb[f][1]);
				if (obs->dPpc!=0.0) y[i+nf+1]= obs->dPpc;
			}
		}
	}
}
/* undifferenced phase/code residuals ----------------------------------------*/
static int zdres(int base, const obsd_t *obs, int n, const double *rs,
				 const double *dts, const int *svh, const nav_t *nav,
				 const double *rr, const prcopt_t *opt, int index, double *y,
				 double *e, double *azel)
{
	double r,rr_[3],pos[3],dant[NFREQ]={0},disp[3];
	double zhd,zazel[]={0.0,90.0*D2R};
	int i,nf=NF(opt),j;
	double sphCoef[ANTEST_NMAX * ANTEST_NMAX];

	trace(3,"zdres   : n=%d\n",n);

	for (i=0;i<n*(nf*2+2*TEST_ADDESTPRM);i++) y[i]=0.0;

	if (norm(rr,3)<=0.0) return 0; /* no receiver position */

	for (i=0;i<3;i++) rr_[i]=rr[i];

	/* earth tide correction */
	if (opt->tidecorr) {
		tidedisp(gpst2utc(obs[0].time),rr_,opt->tidecorr,&nav->erp,
				 opt->odisp[base],disp);
		for (i=0;i<3;i++) rr_[i]+=disp[i];
	}
	ecef2pos(rr_,pos);

	for (i=0;i<n;i++) {
		/* compute geometric-range and azimuth/elevation angle */
		if ((r=geodist(rs+i*6,rr_,e+i*3))<=0.0) continue;
		if (satazel(pos,e+i*3,azel+i*2)<opt->elmin) continue;

		/* excluded satellite? */
		if (satexclude(obs[i].sat,svh[i],opt)) continue;

		/* satellite clock-bias */
		r+=-CLIGHT*dts[i*2];

		/* troposphere delay model (hydrostatic) */
		zhd=tropmodel(obs[0].time,pos,zazel,0.0);
		r+=tropmapf(obs[i].time,pos,azel+i*2,NULL)*zhd;

		/* receiver antenna phase center correction */
#ifndef RTKPOS4PCV
		antmodel(opt->pcvr+index,satsys(obs[i].sat,NULL),opt->antdel[index],azel+i*2,opt->posopt[1],
				 dant);
#else
		if ((opt->antestmode == ANTEST_MODE_PCV) && (index == 0)) {
            /* estimate PCV (for rover) */
            double *xp = rr;
            int num = sphharmonic(azel+i*2, opt, sphCoef);
            dant[0] = -dot(sphCoef, &xp[sph_start], num);
            for (j=1;j<nf;j++) dant[j] = dant[0];
        } else {
            /* use antenna model */
            antmodel(opt->pcvr+index,satsys(obs[i].sat,NULL),opt->antdel[index],azel+i*2,opt->posopt[1],
					dant);
        }
#endif

		/* undifferenced phase/code residual for satellite */
		zdres_sat(base,r,obs+i,nav,azel+i*2,dant,opt,y+i*(nf*2+2*TEST_ADDESTPRM));
    }
    trace(4,"rr_=%.3f %.3f %.3f\n",rr_[0],rr_[1],rr_[2]);
    trace(4,"pos=%.9f %.9f %.3f\n",pos[0]*R2D,pos[1]*R2D,pos[2]);
    for (i=0;i<n;i++) {
        trace(4,"sat=%2d %13.3f %13.3f %13.3f %13.10f %6.1f %5.1f\n",
              obs[i].sat,rs[i*6],rs[1+i*6],rs[2+i*6],dts[i*2],azel[i*2]*R2D,
              azel[1+i*2]*R2D);
    }
    trace(4,"y=\n"); tracemat(4,y,nf*2+2*TEST_ADDESTPRM,n,13,3);
    
    return 1;
}

/* test valid observation data -----------------------------------------------*/
static int validobs(int i, int j, int f, int nf, double *y)
{
	int val = 0;
	return y[f+i*(nf*2+2*TEST_ADDESTPRM)]!=0.0&&y[f+j*(nf*2+2*TEST_ADDESTPRM)]!=0.0&&
           (f<nf||(y[f-nf+i*(nf*2+2*TEST_ADDESTPRM)]!=0.0&&y[f-nf+j*(nf*2+2*TEST_ADDESTPRM)]!=0.0));
}
/* double-differenced measurement error covariance ---------------------------*/
static void ddcov(const int *nb,const int *nb2, int n, int n2, const double *Ri, const double *Rj,
                  int nv, double *R)
{
    int i,j,k=0,b;
    
    trace(3,"ddcov   : n=%d\n",n);
    
    for (i=0;i<nv*nv;i++) R[i]=0.0;
    for (b=0;b<n;k+=nb[b++]) {
        for (i=0;i<nb[b];i++) for (j=0;j<nb[b];j++) {
            R[k+i+(k+j)*nv]=Ri[k+i]+(i==j?Rj[k+i]:0.0);
        }
    }
    for (b=0;b<n2;k+=nb2[b++]) {
        for (i=0;i<nb2[b];i++) {
            R[k+i+(k+i)*nv]=Ri[k+i];
        }
    }
    trace(5,"R=\n"); tracemat(5,R,nv,nv,8,6);
}
/* baseline length constraint ------------------------------------------------*/
static int constbl(rtk_t *rtk, const double *x, const double *P, double *v,
                   double *H, double *Ri, double *Rj, int index)
{
    const double thres=0.1; /* threshold for nonliearity (v.2.3.0) */
    double xb[3],b[3],bb,var=0.0;
    int i;
     
    trace(3,"constbl : \n");
    
    /* no constraint */
    if (rtk->opt.baseline[0]<=0.0) return 0;
    
    /* time-adjusted baseline vector and length */
    for (i=0;i<3;i++) {
        xb[i]=rtk->rb[i]+rtk->rb[i+3]*rtk->sol.age;
        b[i]=x[i]-xb[i];
    }
    bb=norm(b,3);
    
    /* approximate variance of solution */
    if (P) {
        for (i=0;i<3;i++) var+=P[i+i*rtk->nx];
        var/=3.0;
    }
    /* check nonlinearity */
    if (var>thres*thres*bb*bb) {
        trace(3,"constbl : equation nonlinear (bb=%.3f var=%.3f)\n",bb,var);
        return 0;
    }
    /* constraint to baseline length */
    v[index]=rtk->opt.baseline[0]-bb;
    if (H) {
        for (i=0;i<3;i++) H[i+index*rtk->nx]=b[i]/bb;
    }
    Ri[index]=0.0;
    Rj[index]=SQR(rtk->opt.baseline[1]);
    
    trace(4,"baseline len   v=%13.3f R=%8.6f %8.6f\n",v[index],Ri[index],Rj[index]);
    
    return 1;
}
/* precise tropspheric model -------------------------------------------------*/
static double prectrop(gtime_t time, const double *pos, int r,
                       const double *azel, const prcopt_t *opt, const double *x,
                       double *dtdx)
{
    double m_w=0.0,cotz,grad_n,grad_e;
    int i=IT(r,opt);
    
    /* wet mapping function */
    tropmapf(time,pos,azel,&m_w);
    
    if (opt->tropopt>=TROPOPT_ESTG&&azel[1]>0.0) {
        
		/* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz=1.0/tan(azel[1]);
        grad_n=m_w*cotz*cos(azel[0]);
        grad_e=m_w*cotz*sin(azel[0]);
        m_w+=grad_n*x[i+1]+grad_e*x[i+2];
        dtdx[1]=grad_n*x[i];
        dtdx[2]=grad_e*x[i];
    }
    else dtdx[1]=dtdx[2]=0.0;
    dtdx[0]=m_w;
    return m_w*x[i];
}
/* glonass inter-channel bias correction -------------------------------------*/
static double gloicbcorr(int sat1, int sat2, const prcopt_t *opt, double lam1,
                         double lam2, int f)
{
    double dfreq;
    
    if (f>=NFREQGLO||f>=opt->nfreq||!opt->exterr.ena[2]) return 0.0;
    
    dfreq=(CLIGHT/lam1-CLIGHT/lam2)/(f==0?DFRQ1_GLO:DFRQ2_GLO);
    
    return opt->exterr.gloicb[f]*0.01*dfreq; /* (m) */
}
/* double-differenced phase/code residuals -----------------------------------*/
static int ddres(rtk_t *rtk, const nav_t *nav, double dt, const double *x,
                 const double *P, const int *sat, double *y, double *e,
                 double *azel, const int *iu, const int *ir, int ns, double *v,
				 double *H, double *R, int *vflg, const obsd_t *obs, int *ux)
{
    prcopt_t *opt=&rtk->opt;
    double bl,dr[3],posu[3],posr[3],didxi=0.0,didxj=0.0,*im;
	double *tropr,*tropu,*dtdxr,*dtdxu,*Ri,*Rj,s,lami,lamj,fi,fj,df,*Hi=NULL;
  	double freq1,freq2;
	int i,j,k,m,f,ff,nv=0,nb[NFREQ*2*2+2]={0},nb2[2*TEST_ADDESTPRM]={0},b=0,b2=0,sysi,sysj,nf=NF(opt);
    int mmax;
    int isysi,isysj;
    int type;
    int indf;
    int tmp;
    double sphCoefi[ANTEST_NMAX * ANTEST_NMAX];
	double sphCoefj[ANTEST_NMAX * ANTEST_NMAX];

	int ii[2]={-1,-1};
	int il=-1;
	int is[2][2]={{-1,-1},{-1,-1}};
	int i2[2]={-1,-1};
	int ib[2]={-1,-1};

	int nr = NR(opt);

	double gamma,ci1,ci2,cj1,cj2;

	trace(3,"ddres   : dt=%.1f nx=%d ns=%d\n",dt,rtk->nx,ns);

	bl=baseline(x,rtk->rb,dr);
    ecef2pos(x,posu); ecef2pos(rtk->rb,posr);
    
    Ri=mat(ns*(nf*2+2*TEST_ADDESTPRM)+2,1); Rj=mat(ns*(nf*2+2*TEST_ADDESTPRM)+2,1); im=mat(ns,1);
	tropu=mat(ns,1); tropr=mat(ns,1); dtdxu=mat(ns,3); dtdxr=mat(ns,3);
    
	for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
		rtk->ssat[i].resp[j]=rtk->ssat[i].resc[j]=0.0;
    }
#ifdef RTKPOS4PCV
    for (j=0;j<NFREQ;j++) base_satno[j]=0;
#endif

    /* compute factors of ionospheric and tropospheric delay */
	for (i=0;i<ns;i++) {
        if (opt->ionoopt>=IONOOPT_EST) {
            im[i]=(ionmapf(posu,azel+iu[i]*2)+ionmapf(posr,azel+ir[i]*2))/2.0;
		}
		if (opt->tropopt>=TROPOPT_EST) {
			tropu[i]=prectrop(rtk->sol.time,posu,0,azel+iu[i]*2,opt,x,dtdxu+i*3);
			tropr[i]=prectrop(rtk->sol.time,posr,1,azel+ir[i]*2,opt,x,dtdxr+i*3);
		}
	}

	if(opt->diff!=DIFOPT_EXCGLO) mmax = 1;
	else            mmax = 2;

	for (i=0;i<rtk->nx;i++) ux[i]=0;
	for (i=0;i<3;i++) ux[i]=1;
	for (i=0;i<NT(opt);i++) ux[IT(0,opt)+i]=1;


	for (m=0;m<mmax;m++) {/* for each system (0:gps/qzss/sbas,1:glonass) */

		for (f=opt->mode>PMODE_DGPS?0:nf;f<nf*2;f++) {
			if (f<nf) type=0;
			else      type=1;
			ff=f%nf;
			indf=rtk->opt.oprfrq[ff];
			/* search reference satellite with highest elevation */
			for (i=-1,j=0;j<ns;j++) {
				sysi=rtk->ssat[sat[j]-1].sys;
				if ((m==0)&&(sysi==SYS_GLO)&&(opt->diff==DIFOPT_EXCGLO)) continue;
				if ((m==1)&&(sysi!=SYS_GLO)) continue;
				//if ((opt->diff==DIFOPT_EXCGLO) && ((m==0&&sysi==SYS_GLO)||(m==1&&sysi!=SYS_GLO))) {
				//	continue;
				//}
				if (!validobs(iu[j],ir[j],f,nf,y)) {
					continue;
				}
				if (i<0||azel[1+iu[j]*2]>=azel[1+iu[i]*2]) i=j;
			}
			if (i<0) continue;

			sysi=rtk->ssat[sat[i]-1].sys;
			if (sysi!=SYS_GAL) gamma=SQR(nav->lam[sat[i]-1][1])/SQR(nav->lam[sat[i]-1][0]); /* f1^2/f2^2 */
			else               gamma=SQR(nav->lam[sat[i]-1][2])/SQR(nav->lam[sat[i]-1][0]); /* f1^2/f2^2 */
			ci1=gamma/(gamma-1.0);  /*  f1^2/(f1^2-f2^2) */
			ci2=-1.0 /(gamma-1.0);  /* -f2^2/(f1^2-f2^2) */

			/* make double difference */
			for (j=0;j<ns;j++) {
				if (i==j) continue;
				sysj=rtk->ssat[sat[j]-1].sys;
				if (sysj!=SYS_GAL) gamma=SQR(nav->lam[sat[j]-1][1])/SQR(nav->lam[sat[j]-1][0]); /* f1^2/f2^2 */
				else               gamma=SQR(nav->lam[sat[j]-1][2])/SQR(nav->lam[sat[j]-1][0]); /* f1^2/f2^2 */
				cj1=gamma/(gamma-1.0);  /*  f1^2/(f1^2-f2^2) */
				cj2=-1.0 /(gamma-1.0);  /* -f2^2/(f1^2-f2^2) */
				if ((m==0)&&(sysj==SYS_GLO)&&(opt->diff==DIFOPT_EXCGLO)) continue;
				if ((m==1)&&(sysj!=SYS_GLO)) continue;

				if (!validobs(iu[j],ir[j],f,nf,y)) {
					continue;
				}
				//lami=nav->lam[sat[i]-1][ff];
				//lamj=nav->lam[sat[j]-1][ff];
				lami=nav->lam[sat[i]-1][rtk->opt.oprfrq[ff]];
				lamj=nav->lam[sat[j]-1][rtk->opt.oprfrq[ff]];
				if (lami<=0.0||lamj<=0.0) {
					continue;
				}
			   //	indf=rtk->opt.oprfrq[ff];
				if (H) Hi=H+nv*rtk->nx;

				/* double-differenced residual */
				v[nv]=(y[f+iu[i]*(nf*2+2*TEST_ADDESTPRM)]-y[f+ir[i]*(nf*2+2*TEST_ADDESTPRM)])-
					  (y[f+iu[j]*(nf*2+2*TEST_ADDESTPRM)]-y[f+ir[j]*(nf*2+2*TEST_ADDESTPRM)]);

				/* ISB */
//				chk_isb(2, satsys(obs->sat,NULL), rtk->opt, stas, ydifIdb);

				/* partial derivatives by rover position */
				if (H) {
					for (k=0;k<3;k++) {
						Hi[k]=-e[k+iu[i]*3]+e[k+iu[j]*3];
					}
				}
				/* double-differenced ionospheric delay term */
				if (opt->ionoopt==IONOOPT_EST) {
					fi=lami/lam_carr[0]; fj=lamj/lam_carr[0];
					didxi=(f<nf?-1.0:1.0)*fi*fi*im[i];
					didxj=(f<nf?-1.0:1.0)*fj*fj*im[j];
					ii[0]=II(sat[i],opt);
					ii[1]=II(sat[j],opt);
					v[nv]-=didxi*x[ii[0]]-didxj*x[ii[1]];
					if (H) {
						Hi[ii[0]]= didxi;
						Hi[ii[1]]=-didxj;
					}
				}

#ifdef RTKPOS4PCV
				/* estimate PCV */
				if (opt->antestmode == ANTEST_MODE_PCV) {
					if (H) {
						sphharmonic(azel+iu[i]*2, opt, sphCoefi);
						sphharmonic(azel+iu[j]*2, opt, sphCoefj);
						for (k=sph_start; k<rtk->nx; k++) {
							Hi[k] = sphCoefi[k - sph_start] - sphCoefj[k - sph_start];
						}
					}
				}
#endif
				/* double-differenced tropospheric delay term */
				if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
					v[nv]-=(tropu[i]-tropu[j])-(tropr[i]-tropr[j]);
					for (k=0;k<(opt->tropopt<TROPOPT_ESTG?1:3);k++) {
						if (!H) continue;
						Hi[IT(0,opt)+k]= (dtdxu[k+i*3]-dtdxu[k+j*3]);
						Hi[IT(1,opt)+k]=-(dtdxr[k+i*3]-dtdxr[k+j*3]);
					}
				}

				/* double-differenced phase-bias term */
				if (f<nf) {
					/*電離層フリーモード*/
					if(opt->ionoopt==IONOOPT_IFLC){
						/*L1とL2の電離層フリー線形結合の搬送波位相アンビキュイティの二重差をvkとHkに加算*/
						ib[0]=IB(sat[i],f,opt);
						ib[1]=IB(sat[j],f,opt);
						v[nv]-=x[ib[0]]-x[ib[1]];
						if (H) {
							Hi[ib[0]]= 1.0;
							Hi[ib[1]]=-1.0;
						}
					}else{
						/*L1またはL2の搬送波位相アンビキュイティの二重差をvkとHkに加算*/
						ib[0]=IB(sat[i],f,opt);
						ib[1]=IB(sat[j],f,opt);
						v[nv]-=lami*x[ib[0]]-lamj*x[ib[1]];
						if (H) {
							Hi[ib[0]]= lami;
							Hi[ib[1]]=-lamj;
						}
					}
				}
				/*glonass && IFBテーブルを使用する？(3:use IFB table)*/
				if (rtk->opt.glomodear==GLO_ARMODE_IFB&&sysi==SYS_GLO&&sysj==SYS_GLO) {
					/* IFBテーブルからの補正計算 */
					freq1= CLIGHT/lami;
					freq2= CLIGHT/lamj;
					trace(5,"n[%d] v[n] before [%8.6f]lami[%f]lamj[%f]freq1[%f]freq2[%f]\n",nv,v[nv],lami,lamj,freq1,freq2);
					v[nv]-=gloifbcorr(nav,opt,iu[j],ir[j],f,freq1,freq2);
					trace(5,"n[%d] v[n] after  [%8.6f]\n",nv,v[nv]);
				}
				/* glonass receiver h/w bias term */
				else if (rtk->opt.glomodear==GLO_ARMODE_AUTO&&sysi==SYS_GLO&&sysj==SYS_GLO&&ff<NFREQGLO) {
					df=(CLIGHT/lami-CLIGHT/lamj)/1E6; /* freq-difference (MHz) */
					trace(5,"H/W bias n[%d] v[n] before [%8.6f]lami[%f]lamj[%f]\n",nv,v[nv],lami,lamj);
					il=IL(ff,opt);
					v[nv]-=df*x[il];
					trace(5,"H/W bias n[%d] v[n] after  [%8.6f]\n",nv,v[nv]);
					if (H) Hi[il]=df;

					if (H) trace(5,"H/W bias IL[%d] Hi[n]  [%8.6f]\n",il,Hi[il]);
				}
				/* glonass interchannel bias correction */
				else if (sysi==SYS_GLO&&sysj==SYS_GLO) {

					v[nv]-=gloicbcorr(sat[i],sat[j],&rtk->opt,lami,lamj,f);
				}
				/* double-differenced isb */
				if (   (type==0 && ((opt->isb==ISBOPT_EST) || (opt->isb==ISBOPT_EST_L)))
					|| (type==1 && ((opt->isb==ISBOPT_EST) || (opt->isb==ISBOPT_EST_P) || (opt->isb==ISBOPT_EST_0M)))) {
					isysi = sysind(sysi);
					isysj = sysind(sysj);
					if (opt->ionoopt!=IONOOPT_IFLC) {
						if(isysi!=isysj) {
							if(ISYSGPS<isysi) {
								is[0][0]=IS(isysi,type,ff,opt);
								v[nv]-=x[is[0][0]];
								if(H) Hi[is[0][0]]=1;
							}
							if(ISYSGPS<isysj) {
								is[1][0]=IS(isysj,type,ff,opt);
								v[nv]-=-x[is[1][0]];
								if(H) Hi[is[1][0]]=-1;
							}
						}
					}
					else {
						if(ISYSGPS<isysi) {
							is[0][0]=IS(isysi,type,0,opt);
							is[0][1]=IS(isysi,type,1,opt);
							v[nv]-=ci1*x[is[0][0]]+ci2*x[is[0][1]];
							if(H) Hi[is[0][0]]=ci1;
							if(H) Hi[is[0][1]]=ci2;
						}
						if(ISYSGPS<isysj) {
							is[1][0]=IS(isysj,type,0,opt);
							is[1][1]=IS(isysj,type,1,opt);
							v[nv]-=-(cj1*x[is[1][0]]+cj2*x[is[1][1]]);
							if(H) Hi[is[1][0]]=-cj1;
							if(H) Hi[is[1][1]]=-cj2;
						}
					}
				}

				if (opt->gl2bias==GL2OPT_EST) {
					if (f>=nf && indf==ind_L2) {
						if(isL2P(obs[iu[j]].code[1]) && isL2C(obs[iu[i]].code[1])) {
							i2[0] = I2(0,opt);
							v[nv] -= -x[i2[0]];
							if(H) Hi[i2[0]]= -1;
						}
						else if(isL2C(obs[iu[j]].code[1]) && isL2P(obs[iu[i]].code[1])) {
							i2[0] = I2(0,opt);
							v[nv] -= x[i2[0]];
							if(H) Hi[i2[0]]=1;
						}
						if(isL2P(obs[ir[j]].code[1]) && isL2C(obs[ir[i]].code[1])) {
							i2[1] = I2(1,opt);
							v[nv]-= x[i2[1]];
							if(H) Hi[I2(1, opt)]=1;
						}
						else if(isL2C(obs[ir[j]].code[1]) && isL2P(obs[ir[i]].code[1])) {
							i2[1] = I2(1,opt);
							v[nv]-=-x[i2[1]];
							if(H) Hi[i2[1]]= -1;
						}
					}
				}

				if (f<nf) rtk->ssat[sat[j]-1].resc[indf]=v[nv];
				else      rtk->ssat[sat[j]-1].resp[indf]=v[nv];
				/* test innovation */
				if (opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) {
					if (f<nf) {
						rtk->ssat[sat[i]-1].rejc[indf]++;
						rtk->ssat[sat[j]-1].rejc[indf]++;
					}
					errmsg(rtk,"outlier rejected (sat=%3d-%3d %s%d v=%.3f)\n",
						   sat[i],sat[j],f<nf?"L":"P",f%nf+1,v[nv]);
					continue;
				}

				/*一重差の観測誤差分散行列を計算*/
				/* single-differenced measurement error variances */
				Ri[nv]=varerr(sysi,azel[1+iu[i]*2],bl,dt,indf,type,opt,obs,nav,0) + varerr(sysi,azel[1+ir[i]*2],bl,dt,indf,type,opt,obs,nav,1);
				Rj[nv]=varerr(sysj,azel[1+iu[j]*2],bl,dt,indf,type,opt,obs,nav,0) + varerr(sysj,azel[1+ir[j]*2],bl,dt,indf,type,opt,obs,nav,1);

				/* set valid data flags */
				if (opt->mode>PMODE_DGPS) {
					if (f<nf) rtk->ssat[sat[i]-1].vsat[indf]=rtk->ssat[sat[j]-1].vsat[indf]=1;
				}
				else {
					rtk->ssat[sat[i]-1].vsat[indf]=rtk->ssat[sat[j]-1].vsat[indf]=1;
				}
                
				trace(4,"sat=%3d-%3d %s%d v=%13.3f R=%8.6f %8.6f\n",sat[i],
					  sat[j],f<nf?"L":"P",f%nf+1,v[nv],Ri[nv],Rj[nv]);
            
				vflg[nv++]=(sat[i]<<16)|(sat[j]<<8)|((f<nf?0:1)<<4)|(f%nf);
				nb[b]++;

				if(-1!=ii[0]) ux[ii[0]]=1;
				if(-1!=ii[1]) ux[ii[1]]=1;
				if(-1!=il   ) ux[il]   =1;
				if(-1!=is[0][0]) ux[is[0][0]]=1;
				if(-1!=is[1][0]) ux[is[1][0]]=1;
				if(-1!=is[0][1]) ux[is[0][1]]=1;
				if(-1!=is[1][1]) ux[is[1][1]]=1;
				if(-1!=i2[0]) ux[i2[0]]=1;
				if(-1!=i2[1]) ux[i2[1]]=1;
				if(-1!=ib[0]) ux[ib[0]]=1;
				if(-1!=ib[1]) ux[ib[1]]=1;
			}
			/* restore single-differenced residuals assuming sum equal zero */
			if (f<nf) {
				for (j=0,s=0.0;j<ns;j++) s+=rtk->ssat[j].resc[indf];
				s/=nb[b]+1;
				for (j=0;j<MAXSAT;j++) {
					if (j==sat[i]-1||rtk->ssat[j].resc[indf]!=0.0) rtk->ssat[j].resc[indf]-=s;
				}
#ifdef RTKPOS4PCV
                base_satno[f] = sat[i];
#endif
            }
			else {
				for (j=0,s=0.0;j<MAXSAT;j++) s+=rtk->ssat[j].resp[indf];
				s/=nb[b]+1;
				for (j=0;j<MAXSAT;j++) {
					if (j==sat[i]-1||rtk->ssat[j].resp[indf]!=0.0)
						rtk->ssat[j].resp[indf]-=s;
				}
			}
			b++;
		}
	}
    /* end of system loop */

    
    if (opt->gl2bias==GL2OPT_EST) {
        f=1+nf;
        for (i=0;i<2;i++) {
            for (j=0;j<ns;j++) {
				if (H) Hi=H+nv*rtk->nx;
				if     (i==0) k=iu[j];
				else if(i==1) k=ir[j];
				v[nv] = obs[k].dPpc;

				if(obs[k].dPpc!=0.0) {
					i2[0]=I2(i,opt);
					v[nv]-=x[i2[0]];
					if (H) Hi[i2[0]]=1;
					Ri[nv]=varerr(rtk->ssat[sat[j]-1].sys,azel[1+k*2],bl,dt,ind_L2,1,opt,obs,nav,i);
					nv++;
					nb2[b2]++;
					if(-1!=i2[0]) ux[i2[0]]=1;
                }
            }
            for (j=0,s=0.0;j<MAXSAT;j++) s+=rtk->ssat[j].resdpl2;
			s/=nb2[b2]+1;
			for (j=0;j<MAXSAT;j++) {
                if (rtk->ssat[j].resdpl2!=0.0) rtk->ssat[j].resdpl2-=s;     /*0?*/
            }
            b2++;
        }
    }


    
    /* baseline length constraint for moving baseline */
    if (opt->mode==PMODE_MOVEB&&constbl(rtk,x,P,v,H,Ri,Rj,nv)) {
        vflg[nv++]=3<<4;
		nb[b++]++;
	}
	if (H) {trace(5,"H=\n"); tracemat(5,H,rtk->nx,nv,7,4);}
	if (H) {trace(5,"ux=\n"); traceimat(5,ux,rtk->nx,1);}

	/* double-differenced measurement error covariance */
	ddcov(nb,nb2,b,b2,Ri,Rj,nv,R);
    
    free(Ri);
    free(Rj);
    free(im);
    free(tropu);
    free(tropr); free(dtdxu); free(dtdxr);
    
	return nv;
}

/* time-interpolation of residuals (for post-mission) ------------------------*/
static double intpres(gtime_t time, const obsd_t *obs, int n, const nav_t *nav,
                      rtk_t *rtk, double *y)
{
    static obsd_t obsb[MAXOBS];
    static double yb[MAXOBS*NFREQ*2],rs[MAXOBS*6],dts[MAXOBS*2],var[MAXOBS];
    static double e[MAXOBS*3],azel[MAXOBS*2];
    static int nb=0,svh[MAXOBS*2];
    prcopt_t *opt=&rtk->opt;
    double tt=timediff(time,obs[0].time),ttb,*p,*q;
    int i,j,k,nf=NF(opt);
    
    trace(3,"intpres : n=%d tt=%.1f\n",n,tt);
    
    if (nb==0||fabs(tt)<DTTOL) {
        nb=n; for (i=0;i<n;i++) obsb[i]=obs[i];
        return tt;
    }
    ttb=timediff(time,obsb[0].time);
    if (fabs(ttb)>opt->maxtdiff*2.0||ttb==tt) return tt;
    
    satposs(time,obsb,nb,nav,opt->sateph,rs,dts,var,svh);
    
    if (!zdres(1,obsb,nb,rs,dts,svh,nav,rtk->rb,opt,1,yb,e,azel)) {
        return tt;
    }
    for (i=0;i<n;i++) {
        for (j=0;j<nb;j++) if (obsb[j].sat==obs[i].sat) break;
        if (j>=nb) continue;
        for (k=0,p=y+i*(nf*2+2*TEST_ADDESTPRM),q=yb+j*(nf*2+2*TEST_ADDESTPRM);k<(nf*2+2*TEST_ADDESTPRM);k++,p++,q++) {
            if (*p==0.0||*q==0.0) *p=0.0; else *p=(ttb*(*p)-tt*(*q))/(ttb-tt);
        }
    }
    return fabs(ttb)>fabs(tt)?ttb:tt;
}

/* single to double-difference transformation matrix (D') --------------------*/
static int ddmat(rtk_t *rtk, double *D)
{
	int i,j,k,m,n,f,nb=0,nx=rtk->nx,na=rtk->na;
	int nf=NF(&rtk->opt),slip;
    int mmax;

	trace(3,"ddmat   :\n");
    
    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
        rtk->ssat[i].fix[j]=0;
    }
    for (i=0;i<na;i++) D[i+i*nx]=1.0;
    
    if(rtk->opt.diff!=DIFOPT_EXCGLO) mmax = 1;
    else                             mmax = 2;
    for (m=0;m<mmax;m++) { /* 0:gps,1:glonass */
        
        if (m==1&&rtk->opt.glomodear==GLO_ARMODE_OFF) continue;

		for (n=0,k=na;n<nf;n++,k+=MAXSAT) {
			f=rtk->opt.oprfrq[n];
			for (i=k;i<k+MAXSAT;i++) {
				if(rtk->x[i]==0.0) continue;
				if ((m==0)&&(rtk->ssat[i-k].sys==SYS_GLO)&&(rtk->opt.diff==DIFOPT_EXCGLO)) continue;
				if ((m==1)&&(rtk->ssat[i-k].sys!=SYS_GLO)) continue;
				slip=rtk->ssat[i-k].slip[f];
				/*電離層フリーの場合、2波のbit or を取る*/
				if (rtk->opt.ionoopt==IONOOPT_IFLC) {
					slip|=rtk->ssat[i-k].slip[rtk->opt.oprfrq[1]];
				}

				if (rtk->ssat[i-k].lock[f]>0&&!(slip&2)&&
					rtk->ssat[i-k].azel[1]>=rtk->opt.elmaskar) {
					rtk->ssat[i-k].fix[f]=2; /* fix */
					break;
				}
				else rtk->ssat[i-k].fix[f]=1;
			}
			for (j=k;j<k+MAXSAT;j++) {
				if(i==j) continue;
				if(rtk->x[j]==0.0) continue;
				if ((m==0)&&(rtk->ssat[j-k].sys==SYS_GLO)&&(rtk->opt.diff==DIFOPT_EXCGLO)) continue;
				if ((m==1)&&(rtk->ssat[j-k].sys!=SYS_GLO)) continue;
				slip=rtk->ssat[i-k].slip[f];
				/*電離層フリーの場合、2波のbit or を取る*/
				if (rtk->opt.ionoopt==IONOOPT_IFLC) {
					slip|=rtk->ssat[i-k].slip[rtk->opt.oprfrq[1]];
				}

				if (rtk->ssat[j-k].lock[f]>0&&!(slip&2)&&
					rtk->ssat[j-k].azel[1]>=rtk->opt.elmaskar) {
					D[i+(na+nb)*nx]= 1.0;
					D[j+(na+nb)*nx]=-1.0;
					nb++;
					rtk->ssat[j-k].fix[f]=2; /* fix */
				}
				else rtk->ssat[j-k].fix[f]=1;
			}
		}
    }
    trace(5,"D=\n"); tracemat(5,D,nx,na+nb,2,0);
    return nb;
}

/* restore single-differenced ambiguity --------------------------------------*/
static void restamb(rtk_t *rtk, const double *bias, int nb, double *xa)
{
    int i,n,m,f,index[MAXSAT],nv=0;
	int nf=NF(&rtk->opt);
    int mmax;

    trace(3,"restamb :\n");
    
    for (i=0;i<rtk->nx;i++) xa[i]=rtk->x [i];
    for (i=0;i<rtk->na;i++) xa[i]=rtk->xa[i];
    
    if(rtk->opt.diff!=DIFOPT_EXCGLO) mmax = 1;
    else                             mmax = 2;

	for (m=0;m<mmax;m++) for (f=0;f<nf;f++) {
        
        for (n=i=0;i<MAXSAT;i++) {
		   if(rtk->ssat[i].fix[rtk->opt.oprfrq[f]]!=2) continue;
            if((m==0)&&(rtk->ssat[i].sys==SYS_GLO)&&(rtk->opt.diff==DIFOPT_EXCGLO)) continue;
            if((m==1)&&(rtk->ssat[i].sys!=SYS_GLO)) continue;
			index[n++]=IB(i+1,f,&rtk->opt);
        }
        if (n<2) continue;

        xa[index[0]]=rtk->x[index[0]];
        
        for (i=1;i<n;i++) {
            xa[index[i]]=xa[index[0]]-bias[nv++];
        }
    }
}

/* hold integer ambiguity ----------------------------------------------------*/
static void holdamb(rtk_t *rtk, const double *xa)
{
    double *v,*H,*R;
	int i,j,n,m,f,info,index[MAXSAT],nb=rtk->nx-rtk->na,nv=0;
	int nf=NF(&rtk->opt);
	int mmax;

	trace(3,"holdamb :\n");
    
	v=zeros(nb,1); H=zeros(nb,rtk->nx);

	if(rtk->opt.diff!=DIFOPT_EXCGLO) mmax = 1;
	else                             mmax = 2;
	for (m=0;m<mmax;m++) {
		for (j=0;j<nf;j++) {
			f=rtk->opt.oprfrq[j];

			for (n=i=0;i<MAXSAT;i++) {
				if(rtk->ssat[i].azel[1]<rtk->opt.elmaskhold) continue;
				if(rtk->ssat[i].fix[f]!=2) continue;
				if((m==0) && (rtk->ssat[i].sys==SYS_GLO) && (rtk->opt.diff==DIFOPT_EXCGLO)) continue;
				if((m==1) && (rtk->ssat[i].sys!=SYS_GLO)) continue;

				index[n++]=IB(i+1,j,&rtk->opt);
				rtk->ssat[i].fix[f]=3; /* hold */
			}
			/* constraint to fixed ambiguity */
			for (i=1;i<n;i++) {
				v[nv]=(xa[index[0]]-xa[index[i]])-(rtk->x[index[0]]-rtk->x[index[i]]);

				H[index[0]+nv*rtk->nx]= 1.0;
				H[index[i]+nv*rtk->nx]=-1.0;
				nv++;
			}
		}
	}
	if (nv>0) {
		R=zeros(nv,nv);
		for (i=0;i<nv;i++) R[i+i*nv]=VAR_HOLDAMB;

		/* update states with constraints */
		if ((info=filter(rtk->x,rtk->P,H,v,R,rtk->nx,nv, NULL))) {
			errmsg(rtk,"filter error (info=%d)\n",info);
        }
        free(R);
    }
    free(v); free(H);
}

/* resolve integer ambiguity by LAMBDA ---------------------------------------*/
static int resamb_LAMBDA(rtk_t *rtk, double *bias, double *xa)
{
    prcopt_t *opt=&rtk->opt;
	int i,j,ny,nb,info,nx=rtk->nx,na=rtk->na;
    double *D,*DP,*y,*Qy,*b,*db,*Qb,*Qab,*QQ,s[2];
    
    trace(3,"resamb_LAMBDA : nx=%d\n",nx);
    
    rtk->sol.ratio=0.0;
    
	if (rtk->opt.mode<=PMODE_DGPS||rtk->opt.modertkar==RTKAR_OFF||
		rtk->opt.thresar[0]<1.0) {
        return 0;
    }
    /* single to double-difference transformation matrix (D') */
    D=zeros(nx,nx);
	if ((nb=ddmat(rtk,D))<=0) {
        errmsg(rtk,"no valid double-difference\n");
        free(D);
		return 0;
    }
	ny=na+nb; y=mat(ny,1); Qy=mat(ny,ny); DP=mat(ny,nx);
    b=mat(nb,2); db=mat(nb,1); Qb=mat(nb,nb); Qab=mat(na,nb); QQ=mat(na,nb);

    /* transform single to double-differenced phase-bias (y=D'*x, Qy=D'*P*D) */
	matmul("TN",ny, 1,nx,1.0,D ,rtk->x,0.0,y );
	matmul("TN",ny,nx,nx,1.0,D ,rtk->P,0.0,DP);
    matmul("NN",ny,ny,nx,1.0,DP,D     ,0.0,Qy);
    
    /* phase-bias covariance (Qb) and real-parameters to bias covariance (Qab) */
    for (i=0;i<nb;i++) for (j=0;j<nb;j++) Qb [i+j*nb]=Qy[na+i+(na+j)*ny];
    for (i=0;i<na;i++) for (j=0;j<nb;j++) Qab[i+j*na]=Qy[   i+(na+j)*ny];
    
    /* lambda/mlambda integer least-square estimation */
	if (!(info=lambda(nb,2,y+na,Qb,b,s))) {

        trace(4,"N(1)="); tracemat(4,b   ,1,nb,10,3);
        trace(4,"N(2)="); tracemat(4,b+nb,1,nb,10,3);
        
		rtk->sol.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
		if (rtk->sol.ratio>999.9) rtk->sol.ratio=999.9f;
        
        /* validation by popular ratio-test */
        if (s[0]<=0.0||s[1]/s[0]>=opt->thresar[0]) {
            
            /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
			for (i=0;i<na;i++) {
                rtk->xa[i]=rtk->x[i];
				for (j=0;j<na;j++) rtk->Pa[i+j*na]=rtk->P[i+j*nx];
            }
            for (i=0;i<nb;i++) {
                bias[i]=b[i];
                y[na+i]-=b[i];
            }
            if (!matinv(Qb,nb)) {
                matmul("NN",nb,1,nb, 1.0,Qb ,y+na,0.0,db);
                matmul("NN",na,1,nb,-1.0,Qab,db  ,1.0,rtk->xa);
                
                /* covariance of fixed solution (Qa=Qa-Qab*Qb^-1*Qab') */
                matmul("NN",na,nb,nb, 1.0,Qab,Qb ,0.0,QQ);
                matmul("NT",na,na,nb,-1.0,QQ ,Qab,1.0,rtk->Pa);
                
                trace(3,"resamb : validation ok (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                      nb,s[0]==0.0?0.0:s[1]/s[0],s[0],s[1]);
                
				/* restore single-differenced ambiguity */
				restamb(rtk,bias,nb,xa);
            }
            else nb=0;
        }
        else { /* validation failed */
		    errmsg(rtk,"ambiguity validation failed (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                   nb,s[1]/s[0],s[0],s[1]);
			nb=0;
        }
    }
    else {
        errmsg(rtk,"lambda error (info=%d)\n",info);
        nb=0;
    }
    free(D); free(y); free(Qy); free(DP);
    free(b); free(db); free(Qb); free(Qab); free(QQ);
    
    return nb; /* number of ambiguities */
}


/* validation of solution ----------------------------------------------------*/
static int valpos(rtk_t *rtk, const double *v, const double *R, const int *vflg,
                  int nv, double thres)
{
#if 0
    prcopt_t *opt=&rtk->opt;
    double vv=0.0;
#endif
    double fact=thres*thres;
    int i,stat=1,sat1,sat2,type,freq;
    char *stype;
    
    trace(3,"valpos  : nv=%d thres=%.1f\n",nv,thres);
    
    /* post-fit residual test */
    for (i=0;i<nv;i++) {
        if (v[i]*v[i]<=fact*R[i+i*nv]) continue;
        sat1=(vflg[i]>>16)&0xFF;
        sat2=(vflg[i]>> 8)&0xFF;
        type=(vflg[i]>> 4)&0xF;
        freq=vflg[i]&0xF;
        stype=type==0?"L":(type==1?"L":"C");
        errmsg(rtk,"large residual (sat=%2d-%2d %s%d v=%6.3f sig=%.3f)\n",
              sat1,sat2,stype,freq+1,v[i],SQRT(R[i+i*nv]));
    }
#if 0 /* omitted v.2.4.0 */
    if (stat&&nv>NP(opt)) {
        
        /* chi-square validation */
        for (i=0;i<nv;i++) vv+=v[i]*v[i]/R[i+i*nv];
        
        if (vv>chisqr[nv-NP(opt)-1]) {
            errmsg(rtk,"residuals validation failed (nv=%d np=%d vv=%.2f cs=%.2f)\n",
                   nv,NP(opt),vv,chisqr[nv-NP(opt)-1]);
            stat=0;
        }
        else {
            trace(3,"valpos : validation ok (%s nv=%d np=%d vv=%.2f cs=%.2f)\n",
                  rtk->tstr,nv,NP(opt),vv,chisqr[nv-NP(opt)-1]);
        }
    }
#endif
    return stat;
}

static int avepres(rtk_t *rtk, const obsd_t *obs)
{
	ssat_t *ssat;
	int isys;
	int i,j,k,f;
	int v[NSYS][NFREQ]={0};
	int nsys[NFREQ]={0};

	double sum[NSYS][NFREQ]={0};
	int n[NSYS][NFREQ]={0};

	for (i=0,k=0;i<MAXSAT;i++) {
		ssat=rtk->ssat+i;
		if (!ssat->vs) continue;
		isys=sysind(ssat->sys);
		for (j=0;j<NF(&rtk->opt);j++) {
			f=rtk->opt.oprfrq[j];
			if(ssat->resc[f] == 0.0) continue;
			if(v[isys-1][j] == 0) {
				v[isys-1][j] = 1;
				++nsys[j];
			}
		}
	}
	for (j=0;j<NF(&rtk->opt);j++) {
		f=rtk->opt.oprfrq[j];
		if(nsys[j] < 2)        continue;
		if(v[ISYSGPS-1][j] == 0) continue;

		for (i=0,k=0;i<MAXSAT;i++) {
			ssat=rtk->ssat+i;
			if (!ssat->vs) continue;

			isys=sysind(ssat->sys);
			if(v[isys-1][j] == 0) continue;

			++n[isys-1][j];
			sum[isys-1][j] += ssat->resc[f];
			++k;
		}
	}
	for (isys=ISYSGPS;isys<=NSYS;isys++) {
		for (j=0;j<NF(&rtk->opt);j++) {
			if(n[isys-1][j] == 0 ) continue;

			rtk->estisb.n[isys-1][j] += n[isys-1][j];
			rtk->estisb.ave[isys-1][j] += (sum[isys-1][j] - n[isys-1][j]*rtk->estisb.ave[isys-1][j])/rtk->estisb.n[isys-1][j];
		}
	}
	return 0;
}

/* relative positioning ------------------------------------------------------*/
static int relpos(rtk_t *rtk, const obsd_t *obs, int nu, int nr,
				  const nav_t *nav)
{
	prcopt_t *opt=&rtk->opt;
	gtime_t time=obs[0].time;
	double *rs,*dts,*var,*y,*e,*azel,*v,*H,*R,*xp,*Pp,*xa,*bias,dt;
	int i,j,k,f,n=nu+nr,ns,ny,nv,sat[MAXSAT],iu[MAXSAT],ir[MAXSAT],niter;
	int info,vflg[MAXOBS*NFREQ*2+1],svh[MAXOBS*2];
	int stat=rtk->opt.mode<=PMODE_DGPS?SOLQ_DGPS:SOLQ_FLOAT;
	int nf=NF(opt),slip;
	int *ux;
	int t;

	trace(3,"relpos  : nx=%d nu=%d nr=%d\n",rtk->nx,nu,nr);

	dt=timediff(time,obs[nu].time);

	rs=mat(6,n); dts=mat(2,n); var=mat(1,n); y=mat(nf*2+2*TEST_ADDESTPRM,n); e=mat(3,n);
	azel=zeros(2,n);

	for (i=0;i<MAXSAT;i++) {
		rtk->ssat[i].sys=satsys(i+1,NULL);
		for (j=0;j<NFREQ;j++) rtk->ssat[i].vsat[j]=rtk->ssat[i].snr[j]=0;
	}
	/* satellite positions/clocks */
	satposs(time,obs,n,nav,opt->sateph,rs,dts,var,svh);

	/* undifferenced residuals for base station */
	if (!zdres(1,obs+nu,nr,rs+nu*6,dts+nu*2,svh+nu,nav,rtk->rb,opt,1,
			   y+nu*(nf*2+2*TEST_ADDESTPRM),e+nu*3,azel+nu*2)) {
		errmsg(rtk,"initial base station position error\n");

		free(rs); free(dts); free(var); free(y); free(e); free(azel);
		return 0;
	}
	/* time-interpolation of residuals (for post-processing) */
	if (opt->intpref) {
		dt=intpres(time,obs+nu,nr,nav,rtk,y+nu*(nf*2+2*TEST_ADDESTPRM));
	}
	/* select common satellites between rover and base-station */
	if ((ns=selsat(obs,azel,nu,nr,opt,sat,iu,ir))<=0) {
		errmsg(rtk,"no common satellite\n");

		free(rs); free(dts); free(var); free(y); free(e); free(azel);
		return 0;
	}
	/* temporal update of states */
	trace(4,"x(0)="); tracemat(4,rtk->x,1,NR(opt),13,12);
	udstate(rtk,obs,sat,iu,ir,ns,nav);

	trace(4,"x(0)="); tracemat(4,rtk->x,1,NR(opt),13,12);

	xp=zeros(rtk->nx,1); Pp=zeros(rtk->nx,rtk->nx); xa=zeros(rtk->nx,1);
	ux=izeros(rtk->nx,1);
	matcpy(xp,rtk->x,rtk->nx,1);

	ny=ns*(nf*2+2*TEST_ADDESTPRM)+2;
	v=zeros(ny,1); H=zeros(rtk->nx,ny); R=zeros(ny,ny); bias=zeros(rtk->nx,1);

	/* add 2 iterations for baseline-constraint moving-base */
	niter=opt->niter+(opt->mode==PMODE_MOVEB&&opt->baseline[0]>0.0?2:0);

	for (i=0;i<niter;i++) {
		/* undifferenced residuals for rover */
		if (!zdres(0,obs,nu,rs,dts,svh,nav,xp,opt,0,y,e,azel)) {
			errmsg(rtk,"rover initial position error\n");
			stat=SOLQ_NONE;
			break;
		}
		/* double-differenced residuals and partial derivatives */
		if ((nv=ddres(rtk,nav,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,H,R,vflg,obs,ux))<1) {
			errmsg(rtk,"no double-differenced residual\n");
			stat=SOLQ_NONE;
			break;
		}
		trace(5," xp[%15.12f][%15.12f][%15.12f] \n",xp[0],xp[1],xp[2]);
		trace(5," v=\n"); tracemat(5,v,1,ny,13,12);

		/* kalman filter measurement update */
		matcpy(Pp,rtk->P,rtk->nx,rtk->nx);

		trace(5," ny=%d nv=%d nx=%d\n",ny,nv,rtk->nx);
		if ((info=filter(xp,Pp,H,v,R,rtk->nx,nv,ux))) {
			errmsg(rtk,"filter error (info=%d)\n",info);
			stat=SOLQ_NONE;
			break;
		}
		trace(4,"x(%d)=",i+1); tracemat(4,xp,1,NR(opt),13,12);
	}
	if (stat!=SOLQ_NONE&&zdres(0,obs,nu,rs,dts,svh,nav,xp,opt,0,y,e,azel)) {

		/* post-fit residuals for float solution */
		nv=ddres(rtk,nav,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,NULL,R,vflg,obs,ux);

		/* validation of float solution */
		if (valpos(rtk,v,R,vflg,nv,4.0)) {

			/* update state and covariance matrix */
			matcpy(rtk->x,xp,rtk->nx,1);
			matcpy(rtk->P,Pp,rtk->nx,rtk->nx);
			imatcpy(rtk->ux,ux,rtk->nx,1);

			/* update ambiguity control struct */
			rtk->sol.ns=0;
			rtk->sol.sys=0x00;
			for (i=0;i<ns;i++) for (j=0;j<nf;j++) {
				f=rtk->opt.oprfrq[j];
				if (!rtk->ssat[sat[i]-1].vsat[f]) {
					continue;
				}
				rtk->ssat[sat[i]-1].lock[f]++;
				rtk->ssat[sat[i]-1].outc[f]=0;
#ifndef RTKPOS4PCV
				if (j==0) rtk->sol.ns++; /* valid satellite count by oprfrq[0] */
#else
				if (opt->antestmode == ANTEST_MODE_NONE) {
					if (f==0) rtk->sol.ns++; /* valid satellite count by L1 */
				} else {
					if (f==opt->antestfreq) rtk->sol.ns++;
				}
#endif
				rtk->sol.sys |= rtk->ssat[sat[i]-1].sys;
				for(k=0;k<NFREQ;k++) {
					if(rtk->ssat[sat[i]-1].vsat[k]==1) {
						rtk->sol.sys_prd[k] |= rtk->ssat[sat[i]-1].sys;
					}
				}
			}
			/* lack of valid satellites */
			if (rtk->sol.ns<4) stat=SOLQ_NONE;
		}
		else stat=SOLQ_NONE;
	}

	if(stat!=SOLQ_NONE&&opt->ionoopt==IONOOPT_IFLC) {
		if(opt->modertkar==RTKAR_LC) {
			/*WL/NLによるアンビキュイティ決定*/
			if(Amb_WLNL(rtk,obs,nav,sat,ns,iu,ir)){
				stat=SOLQ_FIX;
			}
		}
		else if(opt->modertkar==RTKAR_TCAR) {
			/* resolve integer ambiguity by TCAR */
			if (resamb_TCAR(rtk,obs,sat,iu,ir,ns,nav,azel)) {
				stat=SOLQ_FIX;
			}
		}
	}
	/* resolve integer ambiguity by LAMBDA */
	else if (stat!=SOLQ_NONE&&opt->modertkar==RTKAR_LAMBDA) {
		if(resamb_LAMBDA(rtk,bias,xa)>1) {
			if (zdres(0,obs,nu,rs,dts,svh,nav,xa,opt,0,y,e,azel)) {

				/* post-fit reisiduals for fixed solution */
				nv=ddres(rtk,nav,dt,xa,NULL,sat,y,e,azel,iu,ir,ns,v,NULL,R,vflg,obs,ux);

				/* validation of fixed solution */
				if (valpos(rtk,v,R,vflg,nv,4.0)) {

					/* hold integer ambiguity */
					if (++rtk->nfix>=rtk->opt.minfix&&
						rtk->opt.modear==ARMODE_FIXHOLD) {
						holdamb(rtk,xa);
					}
					stat=SOLQ_FIX;
				}
			}
		}
	}

	/* save solution status */
	if (stat==SOLQ_FIX) {
		for (i=0;i<3;i++) {
			rtk->sol.rr[i]=rtk->xa[i];
			rtk->sol.qr[i]=(float)rtk->Pa[i+i*rtk->na];
			if(opt->dynamics) {
				rtk->sol.rr[i+3]=rtk->xa[i+3];
				rtk->sol.rr[i+6]=rtk->xa[i+6];
				rtk->sol.qr[i+6]=(float)rtk->Pa[(i+3)+(i+3)*rtk->na];
				rtk->sol.qr[i+12]=(float)rtk->Pa[(i+6)+(i+6)*rtk->na];
			}
		}
		rtk->sol.qr[3]=(float)rtk->Pa[1];
		rtk->sol.qr[4]=(float)rtk->Pa[1+2*rtk->na];
		rtk->sol.qr[5]=(float)rtk->Pa[2];
		if(opt->dynamics) {
			rtk->sol.qr[ 9]=(float)rtk->Pa[4+3*rtk->na];
			rtk->sol.qr[10]=(float)rtk->Pa[4+5*rtk->na];
			rtk->sol.qr[11]=(float)rtk->Pa[5+3*rtk->na];
			rtk->sol.qr[15]=(float)rtk->Pa[7+6*rtk->na];
			rtk->sol.qr[16]=(float)rtk->Pa[7+8*rtk->na];
			rtk->sol.qr[17]=(float)rtk->Pa[8+6*rtk->na];
		}
		if ((rtk->opt.isb==ISBOPT_EST) || (rtk->opt.isb==ISBOPT_EST_L)) {
			for(i=ISYSGPS;i<=NSYS;i++) {
				for (j=0;j<nf;j++) {
					rtk->sol.isb[i-ISYSGPS][j][0]=rtk->xa[IS(i,0,j,&rtk->opt)];
				}
			}
		}
		if ((rtk->opt.isb==ISBOPT_EST) || (rtk->opt.isb==ISBOPT_EST_P)) {
			for(i=ISYSGPS;i<=NSYS;i++) {
				for (j=0;j<nf;j++) {
					rtk->sol.isb[i-ISYSGPS][j][1]=rtk->xa[IS(i,1,j,&rtk->opt)];
				}
			}
		}
		if(rtk->opt.isb==ISBOPT_EST_0M) {
			for(i=ISYSGPS;i<=NSYS;i++) {
				for (j=0;j<nf;j++) {
					rtk->sol.isb[i-ISYSGPS][j][1]=rtk->xa[IS(i,1,j,&rtk->opt)];
				}
			}
			avepres(rtk, obs);
			for(i=ISYSGPS;i<=NSYS;i++) {
				for (j=0;j<nf;j++) {
					rtk->sol.isb[i-ISYSGPS][j][0]=-rtk->estisb.ave[i-1][j];
				}
			}
		}

		if(opt->gl2bias==GL2OPT_EST) {
			for(i=0;i<=2;i++) {
				rtk->sol.gl2[i]=rtk->xa[I2(i,&rtk->opt)];
			}
		}
		for(i=0;i<NT(opt);i++) {
			rtk->sol.trop[i]=rtk->xa[IT(0,opt)+i];
		}
		for(i=0;i<MAXSAT;i++) {
			rtk->sol.ion[i]=rtk->xa[II(i+1,opt)];
		}
	}
	else {
		for (i=0;i<3;i++) {
			rtk->sol.rr[i]=rtk->x[i];
			rtk->sol.qr[i]=(float)rtk->P[i+i*rtk->nx];
			if(opt->dynamics) {
				rtk->sol.rr[i+3]=rtk->x[i+3];
				rtk->sol.rr[i+6]=rtk->x[i+6];
				rtk->sol.qr[i+6]=(float)rtk->P[(i+3)+(i+3)*rtk->nx];
				rtk->sol.qr[i+12]=(float)rtk->P[(i+6)+(i+6)*rtk->nx];
			}
		}
		rtk->sol.qr[3]=(float)rtk->P[1];
		rtk->sol.qr[4]=(float)rtk->P[1+2*rtk->nx];
		rtk->sol.qr[5]=(float)rtk->P[2];
		if(opt->dynamics) {
			rtk->sol.qr[ 9]=(float)rtk->P[4+3*rtk->nx];
			rtk->sol.qr[10]=(float)rtk->P[4+5*rtk->nx];
			rtk->sol.qr[11]=(float)rtk->P[5+3*rtk->nx];
			rtk->sol.qr[15]=(float)rtk->P[7+6*rtk->nx];
			rtk->sol.qr[16]=(float)rtk->P[7+8*rtk->nx];
			rtk->sol.qr[17]=(float)rtk->P[8+6*rtk->nx];
		}
		if(rtk->nfix==0) {

			if ((rtk->opt.isb==ISBOPT_EST) || (rtk->opt.isb==ISBOPT_EST_L)) {
				for(i=ISYSGPS;i<=NSYS;i++) {
					for (j=0;j<nf;j++) {
						rtk->sol.isb[i-ISYSGPS][j][0]=rtk->xa[IS(i,0,j,&rtk->opt)];
					}
				}
			}
			if ((rtk->opt.isb==ISBOPT_EST) || (rtk->opt.isb==ISBOPT_EST_P) || (rtk->opt.isb==ISBOPT_EST_0M)) {
				for(i=ISYSGPS;i<=NSYS;i++) {
					for (j=0;j<nf;j++) {
						rtk->sol.isb[i-ISYSGPS][j][1]=rtk->xa[IS(i,1,j,&rtk->opt)];
					}
				}
			}
			if(opt->gl2bias==GL2OPT_EST) {
				for(i=0;i<=2;i++) rtk->sol.gl2[i]=rtk->x[I2(i,&rtk->opt)];
			}
		}
		for(i=0;i<NT(opt);i++) {
			rtk->sol.trop[i]=rtk->x[IT(0,opt)+i];
		}
		for(i=0;i<MAXSAT;i++) {
			rtk->sol.ion[i]=rtk->x[II(i+1,opt)];
		}
		rtk->nfix=0;
	}
	for (i=0;i<n;i++) for (j=0;j<nf;j++) {
		f=rtk->opt.oprfrq[j];
		if (obs[i].L[f]==0.0) continue;
		rtk->ssat[obs[i].sat-1].pt[obs[i].rcv-1][f]=obs[i].time;
		rtk->ssat[obs[i].sat-1].ph[obs[i].rcv-1][f]=obs[i].L[f];
	}
	for (i=0;i<ns;i++) for (j=0;j<nf;j++) {
		f=rtk->opt.oprfrq[j];
		/* output snr of rover receiver */
		rtk->ssat[sat[i]-1].snr[f]=MIN(obs[iu[i]].SNR[f],obs[ir[i]].SNR[f]);
	}
	for (i=0;i<MAXSAT;i++) for (j=0;j<nf;j++) {
		f=rtk->opt.oprfrq[j];
		if (rtk->ssat[i].fix[f]==2&&stat!=SOLQ_FIX) rtk->ssat[i].fix[f]=1;
		slip=rtk->ssat[i].slip[f];
		/*電離層フリーの場合、2波のbit or を取る*/
		if (rtk->opt.ionoopt==IONOOPT_IFLC) {
			slip|=rtk->ssat[i].slip[1];
		}
		if (slip&1) rtk->ssat[i].slipc[f]++;
	}

	if (stat!=SOLQ_NONE) rtk->sol.stat=stat;

	/*単一基線解析結果ファイルに出力する統計データを計算する関数を挿入する。*/
#ifndef RTKPOS4PCV
	stat_stock(rtk,stat,rtk->opt.nfreq,nu,obs,nv,v,vflg);
#else
	/* output antenna file */
	outantest(rtk);
#endif

	free(rs); free(dts); free(var); free(y); free(e); free(azel);
	free(xp); free(Pp);  free(xa);  free(v); free(H); free(R); free(bias);
	free(ux);

	return stat!=SOLQ_NONE;
}

/* precise positioning ---------------------------------------------------------
* input observation data and navigation message, compute rover position by 
* precise positioning
* args   : rtk_t *rtk       IO  rtk control/result struct
*            rtk->sol       IO  solution
*                .time      O   solution time
*                .rr[]      IO  rover position/velocity
*                               (I:fixed mode,O:single mode)
*                .dtr[0]    O   receiver clock bias (s)
*                .dtr[1]    O   receiver glonass-gps time offset (s)
*                .Qr[]      O   rover position covarinace
*                .stat      O   solution status (SOLQ_???)
*                .ns        O   number of valid satellites
*                .age       O   age of differential (s)
*                .ratio     O   ratio factor for ambiguity validation
*            rtk->rb[]      IO  base station position/velocity
*                               (I:relative mode,O:moving-base mode)
*            rtk->nx        I   number of all states
*            rtk->na        I   number of integer states
*            rtk->ns        O   number of valid satellite
*            rtk->tt        O   time difference between current and previous (s)
*            rtk->x[]       IO  float states pre-filter and post-filter
*            rtk->P[]       IO  float covariance pre-filter and post-filter
*            rtk->xa[]      O   fixed states after AR
*            rtk->Pa[]      O   fixed covariance after AR
*            rtk->ssat[s]   IO  sat(s+1) status
*                .sys       O   system (SYS_???)
*                .az   [r]  O   azimuth angle   (rad) (r=0:rover,1:base)
*                .el   [r]  O   elevation angle (rad) (r=0:rover,1:base)
*                .vs   [r]  O   data valid single     (r=0:rover,1:base)
*                .resp [f]  O   freq(f+1) pseudorange residual (m)
*                .resc [f]  O   freq(f+1) carrier-phase residual (m)
*                .vsat [f]  O   freq(f+1) data vaild (0:invalid,1:valid)
*                .fix  [f]  O   freq(f+1) ambiguity flag
*                               (0:nodata,1:float,2:fix,3:hold)
*                .slip [f]  O   freq(f+1) slip flag
*                               (bit8-7:rcv1 LLI, bit6-5:rcv2 LLI,
*                                bit2:parity unknown, bit1:slip)
*                .lock [f]  IO  freq(f+1) carrier lock count
*                .outc [f]  IO  freq(f+1) carrier outage count
*                .slipc[f]  IO  freq(f+1) cycle slip count
*                .rejc [f]  IO  freq(f+1) data reject count
*                .gf        IO  geometry-free phase (L1-L2) (m)
*                .gf2       IO  geometry-free phase (L1-L5) (m)
*            rtk->nfix      IO  number of continuous fixes of ambiguity
*            rtk->neb       IO  bytes of error message buffer
*            rtk->errbuf    IO  error message buffer
*            rtk->tstr      O   time string for debug
*            rtk->opt       I   processing options
*          obsd_t *obs      I   observation data for an epoch
*                               obs[i].rcv=1:rover,2:reference
*                               sorted by receiver and satellte
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation messages
* return : status (0:no solution,1:valid solution)
* notes  : before calling function, base station position rtk->sol.rb[] should
*          be properly set for relative mode except for moving-baseline
*-----------------------------------------------------------------------------*/
extern int rtkpos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	prcopt_t *opt=&rtk->opt;
    sol_t solb={{0}};
    gtime_t time;
    int i,j,nu,nr;
	char msg[128]="";

	trace(3,"rtkpos  : time=%s n=%d\n",time_str(obs[0].time,3),n);
	trace(4,"obs=\n"); traceobs(4,obs,n);

	/* set base staion position */
	if (opt->refpos<=3&&opt->mode!=PMODE_SINGLE&&opt->mode!=PMODE_MOVEB) {
		for (i=0;i<6;i++) rtk->rb[i]=i<3?opt->rb[i]:0.0;
		trace(3,"rtk->rb=%f,%f,%f\n", rtk->rb[0], rtk->rb[1], rtk->rb[2]);
	}
	/* count rover/base station observations */
    for (nu=0;nu   <n&&obs[nu   ].rcv==1;nu++) ;
    for (nr=0;nu+nr<n&&obs[nu+nr].rcv==2;nr++) ;
   
    time=rtk->sol.time; /* previous epoch */
    
    /* rover position by single point positioning */
	if (!pntpos(obs,nu,nav,&rtk->opt,&rtk->sol,NULL,rtk->ssat,msg)) {
		errmsg(rtk,"rover point pos error (%s)\n",msg);
        
		if (!rtk->opt.dynamics) {
			outsolstat(rtk, obs);
            return 0;
        }
    }
    /* base station position by single point positioning */
	if ((opt->mode==PMODE_KINEMA) || (opt->mode==PMODE_STATIC) || (opt->mode==PMODE_FIXED)) {
        /* estimate position/velocity of base station */
		if (!pntpos(obs+nu,nr,nav,&rtk->opt,&solb,NULL,NULL,msg)) {
			errmsg(rtk,"base station position error (%s)\n",msg);
            if (!rtk->opt.dynamics) {
                return 0;
            }
        }
    }

    if (time.time!=0) rtk->tt=timediff(rtk->sol.time,time);
    
    /* single point positioning */
	if (opt->mode==PMODE_SINGLE) {
		pntposoutsolstat(rtk, obs, statlevel, fp_stat);
		return 1;
	}

	/* precise point positioning */
	if ((opt->mode==PMODE_PPP_KINEMA) || (opt->mode==PMODE_PPP_STATIC) || (opt->mode==PMODE_PPP_FIXED)) {
		pppos(rtk,obs,nu,nav);
		pppoutsolstat(rtk,statlevel,fp_stat,obs);
		return 1;
    }
	/* check number of data of base station and age of differential */
	if (nr==0) {
		errmsg(rtk,"no base station observation data for rtk\n");
		outsolstat(rtk, obs);
        return 1;
    }
    if (opt->mode==PMODE_MOVEB) { /*  moving baseline */
        
        /* estimate position/velocity of base station */
		if (!pntpos(obs+nu,nr,nav,&rtk->opt,&solb,NULL,NULL,msg)) {
            errmsg(rtk,"base station position error (%s)\n",msg);
            return 0;
        }
		rtk->sol.age=(float)timediff(rtk->sol.time,solb.time);
        
        if (fabs(rtk->sol.age)>TTOL_MOVEB) {
			errmsg(rtk,"time sync error for moving-base (age=%.1f)\n",rtk->sol.age);
            return 0;
        }
        for (i=0;i<6;i++) rtk->rb[i]=solb.rr[i];
        
        /* time-synchronized position of base station */
        for (i=0;i<3;i++) rtk->rb[i]+=rtk->rb[i+3]*rtk->sol.age;
    }
	else {
		rtk->sol.age=(float)timediff(obs[0].time,obs[nu].time);
		if (fabs(rtk->sol.age)>30) {
		rtk->sol.age=(float)timediff(obs[0].time,obs[nu].time);
		}
		if (fabs(rtk->sol.age)>opt->maxtdiff) {
			errmsg(rtk,"age of differential error (age=%.1f)\n",rtk->sol.age);
            //outsolstat(rtk);
			outsolstat(rtk, obs);
            return 1;
        }
	}
	/* relative potitioning */
	relpos(rtk,obs,nu,nr,nav);

	outsolstat(rtk, obs);
	outrnxobsb_ion(obs, nav, rtk->sol.ion, nu,0);

	return 1;
}
/*  -----------------------------------------------------
* 補正テーブルからGLONASS IFB補正
* args   : nav_t *nav       I  航法データ
*          prcopt_t *opt     I  解析処理オプション
*          int sta1          I  観測点1
*          int sta2          I  観測点2
*          int     f         I  周波数識別（f=0、L1；f=1、L2）
*          double　f1        I  観測点1のGLONASS衛星周波数
*          double　f2        I  観測点2のGLONASS衛星周波数
* return : double  IFB補正値
* note   :
*-------------------------------------------------------*/
double gloifbcorr(const nav_t *nav, prcopt_t *opt,int sta1, int sta2,int f,double f1, double f2){
	const ifb_t *ifb = nav->ifbs;
	int i;
	double bias;

	trace(3,"gloifbcorr: sta1=%d sta2=%d\n",sta1,sta2);
#ifndef RTKPOS4PCV
    if(f>=NFREQGLO || f>=opt->nfreq) return 0.0;
#endif
	/*受信機名比較*/
	for(i=0;i<nav->nifbs;i++){
		if((strcmp(opt->rectype[0],(ifb+i)->rec2)==0)&&(strcmp(opt->rectype[1],(ifb+i)->rec1)==0)){

			/*補正値計算*/
			/*bias=（f1-f2）*gloifb より*/
			bias=(f1-f2)/1E6*((ifb+i)->ifbvar[f]);
			trace(5,"f1[%f] f2[%f]ifvar[%f] bias[%f]\n",f1,f2,(ifb+i)->ifbvar[f],bias);
			return bias;
#ifndef RTKPOS4PCV
		}else if((strcmp(opt->rectype[0],(ifb+i)->rec1)==0)&&(strcmp(opt->rectype[1],(ifb+i)->rec2)==0)){
			bias=(f1-f2)/1E6*((ifb+i)->ifbvar[f]);
			trace(5,"f1[%f] f2[%f]ifvar[%f] bias[%f]\n",f1,f2,(ifb+i)->ifbvar[f],-bias);
			return -bias;
#endif
        }
	}
	return 0;
}

/* carrier-phase LC (m) ------------------------------------------------------*/
static double L_LC(int i, int j, int k, const double *Li, const double *Lj)
{
    const double f1=FREQ1,f2=FREQ2,f5=FREQ5;
    double L1,L2,L5;
    
    if ((i&&(!Li[0]||!Lj[0]))||(j&&(!Li[1]||!Lj[1]))||(k&&(!Li[2]||!Lj[2]))) {
        return 0.0;
    }
	L1=CLIGHT/f1*(Li[0]-Lj[0]);
    L2=CLIGHT/f2*(Li[1]-Lj[1]);
    L5=CLIGHT/f5*(Li[2]-Lj[2]);
	return (i*f1*L1+j*f2*L2+k*f5*L5)/(i*f1+j*f2+k*f5);
}
/* pseudorange LC (m) --------------------------------------------------------*/
static double P_LC(int i, int j, int k, const double *Pi, const double *Pj)
{
    const double f1=FREQ1,f2=FREQ2,f5=FREQ5;
    double P1,P2,P5;
    
    if ((i&&(!Pi[0]||!Pj[0]))||(j&&(!Pi[1]||!Pj[1]))||(k&&(!Pi[2]||!Pj[2]))) {
        return 0.0;
    }
    P1=Pi[0]-Pj[0];
    P2=Pi[1]-Pj[1];
    P5=Pi[2]-Pj[2];
    return (i*f1*P1+j*f2*P2+k*f5*P5)/(i*f1+j*f2+k*f5);
}
/* noise variance of LC (m) --------------------------------------------------*/
static double var_LC(int i, int j, int k, double sig)
{
    const double f1=FREQ1,f2=FREQ2,f5=FREQ5;
    
	return (SQR(i*f1)+SQR(j*f2)+SQR(k*f5))/SQR(i*f1+j*f2+k*f5)*SQR(sig);
}
/* single-difference noise variance ------------------------------------------*/
static double SD_var(double var, double el)
{
    double sinel=sin(el);
    return 2.0*(var+var/sinel/sinel);
}
/* average LC ----------------------------------------------------------------*/
static void average_LC(rtk_t *rtk, const obsd_t *obs, const int *sat,
					   const int *iu, const int *ir, int ns, const nav_t *nav,
                       const double *azel)
{
	ambc_t *amb;
    double LC1,LC2,LC3,var1,var2,var3,err=rtk->opt.err[1]*rtk->opt.eratio[0];
    int i,j,k;
    
	for (i=0;i<ns;i++) {
        if (satsys(sat[i],NULL)!=SYS_GPS) continue;
        j=iu[i]; k=ir[i];

		/* triple-freq carrier and code LC (m) */
		LC1=L_LC(1,-1, 0,obs[j].L,obs[k].L)-P_LC(1,1,0,obs[j].P,obs[k].P);
		LC2=L_LC(0, 1,-1,obs[j].L,obs[k].L)-P_LC(0,1,1,obs[j].P,obs[k].P);
		LC3=L_LC(1,-6, 5,obs[j].L,obs[k].L)-P_LC(1,1,0,obs[j].P,obs[k].P);

		/* measurement noise variance (m) */
		var1=SD_var(var_LC(1,1,0,err),azel[1+2*iu[i]]);
		var2=SD_var(var_LC(0,1,1,err),azel[1+2*iu[i]]);
        var3=SD_var(var_LC(1,1,0,err),azel[1+2*iu[i]]);

		amb=rtk->ambc+sat[i]-1;

		if (rtk->ssat[sat[i]-1].slip[0]||rtk->ssat[sat[i]-1].slip[1]||
			rtk->ssat[sat[i]-1].slip[2]||amb->n[0]==0||
			fabs(timediff(amb->epoch[0],obs[0].time))>MIN_ARC_GAP) {

			amb->fixcnt=0;
			if (LC1) {
				amb->n[0]=1; amb->LC[0]=LC1; amb->LCv[0]=var1;
			}
            if (LC2) {
				amb->n[1]=1; amb->LC[1]=LC2; amb->LCv[1]=var2;
			}
			if (LC3) {
                amb->n[2]=1; amb->LC[2]=LC3; amb->LCv[2]=var3;
            }
        }
        else { /* averaging */
            if (LC1) {
				amb->LC [0]+=(LC1-amb->LC[0])/(++amb->n[0]);
                amb->LCv[0]+=(var1-amb->LCv[0])/amb->n[0];
            }
            if (LC2) {
				amb->LC [1]+=(LC2-amb->LC[1])/(++amb->n[1]);
                amb->LCv[1]+=(var2-amb->LCv[1])/amb->n[1];
            }
            if (LC3) {
				amb->LC [2]+=(LC3-amb->LC[2])/(++amb->n[2]);
                amb->LCv[2]+=(var3-amb->LCv[2])/amb->n[2];
            }
        }
        amb->epoch[0]=obs[0].time;
    }
}
static int fixsol(rtk_t *rtk, const double *v, const double *H, int nv)
{
    double *R;
    int i,j,info;
    
    if (nv<=0) return 0;
    
    R=zeros(nv,nv);
    
    for (i=0;i<nv;i++) R[i+i*nv]=0.0001;
    
    /* update states with constraints */
	if (!(info=filter(rtk->x,rtk->P,H,v,R,rtk->nx,nv,rtk->ux))) {
        
        /* set solution */
        for (i=0;i<rtk->na;i++) {
            rtk->xa[i]=rtk->x[i];
            for (j=0;j<rtk->na;j++) {
                rtk->Pa[i+j*rtk->na]=rtk->Pa[j+i*rtk->na]=rtk->P[i+j*rtk->nx];
            }
        }
    }
    else {
        trace(1,"filter error (info=%d)\n",info);
        nv=0;
    }
    free(R);
    return nv;
}
/* resolve integer ambiguity by WL-NL ----------------------------------------*/
extern int resamb_WLNL(rtk_t *rtk, const obsd_t *obs, const int *sat,
					   const int *iu, const int *ir, int ns, const nav_t *nav,
                       const double *azel)
{
    ambc_t *ambi,*ambj;
    double lam_N,lam_W,C1,C2,BW,vW,BC,vC,B1,v1,*H,*v;
    int i,j,k,l,NW,N1,N2,nv=0,fix;
    int ind[2];
    double lam[2];
    double p0;
    
	if (ns<=0||rtk->opt.ionoopt!=IONOOPT_IFLC||rtk->opt.nfreq<2) return 0;
    
    trace(3,"resamb_WLNL: time=%s ns=%d\n",time_str(obs[0].time,0),ns);
    
    ind[0]=rtk->opt.oprfrq[0];
    ind[1]=rtk->opt.oprfrq[1];
    lam[0]=nav->lam[sat[0]-1][ind[0]];
    lam[1]=nav->lam[sat[0]-1][ind[1]];

    if((ind[0]==0) && (ind[1]==1)) {
        lam_N=lam_LC(1, 1,0);
		lam_W=lam_LC(1,-1,0);
    }
    else if((ind[0]==0) && (ind[1]==2)) {
		lam_N=lam_LC(1,0, 1);
        lam_W=lam_LC(1,0,-1);
    }

    C1= SQR(lam[1])/(SQR(lam[1])-SQR(lam[0]));
    C2=-SQR(lam[0])/(SQR(lam[1])-SQR(lam[0]));
    
    v=zeros(ns,1); H=zeros(rtk->nx,ns);
    
    /* average LC */
	average_LC(rtk,obs,sat,iu,ir,ns,nav,azel);
    
    /* search reference */
    for (i=0,j=1;j<ns;j++) {
        ambi=rtk->ambc+sat[i]-1;
        ambj=rtk->ambc+sat[j]-1;
		if (ambj->n[0]>ambi->n[0]) i=j;
    }
    /* resolve double-difference ambiguity */
    for (j=0;j<ns;j++) {
        if (i==j) continue;
		ambi=rtk->ambc+sat[i]-1;
        ambj=rtk->ambc+sat[j]-1;
		k=IB(sat[i],0,&rtk->opt);
		l=IB(sat[j],0,&rtk->opt);

		if (!ambi->n[0]||!ambj->n[0]) continue;
        
        /* wide lane ambiguity (cycle) */
		BW=(ambi->LC[0]-ambj->LC[0])/lam_W;
        vW=(ambi->LCv[0]/ambi->n[0]+ambj->LCv[0]/ambj->n[0])/SQR(lam_W);
		NW=ROUND(BW);

		/* narrow lane ambiguity (cycle) */
		BC=rtk->x[k]-rtk->x[l];
		vC=rtk->P[k+k*rtk->nx]+rtk->P[l+l*rtk->nx]-2.0*rtk->P[k+l*rtk->nx];
        B1=(BC+C2*lam[1]*NW)/lam_N;
        v1=vC/SQR(lam_N);
        N1=ROUND(B1);
        N2=N1-NW;

		fix=fabs(NW-BW)<=FIX_THRES*1.5&&sqrt(vW)<=FIX_THRES&&
			fabs(N1-B1)<=FIX_THRES*1.5&&sqrt(v1)<=FIX_THRES;

		/* constraint to dd-ambiguity */
		v[nv]=(C1*lam[0]*N1+C2*lam[1]*N2)-(rtk->x[k]-rtk->x[l]);
        H[k+nv*rtk->nx]= 1.0;
        H[l+nv*rtk->nx]=-1.0;
        nv++;
    }
    /* fix solution */
    nv=fixsol(rtk,v,H,nv);
    
    free(v); free(H);
    
    return nv>=ns/2; /* fix if a half ambiguities fixed */
}

/* resolve integer ambiguity by TCAR -----------------------------------------*/
int resamb_TCAR(rtk_t *rtk, const obsd_t *obs, const int *sat,
                       const int *iu, const int *ir, int ns, const nav_t *nav,
                       const double *azel)
{
    ambc_t *ambi,*ambj;
    double lam_E,lam_F,lam_N,C1,C2,BE,BF,vE,vF,BC,vC,B1,v1,*H,*v;
    int i,j,k,l,nv=0,NE,NF,NW,N1,N2,fix;
    int ind[3];
    double lam[3];
	double a;
	double p0;

//	if (ns<=0||rtk->opt.ionoopt!=IONOOPT_IFLC||rtk->opt.nfreq<3) return 0;
	if (ns<=0||rtk->opt.ionoopt!=IONOOPT_IFLC||rtk->opt.modertkar!=RTKAR_TCAR) return 0;

    ind[0]=rtk->opt.oprfrq[0];
    ind[1]=rtk->opt.oprfrq[1];
    ind[2]=rtk->opt.oprfrq[2];
    lam[0]=nav->lam[sat[0]-1][ind[0]];
    lam[1]=nav->lam[sat[0]-1][ind[1]];
    lam[2]=nav->lam[sat[0]-1][ind[2]];

    trace(3,"resamb_TCAR: time=%s ns=%d\n",time_str(obs[0].time,0),ns);
    
    lam_E=lam_LC(0,1,-1);
	lam_F=lam_LC(1,-6,5);
    lam_N=lam_LC(1, 1,0);
	C1= SQR(lam[1])/(SQR(lam[1])-SQR(lam[0]));
    C2=-SQR(lam[0])/(SQR(lam[1])-SQR(lam[0]));
    
    v=zeros(ns,1); H=zeros(rtk->nx,ns);
    
    /* average LC */
	average_LC(rtk,obs,sat,iu,ir,ns,nav,azel);

    /* search reference */
    for (i=0,j=1;j<ns;j++) {
		ambi=rtk->ambc+sat[i]-1;
        ambj=rtk->ambc+sat[j]-1;
        if (ambj->n[0]>ambi->n[0]) i=j;
    }
    /* resolve double-difference bias */
    for (j=0;j<ns;j++) {
		if (i==j) continue;
        ambi=rtk->ambc+sat[i]-1;
		ambj=rtk->ambc+sat[j]-1;
		k=IB(sat[i],0,&rtk->opt);
		l=IB(sat[j],0,&rtk->opt);

        if (!ambi->n[1]||!ambj->n[1]||!ambi->n[2]||!ambj->n[2]) continue;
        
		/* extra wide lane ambiguity (cycle) */
		BE=(ambi->LC[1]-ambj->LC[1])/lam_E;
		BF=(ambi->LC[2]-ambj->LC[2])/lam_F;
		vE=(ambi->LCv[1]/ambi->n[1]+ambj->LCv[1]/ambi->n[1])/SQR(lam_E);
		vF=(ambi->LCv[2]/ambi->n[2]+ambj->LCv[2]/ambi->n[2])/SQR(lam_F);
		NE=ROUND(BE); NF=ROUND(BF); NW=5*NE+NF;


		if(!check_intdd( FIX_THRES_WLNL, BE, (double)NE, vE, &p0)) {
			continue;
		}
		if(!check_intdd( FIX_THRES_WLNL, BF, (double)NF, vF, &p0)) {
			continue;
		}

		/* narrow lane ambiguity (cycle) */
		BC=rtk->x[k]-rtk->x[l];
		vC=rtk->P[k+k*rtk->nx]+rtk->P[l+l*rtk->nx]-2.0*rtk->P[k+l*rtk->nx];
		B1=(BC+C2*lam[1]*NW)/lam_N;
		v1=vC/SQR(lam_N);
		N1=ROUND(B1); N2=N1-NW;

		if(!check_intdd( FIX_THRES_WLNL, B1, (double)N1, v1, &p0)) {
			continue;
		}

		/* constraint to dd-ambiguity */
        v[nv]=(C1*lam[0]*N1+C2*lam[1]*N2)-(rtk->x[k]-rtk->x[l]);
        H[k+nv*rtk->nx]= 1.0;
		H[l+nv*rtk->nx]=-1.0;
        nv++;
    }
    /* fix solution */
	nv=fixsol(rtk,v,H,nv);
    
    free(v); free(H);
    
    return nv>=ns/2; /* fix if a half ambiguities fixed */
}

/*線形結合観測量の時間平均 */
static void waveAvg( const int *satelitelist, int sateliteno,
					 const nav_t *nav,rtk_t *rtk, const obsd_t *obs,
					 const int *indexu, const int *indexr)
{
	double MW,var1,err=rtk->opt.err[1]*rtk->opt.eratio[0];
	int i,j,k;
	double sd1,tsin ;
    double lam[2];
    int ind[2];

	double c0;
	double p0;
	double sig,sigp;
	double c0v,p0v;
	double c1,p1;
	double ydif[2]={0.0,0.0};
	int a,b;

	trace(3,"waveAvg: \n");

	/*衛星数分ループ */
	for (i=0;i<sateliteno;i++) {

		if (satsys(satelitelist[i],NULL)!=SYS_GPS) continue;
		j=indexu[i]; k=indexr[i];
		if(obs[j].azel[1]<=0.0) continue;

		ind[0]=rtk->opt.oprfrq[0];
		ind[1]=rtk->opt.oprfrq[1];
		lam[0]=nav->lam[satelitelist[i]-1][ind[0]];
		lam[1]=nav->lam[satelitelist[i]-1][ind[1]];

		/* [2]観測誤差モデルによる標準偏差 */
		if (obs[j].azel[1]<rtk->opt.elmin) continue;
		sig = sqrt( SQR( rtk->opt.err[1] ) + SQR( rtk->opt.err[2] / sin( obs[j].azel[1] )));
		sigp = sig * rtk->opt.eratio[ind[0]];

		if ((obs[j].L[ind[0]]==0.0||obs[k].L[ind[0]]==0.0)||(obs[j].L[ind[1]]==0.0||obs[k].L[ind[1]]==0.0)) {
			c0 = 0.0;
			c0v = 0.0;
		}else{
			ydif[0] = chk_L2Csft(&obs[j], &rtk->opt, &stas[0]);
			ydif[1] = chk_L2Csft(&obs[k], &rtk->opt, &stas[1]);
			c0 = ( ( obs[j].L[ind[0]]-obs[k].L[ind[0]] ) + (-1)*((obs[j].L[ind[1]]-ydif[0])-(obs[k].L[ind[1]]-ydif[1]) ))
				  /(1/lam[0]+(-1)*1/lam[1]);
			c0v = ( lam[1]*lam[1] + lam[0]*lam[0] ) / ( lam[1] - lam[0] ) / ( lam[1] - lam[0] ) * sig * sig;

		}
		if ((obs[j].P[ind[0]]==0.0||obs[k].P[ind[0]]==0.0)||(obs[j].P[ind[1]]==0.0||obs[k].P[ind[1]]==0.0)) {
			p0 = 0.0;
			p0v=0.0;
		}else{
			p0 = ( lam[1]*(obs[j].P[0]-obs[k].P[0]) + lam[0]*(obs[j].P[1]-obs[k].P[1]) )
				 /(lam[1]+lam[0]);
			p0v = ( lam[1]*lam[1] + lam[0]*lam[0] ) / ( lam[1] + lam[0] ) / ( lam[1] + lam[0] ) * sigp * sigp;
		}
		MW=c0-p0;
		var1 = c0v+ p0v;

		if (rtk->ssat[satelitelist[i]-1].slip[0]||rtk->ssat[satelitelist[i]-1].slip[1]||
			rtk->ssat[satelitelist[i]-1].slip[2]||rtk->pass[satelitelist[i]-1].nd==0||
			fabs(timediff(rtk->pass[satelitelist[i]-1].te,obs[0].time))>MIN_ARC_GAP) {

				rtk->pass[satelitelist[i]-1].nd=0;
				rtk->pass[satelitelist[i]-1].mw = 0.0;
				rtk->pass[satelitelist[i]-1].mwv= 0.0;
		}
		if (MW) {
			rtk->pass[satelitelist[i]-1].nd++;
			rtk->pass[satelitelist[i]-1].mw+= (MW-rtk->pass[satelitelist[i]-1].mw)/(rtk->pass[satelitelist[i]-1].nd);
			rtk->pass[satelitelist[i]-1].mwv +=(var1-rtk->pass[satelitelist[i]-1].mwv)/rtk->pass[satelitelist[i]-1].nd;
		}

		rtk->pass[satelitelist[i]-1].te =obs[0].time;
	}
}

/* Amb_WLNL---------------------------------------------------
* WL/NLによるアンビギュイティ決定
* args   : rtk_t *rtk        IO  RTKデータ構造体
*          const obsd_t *obs I   観測データ
*          const nav_t *nav  I   航法データ
*          const int *sat    I   衛星番号リスト
*          int num           I   衛星の数
*          const int *idx1   I   移動局の観測データのインデックス
*          const int *idx2   I   基準局の観測データのインデックス
*--------------------------------------------------------------*/
int Amb_WLNL(rtk_t *rtk, const obsd_t *obs,
			 const nav_t *nav, const int *sat ,
			 const int num ,
			  const int *idx1, const int *idx2)
{
	double lam_N,lam_W,C1,C2,BW,vW,BC,vC,B1,v1,p0,*H,*v,*R;
	int i,j,k,l,NW,N1,nv=0,info;
	double elmask;
    double lam[2];
	trace(3,"Amb_WLNL: time=%s num=%d\n",time_str(obs[0].time,0),num);

	if (num<=0||rtk->opt.ionoopt!=IONOOPT_IFLC||rtk->opt.modertkar!=RTKAR_LC) return 0;

	/*初期設定*/

    v=zeros(num,1); H=zeros(rtk->nx,num);
	/*MW線形結合観測量の算出と時間平均*/
	waveAvg(sat,num,nav,rtk,obs,idx1,idx2);

	elmask=rtk->opt.elmaskar>0.0?rtk->opt.elmaskar:rtk->opt.elmin;

	/*基準衛星の決定*/
	for (i=0,j=1;j<num;j++) {
		lam[0] = nav->lam[sat[i]-1][rtk->opt.oprfrq[0]];
		lam[1] = nav->lam[sat[i]-1][rtk->opt.oprfrq[1]];
		if((j==1) || (rtk->pass[sat[j]-1].nd > rtk->pass[sat[i]-1].nd)) {
			if(rtk->pass[sat[j]-1].nd > rtk->pass[sat[i]-1].nd){
				i=j;
			}
			lam[0] = nav->lam[sat[i]-1][rtk->opt.oprfrq[0]];
			lam[1] = nav->lam[sat[i]-1][rtk->opt.oprfrq[1]];
			lam_N=1.0/(1.0/lam[0]+1.0/lam[1]);
			lam_W=1.0/(1.0/lam[0]-1.0/lam[1]);
			C1= SQR(lam[1])/(SQR(lam[1])-SQR(lam[0]));
			C2=-SQR(lam[0])/(SQR(lam[1])-SQR(lam[0]));
		}
	}
	/*衛星数分ループ*/
	for (j=0;j<num;j++) {
		if (i==j) continue;

		k=IB(sat[i],0,&rtk->opt);
		l=IB(sat[j],0,&rtk->opt);

		if((lam[0]!=nav->lam[sat[j]-1][rtk->opt.oprfrq[0]]) || (lam[1]!=nav->lam[sat[j]-1][rtk->opt.oprfrq[1]])) {
            continue;
        }

		if (!rtk->pass[sat[i]-1].nd||!rtk->pass[sat[j]-1].nd) {
            continue;
        }

		if (!rtk->ssat[sat[i]-1].vsat[0]||
			!rtk->ssat[sat[j]-1].vsat[0]||
			obs[i].azel[1]<elmask||obs[j].azel[1]<elmask) {
                continue;
        }
		/*MW線形結合の二重差より四捨五入でWLアンビギュイティ決定*/
		BW=(rtk->pass[sat[i]-1].mw-rtk->pass[sat[j]-1].mw)/lam_W;
		vW=(rtk->pass[sat[i]-1].mwv/rtk->pass[sat[i]-1].nd + rtk->pass[sat[j]-1].mwv/rtk->pass[sat[j]-1].nd)/SQR(lam_W);
		NW=ROUND(BW);

		trace(5,"i[%d]j[%d]mw1[%f] mw2[%f]sat[i]-1[%d] sat[j]-1[%d]\n",i,j,rtk->pass[sat[i]-1].mw,rtk->pass[sat[j]-1].mw,sat[i]-1,sat[j]-1);
		trace(5,"mwv[%f] nd[%d] mwv[%f] nd[%d] lam_W[%f] \n",rtk->pass[sat[i]-1].mwv,rtk->pass[sat[i]-1].nd,rtk->pass[sat[j]-1].mwv,rtk->pass[sat[j]-1].nd,lam_W);

		if( vW < 1e-30 ) {
			showmsg( "vW sigma = 0" );
			trace( 2, "vW sigma = 0\n" );
			continue;
		}
		if(!check_intdd( FIX_THRES_WLNL, BW, (double)NW, vW, &p0)) {
			trace(2,"BW, %d, %d, %f, %f, %f,false\n",sat[i], sat[j], BW,(double)NW,vW);
			continue;
		}
		trace(2,"BW, %d, %d, %f, %f, %f,true\n",sat[i], sat[j], BW,(double)NW,vW);
		/*電離層フリーより四捨五入でNLアンビギュイティ決定*/
		BC=rtk->x[k]-rtk->x[l];
		vC=rtk->P[k+k*rtk->nx]+rtk->P[l+l*rtk->nx]-2.0*rtk->P[k+l*rtk->nx];
		B1=(BC+C2*lam[1]*NW)/lam_N;
		v1=vC/SQR(lam_N);
		N1=ROUND(B1);

		if( v1 < 1e-30 ) {
			showmsg( "v1 sigma = 0" );
			trace( 2, "v1 sigma = 0\n" );
			continue;
		}
		trace(1,"Amb_WLNL %d %d %f %f %d %f %f %f %f %d \n",sat[i], sat[j], BW,vW,NW,BC,vC,B1,v1,N1);
		if(!check_intdd( FIX_THRES_WLNL, B1, (double)N1, v1, &p0)) {
			continue;
		}

		/*FIX検定OK?*/
		/*残差行列vと偏微分行列Hの更新*/
		v[nv]=(lam_N*N1-C2*lam[1]*NW)-(rtk->x[k]-rtk->x[l]);
		trace(2,"nv[%d] v[nv][%f] \n",nv,v[nv]);

		H[k+nv*rtk->nx]= 1.0;
		H[l+nv*rtk->nx]=-1.0;
		nv++;
    }
    /*FIX解の計算*/
    if (nv<=0) {
        free(v); free(H);
        return 0;
    }

    R=zeros(nv,nv);
	for (i=0;i<nv;i++) R[i+i*nv]=0.0001;
	info=filter(rtk->x,rtk->P,H,v,R,rtk->nx,nv,rtk->ux) ;
	if (info!=0) {
		trace(1,"filter error (info=%d)\n",info);
		nv=0;
	}else{

		for (i=0;i<rtk->na;i++) {
			rtk->xa[i]=rtk->x[i];
			for (j=0;j<rtk->na;j++) {
				rtk->Pa[i+j*rtk->na]=rtk->Pa[j+i*rtk->na]=rtk->P[i+j*rtk->nx];
			}
		}
	}

	free(R);free(v); free(H);
	return nv>=num/2;
}

extern int rtkit(const int r,const prcopt_t *opt) {
	return IT(r,opt);
}
extern int rtknt(const prcopt_t *opt) {
	return NT(opt);
}
extern int rtknr(const prcopt_t *opt)
{
	return NR(opt);
}
extern int rtknx(const prcopt_t *opt)
{
	return NX(opt);
}
