/*------------------------------------------------------------------------------
* pntpos.c : standard positioning
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

static const char rcsid[]="$Id:$";

/* constants -----------------------------------------------------------------*/
#define VAR_ISBP    SQR(0.001)      /*   ISB(pseudorange) (m^2) */

#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay std (m) */
#define ERR_TROP    3.0         /* tropspheric delay std (m) */
#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5         /* broadcast iono model error factor */
#define ERR_CBIAS   0.3         /* code bias error std (m) */
#define REL_HUMI    0.7         /* relative humidity for saastamoinen model */

#define NP(opt)     (3)
#define NC(opt)     (1+(    ((opt)->tsyscorr>1)  &&  ( ((opt)->isb!=ISBOPT_EST)&&((opt)->isb!=ISBOPT_EST_P)&&((opt)->isb!=ISBOPT_EST_0M) )?1:0))

#define NS(opt)     (((opt)->isb==ISBOPT_EST)||((opt)->isb==ISBOPT_EST_P)||((opt)->isb==ISBOPT_EST_0M)?(NSYS):0)
#define N2(opt)     ((((opt)->gl2bias!=GL2OPT_EST) || ((opt)->ionoopt!=IONOOPT_IFLC))?0:1)
#define NX(opt)     (NP(opt)+NC(opt)+NS(opt)+N2(opt))

#define IC(s,opt)   (NP(opt)+(NC(opt)==2?(s):0)) /* state index of clocks (s=0:gps,1:glo) */
#define IS(sys,opt) (NP(opt)+NC(opt)+(sys-NSYSGPS))   /* isb */
#define I2(opt)     (NP(opt)+NS(opt))   /* gl2bias */



/* pseudorange measurement error variance ------------------------------------*/
static double varerr(const prcopt_t *opt, double el, int sys, const obsd_t *obs, const nav_t *nav)
{
	double fact;
	double a=opt->err[1], b=opt->err[2];

	/*衛星システムから分散荷重係数fの決定*/
	fact=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);

	trace(4,"errtbl_before: a = %f,b = %f\n",a,b);
	/*観測誤差テーブルを使用する？*/
	if(opt->errmodel == ERRMODEL_TABLE){
		/* roverで検索 */
		int ret = search_errtbl(opt->rectype[0], opt->anttype[0], sys, obs->code[0], 0, nav, &a, &b);
		if (ret == 0) {
			/* baseで検索 */
			search_errtbl(opt->rectype[1], opt->anttype[1], sys, obs->code[0], 0, nav, &a, &b);
		}
	}
	trace(4,"errtbl_after: a = %f,b = %f\n",a,b);

	/*DCBを補正できない？*/
	if (has_dcb(obs->sat, obs->code[0], nav, NULL) == 0) {
		/*設定値Code Error Ratioでfを荷重する */
		fact *= opt->dcbratio;
	}

	/*電離層フリーモード*/
	if (opt->ionoopt == IONOOPT_IFLC) {
		/* 分散荷重係数 f=3.0*f */
		fact *= 3.0;
	}

	/*疑似距離観測誤差の計算*/
#if TEST_CMPMODEL == 0		//	プロトタイプと一致
	return SQR(fact)*(SQR(opt->err[0])*SQR(a) + SQR(b)/SQR(sin(el)));		
#elif TEST_CMPMODEL == 1	//	RTKLIBと同一
	return SQR(fact)*SQR(opt->err[0])*(SQR(a) + SQR(b)/sin(el));
#else						// 正解式
	return SQR(fact)*SQR(opt->err[0])*(SQR(a) + SQR(b)/SQR(sin(el)));
#endif
}

/* get tgd parameter (m) -----------------------------------------------------*/
static double gettgd(int sat, const nav_t *nav)
{
    int i;
    for (i=0;i<nav->n;i++) {
        if (nav->eph[i].sat!=sat) continue;
        return CLIGHT*nav->eph[i].tgd[0];
    }
    return 0.0;
}
/* psendorange with code bias correction -------------------------------------*/
static double prange(const obsd_t *obs, const nav_t *nav, const double *azel,
                     int iter, const prcopt_t *opt, double *var)
{
    const double *lam=nav->lam[obs->sat-1];
    double PC,P1,P2,P1_P2,P1_C1,P2_C2,gamma,rP1_C1=0.0,rP2_C2=0.0;
    int i=0,j=1,k,l,sys,res;
    int s_code;
    double c1,c2;
    double ydifIdb[NFREQ][2] = {0};
    *var=0.0;
    
    if (!(sys=satsys(obs->sat,NULL))) return 0.0;
    
    /* L1-L2 for GPS/GLO/QZS, L1-L5 for GAL/SBS */
    if (NFREQ>=3&&(sys&(SYS_GAL|SYS_SBS))) j=2;
    
    if (NFREQ<2||lam[i]==0.0||lam[j]==0.0) return 0.0;
    
    /* test snr mask */
    if (iter>0) {
        if (testsnr(0,i,azel[1],obs->SNR[i]*0.25,&opt->snrmask)) {
            trace(4,"snr mask: %s sat=%2d el=%.1f snr=%.1f\n",
                  time_str(obs->time,0),obs->sat,azel[1]*R2D,obs->SNR[i]*0.25);
            return 0.0;
        }
        if (opt->ionoopt==IONOOPT_IFLC) {
            if (testsnr(0,j,azel[1],obs->SNR[j]*0.25,&opt->snrmask)) return 0.0;
        }
    }
    gamma=SQR(lam[j])/SQR(lam[i]); /* f1^2/f2^2 */
    c1=gamma/(gamma-1.0);  /*  f1^2/(f1^2-f2^2) */
    c2=-1.0 /(gamma-1.0);  /* -f2^2/(f1^2-f2^2) */
    P1=obs->P[i];
    P2=obs->P[j];
    P1_P2=nav->cbias[obs->sat-1][0];
    P1_C1=nav->cbias[obs->sat-1][1];
    P2_C2=nav->cbias[obs->sat-1][2];

	s_code = sysind((char)sys);
	res=0;

    
    rP1_C1=stas[0].dcb[s_code-1][1];
    if(GL2OPT_TABLE==opt->gl2bias) rP2_C2=stas[0].dcb[s_code-1][2];

    /* if no P1-P2 DCB, use TGD instead */
    if (P1_P2==0.0&&(sys&(SYS_GPS|SYS_GAL|SYS_QZS))) {
        P1_P2=(1.0-gamma)*gettgd(obs->sat,nav);
    }

	/* isb correction */
	chk_isb(sys, opt, stas, ydifIdb);

	if (opt->ionoopt==IONOOPT_IFLC) { /* dual-frequency */
        
        if (P1==0.0||P2==0.0) return 0.0;
        if (obs->code[i]==CODE_L1C) P1+=P1_C1+rP1_C1; /* C1->P1 + 地上局C1->P1*/
        /*L2C/A?CODE_L2L?CODE_L2X?*/
        /*L2疑似距離の修正*/
        if ((obs->code[j]==CODE_L2C)||(obs->code[j]==CODE_L2L)||(obs->code[j]==CODE_L2X)) P2+=P2_C2+rP2_C2; /* C2->P2 + 地上局C2->P2*/
        
        /* iono-free combination */
		PC=(gamma*P1-P2)/(gamma-1.0);
        PC-=c1*ydifIdb[i][1] + c2*ydifIdb[j][1];
    }
    else { /* single-frequency */
        
        if (P1==0.0) return 0.0;
        if (obs->code[i]==CODE_L1C) P1+=P1_C1+rP1_C1; /* C1->P1 + 地上局C1->P1*/
		PC=P1-P1_P2/(1.0-gamma);
        PC-=ydifIdb[i][1];
    }
    if (opt->sateph==EPHOPT_SBAS) PC-=P1_C1; /* sbas clock based C1 */

    
    *var=SQR(ERR_CBIAS);
    
    return PC;
}

/* ionospheric correction ------------------------------------------------------
* compute ionospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          int    sat       I   satellite number
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    ionoopt   I   ionospheric correction option (IONOOPT_???)
*          double *ion      O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric delay (L1) variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var)
{
    trace(4,"ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
          time_str(time,3),ionoopt,sat,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* broadcast model */
    if (ionoopt==IONOOPT_BRDC) {
        *ion=ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* sbas ionosphere model */
    if (ionoopt==IONOOPT_SBAS) {
        return sbsioncorr(time,nav,pos,azel,ion,var);
    }
    /* ionex tec model */
    if (ionoopt==IONOOPT_TEC) {
        return iontec(time,nav,pos,azel,1,ion,var);
    }
    /* qzss broadcast model */
    if (ionoopt==IONOOPT_QZS&&norm(nav->ion_qzs,8)>0.0) {
        *ion=ionmodel(time,nav->ion_qzs,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* lex ionosphere model */
    if (ionoopt==IONOOPT_LEX) {
        return lexioncorr(time,nav,pos,azel,ion,var);
    }
    *ion=0.0;
    *var=ionoopt==IONOOPT_OFF?SQR(ERR_ION):0.0;
    return 1;
}
/* tropospheric correction -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    tropopt   I   tropospheric correction option (TROPOPT_???)
*          double *trp      O   tropospheric delay (m)
*          double *var      O   tropospheric delay variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var)
{
    trace(4,"tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f\n",
          time_str(time,3),tropopt,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* saastamoinen model */
    if (tropopt==TROPOPT_SAAS||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=tropmodel(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        return 1;
    }
    /* sbas troposphere model */
    if (tropopt==TROPOPT_SBAS) {
        *trp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    /* no correction */
    *trp=0.0;
    *var=tropopt==TROPOPT_OFF?SQR(ERR_TROP):0.0;
    return 1;
}
/* pseudorange residuals -----------------------------------------------------*/
//static int rescode(int iter, const obsd_t *obs, int n, const double *rs,
static int rescode(int iter, obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   double *v, double *H, double *var, double *azel, int *vsat,
				   double *resp, int *nx, int* ic, int* is, int* i2, int* ix)
{
	double r,dion,dtrp,vmeas,vion,vtrp,rr[3],pos[3],dtr,e[3],P;
	int i,j,nv=0,ns[NSYS]={0},sys,k;
    int isys,jsys;
	int sysa[NSYS];
	int nsys=0;
	int nxt;
    int jrm;
    double gamma;
    double c2;
	double tmp;
	int ncs[2]={0};
   
	trace(3,"resprng : n=%d\n",n);

    *nx=NX(opt);
	nxt = NX(opt);
    
	for (i=0;i<3;i++) rr[i]=x[i];
	dtr=x[IC(0,opt)];

    ic[0]=IC(0,opt);
	ic[1]=IC(1,opt);
	for (i=ISYSGPS;i<=NSYS;i++) is[i-ISYSGPS]=IS(i,opt);
	i2[0]=I2(opt);
	for (i=0;i<nxt;i++) ix[i] = 1;

	ecef2pos(rr,pos);

	for (i=0;i<n&&i<MAXOBS;i++) {
		vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;

		if (!(sys=satsys(obs[i].sat,NULL))) continue;

		isys = sysind(sys);

		/* reject duplicated observation data */
        if (i<n-1&&i<MAXOBS-1&&obs[i].sat==obs[i+1].sat) {
			trace(2,"duplicated observation data %s sat=%2d\n",
				  time_str(obs[i].time,3),obs[i].sat);
			i++;
			continue;
		}
		trace(3,"rr:%f,%f,%f\n", rr[0], rr[1], rr[2]);
		/* geometric distance/azimuth/elevation angle */
		r=geodist(rs+i*6,rr,e);
		if (r<=0.0){
			continue;
		}
		satazel(pos,e,azel+i*2);
		memcpy(obs[i].azel,azel+i*2,sizeof(double)*2);
		if(azel[1+i*2]<opt->elmin){
			 continue;
		}

		/* psudorange with code bias correction */
		if ((P=prange(obs+i,nav,azel+i*2,iter,opt,&vmeas))==0.0) continue;

		/* excluded satellite? */
		if (satexclude(obs[i].sat,svh[i],opt)) continue;

		/* ionospheric corrections */
		if (!ionocorr(obs[i].time,nav,obs[i].sat,pos,azel+i*2,
					  iter>0?opt->ionoopt:IONOOPT_BRDC,&dion,&vion)) continue;

		/* tropospheric corrections */
		if (!tropcorr(obs[i].time,nav,pos,azel+i*2,
					  iter>0?opt->tropopt:TROPOPT_SAAS,&dtrp,&vtrp)) {
			continue;
		}
		/* pseudorange residual */
		v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);

		/* design matrix */
		for (j=0;j<3;j++) H[j+nv*nxt]=-e[j];
		for (j=3;j<nxt;j++) H[j+nv*nxt]=0.0;

		/* time system and receiver bias offset */

		if(ns[isys-1]==0) nsys++;
		ns[isys-1]++;

		H[IC(0,opt)+nv*(*nx)]=1.0;
		if ((sys==SYS_GLO) && (NC(opt)==2)) {
			v[nv]-=x[IC(1,opt)];
			H[IC(1,opt)+nv*(*nx)]=1.0;
			ncs[1]++;
		}
		else {
			ncs[0]++;
		}

		if((opt->isb==ISBOPT_EST) || (opt->isb==ISBOPT_EST_P) || (opt->isb==ISBOPT_EST_0M)){
			if (ISYSGPS<isys) {
				j=IS(isys,opt);
				v[nv]-=x[IS(isys,opt)];
				H[IS(isys,opt)+nv*nxt]=1.0;
			}
		}

		if((opt->gl2bias==GL2OPT_EST) && (opt->ionoopt==IONOOPT_IFLC)) {
			if(isL2C(obs[i].code[1]) && (isys==ISYSGPS)) {
				gamma=SQR(nav->lam[obs[i].sat-1][1])/SQR(nav->lam[obs[i].sat-1][0]);
				c2=-1.0 /(gamma-1.0);
				v[nv]-=c2*x[I2(opt)];
				H[I2(opt)+nv*nxt]=c2;
			}
		}

        vsat[i]=1; resp[i]=v[nv];
        
        /* error variance */
		var[nv++]=varerr(opt,azel[1+i*2],sys,obs+i,nav)+vare[i]+vmeas+vion+vtrp;
		trace(3,"var:%f,%f,%f,%f,%f\n", varerr(opt,azel[1+i*2],sys,obs+i,nav),vare[i],vmeas,vion,vtrp);

        trace(4,"sat=%2d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n",obs[i].sat,
              azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i],sqrt(var[nv-1]));

        if((opt->gl2bias==GL2OPT_EST) && (opt->ionoopt==IONOOPT_IFLC)) {
            for (j=0;j<nxt;j++) H[j+nv*nxt]=0.0;

			if((isys==ISYSGPS) && (obs[i].dPpc!=0.0)) {
				v[nv]=obs[i].dPpc - (-x[I2(opt)]);
				H[I2(opt)+nv*nxt]=-1;
				var[nv++]=sqrt(2.0)*(varerr(opt,azel[1+i*2],sys,obs+i,nav)+vare[i]+vmeas+vion+vtrp);
            }
        }

	}
	/* shrink design matrix */
	for(isys=0;isys<NC(opt);isys++) {
		if (ncs[isys]>0) continue;
		ic[isys]=-1;
        jrm=IC(isys,opt);
		ix[jrm] = 0;
		for (i=0;i<nv;i++) {
			for (j=0    ;j<jrm;j++) {
                H[j-(i  ) + i*nxt]=H[j+i*nxt];
            }
            for (j=jrm+1;j<nxt;j++) {
                H[j-(i+1) + i*nxt]=H[j+i*nxt];
            }
        }
		nxt--;

		for(i=0;i<NSYS;i++) {
			if(is[i]>=0) is[i]-=1;
		}
		for(i=0;i<1;i++) {
			if(i2[i]>=0) i2[i]-=1;
		}
	}
	trace(2,"H=\n"); tracemat(2,H,*nx,nv,8,3);
	if((opt->isb==ISBOPT_EST) || (opt->isb==ISBOPT_EST_P) || (opt->isb==ISBOPT_EST_0M)) {
		for (isys=ISYSGPS;isys<=NSYS;isys++) {
			if ((ns[isys-1]==0) || (nsys<=1) || (ISYSGPS==isys)) {

				jrm=is[isys-1];
				is[isys-1]=-1;
				for (jsys=isys+1;jsys<=NSYS;jsys++) is[jsys-1]-=1;
				k = IS(isys,opt);
				ix[IS(isys,opt)] = 0;
				for (i=0;i<nv;i++) {
					for (j=    0;j<jrm;j++) {
						H[j-(i  ) + i*nxt]=H[j+i*nxt];
					}
					for (j=jrm+1;j<nxt;j++) {
						H[j-(i+1) + i*nxt]=H[j+i*nxt];
					}
				}
				nxt--;
				for(i=0;i<1;i++) {
					if(i2[i]>=0) i2[i]-=1;
				}
			}
		}
	}

	*nx=nxt;

	return nv;
}

/* validate solution ---------------------------------------------------------*/
static int valsol(const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, int nv, int nx,
                  char *msg)
{
    double azels[MAXOBS*2],dop[4],vv;
    int i,ns;
    
	trace(3,"valsol  : n=%d nv=%d\n",n,nv);

    /* chi-square validation of residuals */
	if(opt->chisqr==CHIOPT_ON) {
        vv=dot(v,v,nv);
        if (nv>nx&&vv>chisqr[nv-nx-1]) {
            sprintf(msg,"chi-square error nv=%d vv=%.1f cs=%.1f",nv,vv,chisqr[nv-nx-1]);
            return 0;
        }
    }
    /* large gdop check */
    for (i=ns=0;i<n;i++) {
        if (!vsat[i]) continue;
        azels[  ns*2]=azel[  i*2];
        azels[1+ns*2]=azel[1+i*2];
        ns++;
    }
    dops(ns,azels,opt->elmin,dop);
    if (dop[0]<=0.0||dop[0]>opt->maxgdop) {
        sprintf(msg,"gdop error nv=%d gdop=%.1f",nv,dop[0]);
        return 0;
	}
    return 1;
}
/* estimate receiver position ------------------------------------------------*/
static int estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
				  const double *vare, const int *svh, const nav_t *nav,
				  const prcopt_t *opt, sol_t *sol, double *azel, int *vsat,
                  double *resp, char *msg)
{
	double *x,*dx,*Q,*v,*H,*var,sig;
	int i,j,k,info,stat,nx0,nx,nv;
	int is[NSYS],ic[2],i2[1];
	int *ix;
	int isys;

	trace(3,"estpos  : n=%d\n",n);
	nx0=NX(opt);

	if((opt->gl2bias!=GL2OPT_EST) || (opt->ionoopt!=IONOOPT_IFLC)) {
		v=mat(n,1); H=mat(nx0,n); var=mat(n,1);
	}
	else {
		v=mat(n*2,1); H=mat(nx0,n*2); var=mat(n*2,1);
	}


	x=zeros(nx0,1); dx=zeros(nx0,1); Q=zeros(nx0,nx0);
	ix=izeros(nx0,1);

	for (i=0;i<3;i++) x[i]=sol->rr[i];

	for (i=0;i<MAXITR;i++) {

		/* pseudorange residuals */
		nv=rescode(i,obs,n,rs,dts,vare,svh,nav,x,opt,v,H,var,azel,vsat,resp,&nx,ic,is,i2,ix);

		if (nv<nx) {
			sprintf(msg,"lack of valid sats ns=%d",nv);
			break;
		}
		/* weight by variance */
		for (j=0;j<nv;j++) {
			sig=sqrt(var[j]);
			v[j]/=sig;
			for (k=0;k<nx;k++) H[k+j*nx]/=sig;
		}
		/* least square estimation */
		if ((info=lsq(H,v,nx,nv,dx,Q))) {
			sprintf(msg,"lsq error info=%d",info);
			break;
		}
		for (j=0,k=0;j<nx0;j++) {
			if(!ix[j]) continue;
			x[j]+=dx[k];
			++k;
		}
		if (norm(dx,nx)<1E-4) {
			if(sol) {
				sol->type=0;
				sol->time=timeadd(obs[0].time,-x[3]/CLIGHT);

				if(ic[0]>0) sol->dtr[0]=x[ic[0]]/CLIGHT; /* receiver clock bias (s) */
				else        sol->dtr[0]=0.0;
				if(ic[1]>0) sol->dtr[1]=sol->dtr[0]-x[ic[1]]/CLIGHT; /* glonass-gps time offset (s) */
				else        sol->dtr[1]=0.0;

				for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
				for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*nx];
				sol->qr[3]=(float)Q[1];    /* cov xy */
				sol->qr[4]=(float)Q[2+nx]; /* cov yz */
				sol->qr[5]=(float)Q[2];    /* cov zx */
				sol->ns=(unsigned char)nv;
				sol->age=sol->ratio=0.0;
				if((opt->isb==ISBOPT_EST) || (opt->isb==ISBOPT_EST_P) || (opt->isb==ISBOPT_EST_0M)) {
					if (((opt->mode<PMODE_KINEMA) || (opt->mode>PMODE_FIXED)) || (opt->isb!=ISBOPT_EST_0M)) {

						for(j=ISYSGPS;j<=NSYS;j++) {
							if(is[j-1]==-1) {
								sol->isb[j-ISYSGPS][0][1]=0.0;
							}
							else {
								sol->isb[j-ISYSGPS][0][1]=x[IS(j,opt)];
							}
						}
					}
				}
				if((opt->gl2bias==GL2OPT_EST) && (opt->ionoopt==IONOOPT_IFLC)) {
					if(i2[0]>=0) sol->gl2[0]=x[I2(opt)];
					else         sol->gl2[0]=0.0;
				}

				/* validate solution */
				if ((stat=valsol(azel,vsat,n,opt,v,nv,nx,msg))) {
					sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
				}
			}
			free(v); free(H); free(var);
			free(x); free(dx); free(Q);
			free(ix);

			return stat;
		}
	}
	if (i>=MAXITR) sprintf(msg,"iteration divergent i=%d, norm=%.4f",i,norm(dx,nx));

	free(v); free(H); free(var);
	free(x); free(dx); free(Q);
	free(ix);
    
	return 0;
}
/* raim fde (failure detection and exclution) -------------------------------*/
static int raim_fde(const obsd_t *obs, int n, const double *rs,
                    const double *dts, const double *vare, const int *svh,
                    const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                    double *azel, int *vsat, double *resp, char *msg)
{
    obsd_t *obs_e;
    sol_t sol_e={{0}};
    char tstr[32],name[16],msg_e[128];
    double *rs_e,*dts_e,*vare_e,*azel_e,*resp_e,rms_e,rms=100.0;
    int i,j,k,nvsat,stat=0,*svh_e,*vsat_e,sat=0;
    
    trace(3,"raim_fde: %s n=%2d\n",time_str(obs[0].time,0),n);
    
    if (!(obs_e=(obsd_t *)malloc(sizeof(obsd_t)*n))) return 0;
    rs_e = mat(6,n); dts_e = mat(2,n); vare_e=mat(1,n); azel_e=zeros(2,n);
    svh_e=imat(1,n); vsat_e=imat(1,n); resp_e=mat(1,n); 
    
    for (i=0;i<n;i++) {
        
        /* satellite exclution */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            obs_e[k]=obs[j];
            matcpy(rs_e +6*k,rs +6*j,6,1);
            matcpy(dts_e+2*k,dts+2*j,2,1);
            vare_e[k]=vare[j];
            svh_e[k++]=svh[j];
        }
        /* estimate receiver position without a satellite */
		if (!estpos(obs_e,n-1,rs_e,dts_e,vare_e,svh_e,nav,opt,&sol_e,azel_e,
                    vsat_e,resp_e,msg_e)) {
            trace(3,"raim_fde: exsat=%2d (%s)\n",obs[i].sat,msg);
            continue;
        }
        for (j=nvsat=0,rms_e=0.0;j<n-1;j++) {
            if (!vsat_e[j]) continue;
            rms_e+=SQR(resp_e[j]);
            nvsat++;
        }
        if (nvsat<5) {
            trace(3,"raim_fde: exsat=%2d lack of satellites nvsat=%2d\n",
                  obs[i].sat,nvsat);
            continue;
        }
        rms_e=sqrt(rms_e/nvsat);
        
        trace(3,"raim_fde: exsat=%2d rms=%8.3f\n",obs[i].sat,rms_e);
        
        if (rms_e>rms) continue;
        
        /* save result */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            matcpy(azel+2*j,azel_e+2*k,2,1);
            vsat[j]=vsat_e[k];
            resp[j]=resp_e[k++];
        }
        stat=1;
        if(sol) *sol=sol_e;
        sat=obs[i].sat;
        rms=rms_e;
        vsat[i]=0;
        strcpy(msg,msg_e);
    }
    if (stat) {
        time2str(obs[0].time,tstr,2); satno2id(sat,name);
        trace(2,"%s: %s excluded by raim\n",tstr+11,name);
    }
    free(obs_e);
    free(rs_e ); free(dts_e ); free(vare_e); free(azel_e);
    free(svh_e); free(vsat_e); free(resp_e);
    return stat;
}
/* doppler residuals ---------------------------------------------------------*/
static int resdop(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const nav_t *nav, const double *rr, const double *x,
                  const double *azel, const int *vsat, double *v, double *H)
{
    double lam,rate,pos[3],E[9],a[3],e[3],vs[3],cosel;
    int i,j,nv=0;

    trace(3,"resdop  : n=%d\n",n);
    
    ecef2pos(rr,pos); xyz2enu(pos,E);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        lam=nav->lam[obs[i].sat-1][0];
        
        if (obs[i].D[0]==0.0||lam==0.0||!vsat[i]||norm(rs+3+i*6,3)<=0.0) {
            continue;
        }
        /* line-of-sight vector in ecef */
        cosel=cos(azel[1+i*2]);
        a[0]=sin(azel[i*2])*cosel;
        a[1]=cos(azel[i*2])*cosel;
        a[2]=sin(azel[1+i*2]);
        matmul("TN",3,1,3,1.0,E,a,0.0,e);
        
        /* satellite velocity relative to receiver in ecef */
        for (j=0;j<3;j++) vs[j]=rs[j+3+i*6]-x[j];
        
        /* range rate with earth rotation correction */
        rate=dot(vs,e,3)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                      rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);
        
        /* doppler residual */
        v[nv]=-lam*obs[i].D[0]-(rate+x[3]-CLIGHT*dts[1+i*2]);
        
        /* design matrix */
        for (j=0;j<4;j++) H[j+nv*4]=j<3?-e[j]:1.0;
        
        nv++;
    }
    return nv;
}
/* estimate receiver velocity ------------------------------------------------*/
static void estvel(const obsd_t *obs, int n, const double *rs, const double *dts,
                   const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                   const double *azel, const int *vsat)
{
    double x[4]={0},dx[4],Q[16],*v,*H;
	int i,j,k,nv;
    
    trace(3,"estvel  : n=%d\n",n);
    
    v=mat(n,1); H=mat(4,n);
    
    for (i=0;i<MAXITR;i++) {
        
		/* doppler residuals */
        if ((nv=resdop(obs,n,rs,dts,nav,opt->rb,x,azel,vsat,v,H))<4) {
            break;
        }
        /* least square estimation */
        if (lsq(H,v,4,nv,dx,Q)) break;
        
        for (j=0;j<4;j++) x[j]+=dx[j];
        
        if (norm(dx,4)<1E-6) {
			if(sol) for (k=0;k<3;k++) sol->rr[k+3]=x[k];
            break;
        }
    }
    free(v); free(H);
}

/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          prcopt_t *opt    I   processing options
*          sol_t  *sol      IO  solution
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
*          ssat_t *ssat     IO  satellite status              (NULL: no output)
*          char   *msg      O   error message for error exit
* return : status(1:ok,0:error)
* notes  : assuming sbas-gps, galileo-gps, qzss-gps, compass-gps time offset and
*          receiver bias are negligible (only involving glonass-gps time offset
*          and receiver bias)
*-----------------------------------------------------------------------------*/
extern int pntpos(const obsd_t *obs, int n, const nav_t *nav,
				  const prcopt_t *opt, sol_t *sol, double *azel, ssat_t *ssat,
                  char *msg)
{
	prcopt_t opt_=*opt;
    double *rs,*dts,*var,*azel_,*resp;
	int i,j,stat,vsat[MAXOBS]={0},svh[MAXOBS];
	int ci=1;

    if(sol) {
        sol->stat=SOLQ_NONE;
    }
    if (n<=0) {strcpy(msg,"no observation data"); return 0;}
    
    trace(3,"pntpos  : tobs=%s n=%d\n",time_str(obs[0].time,3),n);

    if(sol) sol->time=obs[0].time;
	msg[0]='\0';

	if (opt_.mode!=PMODE_SINGLE) { /* for precise positioning */
#if TEST_CMPMODEL == 0
		opt_.sateph =EPHOPT_BRDC;
#else
#if 0
		opt_.sateph =EPHOPT_BRDC;
#endif
#endif
		opt_.ionoopt=IONOOPT_BRDC;
		opt_.tropopt=TROPOPT_SAAS;
//		if(opt_.isb!=ISBOPT_EST_0M) {
//		//	opt_.isb=ISBOPT_EST;
//		}
		ci=2;
	}

	for(j=0;j<ci;++j) {

		if(opt_.mode!=PMODE_SINGLE) {
			if(j==0) opt_.isb=ISBOPT_EST;
			else     opt_.isb=ISBOPT_OFF;
		}

		rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel_=zeros(2,n); resp=mat(1,n);

		/* satellite positons, velocities and clocks */
		satposs(obs[0].time,obs,n,nav,opt_.sateph,rs,dts,var,svh);

		/* estimate receiver position with pseudorange */
		stat=estpos(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);

		/* raim fde */
		if (!stat&&n>=6&&opt->posopt[4]) {
			stat=raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);
		}
		/* estimate receiver velocity with doppler */
		if (stat) estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);

		if (azel) {
			for (i=0;i<n*2;i++) azel[i]=azel_[i];
		}
		if (ssat) {

			for (i=0;i<MAXSAT;i++) {
				ssat[i].vs=0;
				ssat[i].azel[0]=ssat[i].azel[1]=0.0;
				ssat[i].resp[0]=ssat[i].resc[0]=0.0;
				ssat[i].snr[0]=0;
			}
			for (i=0;i<n;i++) {
				ssat[obs[i].sat-1].azel[0]=azel_[  i*2];
				ssat[obs[i].sat-1].azel[1]=azel_[1+i*2];
				ssat[obs[i].sat-1].snr[0]=obs[i].SNR[0];
				if (!vsat[i]) continue;
				ssat[obs[i].sat-1].vs=1;
				ssat[obs[i].sat-1].resp[0]=resp[i];
			}
		}
		free(rs); free(dts); free(var); free(azel_); free(resp);

		if (stat) break;
	}
	return stat;
}

/* output solution status ----------------------------------------------------*/
extern void pntposoutsolstat(rtk_t *rtk, const obsd_t *obs, int level, FILE *fp)
{
    ssat_t *ssat;
    double tow,pos[3],vel[3],acc[3],vela[3]={0},acca[3]={0};
    int i,j,k,week,est,nfreq;
    char id[32];
    int type,sys;
    
    if (level<=0||!fp) return;
    
    trace(3,"outsolstat:\n");
    
    est=rtk->opt.mode>=PMODE_DGPS;

    tow=time2gpst(rtk->sol.time,&week);
    
    /* receiver position */
    fprintf(fp,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
            rtk->sol.stat,rtk->sol.rr[0],rtk->sol.rr[1],rtk->sol.rr[2],
            0.0,0.0,0.0);

    /* receiver velocity and acceleration */
    ecef2pos(rtk->sol.rr,pos);
    ecef2enu(pos,rtk->sol.rr+3,vel);
    fprintf(fp,"$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
            week,tow,rtk->sol.stat,vel[0],vel[1],vel[2],
            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);

    /* receiver clocks */
    fprintf(fp,"$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
			week,tow,rtk->sol.stat,1,rtk->sol.dtr[0]*1E9,
			(rtk->sol.dtr[0]+rtk->sol.dtr[1])*1E9,0.0,0.0);
    
	/* ISB */
	if ((rtk->opt.isb==ISBOPT_EST) || (rtk->opt.isb==ISBOPT_EST_P) || (rtk->opt.isb==ISBOPT_EST_0M)) {
		for (sys=SYS_GPS,i=1;sys<=rtk->opt.navsys;sys<<=1,i++) {
			if(!(rtk->opt.navsys & sys)) continue;
            j=IS(sysind(sys),&rtk->opt);
            fprintf(fp,"$ISB,%d,%.3f,%d,%d,%d,%.4f,%.4f\n",week,tow,rtk->sol.stat,
				  sys,1,rtk->x[j],0.0);
		}
	}
	/* receiver L2P-L2C */
	if (rtk->opt.gl2bias==GL2OPT_EST) {
		j=I2(&rtk->opt);
		fprintf(fp,"$RCVL2B,%d,%.3f,%d,%d,%.4f,%.4f\n",week,tow,rtk->sol.stat,
				1,rtk->sol.gl2[0],0.0);
    }
    if (rtk->sol.stat==SOLQ_NONE||level<=1) return;
    
    /* residuals and status */
    for (i=0,k=0;i<MAXSAT;i++) {
        ssat=rtk->ssat+i;
        if (!ssat->vs) continue;
        while((obs[k].sat!=i+1)) {
            ++k;
        }
		satno2id(i+1,id);
		j=0;
		fprintf(fp,"$SAT,%d,%.3f,%s,%d,%d,%.1f,%.1f,%.4f,%.4f,%d,%.0f,%d,%d,%d,%d,%d,%d\n",
                week,tow,id,j+1,obs[k].code[j],ssat->azel[0]*R2D,ssat->azel[1]*R2D,
				ssat->resp [j],ssat->resc[j],  ssat->vsat[j],ssat->snr[j]*0.25,
                ssat->fix  [j],ssat->slip[j]&3,ssat->lock[j],ssat->outc[j],
                ssat->slipc[j],ssat->rejc[j]);
        ++k;
    }
}