/*------------------------------------------------------------------------------
* ppp.c : precise point positioning
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
* options : -DIERS_MODEL use IERS tide model
*
* references :
*     [1] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*     [2] D.D.McCarthy and G.Petit, IERS Technical Note 32, IERS Conventions
*         2003, November 2003
*     [3] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*         Space Technology Library, 2004
*     [4] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*         May 2009
*     [5] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*         Code Biases, URA
*     [6] MacMillan et al., Atmospheric gradients and the VLBI terrestrial and
*         celestial reference frames, Geophys. Res. Let., 1997
*     [7] G.Petit and B.Luzum (eds), IERS Technical Note No. 36, IERS
*         Conventions (2010), 2010
*
* history : 2015/01/08  1.0  new
*-----------------------------------------------------------------------------*/

#include "rtklib.h"

static const char rcsid[]="$Id:$";

#define AS2R        (D2R/3600.0)    /* arc sec to radian */
#define GME         3.986004415E+14 /* earth gravitational constant */
#define GMS         1.327124E+20    /* sun gravitational constant */
#define GMM         4.902801E+12    /* moon gravitational constant */

                                    /* initial variances */
#define VAR_POS     SQR(100.0)      /*   receiver position (m^2) */
#define VAR_CLK     SQR(100.0)      /*   receiver clock (m^2) */
#define VAR_ZTD     SQR(  0.3)      /*   ztd (m^2) */
#define VAR_GRA     SQR(0.001)      /*   gradient (m^2) */
#define VAR_BIAS    SQR(100.0)      /*   phase-bias (m^2) */
#define VAR_ISBL    SQR(1)/* initial variance of isb(phase) (m^2) */
#define VAR_ISBP    SQR(100.0)/* initial variance of isb(pseudorange) (m^2) */

//#define VAR_GL2     SQR(0.001)      /*   GPS L2P-L2C (m^2) */
#define VAR_GL2     SQR(1)      /*   GPS L2P-L2C (m^2) */

#define VAR_IONO_OFF SQR(10.0)      /* variance of iono-model-off */

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity for saastamoinen model */

#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nfreq)
#define NP(opt)     ((opt)->dynamics?9:3)               /* number of pos solution */
#define NC(opt)     (1+((opt)->tsyscorr>1?1:0))
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))

#define NS0(opt)    ((opt)->isb<ISBOPT_EST?0:(opt)->isb==ISBOPT_EST?(NSYS)*2:(NSYS))
#define NS(opt)     (NS0(opt)*((opt)->ionoopt==IONOOPT_IFLC?2:(opt)->nfreq))

#define N2(opt)     ((opt)->gl2bias!=GL2OPT_EST?0:1)
#define NB(opt)     (MAXSAT*NF(opt))

#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NS(opt)+N2(opt))
#define NX(opt)     (NR(opt)+NB(opt))

#define IC(s,opt)        (NP(opt)+((opt)->tsyscorr>1?(s):0)) /* state index of clocks (s=0:gps,1:glo) */
#define IT(opt)          (NP(opt)+NC(opt))                    /* state index of tropos */
#define IS(sys,type,f,opt) (NP(opt)+NC(opt)+NT(opt)+NS0(opt)*(f)+((opt)->isb==ISBOPT_EST?(NSYS)*(type):0)+(sys-NSYSGPS))   /* isb */
#define I2(opt)          (NP(opt)+NC(opt)+NT(opt)+NS(opt))   /* gl2bias (r:0=rov)         */
#define IB(s,f,opt)      (NR(opt)+MAXSAT*(f)+(s)-1)          /* phase bias (s:satno,f:freq) */


/* function prototypes -------------------------------------------------------*/
#ifdef IERS_MODEL
extern int dehanttideinel_(double *xsta, int *year, int *mon, int *day,
                           double *fhr, double *xsun, double *xmon,
                           double *dxtide);
#endif

static int corr_ion(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var,
                    int *brk);

extern void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
    double rsun[3],esun[3],r,ang,erpv[5]={0},cosa;
	int i,j;
    const char *type;
    
    trace(3,"testeclipse:\n");
    
    /* unit vector of sun direction (ecef) */
    sunmoonpos(gpst2utc(obs[0].time),erpv,rsun,NULL,NULL);
    normv3(rsun,esun);
    
    for (i=0;i<n;i++) {
        type=nav->pcvs[obs[i].sat-1].type;
        
        if ((r=norm(rs+i*6,3))<=0.0) continue;
#if 1
        /* only block IIA */
        if (*type&&!strstr(type,"BLOCK IIA")) continue;
#endif
        /* sun-earth-satellite angle */
        cosa=dot(rs+i*6,esun,3)/r;
        cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
        ang=acos(cosa);
        
        /* test eclipse */
        if (ang<PI/2.0||r*sin(ang)>RE_WGS84) continue;
        
        trace(2,"eclipsing sat excluded %s sat=%2d\n",time_str(obs[0].time,0),
              obs[i].sat);
        
        for (j=0;j<3;j++) rs[j+i*6]=0.0;
    }
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

/*static*/extern int corrmeas(const obsd_t *obs, const nav_t *nav, const double *pos,
                    const double *azel, const prcopt_t *opt,
					const double *dantr, const double *dants, double phw,
					double *meas, double *var, int *brk)
{
	const double *lam=nav->lam[obs->sat-1];
	double ion=0.0,L1,P1,PC,P1_P2,P1_C1,vari,gamma;
	int i,k,f;
	double ydifIdb[NFREQ][2] = {0};
	int nf;
	double c;
	int testflag=1;
	int valnf=0;
	int valf[MAXFREQ]={0};
	nf=NF(opt);

	trace(4,"corrmeas:\n");

	for(i=0;i<nf;++i) {
		meas[i   ]=0.0;
		meas[i+nf]=0.0;
		var[i   ]=0.0;
		var[i+nf]=0.0;
	}

	/* iono-free LC */
	if (opt->ionoopt==IONOOPT_IFLC) {
		return ifmeas(obs,nav,azel,opt,dantr,dants,phw,meas,var,&stas[0]);
	}

	for(k=0;k<nf;k++) {
		f=opt->oprfrq[k];

		if (testflag==1&&obs->L[f]==0.0&&obs->P[f]==0.0) continue;

		if (lam[f]==0.0||obs->L[f]==0.0||obs->P[f]==0.0) return 0;
		if (testsnr(0,0,azel[1],obs->SNR[f]*0.25,&opt->snrmask)) return 0;
		valnf++;
		valf[k]=1;
	}
	if(valnf==0) return 0;

	L1=obs->L[0]*lam[0];
	P1=obs->P[0];

	/* dcb correction */
	gamma=SQR(lam[1]/lam[0]); /* f1^2/f2^2 */
	P1_P2=nav->cbias[obs->sat-1][0];
	P1_C1=nav->cbias[obs->sat-1][1];
	if (P1_P2==0.0&&(satsys(obs->sat,NULL)&(SYS_GPS|SYS_GAL|SYS_QZS))) {
		P1_P2=(1.0-gamma)*gettgd(obs->sat,nav);
	}

	if (obs->code[0]==CODE_L1C) P1+=P1_C1; /* C1->P1 */
	PC=P1-P1_P2/(1.0-gamma);               /* P1->PC */

	/* isb correction */
	chk_isb(satsys(obs->sat,NULL), opt, stas, ydifIdb);

	/* slant ionospheric delay L1 (m) */
	if (!corr_ion(obs->time,nav,obs->sat,pos,azel,opt->ionoopt,&ion,&vari,brk)) {

		trace(2,"iono correction error: time=%s sat=%2d ionoopt=%d\n",
			  time_str(obs->time,2),obs->sat,opt->ionoopt);
		return 0;
	}

	for(k=0;k<nf;k++){
		if(valf[k]==0) continue;
		f=opt->oprfrq[k];
		c=lam[f]*lam[f]/lam[0]/lam[0];
		/* ionosphere and windup corrected phase and code */
		meas[k   ]=obs->L[f]*lam[f]+ion*c-lam[f]*phw;
		meas[k+nf]=obs->P[f]-ion*c;
		if(f==0) {
			meas[k+nf]-=P1_P2/(1.0-gamma);
			if(obs->code[0]==CODE_L1C)  meas[k+nf]+=P1_C1;
		}
		else if(f==1) {
			meas[k+nf]-=P1_P2/(1.0-gamma)*gamma;
		}

		var[k   ]+=vari;
		var[k+nf]+=vari+SQR(ERR_CBIAS);

		/* isb corrected phase and code */
		meas[k   ]-=ydifIdb[f][0];
		meas[k+nf]-=ydifIdb[f][1];

		/* antenna phase center variation correction */
		if(dants) {
			meas[k   ]-=dants[f];
			meas[k+nf]-=dants[f];
		}
		if(dantr) {
			meas[k   ]-=dantr[f];
			meas[k+nf]-=dantr[f];
		}
	}
	trace(2,"meas=%.5f,%.5f,%.5f,%.5f\n",meas[0],meas[0+nf],meas[1],meas[1+nf]);
	return 1;
}

/* output solution status for PPP --------------------------------------------*/
extern void pppoutsolstat(rtk_t *rtk, int level, FILE *fp, const obsd_t *obs)
{
	ssat_t *ssat;
    double tow,pos[3],vel[3],acc[3];
	int i,j,week,nfreq=NF(&rtk->opt),k,f,fr,is;
	char id[32];
    int sys;
    int type;
    
    if (level<=0||!fp) return;
    
    trace(3,"pppoutsolstat:\n");
    
    tow=time2gpst(rtk->sol.time,&week);
    
    /* receiver position */
    fprintf(fp,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
            rtk->sol.stat,rtk->x[0],rtk->x[1],rtk->x[2],0.0,0.0,0.0);
    
    /* receiver velocity and acceleration */
    if (rtk->opt.dynamics) {
        ecef2pos(rtk->sol.rr,pos);
		ecef2enu(pos,rtk->x+3,vel);
        ecef2enu(pos,rtk->x+6,acc);
        fprintf(fp,"$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
                week,tow,rtk->sol.stat,vel[0],vel[1],vel[2],acc[0],acc[1],acc[2],
                0.0,0.0,0.0,0.0,0.0,0.0);
    }
    /* receiver clocks */
	i=IC(0,&rtk->opt);
	if(rtk->opt.tsyscorr>1) {
		fprintf(fp,"$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
				week,tow,rtk->sol.stat,1,rtk->x[i]*1E9/CLIGHT,rtk->x[i+1]*1E9/CLIGHT,
				0.0,0.0);
	}
	else {
		fprintf(fp,"$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
				week,tow,rtk->sol.stat,1,rtk->x[i]*1E9/CLIGHT,0.0,
				0.0,0.0);
	}

    /* tropospheric parameters */
    if (rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG) {
		i=IT(&rtk->opt);
        fprintf(fp,"$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n",week,tow,rtk->sol.stat,
                1,rtk->x[i],0.0);
    }
    /* ISB */
	if ((rtk->opt.isb==ISBOPT_EST) || (rtk->opt.isb==ISBOPT_EST_P) || (rtk->opt.isb==ISBOPT_EST_L) || (rtk->opt.isb==ISBOPT_EST_0M)) {
		for (sys=SYS_GPS,i=1;sys<=rtk->opt.navsys;sys<<=1,i++) {
			if(!(rtk->opt.navsys & sys)) continue;
			for (type=0;type<2;type++) {
				if (type==0 && rtk->opt.isb==ISBOPT_EST_P ) continue;
				if (type==0 && rtk->opt.isb==ISBOPT_EST_0M) continue;
				if (type==1 && rtk->opt.isb==ISBOPT_EST_L ) continue;
				for (j=0;j<nfreq;j++) {
					f=rtk->opt.oprfrq[j];
					if(f<0) continue;
					if(f==2) fr=5;
					else fr=f+1;
					is=IS(i,type,j,&rtk->opt);
					fprintf(fp,"$ISB,%d,%.3f,%d,%d,%d,%d,%.4f\n",week,tow,rtk->sol.stat,
							sys,type,fr,rtk->x[is]);
				}
			}
        }
    }
    /* receiver L2P-L2C */
    if (rtk->opt.gl2bias==GL2OPT_EST) {
        j=I2(&rtk->opt);
        fprintf(fp,"$RCVL2B,%d,%.3f,%d,%.4f\n",week,tow,rtk->sol.stat,rtk->x[j]);
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
		for (j=0;j<nfreq;j++) {
			f=rtk->opt.oprfrq[j];
			if(f<0) continue;
			if(f==2) fr=5;
			else fr=f+1;
			fprintf(fp,"$SAT,%d,%.3f,%s,%d,%d,%.1f,%.1f,%.15f,%.15f,%d,%.0f,%d,%d,%d,%d,%d,%d\n",
					week,tow,id,fr,obs[k].code[f],ssat->azel[0]*R2D,ssat->azel[1]*R2D,
					ssat->resp[f],ssat->resc[f],ssat->vsat[f],ssat->snr[f]*0.25,
					ssat->fix[f],ssat->slip[f]&3,ssat->lock[f],ssat->outc[f],
					ssat->slipc[f],ssat->rejc[f]);
		}
		++k;
	}
}
/* solar/lunar tides (ref [2] 7) ---------------------------------------------*/
static void tide_pl(const double *eu, const double *rp, double GMp,
                    const double *pos, double *dr)
{
    const double H3=0.292,L3=0.015;
    double r,ep[3],latp,lonp,p,K2,K3,a,H2,L2,dp,du,cosp,sinl,cosl;
    int i;
    
    trace(4,"tide_pl : pos=%.3f %.3f\n",pos[0]*R2D,pos[1]*R2D);
    
    if ((r=norm(rp,3))<=0.0) return;
    
    for (i=0;i<3;i++) ep[i]=rp[i]/r;
    
    K2=GMp/GME*SQR(RE_WGS84)*SQR(RE_WGS84)/(r*r*r);
    K3=K2*RE_WGS84/r;
    latp=asin(ep[2]); lonp=atan2(ep[1],ep[0]);
    cosp=cos(latp); sinl=sin(pos[0]); cosl=cos(pos[0]);
    
    /* step1 in phase (degree 2) */
    p=(3.0*sinl*sinl-1.0)/2.0;
    H2=0.6078-0.0006*p;
    L2=0.0847+0.0002*p;
    a=dot(ep,eu,3);
    dp=K2*3.0*L2*a;
    du=K2*(H2*(1.5*a*a-0.5)-3.0*L2*a*a);
    
    /* step1 in phase (degree 3) */
    dp+=K3*L3*(7.5*a*a-1.5);
    du+=K3*(H3*(2.5*a*a*a-1.5*a)-L3*(7.5*a*a-1.5)*a);
    
    /* step1 out-of-phase (only radial) */
    du+=3.0/4.0*0.0025*K2*sin(2.0*latp)*sin(2.0*pos[0])*sin(pos[1]-lonp);
    du+=3.0/4.0*0.0022*K2*cosp*cosp*cosl*cosl*sin(2.0*(pos[1]-lonp));
    
    dr[0]=dp*ep[0]+du*eu[0];
    dr[1]=dp*ep[1]+du*eu[1];
    dr[2]=dp*ep[2]+du*eu[2];
    
    trace(5,"tide_pl : dr=%.3f %.3f %.3f\n",dr[0],dr[1],dr[2]);
}
/* displacement by solid earth tide (ref [2] 7) ------------------------------*/
static void tide_solid(const double *rsun, const double *rmoon,
                       const double *pos, const double *E, double gmst, int opt,
                       double *dr)
{
    double dr1[3],dr2[3],eu[3],du,dn,sinl,sin2l;
    
    trace(3,"tide_solid: pos=%.3f %.3f opt=%d\n",pos[0]*R2D,pos[1]*R2D,opt);
    
    /* step1: time domain */
    eu[0]=E[2]; eu[1]=E[5]; eu[2]=E[8];
    tide_pl(eu,rsun, GMS,pos,dr1);
    tide_pl(eu,rmoon,GMM,pos,dr2);
    
    /* step2: frequency domain, only K1 radial */
    sin2l=sin(2.0*pos[0]);
    du=-0.012*sin2l*sin(gmst+pos[1]);
    
    dr[0]=dr1[0]+dr2[0]+du*E[2];
    dr[1]=dr1[1]+dr2[1]+du*E[5];
    dr[2]=dr1[2]+dr2[2]+du*E[8];
    
    /* eliminate permanent deformation */
    if (opt&8) {
        sinl=sin(pos[0]); 
        du=0.1196*(1.5*sinl*sinl-0.5);
        dn=0.0247*sin2l;
        dr[0]+=du*E[2]+dn*E[1];
        dr[1]+=du*E[5]+dn*E[4];
        dr[2]+=du*E[8]+dn*E[7];
    }
    trace(5,"tide_solid: dr=%.3f %.3f %.3f\n",dr[0],dr[1],dr[2]);
}
/* displacement by ocean tide loading (ref [2] 7) ----------------------------*/
static void tide_oload(gtime_t tut, const double *odisp, double *denu)
{
    const double args[][5]={
        {1.40519E-4, 2.0,-2.0, 0.0, 0.00},  /* M2 */
        {1.45444E-4, 0.0, 0.0, 0.0, 0.00},  /* S2 */
        {1.37880E-4, 2.0,-3.0, 1.0, 0.00},  /* N2 */
        {1.45842E-4, 2.0, 0.0, 0.0, 0.00},  /* K2 */
        {0.72921E-4, 1.0, 0.0, 0.0, 0.25},  /* K1 */
        {0.67598E-4, 1.0,-2.0, 0.0,-0.25},  /* O1 */
        {0.72523E-4,-1.0, 0.0, 0.0,-0.25},  /* P1 */
        {0.64959E-4, 1.0,-3.0, 1.0,-0.25},  /* Q1 */
        {0.53234E-5, 0.0, 2.0, 0.0, 0.00},  /* Mf */
        {0.26392E-5, 0.0, 1.0,-1.0, 0.00},  /* Mm */
        {0.03982E-5, 2.0, 0.0, 0.0, 0.00}   /* Ssa */
    };
    const double ep1975[]={1975,1,1,0,0,0};
    double ep[6],fday,days,t,t2,t3,a[5],ang,dp[3]={0};
    int i,j;

	trace(3,"tide_oload:\n");
    
    /* angular argument: see subroutine arg.f for reference [1] */
    time2epoch(tut,ep);
    fday=ep[3]*3600.0+ep[4]*60.0+ep[5];
    ep[3]=ep[4]=ep[5]=0.0;
    days=timediff(epoch2time(ep),epoch2time(ep1975))/86400.0;
    t=(27392.500528+1.000000035*days)/36525.0;
    t2=t*t; t3=t2*t;
    
    a[0]=fday;
    a[1]=(279.69668+36000.768930485*t+3.03E-4*t2)*D2R; /* H0 */
    a[2]=(270.434358+481267.88314137*t-0.001133*t2+1.9E-6*t3)*D2R; /* S0 */
    a[3]=(334.329653+4069.0340329577*t+0.010325*t2-1.2E-5*t3)*D2R; /* P0 */
    a[4]=2.0*PI;
    
    /* displacements by 11 constituents */
    for (i=0;i<11;i++) {
        ang=0.0;
        for (j=0;j<5;j++) ang+=a[j]*args[i][j];
        for (j=0;j<3;j++) dp[j]+=odisp[j+i*6]*cos(ang-odisp[j+3+i*6]*D2R);
    }
    denu[0]=-dp[1];
    denu[1]=-dp[2];
    denu[2]= dp[0];
    
    trace(5,"tide_oload: denu=%.3f %.3f %.3f\n",denu[0],denu[1],denu[2]);
}
/* iers mean pole (ref [7] eq.7.25) ------------------------------------------*/
static void iers_mean_pole(gtime_t tut, double *xp_bar, double *yp_bar)
{
    const double ep2000[]={2000,1,1,0,0,0};
    double y,y2,y3;
    
    y=timediff(tut,epoch2time(ep2000))/86400.0/365.25;
    
    if (y<3653.0/365.25) { /* until 2010.0 */
        y2=y*y; y3=y2*y;
        *xp_bar= 55.974+1.8243*y+0.18413*y2+0.007024*y3; /* (mas) */
        *yp_bar=346.346+1.7896*y-0.10729*y2-0.000908*y3;
    }
    else { /* after 2010.0 */
        *xp_bar= 23.513+7.6141*y; /* (mas) */
        *yp_bar=358.891-0.6287*y;
    }
}
/* displacement by pole tide (ref [7] eq.7.26) --------------------------------*/
static void tide_pole(gtime_t tut, const double *pos, const double *erpv,
                      double *denu)
{
    double xp_bar,yp_bar,m1,m2,cosl,sinl;
    
    trace(3,"tide_pole: pos=%.3f %.3f\n",pos[0]*R2D,pos[1]*R2D);
    
    /* iers mean pole (mas) */
    iers_mean_pole(tut,&xp_bar,&yp_bar);
    
    m1= erpv[0]/AS2R-xp_bar*1E-3; /* (as) */
    m2=-erpv[1]/AS2R-yp_bar*1E-3;
    
    /* sin(2*theta) = sin(2*phi), cos(2*theta)=-cos(2*phi) */
    cosl=cos(pos[1]);
    sinl=sin(pos[1]);
    denu[0]=  9E-3*sin(pos[0])    *(m1*sinl-m2*cosl); /* de= Slambda (m) */
    denu[1]= -9E-3*cos(2.0*pos[0])*(m1*cosl+m2*sinl); /* dn=-Stheta  (m) */
    denu[2]=-33E-3*sin(2.0*pos[0])*(m1*cosl+m2*sinl); /* du= Sr      (m) */

	trace(5,"tide_pole : denu=%.3f %.3f %.3f\n",denu[0],denu[1],denu[2]);
}
/* tidal displacement ----------------------------------------------------------
* displacements by earth tides
* args   : gtime_t tutc     I   time in utc
*          double *rr       I   site position (ecef) (m)
*          int    opt       I   options (or of the followings)
*                                 1: solid earth tide
*                                 2: ocean tide loading
*                                 4: pole tide
*                                 8: elimate permanent deformation
*          double *erp      I   earth rotation parameters (NULL: not used)
*          double *odisp    I   ocean loading parameters  (NULL: not used)
*                                 odisp[0+i*6]: consituent i amplitude radial(m)
*                                 odisp[1+i*6]: consituent i amplitude west  (m)
*                                 odisp[2+i*6]: consituent i amplitude south (m)
*                                 odisp[3+i*6]: consituent i phase radial  (deg)
*                                 odisp[4+i*6]: consituent i phase west    (deg)
*                                 odisp[5+i*6]: consituent i phase south   (deg)
*                                (i=0:M2,1:S2,2:N2,3:K2,4:K1,5:O1,6:P1,7:Q1,
*                                   8:Mf,9:Mm,10:Ssa)
*          double *dr       O   displacement by earth tides (ecef) (m)
* return : none
* notes  : see ref [1], [2] chap 7
*          see ref [4] 5.2.1, 5.2.2, 5.2.3
*          ver.2.4.0 does not use ocean loading and pole tide corrections
*-----------------------------------------------------------------------------*/
extern void tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp,
                     const double *odisp, double *dr)
{
    gtime_t tut;
    double pos[2],E[9],drt[3],denu[3],rs[3],rm[3],gmst,erpv[5]={0};
    int i;
#ifdef IERS_MODEL
    double ep[6],fhr;
    int year,mon,day;
#endif
    
    trace(3,"tidedisp: tutc=%s\n",time_str(tutc,0));
    
    if (erp) geterp(erp,tutc,erpv);
    
    tut=timeadd(tutc,erpv[2]);
    
    dr[0]=dr[1]=dr[2]=0.0;
    
    if (norm(rr,3)<=0.0) return;
    
    pos[0]=asin(rr[2]/norm(rr,3));
    pos[1]=atan2(rr[1],rr[0]);
    xyz2enu(pos,E);
    
    if (opt&1) { /* solid earth tides */
        
        /* sun and moon position in ecef */
        sunmoonpos(tutc,erpv,rs,rm,&gmst);
        
#ifdef IERS_MODEL
        time2epoch(tutc,ep);
        year=(int)ep[0];
        mon =(int)ep[1];
        day =(int)ep[2];
        fhr =ep[3]+ep[4]/60.0+ep[5]/3600.0;
        
        /* call DEHANTTIDEINEL */
        dehanttideinel_((double *)rr,&year,&mon,&day,&fhr,rs,rm,drt);
#else
        tide_solid(rs,rm,pos,E,gmst,opt,drt);
#endif
        for (i=0;i<3;i++) dr[i]+=drt[i];
    }
    if ((opt&2)&&odisp) { /* ocean tide loading */
        tide_oload(tut,odisp,denu);
        matmul("TN",3,1,3,1.0,E,denu,0.0,drt);
        for (i=0;i<3;i++) dr[i]+=drt[i];
    }
    if ((opt&4)&&erp) { /* pole tide */
        tide_pole(tut,pos,erpv,denu);
        matmul("TN",3,1,3,1.0,E,denu,0.0,drt);
        for (i=0;i<3;i++) dr[i]+=drt[i];
    }
    trace(5,"tidedisp: dr=%.3f %.3f %.3f\n",dr[0],dr[1],dr[2]);
}
/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/

static double varerr(const prcopt_t *opt, double el, int sys, const obsd_t *obs, const int type, const int freq, const nav_t *nav) {
	double fact;
    double a=opt->err[1], b=opt->err[2];
    double c=1.0;

    /*衛星システムから分散荷重係数fの決定*/
    fact = (sys == SYS_GLO)? EFACT_GLO: ((sys == SYS_SBS)? EFACT_SBS: EFACT_GPS);

	trace(4,"errtbl_before(ppp): a = %f,b = %f\n",a,b);
	/*観測誤差テーブルを使用する？*/
	if(opt->errmodel == ERRMODEL_TABLE){
		/* roverで検索 */
		int ret = search_errtbl(opt->rectype[0], opt->anttype[0], sys, obs->code[0], type, nav, &a, &b);
		if (ret == 0) {
			/* baseで検索 */
			search_errtbl(opt->rectype[1], opt->anttype[1], sys, obs->code[0], type, nav, &a, &b);
		}
	} else {
		if (type == 1) {
			/* タイプが疑似距離の場合 */
			c = opt->eratio[0];
		}
	}
	trace(4,"errtbl_after(ppp): a = %f,b = %f\n",a,b);

	/*疑似距離 かつDCBを補正できない？*/
	if ((type == 1)&&(has_dcb(obs->sat, obs->code[0], nav, NULL) == 0)) {
		/*設定値Code Error Ratioでfを荷重する */
		fact *= opt->dcbratio;
	}
	/*電離層フリーモード*/
	if (opt->ionoopt == IONOOPT_IFLC) {
		/* 分散荷重係数 f=3.0*f */
		fact *= 3.0;
	}

	/*疑似距離観測誤差の計算*/
	return SQR(fact)*SQR(c)*(SQR(a) + SQR(b)/SQR(sin(el)));
}

/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    rtk->x[i]=xi;
    for (j=0;j<rtk->nx;j++) {
        rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=i==j?var:0.0;
    }
}
/* dual-frequency iono-free measurements -------------------------------------*/
int ifmeas(const obsd_t *obs, const nav_t *nav, const double *azel,
                  const prcopt_t *opt, const double *dantr, const double *dants,
                  double phw, double *meas, double *var, sta_t* sta)
{
	const double *lam=nav->lam[obs->sat-1];
    double c1,c2,L1,L2,P1,P2,P1_C1,P2_C2,gamma,rP1_C1=0.0,rP2_C2=0.0;
    int i=0,j=1,k;
    int s_code=0,sys,l,res;
    double ydifIdb[NFREQ][2] = {0};
    trace(4,"ifmeas  :\n");
    
    /* L1-L2 for GPS/GLO/QZS, L1-L5 for GAL/SBS */
	if (NFREQ>=3&&(satsys(obs->sat,NULL)&(SYS_GAL|SYS_SBS))) j=2;
    
	if (NFREQ<2||lam[i]==0.0||lam[j]==0.0) return 0;
    
    /* test snr mask */
    if (testsnr(0,i,azel[1],obs->SNR[i]*0.25,&opt->snrmask)||
        testsnr(0,j,azel[1],obs->SNR[j]*0.25,&opt->snrmask)) {
		return 0;
	}
	gamma=SQR(lam[j])/SQR(lam[i]); /* f1^2/f2^2 */
	c1=gamma/(gamma-1.0);  /*  f1^2/(f1^2-f2^2) */
	c2=-1.0 /(gamma-1.0);  /* -f2^2/(f1^2-f2^2) */
    
    L1=obs->L[i]*lam[i]; /* cycle -> m */
	L2=obs->L[j]*lam[j];
	if (isL2C(obs->code[j]) && L2 != 0.0) {
		/* L2C位相波の場合、1/4波長補正チェック
		 *  in case of GPS/GLO/QZS
		 */
		L2 += chk_L2Csft(obs, opt, sta) * lam[j];

    }
	P1=obs->P[i];
    P2=obs->P[j];
	P1_C1=nav->cbias[obs->sat-1][1];
	P2_C2=nav->cbias[obs->sat-1][2];
    
    if (!(sys=satsys(obs->sat,NULL))) return 0;
    switch(sys){
        case SYS_GPS:
            s_code = NSYSGPS;
            break;
        case SYS_GLO:
            s_code = NSYSGPS+NSYSGLO;
			break;
        case SYS_GAL:
            s_code = NSYSGPS+NSYSGLO+NSYSGAL;
            break;
        case SYS_QZS:
            s_code = NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS;
            break;
        case SYS_CMP:
            s_code = NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS+NSYSCMP;
            break;
        default:
            break;
    }
    res=0;

	rP1_C1=sta->dcb[s_code-1][1];
	if(GL2OPT_TABLE==opt->gl2bias) rP2_C2=sta->dcb[s_code-1][2];
    
    if (opt->sateph==EPHOPT_LEX) {
        P1_C1=nav->lexeph[obs->sat-1].isc[0]*CLIGHT; /* ISC_L1C/A */
    }
    if (L1==0.0||L2==0.0||P1==0.0||P2==0.0) return 0;

	/* isb correction */
	chk_isb(sys, opt, sta, ydifIdb);

	/* iono-free phase with windup correction */
	meas[0]=c1*L1+c2*L2-(c1*lam[i]+c2*lam[j])*phw;

	trace(2,"meas[0],%d,%.15f,%.15f,%.15f\n",obs->sat,meas[0],obs->L[j],chk_L2Csft(obs, opt, sta));

    /*電離層フリー結合のDCB補正*/
    /* iono-free code with dcb correction */
    if (obs->code[i]==CODE_L1C) P1+=P1_C1+rP1_C1; /* C1->P1 + 地上局C1->P1*/
    if ((obs->code[j]==CODE_L2C)||(obs->code[j]==CODE_L2S)||(obs->code[j]==CODE_L2L)||(obs->code[j]==CODE_L2X)) P2+=P2_C2+rP2_C2;/* C2->P2 + 地上局C2->P2*/
    meas[1]=c1*P1+c2*P2;
    var[1]=SQR(ERR_CBIAS);
    
    if (opt->sateph==EPHOPT_SBAS) meas[1]-=P1_C1; /* sbas clock based C1 */

    /* isb corrected phase and code */
    meas[0]-=c1*ydifIdb[i][0] + c2*ydifIdb[j][0];
    meas[1]-=c1*ydifIdb[i][1] + c2*ydifIdb[j][1];
    
    /* gps-glonass h/w bias correction for code */
    if (opt->exterr.ena[3]&&satsys(obs->sat,NULL)==SYS_GLO) {
        meas[1]+=c1*opt->exterr.gpsglob[0]+c2*opt->exterr.gpsglob[1];
    }
    /* antenna phase center variation correction */
	for (k=0;k<2;k++) {
		if (dants) meas[k]-=c1*dants[i]+c2*dants[j];
        if (dantr) meas[k]-=c1*dantr[i]+c2*dantr[j];
	}
	return 1;
}

/* slant ionospheric delay ---------------------------------------------------*/
static int corr_ion(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var,
                    int *brk)
{
#ifdef EXTSTEC
    double rate;
#endif
    /* sbas ionosphere model */
    if (ionoopt==IONOOPT_SBAS) {
        return sbsioncorr(time,nav,pos,azel,ion,var);
    }
    /* ionex tec model */
    if (ionoopt==IONOOPT_TEC) {
        return iontec(time,nav,pos,azel,1,ion,var);
    }
#ifdef EXTSTEC
    /* slant tec model */
    if (ionoopt==IONOOPT_STEC) {
        return stec_ion(time,nav,sat,pos,azel,ion,&rate,var,brk);
    }
#endif
    /* broadcast model */
    if (ionoopt==IONOOPT_BRDC) {
	*ion=ionmodel(time,nav->ion_gps,pos,azel);
    *var=SQR(*ion*ERR_BRDCI);
    return 1;
    }
    /* ionosphere model off */
    *ion=0.0;
    *var=VAR_IONO_OFF;
    return 1;
}
/* ionosphere and antenna corrected measurements -----------------------------*/


/* L1/L2 geometry-free phase measurement -------------------------------------*/


/* temporal update of position -----------------------------------------------*/
static void udpos_ppp(rtk_t *rtk)
{
    int i;
    
    trace(3,"udpos_ppp:\n");
    
    /* fixed mode */
    if (rtk->opt.mode==PMODE_PPP_FIXED) {
        for (i=0;i<3;i++) initx(rtk,rtk->opt.ru[i],1E-8,i);
        return;
    }
    /* initialize position for first epoch */
    if (norm(rtk->x,3)<=0.0) {
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
    }
    /* static ppp mode */
    if (rtk->opt.mode==PMODE_PPP_STATIC) return;
    
    /* kinmatic mode without dynamics */
    for (i=0;i<3;i++) {
        initx(rtk,rtk->sol.rr[i],VAR_POS,i);
    }
}

/* temporal update of clock --------------------------------------------------*/
static void udclk_ppp(rtk_t *rtk)
{
	double dtr;
    int i;
    
    trace(3,"udclk_ppp:\n");
    
    /* initialize every epoch for clock (white noise) */
    for (i=0;i<NSYS;i++) {
        if (rtk->opt.sateph==EPHOPT_PREC) {
            /* time of prec ephemeris is based gpst */
            /* negelect receiver inter-system bias  */
			dtr=rtk->sol.dtr[0];
		}
        else {
			dtr=i==0?rtk->sol.dtr[0]:rtk->sol.dtr[0]+rtk->sol.dtr[i];
        }
		initx(rtk,CLIGHT*dtr,VAR_CLK,IC(i,&rtk->opt));
    }
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop_ppp(rtk_t *rtk)
{
    double pos[3],azel[]={0.0,PI/2.0},ztd,var;
    int i=IT(&rtk->opt),j;
    
    trace(3,"udtrop_ppp:\n");
    
    if (rtk->x[i]==0.0) {
        ecef2pos(rtk->sol.rr,pos);
        ztd=sbstropcorr(rtk->sol.time,pos,azel,&var);
        initx(rtk,ztd,var,i);
        
        if (rtk->opt.tropopt>=TROPOPT_ESTG) {
            for (j=0;j<2;j++) initx(rtk,1E-6,VAR_GRA,++i);
        }
    }
    else {
        rtk->P[i*(1+rtk->nx)]+=SQR(rtk->opt.prn[2])*fabs(rtk->tt);
        
        if (rtk->opt.tropopt>=TROPOPT_ESTG) {
            for (j=0;j<2;j++) {
                rtk->P[++i*(1+rtk->nx)]+=SQR(rtk->opt.prn[2]*0.1)*fabs(rtk->tt);
            }
        }
    }
}
/* temporal update of ISB parameters --------------------------------*/
static void udisb_ppp(rtk_t *rtk)
{
	int i,j,f,type;
	int nf;

	trace(3,"udisb_ppp\n");
	nf=NF(&(rtk->opt));
	if (rtk->opt.ionoopt==IONOOPT_IFLC) nf=2;

	for (i=ISYSGPS;i<=NSYS;i++) {
		for (f=0;f<nf;f++) {
			for (type=0;type<2;type++) {
				if (type==0 && rtk->opt.isb==ISBOPT_EST_P ) continue;
				if (type==0 && rtk->opt.isb==ISBOPT_EST_0M) continue;
				if (type==1 && rtk->opt.isb==ISBOPT_EST_L ) continue;
				j=IS(i,type,f,&rtk->opt);

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
/* temporal update of receiver GPS L2P-L2C ------------------------------------*/
static void udgl2bias_ppp(rtk_t *rtk)
{
    int j;
    
    trace(3,"udgl2:\n");

    j=I2(&rtk->opt);
    if (rtk->x[j]==0.0) {
        initx(rtk,1E-6,VAR_GL2,j);
    }
    else {
        rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[7])*rtk->tt;
    }
}
/* detect cycle slip by LLI --------------------------------------------------*/
static void detslp_ll(rtk_t *rtk, const obsd_t *obs, int n)
{
    int i,j;
    
    trace(3,"detslp_ll: n=%d\n",n);
    
    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<rtk->opt.nfreq;j++) {
        if (obs[i].L[j]==0.0||!(obs[i].LLI[j]&3)) continue;
        
        trace(3,"detslp_ll: slip detected sat=%2d f=%d\n",obs[i].sat,j+1);
        
        rtk->ssat[obs[i].sat-1].slip[j]=1;
    }
}
/* detect cycle slip by geometry free phase jump -----------------------------*/
static void detslp_gf(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    double g0,g1;
    int i,j;
    
    trace(3,"detslp_gf: n=%d\n",n);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        if ((g1=gfmeas(obs+i,nav))==0.0) continue;
        
        g0=rtk->ssat[obs[i].sat-1].gf;
        rtk->ssat[obs[i].sat-1].gf=g1;

        trace(4,"detslip_gf: sat=%2d gf0=%8.3f gf1=%8.3f\n",obs[i].sat,g0,g1);
        
        if (g0!=0.0&&fabs(g1-g0)>rtk->opt.thresslip) {
            trace(3,"detslip_gf: slip detected sat=%2d gf=%8.3f->%8.3f\n",
                  obs[i].sat,g0,g1);
            
			for (j=0;j<rtk->opt.nfreq;j++) rtk->ssat[obs[i].sat-1].slip[j]|=1;
		}
	}
}
/* temporal update of phase biases -------------------------------------------*/
static void udbias_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	double offset=0.0,pos[3]={0};
	int i,j,k,f,m,sat,brk=0;
	double *meas;
	double *var;
	double *bias;
	int nf;

	trace(3,"udbias  : n=%d\n",n);

	nf=NF(&rtk->opt);
	meas=zeros(nf*2*n,1);
	var=zeros(nf*2*n,1);

	for (i=0;i<MAXSAT;i++) for (j=0;j<rtk->opt.nfreq;j++) {
		rtk->ssat[i].slip[j]=0;
	}
	/* detect cycle slip by LLI */
	detslp_ll(rtk,obs,n);

	/* detect cycle slip by geometry-free phase jump */
	detslp_gf(rtk,obs,n,nav);


	for (f=0;f<nf;f++) {
		m=rtk->opt.oprfrq[f];


		/* reset phase-bias if expire obs outage counter */
		for (i=0;i<MAXSAT;i++) {
			if (++rtk->ssat[i].outc[m]>(unsigned int)rtk->opt.maxout) {
				initx(rtk,0.0,0.0,IB(i+1,f,&rtk->opt));
			}
		}
	}
		ecef2pos(rtk->sol.rr,pos);

	for (i=k=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat;
		if (!corrmeas(obs+i,nav,pos,rtk->ssat[sat-1].azel,&rtk->opt,NULL,NULL,
					  0.0,meas+nf*2*i,var+nf*2*i,&brk)) {
			continue;
		}
		if (brk) {
			rtk->ssat[sat-1].slip[m]=1;
			trace(2,"%s: sat=%2d correction break\n",time_str(obs[i].time,0),sat);
		}
	}
	for (f=0;f<nf;f++) {
		bias=zeros(n,1);
		m=rtk->opt.oprfrq[f];
		for (i=k=0;i<n&&i<MAXOBS;i++) {
			sat=obs[i].sat;
			j=IB(sat,f,&rtk->opt);
			bias[i]=meas[f+nf*2*i]-meas[f+nf+nf*2*i];
			if (rtk->x[j]==0.0||
				rtk->ssat[sat-1].slip[m]||rtk->ssat[sat-1].slip[1]) continue;
			offset+=bias[i]-rtk->x[j];
			k++;
		}
		/* correct phase-code jump to enssure phase-code coherency */
		if (k>=2&&fabs(offset/k)>0.0005*CLIGHT) {
			for (i=0;i<MAXSAT;i++) {
				j=IB(i+1,f,&rtk->opt);
				if (rtk->x[j]!=0.0) rtk->x[j]+=offset/k;
			}
			trace(2,"phase-code jump corrected: %s n=%2d dt=%12.9fs\n",
				  time_str(rtk->sol.time,0),k,offset/k/CLIGHT);
		}
		for (i=0;i<n&&i<MAXOBS;i++) {
			sat=obs[i].sat;
			j=IB(sat,f,&rtk->opt);

			rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[0])*fabs(rtk->tt);

			if (rtk->x[j]!=0.0&&
				!rtk->ssat[sat-1].slip[m]&&!rtk->ssat[sat-1].slip[1]) continue;

			if (bias[i]==0.0) continue;

			/* reinitialize phase-bias if detecting cycle slip */
			initx(rtk,bias[i],VAR_BIAS,IB(sat,f,&rtk->opt));

		}
		free(bias);
	}

	free(var);
	free(meas);
}

/* temporal update of states --------------------------------------------------*/
static void udstate_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    trace(3,"udstate_ppp: n=%d\n",n);
    
	/* temporal update of position */
    udpos_ppp(rtk);
    
    /* temporal update of clock */
    udclk_ppp(rtk);
    
    /* temporal update of tropospheric parameters */
	if (rtk->opt.tropopt>=TROPOPT_EST) {
        udtrop_ppp(rtk);
    }
    /* temporal update of isb */
	if ((rtk->opt.isb==ISBOPT_EST) || (rtk->opt.isb==ISBOPT_EST_P) || (rtk->opt.isb==ISBOPT_EST_L) || (rtk->opt.isb==ISBOPT_EST_0M)) {
		udisb_ppp(rtk);
    }
    /* temporal update of GPS L2P-L2C bias */
    if (rtk->opt.gl2bias==ISBOPT_EST) {
        udgl2bias_ppp(rtk);
    }
    /* temporal update of phase-bias */
	udbias_ppp(rtk,obs,n,nav);
}
/* satellite antenna phase center variation ----------------------------------*/
/*static*/extern void satantpcv(const double *rs, const double *rr, const pcv_t *pcv,
					  double *dant)
{
    double ru[3],rz[3],eu[3],ez[3],nadir,cosa;
    int i;

    for (i=0;i<3;i++) {
		ru[i]=rr[i]-rs[i];
        rz[i]=-rs[i];
    }
	if (!normv3(ru,eu)||!normv3(rz,ez)) return;

    cosa=dot(eu,ez,3);
    cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
    nadir=acos(cosa);

	antmodel_s(pcv,nadir,dant);
}

/* precise tropospheric model ------------------------------------------------*/
static double prectrop(gtime_t time, const double *pos, const double *azel,
                       const prcopt_t *opt, const double *x, double *dtdx,
                       double *var)
{
    const double zazel[]={0.0,PI/2.0};
    double zhd,m_h,m_w,cotz,grad_n,grad_e;
    
    /* zenith hydrostatic delay */
    zhd=tropmodel(time,pos,zazel,0.0);
    
    /* mapping function */
    m_h=tropmapf(time,pos,azel,&m_w);
    
    if ((opt->tropopt==TROPOPT_ESTG||opt->tropopt==TROPOPT_CORG)&&azel[1]>0.0) {
        
        /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz=1.0/tan(azel[1]);
        grad_n=m_w*cotz*cos(azel[0]);
        grad_e=m_w*cotz*sin(azel[0]);
        m_w+=grad_n*x[1]+grad_e*x[2];
        dtdx[1]=grad_n*(x[0]-zhd);
        dtdx[2]=grad_e*(x[0]-zhd);
    }
    dtdx[0]=m_w;
    *var=SQR(0.01);
    return m_h*zhd+m_w*(x[0]-zhd);
}
/* phase and code residuals --------------------------------------------------*/
static int res_ppp(int iter, obsd_t *obs, int n, const double *rs,
				   const double *dts, const double *vare, const int *svh,
				   const nav_t *nav, const double *x, rtk_t *rtk, double *v,
				   double *H, double *R, double *azel,int *ix)
{
	prcopt_t *opt=&rtk->opt;
    double r,rr[3],disp[3],pos[3],e[3],dtdx[3],dantr[NFREQ]={0};
	double dants[NFREQ]={0},var[MAXOBS*2],dtrp=0.0,vart=0.0/*,varm[2]={0}*/;
	int c,i,j,k,f,l,sat,sys,nv=0,nx=rtk->nx,brk,tideopt;
	int nf=NF(opt);
	gtime_t tsys_epoch;
	gtime_t gtm;
	double	dsec,ep[6];
	char* month;
	int date;
	int isys;
	double t;
	int nsys=NSYS;
	int ic=-1,it[3]={-1,-1,-1},is[2]={-1,-1},i2=-1,ib=-1;
	double* meas;
	double* varm;
	int ff;
	double* lam;
	double gamma;
	double c1,c2;

	meas=zeros(2*nf,1);
	varm=zeros(2*nf,1);
	trace(3,"res_ppp : n=%d nx=%d dts=%f\n",n,nx,dts);

	for (j=0;j<rtk->nx;j++) ix[j]=0;

	for (i=0;i<MAXSAT;i++){
		rtk->ssat[i].vsat[0]=0;
		rtk->ssat[i].sys=satsys(i+1,NULL);
		/*先頭時間をUTCに変換）*/
		if(rtk->ssat[i].sys==SYS_GLO){
			if(opt->tsyscorr==TSYSCORR_CORR){
				gtm = gpst2utc(obs[0].time);
				time2epoch(gtm,ep);
				for(l=0;l<nav->nbipm;l++){
					if(nav->bipm[l].month>0){
						date=nav->bipm[l].month;
						if((ep[0]==nav->bipm[l].year) && (ep[1]==date)&& ((int)ep[2]==nav->bipm[l].day)){
							memcpy(&nav->gps_glo[0],&nav->bipm[l].c0,sizeof(double));
							memcpy(&nav->gps_glo[1],&nav->bipm[l].c1,sizeof(double));
							break;
						}
					}
				}
			}
		}
	}
	for (i=0;i<3;i++) rr[i]=x[i];

    /* earth tides correction */
    if (opt->tidecorr) {
        tideopt=opt->tidecorr==1?1:7; /* 1:solid, 2:solid+otl+pole */
        
        tidedisp(gpst2utc(obs[0].time),rr,tideopt,&nav->erp,opt->odisp[0],
                 disp);
		for (i=0;i<3;i++) rr[i]+=disp[i];
	}
	ecef2pos(rr,pos);

	for (i=0;i<3;i++) ix[i]=1;
	for (i=0;i<NT(opt);i++) ix[IT(opt)]=1;

	for (i=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat;
        if (!(sys=satsys(sat,NULL))||!rtk->ssat[sat-1].vs) continue;
        
        /* geometric distance/azimuth/elevation angle */
        if ((r=geodist(rs+i*6,rr,e))<=0.0||
            satazel(pos,e,azel+i*2)<opt->elmin) continue;

		memcpy(obs[i].azel,azel+i*2,sizeof(double)*2);
        
        /* excluded satellite? */
        if (satexclude(obs[i].sat,svh[i],opt)) continue;
        
        /* tropospheric delay correction */
        if (opt->tropopt==TROPOPT_SAAS) {
            dtrp=tropmodel(obs[i].time,pos,azel+i*2,REL_HUMI);
            vart=SQR(ERR_SAAS);
        }
        else if (opt->tropopt==TROPOPT_SBAS) {
            dtrp=sbstropcorr(obs[i].time,pos,azel+i*2,&vart);
        }
        else if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
            dtrp=prectrop(obs[i].time,pos,azel+i*2,opt,x+IT(opt),dtdx,&vart);
        }
        else if (opt->tropopt==TROPOPT_COR||opt->tropopt==TROPOPT_CORG) {
			dtrp=prectrop(obs[i].time,pos,azel+i*2,opt,x,dtdx,&vart);
        }
		/* satellite antenna model */
		if (opt->posopt[0]) {
			satantpcv(rs+i*6,rr,nav->pcvs+sat-1,dants);
        }
        /* receiver antenna model */
		antmodel(opt->pcvr,satsys(obs[i].sat,NULL),opt->antdel[0],azel+i*2,opt->posopt[1],dantr);

		/* phase windup correction */
		if (opt->posopt[2]) {
			windupcorr(rtk->sol.time,rs+i*6,rr,&rtk->ssat[sat-1].phw);
        }
        /* ionosphere and antenna phase corrected measurements */
		if (!corrmeas(obs+i,nav,pos,azel+i*2,&rtk->opt,dantr,dants,
					  rtk->ssat[sat-1].phw,meas,varm,&brk)) {
            continue;
        }


		if(opt->tsyscorr != TSYSCORR_OFF){

			/*放送暦時計変換パラメータよりGPS時系へ変換*/
			if(sys==SYS_QZS){
				tsys_epoch = gpst2time((int)nav->qzs_gps[2],nav->qzs_gps[3]);
				/* ΔTは先頭の観測データの時刻を用いる
				 * QZSSは衛星クロックパラメータがであり、qzs_gpsゼロが設定されるが、
				 * とりあえず処理を行う。
				 */
				dsec = nav->qzs_gps[0] + nav->qzs_gps[1]*timediff(obs[0].time,tsys_epoch);
				/* dsecを足す			*/
				/* GPST = QZSST + ΔT	*/
				obs[i].time = timeadd(obs[i].time,dsec);
			}
			else if(sys==SYS_GAL){
				tsys_epoch = gst2time((int)nav->gps_gal[2],nav->gps_gal[3]);
				dsec = nav->gps_gal[0] + nav->gps_gal[1]*timediff(obs[0].time,tsys_epoch);
				/* dsecを引く		*/
				/* GPST=GST－ΔT	*/
				obs[i].time = timeadd(obs[i].time,-dsec);
			}
			else if(sys==SYS_GLO){
				if(opt->tsyscorr==TSYSCORR_CORR){
					/*BIPM CircularのパラメータよりGPS時系へ変換*/
					/* gpst-utcについては、rinex.cで補正済み。		*/

					trace(5,"res_ppp bipm before obs n=%d time=%s\n",i,time_str(obs[i].time,9));

					obs[i].time = timeadd(utc2gpst(obs[i].time),
					(-nav->gps_glo[0]+nav->gps_glo[1])*1E-9);
					trace(5,"res_ppp bipm after  obs n=%d time=%s \n",i,time_str(obs[i].time,9));

				}
			}
		}


        /* satellite clock and tropospheric delay */
		r+=-CLIGHT*dts[i*2]+dtrp;
        
		trace(5,"sat=%2d azel=%6.1f %5.1f dtrp=%.3f dantr=%6.3f %6.3f dants=%6.3f %6.3f phw=%6.3f\n",
              sat,azel[i*2]*R2D,azel[1+i*2]*R2D,dtrp,dantr[0],dantr[1],dants[0],
              dants[1],rtk->ssat[sat-1].phw);


		trace(5,"dts[i*2]=%f dtrp=%f\n",dts[i*2],dtrp);


		if (opt->ionoopt==IONOOPT_IFLC) {
			lam=nav->lam[sat-1];
			if (sys!=SYS_GAL) gamma=SQR(lam[1])/SQR(lam[0]); /* f1^2/f2^2 */
			else              gamma=SQR(lam[2])/SQR(lam[1]); /* f1^2/f2^2 */
			c1=gamma/(gamma-1.0);  /*  f1^2/(f1^2-f2^2) */
			c2=-1.0 /(gamma-1.0);  /* -f2^2/(f1^2-f2^2) */
		}

		for (c=0;c<2;c++) { /* for phase and code */
			for (f=0;f<nf;f++) {
				ic=-1;
				it[0]=-1;
				it[1]=-1;
				it[2]=-1;
				is[0]=-1;
				is[1]=-1;
				i2=-1;
				ib=-1;
				ff=opt->oprfrq[f];

				j=c*nf+f;
				if (meas[j]==0.0) continue;
            
				for (k=0;k<nx;k++) H[k+nx*nv]=0.0;

				v[nv]=meas[j]-r;

				for (k=0;k<3;k++) H[k+nx*nv]=-e[k];

				if (sys!=SYS_GLO) ic=IC(0,opt);
				else              ic=IC(1,opt);
				v[nv]-=x[ic];
				H[ic+nx*nv]=1.0;

				if (opt->tropopt>=TROPOPT_EST) {
					for (k=0;k<(opt->tropopt>=TROPOPT_ESTG?3:1);k++) {
						it[k]=IT(opt)+k;
						H[it[k]+nx*nv]=dtdx[k];
					}
				}
				if((c==0) && ((opt->isb==ISBOPT_EST) || (opt->isb==ISBOPT_EST_L))) {
					isys = sysind(sys);

					if (opt->ionoopt!=IONOOPT_IFLC) {
						if(f!=0 || ISYSGPS!=isys) {
							is[0]=IS(isys,1,f,opt);
							v[nv]-= x[is[0]];
							H[is[0]+nx*nv]=1.0;
						}
					}
					else {
						is[0]=IS(isys,1,0,opt);
						is[1]=IS(isys,1,1,opt);
						v[nv]-= c1*x[is[0]]+c2*x[is[1]];
						H[is[0]+nx*nv]=c1;
						H[is[1]+nx*nv]=c2;
					}
				}
				else if((c==1) && ((opt->isb==ISBOPT_EST) || (opt->isb==ISBOPT_EST_P) || (opt->isb==ISBOPT_EST_0M))) {
					isys = sysind(sys);
					if (opt->ionoopt!=IONOOPT_IFLC) {
						if(f!=0 || ISYSGPS!=isys) {
							is[0]=IS(isys,1,f,opt);
							v[nv]-= x[is[0]];
							H[is[0]+nx*nv]=1.0;
						}
					}
					else {
						is[0]=IS(isys,1,0,opt);
						is[1]=IS(isys,1,1,opt);
						v[nv]-= c1*x[is[0]]+c2*x[is[1]];
						H[is[0]+nx*nv]=c1;
						H[is[1]+nx*nv]=c2;
					}
				}
				if ((c==1) && (opt->gl2bias==GL2OPT_EST) && (isL2C(obs[i].code[1]))) {
					i2=I2(opt);
					t = x[i2];
					v[nv]-=x[i2];
					H[i2+nx*nv]=1.0;
				}
				if (c==0) {
					ib=IB(obs[i].sat,f,opt);
					v[nv]-=x[ib];
					H[ib+nx*nv]=1.0;
				}
				/*計画・残差・誤差分散行列の生成*/
				var[nv]=varerr(opt,azel[1+i*2],sys,obs+i,c,ff,nav)+varm[j]+vare[i]+vart;

				if (c==0) rtk->ssat[sat-1].resc[ff]=v[nv]; // carrier phase
				else      rtk->ssat[sat-1].resp[ff]=v[nv]; // pseudorange

				/* test innovation */
#if 0
				if (opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) {
#else
				if (opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno&&sys!=SYS_GLO) {
#endif
					trace(2,"ppp outlier rejected %s sat=%2d type=%d v=%.3f\n",
						  time_str(obs[i].time,0),sat,j,v[nv]);
					rtk->ssat[sat-1].rejc[ff]++;
					continue;
				}
				if (c==0) rtk->ssat[sat-1].vsat[ff]=1;
				if(-1!=ic)    ix[ic]=1;
				if(-1!=it[0])    ix[it[0]]=1;
				if(-1!=it[1])    ix[it[1]]=1;
				if(-1!=it[2])    ix[it[2]]=1;
				if(-1!=is[0]) ix[is[0]]=1;
				if(-1!=is[1]) ix[is[1]]=1;
				if(-1!=i2)    ix[i2]=1;
				if(-1!=ib)    ix[ib]=1;
				nv++;
				if(c==0) trace(2,"v[nv],%d,%.15f,%.15f\n",sat,v[nv],x[ib]);
			}
		}
		if ((opt->gl2bias==GL2OPT_EST) && (obs[i].dPpc!=0.0)) {
			for (k=0;k<nx;k++) H[k+nx*nv]=0.0;
			v[nv]=obs[i].dPpc;
			v[nv]-=-x[I2(opt)];
			H[I2(opt)+nx*nv]=-1;
			var[nv]=varerr(opt,azel[1+i*2],sys,obs+i,1,1,nav);
			rtk->ssat[sat-1].resdpl2=v[nv];
			if(-1!=i2) ix[i2]=1;
			nv++;
		}
    }


    for (i=0;i<nv;i++) for (j=0;j<nv;j++) {
		R[i+j*nv]=i==j?var[i]:0.0;
    }
    trace(5,"x=\n"); tracemat(5,x, 1,nx,8,3);
	trace(5,"v=\n"); tracemat(5,v, 1,nv,8,3);
	trace(5,"H=\n"); tracemat(5,H,nx,nv,8,3);
	trace(5,"R=\n"); tracemat(5,R,nv,nv,8,5);

	trace(2,"Hppp=\n"); tracemat(2,H,nx,nv,8,3);
	trace(2,"Rppp=\n"); tracemat(2,R,nv,nv,8,5);

	free(meas);
	free(varm);

	return nv;
}

extern int pppnf(const prcopt_t *opt)
{
	return NF(opt);
}

extern int pppnx(const prcopt_t *opt)
{
	return NX(opt);
}

extern int pppit(const prcopt_t *opt)
{
	return IT(opt);
}

extern int pppnt(const prcopt_t *opt)
{
	return NT(opt);
}

/* precise point positioning -------------------------------------------------*/
extern void pppos(rtk_t *rtk, obsd_t *obs, int n, const nav_t *nav)
{
    const prcopt_t *opt=&rtk->opt;
	double *rs,*dts,*var,*v,*H,*R,*azel,*xp,*Pp;
	int i,j,f,nv,info,svh[MAXOBS],stat=SOLQ_SINGLE;
	int is0 = IS(1,0,0,&rtk->opt);
	int is1 = IS(1,1,0,&rtk->opt);
	int nc=NC(&rtk->opt);
	int np=NP(&rtk->opt);
	int nf=NF(&rtk->opt);
	int *ux;
	int ff;

	trace(3,"pppos   : nx=%d n=%d\n",rtk->nx,n);

    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel=zeros(2,n);

    for (i=0;i<MAXSAT;i++) rtk->ssat[i].fix[0]=0;
    
    /* temporal update of states */
	udstate_ppp(rtk,obs,n,nav);

    trace(4,"x(0)="); tracemat(4,rtk->x,1,NR(opt),13,4);
    
    /* satellite positions and clocks */
    satposs(obs[0].time,obs,n,nav,rtk->opt.sateph,rs,dts,var,svh);
    
	/* exclude measurements of eclipsing satellite */
    if (rtk->opt.posopt[3]) {
		testeclipse(obs,n,nav,rs);
	}
	xp=mat(rtk->nx,1); Pp=zeros(rtk->nx,rtk->nx);
	ux=izeros(rtk->nx,1);
	matcpy(xp,rtk->x,rtk->nx,1);
	nv=n*(rtk->opt.nfreq*2+TEST_ADDESTPRM); v=mat(nv,1); H=mat(rtk->nx,nv); R=mat(nv,nv);

	for (i=0;i<rtk->opt.niter;i++) {
		/* phase and code residuals */
		if ((nv=res_ppp(i,obs,n,rs,dts,var,svh,nav,xp,rtk,v,H,R,azel,ux))<=0) break;

		/* measurement update */
		matcpy(Pp,rtk->P,rtk->nx,rtk->nx);

		if ((info=filter(xp,Pp,H,v,R,rtk->nx,nv,ux))) {
			trace(2,"ppp filter error %s info=%d\n",time_str(rtk->sol.time,0),
				  info);
			break;
		}
		trace(4,"x(%d)=",i+1); tracemat(4,xp,1,NR(opt),13,4);

		stat=SOLQ_PPP;
	}
	if (stat==SOLQ_PPP) {
		/* postfit residuals */
		res_ppp(1,obs,n,rs,dts,var,svh,nav,xp,rtk,v,H,R,azel,ux);

		matcpy(rtk->x,xp,rtk->nx,1);
		matcpy(rtk->P,Pp,rtk->nx,rtk->nx);
		imatcpy(rtk->ux,ux,rtk->nx,1);
		/* update state and covariance matrix */

        /* ambiguity resolution in ppp */
		if (opt->modepppar==PPPAR_CNES||opt->modepppar==PPPAR_CNES_ILS) {
			if (pppamb(rtk,obs,n,nav,azel)) stat=SOLQ_FIX;
		}
        if (opt->modepppar==PPPAR_FCB) {
			if (pppar(rtk,obs,n,nav)) stat=SOLQ_FIX;
        }
        
        /* update solution status */
        rtk->sol.ns=0;
        rtk->sol.sys=0x00;
        for (i=0;i<n&&i<MAXOBS;i++) {
            if (!rtk->ssat[obs[i].sat-1].vsat[0]) continue;
            rtk->ssat[obs[i].sat-1].lock[0]++;
            rtk->ssat[obs[i].sat-1].outc[0]=0;
			rtk->ssat[obs[i].sat-1].fix [0]=4;
            rtk->sol.ns++;
			rtk->sol.sys |= rtk->ssat[obs[i].sat-1].sys;
			for(f=0;f<NFREQ;f++) {
				if(rtk->ssat[obs[i].sat-1].vsat[f] ==1){
					rtk->sol.sys_prd[f] |= rtk->ssat[obs[i].sat-1].sys;
				}
			}
		}
		rtk->sol.stat=stat;
        
		for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->x[i];
			rtk->sol.qr[i]=(float)rtk->P[i+i*rtk->nx];
        }
        rtk->sol.qr[3]=(float)rtk->P[1];
        rtk->sol.qr[4]=(float)rtk->P[2+rtk->nx];
		rtk->sol.qr[5]=(float)rtk->P[2];
		rtk->sol.dtr[0]=rtk->x[IC(0,opt)]/CLIGHT;
		if(rtk->opt.tsyscorr>1) {
			rtk->sol.dtr[1]=0.0;
		}
		else {
			rtk->sol.dtr[1]=(rtk->x[IC(1,opt)]-rtk->x[IC(0,opt)])/CLIGHT;
		}
		for (i=0;i<n&&i<MAXOBS;i++) {
			rtk->ssat[obs[i].sat-1].snr[0]=MIN(obs[i].SNR[0],obs[i].SNR[1]);
		}
		for (i=0;i<MAXSAT;i++) {
			if (rtk->ssat[i].slip[0]&3) rtk->ssat[i].slipc[0]++;
		}
		for(i=ISYSGPS;i<=NSYS;i++) {
			if((rtk->opt.mode!=PMODE_SINGLE) && (rtk->opt.mode!=PMODE_DGPS)) {
				if((rtk->opt.isb==ISBOPT_EST) || (rtk->opt.isb==ISBOPT_EST_L)) {
					if (rtk->opt.ionoopt!=IONOOPT_IFLC) {
						for(f=0;f<nf;f++) {
							rtk->sol.isb[i-ISYSGPS][f][0]=rtk->x[IS(i,1,f,&rtk->opt)];
						}
					}
					else {
						rtk->sol.isb[i-ISYSGPS][0][0]=rtk->x[IS(i,1,0,&rtk->opt)];
						rtk->sol.isb[i-ISYSGPS][1][0]=rtk->x[IS(i,1,1,&rtk->opt)];
					}
				}
			}
			if((rtk->opt.isb==ISBOPT_EST) || (rtk->opt.isb==ISBOPT_EST_P) || (rtk->opt.isb==ISBOPT_EST_0M)) {
				if (rtk->opt.ionoopt!=IONOOPT_IFLC) {
					for(f=0;f<nf;f++) {
						rtk->sol.isb[i-ISYSGPS][f][1]=rtk->x[IS(i,1,f,&rtk->opt)];
					}
				}
				else {
					rtk->sol.isb[i-ISYSGPS][0][1]=rtk->x[IS(i,1,0,&rtk->opt)];
					rtk->sol.isb[i-ISYSGPS][1][1]=rtk->x[IS(i,1,1,&rtk->opt)];
				}
			}
		}
		if(rtk->opt.gl2bias==GL2OPT_EST) {
			rtk->sol.gl2[0]=rtk->x[I2(&rtk->opt)];
		}
		for(i=0;i<NT(opt);i++) {
			rtk->sol.trop[i]=rtk->x[IT(opt)+i];
		}
	}
	free(rs); free(dts); free(var); free(azel);
	free(xp); free(Pp); free(v); free(H); free(R);
	free(ux);
}

