/*------------------------------------------------------------------------------
* pppar.h : ppp ambiguity resolution (FCB)
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

#include <stdio.h>
#include "rtklib.h"

#define NINCFCB     1024

#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nfreq)
#define NP(opt)     ((opt)->dynamics?9:3)               /* number of pos solution */
#define NC(opt)     (1+((opt)->tsyscorr>1?1:0))
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))

#define NS0(opt)    ((opt)->isb<ISBOPT_EST?0:(opt)->isb==ISBOPT_EST?(NSYS-NSYSGPS)*2:(NSYS-NSYSGPS))
#define NS(opt)     (NS0(opt)*((opt)->ionoopt==IONOOPT_IFLC?2:(opt)->nfreq))

#define N2(opt)     ((opt)->gl2bias!=GL2OPT_EST?0:1)
#define NB(opt)     (MAXSAT*NF(opt))

#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NS(opt)+N2(opt))
#define NX(opt)     (NR(opt)+NB(opt))

#define IC(s,opt)        (NP(opt)+((opt)->tsyscorr>1?(s):0)) /* state index of clocks (s=0:gps,1:glo) */
#define IT(opt)          (NP(opt)+NC(opt))                    /* state index of tropos */
#define IS(sys,type,f,opt) (NP(opt)+NC(opt)+NT(opt)+NS0(opt)*(f)+((opt)->isb==ISBOPT_EST?(NSYS-NSYSGPS)*(type):0)+(sys-NSYSGPS-1))   /* isb */
#define I2(opt)          (NP(opt)+NC(opt)+NT(opt)+NS(opt))   /* gl2bias (r:0=rov)         */
#define IB(s,f,opt)      (NR(opt)+MAXSAT*(f)+(s)-1)          /* phase bias (s:satno,f:freq) */

#define VAR_MW		0.3
#define VAR_FIXSD	0.005

/* get SD variance from rtk->P ------------------------------------------------
*	args   : rtk	I	RTK structure
*	         i		I	bias-1 index
*	         j		I	bias-2 index
*	return : SD variance
*	note   : assuming upper triangular portion of P is filled
*	         SD variance = P(i,i) + P(j,j) - 2*P(i,j)
*----------------------------------------------------------------------------*/
static double sdvar(const rtk_t *rtk, int i, int j)
{
	return rtk->P[i+i*rtk->nx] + rtk->P[j+j*rtk->nx]
		-2.0 * (i<j ? rtk->P[i+j*rtk->nx] : rtk->P[j+i*rtk->nx]);
}

/* FIX解の計算 */
static int fixsol( rtk_t *rtk, nav_t *nav ) {
	int i, j, k, info, ic=0, nf;

	/* [0]変数定義 */
	double *v, *H, *R;
	trace( 3, "fixsol:\n" );

	/* count fixed SD ambiguities */
	for (nf=0,i=0;i<nav->fcb->n;i++) {
		if (nav->fcb->sdd[i].mwfix && nav->fcb->sdd[i].b1fix) nf++;
	}

	/* [1]EKF用変数確保 */
	v=zeros(nf,1); H=zeros(rtk->nx,nf); R=zeros(nf,nf);

	/* [2]FIXEDアンビギュイティによる拘束 */
	for( i = 0; i < nav->fcb->n; i++ ) {
		if ( nav->fcb->sdd[i].mwfix != 1 || nav->fcb->sdd[i].b1fix != 1 ) continue;
		j = IB( nav->fcb->sdd[i].s1, 0, &rtk->opt );
		k = IB( nav->fcb->sdd[i].s2, 0, &rtk->opt );
		v[ic] = nav->fcb->sdd[i].flc - ( rtk->x[j] - rtk->x[k] );
		H[j+ic*nf] = 1.0;
		H[k+ic*nf] = -1.0;
		R[ic+ic*nf] = VAR_FIXSD;
        trace(2,"FX,%6.0f,%2d,%2d,%7.4f,%.4f\n",time2gpst(nav->fcb->sdd[i].te,NULL),
            nav->fcb->sdd[i].s1,nav->fcb->sdd[i].s2,v[ic],sqrt(R[ic+ic*nf]));
		ic ++;
	}

	if ( ic == 0 ) return(0);

    trace(5,"x=\n"); tracemat(5,rtk->x,1,rtk->nx,8,3);
    trace(5,"v=\n"); tracemat(5,v,1,nf,8,3);
    trace(5,"H=\n"); tracemat(5,H,rtk->nx,nf,8,3);
    trace(5,"R=\n"); tracemat(5,R,nf,nf,8,5);

	/* [3]一重差による拘束条件で解を更新 */
	if(( info = filter( rtk->x, rtk->P, H, v, R, rtk->nx, nf, rtk->ux))) {
		trace(1, "filter error( info = %d , ic = %d )\n", info, ic );
		free( v ); free( H ); free( R );
		return( 0 );
	}
	trace( 5, "info=%d", info );

	/* [5]終了処理 */
	free( v ); free( H ); free( R );
	return( 1 );
}

/* FCB入力 */
fcb_t *loadfcb( prcopt_t *popt ) {
	FILE *fp;
	char buf[1024];

	gtime_t ts, te;
	int nwl, nnl, itv;
	int wlcnt = 0, nlcnt = 0;
	double std;
	char cts1[256], cts2[256], cte1[256], cte2[256], cwrk[256];

	/* [0]変数定義 */
	fcb_t *fcb = NULL;
	trace( 3, "loadfcb:\n" );

	/* [1]ファイルオープン */
	if(( fp = fopen( popt->mopt.iffcb, "r" )) == (FILE *)NULL ) {
		showmsg( "error: fcb file open" );
		trace( 1, "error: fcb file open( file = %s )\n", popt->mopt.iffcb );
		return( fcb );
	}

	/* 終端まで読み込み */
	while( fgets( buf, sizeof( buf ), fp )) {

		if( buf[0] == '#' || buf[0] == '\n' ) continue; /* コメント行/空行 */

		/* [2]ヘッダ読み込み */
		if( !fcb ) {
			sscanf( buf, "%s %s %s %s %d %d %d", cts1, cts2, cte1, cte2,
					&nwl, &nnl, &itv );
			cts1[4] = cts1[7] = cte1[4] = cte1[7] = ' ';
			cts2[2] = cts2[5] = cte2[2] = cte2[5] = ' ';
			sprintf( cwrk, "%s %s", cts1, cts2 );
			str2time( cwrk, 0, 26, &ts );
			sprintf( cwrk, "%s %s", cte1, cte2 );
			str2time( cwrk, 0, 26, &te );

			/* [3]FCB構造体の構築 */
			fcb = fcbnew( nnl );
			fcb->ts = ts;
			fcb->te = te;
		}
		else {
			/* [4]WL-FCB読み込み */
			if( wlcnt < nwl ) {
				sscanf( buf, "%d %lf %lf", &fcb->wlfcb[wlcnt].n,
						&fcb->wlfcb[wlcnt].ave, &std );
				fcb->wlfcb[wlcnt].var=SQR(std);
				wlcnt++;
			}

			/* [5]NL-FCB読み込み */
			else if( nlcnt < nnl * nwl ) {
				sscanf( buf, "%s %s %lf %lf", cts1, cts2,
						&fcb->nlfcb[nlcnt].ave, &std );
				fcb->nlfcb[nlcnt].var=SQR(std);
				cts1[4] = cts1[7] = ' ';
				cts2[2] = cts2[5] = ' ';
				sprintf( cwrk, "%s %s", cts1, cts2 );
				str2time( cwrk, 0, 26, &fcb->nlfcb[nlcnt].ts );
				fcb->nlfcb[nlcnt].te = timeadd( fcb->nlfcb[nlcnt].ts, itv );
				nlcnt++;
			}
		}
	}

	/* [6]ファイルクローズ */
	fclose( fp );

	/* [7]終了処理 */
	if( nnl == 0 || nlcnt != nnl * nwl ) {
		showmsg( "error: fcb file read" );
		trace( 1, "error: fcb file read( file = %s, nnl = %d, nlcnt = %d )\n", popt->mopt.iffcb, nnl, nlcnt );
		if( !fcb ) {
			free( fcb->nlfcb );
			free( fcb );
		}
		return( NULL );
	}
	return( fcb );
}

/* FCB推定構造体の初期化 */
fcb_t *fcbnew( int nt ) {
	fcb_t *fcb;
	int i;

	/* [0]変数定義 */
	int m;
	trace( 3, "fcbnew:\n" );

	/* [1]FBC推定構造体の設定 */
	fcb = (fcb_t *)malloc( sizeof( fcb_t ));
	fcb->n = fcb->nmax = 0;
	fcb->sdd = NULL;
	m = NSATGPS * ( NSATGPS - 1 ) / 2;
	for( i = 0; i < m; i++ ) {
		fcb->wlfcb[i].ave = fcb->wlfcb[i].var = 0.0;
	}

	/* [2]NL-FCBプロダクト構造体の設定 */
	fcb->ntfcb = nt;
	fcb->nlfcb = (nlfcb_t *)malloc( sizeof( nlfcb_t ) * nt * m );
	for( i = 0; i < nt * m; i++ ) {
		fcb->nlfcb[i].ave = fcb->nlfcb[i].var = 0.0;
	}

	return( fcb );
}

/* PPP-AR測位 */
int pppar( rtk_t *rtk, obsd_t *obs, int n, nav_t *nav ) {
	int i,k, rc = 0;
	char cwrk1[256], cwrk2[256];
	int nf;
	trace( 3, "pppar: nobs=%d\n", n );

	for( i = 0; i < MAXSAT; i++ ) {
		k=IB(i+1,0,&rtk->opt);
		rtk->pass[i].lcamb = rtk->x[k];
		rtk->pass[i].lcambv = rtk->P[k+k*rtk->nx];
		trace( 4, "i=%d,lcamb=%e,lcambv=%e\n", i, rtk->pass[i].lcamb, rtk->pass[i].lcambv );
	}

	/* [1]パス情報の更新 */
	updpass( rtk, obs, n, nav );
	for( i = 0; i < MAXSAT; i++ ) {
		time2str( rtk->pass[i].ts, cwrk1, 1 );
		time2str( rtk->pass[i].te, cwrk2, 1 );
		trace( 5, "pass[%d]:ts=%s,te=%s,nd=%d,mv=%e,mwv=%e\n", i, cwrk1, cwrk2,
				rtk->pass[i].nd, rtk->pass[i].mw, rtk->pass[i].mwv );
	}

	/* [2]一重差パス構造体の作成 */
	nav->fcb->n = 0;
	addsdd( rtk, obs, n, nav );
	trace(4,"SD sow=%6.0f nobs=%3d nsd=%2d\n",time2gpst(obs[0].time,NULL),n,nav->fcb->n);

	/* [3]ワイドレーンFIX */
	fixwfcb( nav->fcb, &rtk->opt );
	for (nf=i=0;i<nav->fcb->n;i++) if (nav->fcb->sdd[i].mwfix) nf++;
	trace(4,"WF sow=%6.0f fix=%3d/%3d\n",time2gpst(obs[0].time,NULL),nf,nav->fcb->n);

	/* [4]ナローレーンFIX */
	nf=0;
	for( i = 0; i < nav->fcb->n; i++ ) {
		trace( 5, "sdd[%d].exc=%d,mwfix=%d\n", i, nav->fcb->sdd[i].exc, nav->fcb->sdd[i].mwfix );
		if( nav->fcb->sdd[i].exc || !nav->fcb->sdd[i].mwfix ) continue;
        fixnfcb( nav->fcb, &rtk->opt, nav, i );
        if (nav->fcb->sdd[i].b1fix) nf++;
    }
	trace(4,"NF sow=%6.0f fix=%3d/%3d\n",time2gpst(obs[0].time,NULL),nf,nav->fcb->n);

	/* [5]FIX解の計算 */
	rc = fixsol( rtk, nav );

	return( rc );
}

/* パス情報の更新 */
void updpass( rtk_t *rtk, obsd_t *obs, int n, nav_t *nav ) {

	/* [0]変数定義 */
	int i, sat;
	passd_t *p, tmp;
	ssat_t *s;
	double wn, mw, mwv, lcamb;
	unsigned char code[NFREQ];
	trace( 3, "updpass:\n" );

	/* [1]OBSループ */
	for( i = 0; i < n; i++ ) {
		sat = obs[i].sat;
		p = &rtk->pass[sat-1];
		s = &rtk->ssat[sat-1];
		if (!s->vsat[0]) continue;
		if( satsys( sat, NULL ) != SYS_GPS ) continue;
		if( !obs[i].azel[0] || !obs[i].azel[1] ) continue;
		if( obs[i].P[0] < 1.0 || obs[i].P[1] < 1.0 ) continue;

		/* [11]MWデータの計算 */
		trace( 5, "[%d] sat=%d, mw=%e, mwv=%e, P1=%f, P2=%f, C1=%f, C2=%f\n", i, sat, mw, mwv, obs[i].P[0],obs[i].P[1],obs[i].L[0],obs[i].L[1] );

		/* [12]新パスの処理 */
		if( s->slip[0] || s->slip[1] || p->nd == 0 ||
			fabs( timediff( obs[i].time, p->te )) > MIN_ARC_GAP ) {
			p->nd = 0;
			p->wsum = p->mw = p->mwv = 0.0;
			p->ts = obs[i].time;
			trace( 2, "pass reset! pass no=%d, sat=%d, slip[0]=%d, slip[1]=%d\n", i, sat, s->slip[0], s->slip[1] );		}

		/* [13]平均と分散の計算 */
		tmp = *p;
		wn = obs[i].azel[1] >= 30.0 * D2R ? 1.0 : 2.0 * sin( obs[i].azel[1] );
		p->nd++;
		p->wsum += wn;
		p->mw = ( tmp.wsum * tmp.mw + wn * mw) / p->wsum;
		if( p->nd > 1 ) {
            p->mwv=(tmp.wsum*tmp.mwv+wn*(mw-p->mw)*(mw-tmp.mw))/p->wsum;
		}
		else {
			p->mwv = VAR_MW;
		}
		trace(5,"MW,%6.0f,%2d,%.2f,%.2f,%.4f,%.2f,%.2f,%.2f\n",
				time2gpst(obs[i].time,NULL),obs[i].sat,mw,mwv,lcamb,wn,p->mw,p->mwv);

		p->te = obs[i].time;
	}

	return;
}

/* 一重差パス構造体の作成 */
void addsdd( rtk_t *rtk, obsd_t *obs, int n, nav_t *nav ) {

	/* [0]変数定義 */
	int i, j;
	double d;
	gtime_t tsi, tsj;
	sdd_t sdd = {{0}};
	passd_t *p1,*p2;
	int b1,b2;
	trace( 3, "addsdd:\n" );

	for( i = 0; i < n - 1; i++ ) {
		if (!rtk->ssat[obs[i].sat-1].vsat[0]) continue;
		if( satsys( obs[i].sat, NULL ) != SYS_GPS ) continue;
		sdd.s1 = obs[i].sat;
		p1 = rtk->pass + sdd.s1 - 1;
		b1 = IB(sdd.s1,0,&rtk->opt);
		for( j = i + 1; j < n; j++ ) {
			if (!rtk->ssat[obs[j].sat-1].vsat[0]) continue;
			if( satsys( obs[j].sat, NULL ) != SYS_GPS ) continue;

			/* [11]一重差諸量の計算 */
			sdd.s2 = obs[j].sat;
			p2 = rtk->pass + sdd.s2 - 1;
			b2 = IB(sdd.s2,0,&rtk->opt);
			sdd.r = obs[i].rcv;
			d = p1->mw - p2->mw;
			sdd.mwv = p1->mwv/p1->wsum + p2->mwv/p2->wsum;
			sdd.mwi = floor( d );
			sdd.mwf = d - sdd.mwi;
			tsi = p1->ts;
			tsj = p2->ts;
			sdd.ts = timediff( tsi, tsj ) > 0.0 ? tsi : tsj;
			sdd.te = obs[i].time;
			sdd.lc = rtk->x[b1] - rtk->x[b2];
			sdd.lcv = sdvar(rtk,b1,b2);
			trace(5,"S0,%6.0f,%2d,%2d,%7.2f,%5.2f,%7.2f,%5.2f\n",
					time2gpst(obs[0].time,NULL),sdd.s1,sdd.s2,
					d,sqrt(sdd.mwv),sdd.lc,sqrt(sdd.lcv));

			/* [12]一重差データの追加 */
			addsdpass( nav->fcb, &sdd );
		}
	}

	return;
}

/* 一重差パスの追加 */
void addsdpass( fcb_t *fcb, sdd_t *sdd ) {
	sdd_t *sddw;
	/* trace( 4, "addsdpass()\n" ); */

	/* [1]バッファの拡張 */
	if( fcb->nmax <= fcb->n ) {
		if( fcb->nmax <= 0 ) {
			fcb->nmax = NINCFCB;
		}
		else {
			fcb->nmax *= 2;
		}
		if (( sddw = (sdd_t *)realloc( fcb->sdd, sizeof( sdd_t ) * fcb->nmax )) == (sdd_t *)NULL ) {
			free( fcb->sdd );
			fcb->sdd = (sdd_t *)NULL;
			fcb->n = fcb->nmax = 0;
			return;
		}
		fcb->sdd = sddw;
	}

	/* [2]メンバ変数の設定 */
	fcb->sdd[fcb->n].ts = sdd->ts;
	fcb->sdd[fcb->n].te = sdd->te;
	fcb->sdd[fcb->n].s1 = sdd->s1;
	fcb->sdd[fcb->n].s2 = sdd->s2;
	fcb->sdd[fcb->n].r = sdd->r;
	fcb->sdd[fcb->n].ts = sdd->ts;
	fcb->sdd[fcb->n].te = sdd->te;
	fcb->sdd[fcb->n].mwi = sdd->mwi;
	fcb->sdd[fcb->n].mwf = sdd->mwf;
	fcb->sdd[fcb->n].mwv = sdd->mwv;
	fcb->sdd[fcb->n].mwfix = fcb->sdd[fcb->n].b1fix = fcb->sdd[fcb->n].exc = 0;
	fcb->sdd[fcb->n].b1i = fcb->sdd[fcb->n].b1f = fcb->sdd[fcb->n].b1v = 0.0;
	fcb->sdd[fcb->n].lc = sdd->lc;
	fcb->sdd[fcb->n].lcv = sdd->lcv;
	fcb->sdd[fcb->n].flc = fcb->sdd[fcb->n].flcv = 0.0;
	fcb->n++;

	return;
}

/* ワイドレーンアンビギュイティFIX処理 */
void fixwfcb( fcb_t *fcb, prcopt_t *popt ) {

	 /* [0]変数定義 */
	int i, adr;
	double nw, nwfix, nwv, p0;
	trace( 3, "fixwfcb:\n" );

	/* [1]一重差データ整数化判定 */
	for( i = 0; i < fcb->n; i++ ) {

		adr = ( fcb->sdd[i].s1 - 1 ) * NSATGPS - fcb->sdd[i].s1
				* ( fcb->sdd[i].s1 + 1 ) / 2 + fcb->sdd[i].s2 - 1;
		nw = fcb->sdd[i].mwi + fcb->sdd[i].mwf - fcb->wlfcb[adr].ave;
		trace(5,"W0,%6.0f,%2d,%2d,%7.2f,%5.2f,%7.2f\n",time2gpst(fcb->sdd[i].te,NULL),
				fcb->sdd[i].s1,fcb->sdd[i].s2,
				fcb->sdd[i].mwi+fcb->sdd[i].mwf,fcb->wlfcb[adr].ave,nw);
		trace( 5, "i=%d,adr=%d,sdd[%d].s1=%d,sdd[%d].s2=%d,sdd[%d].mwi=%e,sdd[%d].mwf=%e,wlfcb[%d].ave=%e\n",
				i, adr, i, fcb->sdd[i].s1, i, fcb->sdd[i].s2,
				i, fcb->sdd[i].mwi, i, fcb->sdd[i].mwf, adr, fcb->wlfcb[adr].ave );
		nwfix = ROUND( nw );
		nwv = fcb->sdd[i].mwv + fcb->wlfcb[adr].var;
		trace( 5, "nw=%e,nwfix=%e,sdd[%d].mwv=%e,wlfcb[%d].var=%e,nwv=%e\n",
				nw, nwfix, i, fcb->sdd[i].mwv, adr, fcb->wlfcb[adr].var, nwv );

		if( nwv < 1.0e-5 || nwv > 2.0 || fabs(nw-nwfix) > 3.0*sqrt(nwv) ) {
			trace( 2, "warning : nwv=%f i=%d\n", nwv, i );
			fcb->sdd[i].exc = 1;
			continue;
		}
		fcb->sdd[i].mwfix = check_intdd( popt->mopt.minconfw, nw, nwfix, nwv, &p0 );
		trace( 5, "sdd[%d].mwfix=%d, p0=%f\n", i, fcb->sdd[i].mwfix, p0 );
		trace(2,"W1,%6.0f,%2d,%2d,%7.2f,%7.0f,%5.2f,%8.5f,%d\n",
            time2gpst(fcb->sdd[i].te,NULL),fcb->sdd[i].s1,fcb->sdd[i].s2,
            nw-nwfix,nwfix,sqrt(nwv),p0,fcb->sdd[i].mwfix);
	}

	return;
}

/* 一重差のFIXED LCアンビギュイティの算出 */
void fixnfcb( fcb_t *fcb, prcopt_t *popt, nav_t *nav, int i ) {

	/* [0]変数定義 */
	double b1k;
	double fai1, fai1v;
	double n1k, n1ka, n1kv;
	double C1,C2,p0,*lam=nav->lam[fcb->sdd[i].s1-1];
	int j, fain, adr;
	char cwrk1[256], cwrk2[256];
	trace( 4, "fixnfcb:\n" );

	C1=lam[0]*lam[1]/(lam[0]+lam[1]);
	C2=lam[0]/(lam[1]-lam[0]);

	/* [1]b1kの算出 */
	b1k = fcb->sdd[i].lc / C1 - C2 * fcb->sdd[i].mwi;
	trace( 5, "i=%d,sdd[%d].lc=%e,sdd[%d].mwi=%e\n",
			i, i, fcb->sdd[i].lc, i, fcb->sdd[i].mwi );

	/* [2]b1kの整数・小数部分離 */
	fcb->sdd[i].b1i = ROUND( b1k );
	fcb->sdd[i].b1f = b1k - fcb->sdd[i].b1i;

	/* [3]NL-SD-FCBの算出 */
	fai1 = fai1v = 0.0;
	fain = 0;
	time2str( fcb->sdd[i].ts, cwrk1, 1 ); time2str( fcb->sdd[i].te, cwrk2, 1 );
	trace( 5, "sdd[%d].ts=%s,te=%s,b1i=%e,b1f=%e\n",
			i, cwrk1, cwrk2, fcb->sdd[i].b1i, fcb->sdd[i].b1f );
	for( j = 0; j < fcb->ntfcb; j++ ) {
		adr = ( fcb->sdd[i].s1 - 1 ) * NSATGPS - fcb->sdd[i].s1
				* ( fcb->sdd[i].s1 + 1 ) / 2 + fcb->sdd[i].s2 - 1
				+ NSATGPS * ( NSATGPS - 1 ) / 2 * j;
		time2str( fcb->nlfcb[adr].ts, cwrk1, 1 ); time2str( fcb->nlfcb[adr].te, cwrk2, 1 );
		trace( 5, "i=%d,j=%d,adr=%d,sdd[%d].s1=%d,sdd[%d].s2=%d\n",
				i, j, adr, i, fcb->sdd[i].s1, i, fcb->sdd[i].s2 );
		trace( 5, "nlfcb[%d].ts=%s,te=%s,ave=%e,var=%e\n",
				adr, cwrk1, cwrk2, fcb->nlfcb[adr].ave, fcb->nlfcb[adr].var );

		/* テーブル期間がパス期間を過ぎたらbreak */
		if( timediff( fcb->nlfcb[adr].ts, fcb->sdd[i].te ) > 0.0 ||
			 timediff( fcb->nlfcb[adr].te, fcb->sdd[i].ts ) < 0.0 ) continue;

		/* テーブル期間とパス期間の重複期間がテーブル期間/2より長ければ平均処理 */
		if( timediff( timediff( fcb->sdd[i].te, fcb->nlfcb[adr].te ) < 0.0 ?
				fcb->sdd[i].te : fcb->nlfcb[adr].te,
				timediff( fcb->sdd[i].ts, fcb->nlfcb[adr].ts) > 0.0 ?
				fcb->sdd[i].ts : fcb->nlfcb[adr].ts )
				> timediff( fcb->nlfcb[adr].te, fcb->nlfcb[adr].ts ) / 2.0 ) {

			/* 分散値が0なら棄却 */
			if( fabs( fcb->nlfcb[adr].var ) < 1.0e-5 || fabs( fcb->nlfcb[adr].var ) > 2.0 ) continue;

			fai1 += fcb->nlfcb[adr].ave;
			fai1v += fcb->nlfcb[adr].var;
			fain++;
		}
	}

	/* 平均点数が0なら異常で終了 */
	if( !fain ) {
		trace( 2, "warning : fain=%d i=%d\n",fain, i );
		fcb->sdd[i].exc = 1;
		return;
	}
	fai1 /= (double)fain;
	fai1v /= (double)fain;

	/* 整数化判定 */
	n1k = b1k - fai1;
	n1ka = ROUND( n1k );
	n1kv = fcb->sdd[i].lcv / pow( C1, 2.0 ) + fai1v;
	trace(4,"N0,%6.0f,%2d,%2d,%7.2f,%5.2f,%7.2f\n",time2gpst(fcb->sdd[i].te,NULL),
			fcb->sdd[i].s1,fcb->sdd[i].s2,b1k,fai1,n1k);
	trace( 5, "b1k=%e, fai1=%e, n1ka=%e, sdd[%d].lcv=%e, n1kv=%e, fai1v=%e\n",
			b1k, fai1, n1k, n1ka, i, fcb->sdd[i].lcv, n1kv, fai1v );

	if( n1kv < 1.0e-5 || n1kv > 2.0 || fabs(n1k-n1ka) > 1.0*sqrt(n1kv) ) {
		trace( 2, "warning : n1kv=%e i=%d\n", n1kv, i );
		fcb->sdd[i].exc = 1;
		return;
	}

	fcb->sdd[i].b1fix = check_intdd( popt->mopt.minconf1, n1k, n1ka, n1kv, &p0 );
	if( fcb->sdd[i].b1fix ) {
		fcb->sdd[i].flc = C1 * ( n1ka + fai1 + C2 * fcb->sdd[i].mwi);
		fcb->sdd[i].flcv = pow(C1,2.0)*fai1v;
	}
	trace( 5, "sdd[%d].b1fix=%d, p0=%f, flc=%e, flcv=%e\n",
			i, fcb->sdd[i].b1fix, p0, fcb->sdd[i].flc, fcb->sdd[i].flcv );

	return;
}

