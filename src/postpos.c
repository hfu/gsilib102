/*------------------------------------------------------------------------------
* postpos.c : post-processing positioning
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

static const char rcsid[]="$Id: postpos.c,v 1.1 2008/07/17 21:48:06 ttaka Exp $";

/* constants/global variables ------------------------------------------------*/
static sol_t *solf;             /* forward solutions */
static sol_t *solb;             /* backward solutions */
static double *rbf;             /* forward base positions */
static double *rbb;             /* backward base positions */
static int isolf=0;             /* current forward solutions index */
static int isolb=0;             /* current backward solutions index */
static char fsb[MAXSTRPATH]=""; /* 単一基線結果ファイル出力先 */

/*** 変数定義 ***/
static const char outfbase[]="staticfile"; /* 単一基線解析結果ファイル名 */


/* open procssing session ----------------------------------------------------*/
static int openses(const prcopt_t *popt, const solopt_t *sopt,
                   const filopt_t *fopt, nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr)
{
    char *ext;
    int i,j;

    trace(3,"openses :\n");

    for (i=ISYSGPS;i<=NSYS;i++) {
        for (j=1;j<=MAXFREQ;j++) {
            setcodepri( sysno(i), j, popt->codepri[i][j-1]);
        }
    }


	/*L2Cのプライオリティ設定*/
	//setrnxcodepri(1,popt->l2cprior ==0?"CPYWMNDSLX":"SLXCPYWMND");
	if (popt->l2cprior ==1) {
	  setcodepri(SYS_GPS, 2, "SLXCPYWMND");
	}
	else {
	  setcodepri(SYS_GPS, 2, "CPYWMNDSLX");
	}
    /* read satellite antenna parameters */
    if (*fopt->satantp&&!(readpcv(fopt->satantp,pcvs))) {
		showmsg("error : no sat ant pcv in %s",fopt->satantp);
		trace(1,"sat antenna pcv read error: %s\n",fopt->satantp);
        return 0;
    }
    /* read receiver antenna parameters */
    if (*fopt->rcvantp&&!(readpcv(fopt->rcvantp,pcvr))) {
		showmsg("error : no rec ant pcv in %s",fopt->rcvantp);
        trace(1,"rec antenna pcv read error: %s\n",fopt->rcvantp);
        return 0;
    }
    /* read dcb parameters */
    if (*fopt->dcb) {
        readdcb(fopt->dcb,nav);
    }

    /* read isb parameters */
    if (*fopt->isb) {
        readisb(fopt->isb,nav);
    }

	/*1/4波長シフトテーブルを読み込み*/
	if(*popt->mopt.ifpcs){
		if((readL2C(popt->mopt.ifpcs,nav))!=0){
			showmsg("error : no 1/4cycle phase correction file");
			trace(1,"no 1/4cycle phase correction file\n");
			return 0;
		}
	}

	/*GLONASS IFBテーブルを使用？*/
	if(popt->glomodear==GLO_ARMODE_IFB && popt->navsys & SYS_GLO){
		/*IFBテーブルの読み込み*/
		if((readifb(popt->mopt.ififb,nav))!=0){
			showmsg("error : no GLONASS IFB table file");
			trace(1,"no GLONASS IFB table file\n");
			return 0;
		}
	}
	/*観測誤差モデルを使用？*/
	if(popt->errmodel==ERRMODEL_TABLE){
		/*観測誤差モデルの読み込み*/
		if((readerr(popt->mopt.iferr,nav))!=0){
			showmsg("error : no error model file");
			trace(1,"no error model file\n");
        	return 0;
		}
	}

	/*時系補正「Correction」モード */
	if(popt->tsyscorr==TSYSCORR_CORR){
		/*GLONASS時系変換パラメータの読み込み*/
		if((readcirt(fopt->cirtfile,nav))!=0){
			showmsg("error : no BIPM Circular T file");
			trace(1,"no BIPM Circular T file\n");
			return 0;
		}
	}

    /* read ionosphere data file */
    if (*fopt->iono&&(ext=strrchr(fopt->iono,'.'))) {
        if (strlen(ext)==4&&(ext[3]=='i'||ext[3]=='I')) {
            readtec(fopt->iono,nav,0);
        }
#ifdef EXTSTEC
        else if (!strcmp(ext,".stec")||!strcmp(ext,".STEC")) {
            stec_read(fopt->iono,nav);
        }
#endif
    }
    /* open geoid data */
    if (sopt->geoid>0&&*fopt->geoid) {
        if (!opengeoid(sopt->geoid,fopt->geoid)) {
            showmsg("error : no geoid data %s",fopt->geoid);
            trace(2,"no geoid data %s\n",fopt->geoid);
        }
    }
    /* read erp data */
    if (*fopt->eop) {
        if (!readerp(fopt->eop,&nav->erp)) {
            showmsg("error : no erp data %s",fopt->eop);
            trace(2,"no erp data %s\n",fopt->eop);
        }
    }
    return 1;
}

/* close procssing session ---------------------------------------------------*/
static void closeses(nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr)
{
    trace(3,"closeses:\n");

    /* free antenna parameters */
    free(pcvs->pcv); pcvs->pcv=NULL; pcvs->n=pcvs->nmax=0;
    free(pcvr->pcv); pcvr->pcv=NULL; pcvr->n=pcvr->nmax=0;

    /* close geoid data */
    closegeoid();

    /* free erp data */
    free(nav->erp.data); nav->erp.data=NULL; nav->erp.n=nav->erp.nmax=0;

    /* close solution statistics and debug trace */
    rtkclosestat();
//	rtkclosecrinex();
	closernxobsh_ion();
//    rtkclosesnx(NULL);
    traceclose();
}


static void degtodms(double deg, double *dms)
{
	double sgn=1.0;
	if (deg<0.0) {deg=-deg; sgn=-1.0;}
	dms[0]=floor(deg);
	dms[1]=floor((deg-dms[0])*60.0);
	dms[2]=(deg-dms[0]-dms[1]/60.0)*3600;
	dms[0]*=sgn;
}

int addpasstime(passtime_t *pass, gtime_t ts){

	gtime_t *st;
	gtime_t *et;


	trace( 4, "addpasstime()\n" );

	/* passd構造体の増設 */
	if( pass->nmax <= pass->n ) {
		if( pass->nmax <= 0 ) {
			pass->nmax = NUM_PASSD_ADD;
		}
		else {
			pass->nmax *= 2;
		}
		if (( st = (gtime_t *)realloc( pass->st, sizeof( gtime_t ) * pass->nmax ))
			 == (gtime_t *)NULL ) {
			free( pass->st );
			pass->st = (gtime_t *)NULL;
			pass->n = pass->nmax = 0;
			return( -1 );
		}
		pass->st = st;

		if (( et = (gtime_t *)realloc( pass->et, sizeof( gtime_t ) * pass->nmax ))
			 == (gtime_t *)NULL ) {
			free( pass->et );
			pass->et = (gtime_t *)NULL;
			free( pass->st );
			pass->st = (gtime_t *)NULL;
			pass->n = pass->nmax = 0;
			return( -1 );
		}
		pass->et = et;

	}
	pass->st[pass->n]=ts;
	pass->et[pass->n++]=ts;

	return 0;

}

static int setPassData(rtk_t *rtk){

	int i,j;
	obsd_t obs;
	double difftime;

	for (i=0;i<obss.n;i++){
		for(j=0;j<NFREQ;j++){
			obs=obss.data[i];

			if(obs.code[j]==0){
				continue;
			}
            if(obs.rcv>2) continue;

			/* パス情報格納領域確認・確保 */
			if(rtk->passtime[obs.rcv-1][obs.sat-1][j].n==0){
				/* 領域を確保(NUM_PASS_ADD分)、初期化 */
				if(addpasstime(&rtk->passtime[obs.rcv-1][obs.sat-1][j],obs.time)!=0){
					trace(1,"setPassData : malloc error");
					return -1;
				}
				rtk->gsi.ssat[obs.rcv-1][obs.sat-1]=1;
				rtk->gsi.sys|=satsys(obs.sat-1,NULL);

				continue;
			}

			/* 最後のパスデータの取得 */
			/* 最後のパスの終了時刻と、データ時刻の差を計算 */
			difftime = fabs(timediff( obs.time, rtk->passtime[obs.rcv-1][obs.sat-1][j].et[rtk->passtime[obs.rcv-1][obs.sat-1][j].n-1] ));

			if( difftime < MIN_ARC_GAP ){
				/* パス情報保存 */
				/* 最後のパスの終了時刻と、データ時刻の差がパス分割時間間隔(MIN_ARC_GAP)より小さい場合、
				最後のパスの終了時刻＝データ時刻とする */
				rtk->passtime[obs.rcv-1][obs.sat-1][j].et[rtk->passtime[obs.rcv-1][obs.sat-1][j].n-1] = obs.time;

			}else{
				/* 最後のパスの終了時刻と、データ時刻の差がパス分割時間間隔(MIN_ARC_GAP)以上の場合、
				次の配列に開始終了時刻を設定する。*/
				if(addpasstime(&rtk->passtime[obs.rcv-1][obs.sat-1][j],obs.time)!=0){
					trace(1,"setPassData : malloc error");
					return -1;
				}
				continue;
			}
		}
	}

	return 0;
}

static int makeTrackingStatus(int sta,char* path_str,const prcopt_t *popt, rtk_t *rtk){
		/* パス情報出力 */
		/* csvファイルオープン */
	char path[1024]={0};
	FILE *fp;
	int i,j,k,prn;
	struct tm *time_inf;
	double durationtime;

		strcpy( path, path_str);
		if(sta==0){
			strcat( path, ".rover_FieldbookStatic.out.TrackingStatus.csv" );
		}else{
			strcat( path, ".base_FieldbookStatic.out.TrackingStatus.csv" );
		}

		if(( fp = fopen( path, "w" )) == NULL ) {
	/*
			showmsg("error : cannot open solution output file");
			trace(1,"cannot open solution output file\n");
	*/
			return( 0 );
		}

		/* パス情報ソート */
		for(i=0;i<MAXSAT;i++){
			for(j=0;j<NFREQ;j++){
				/*if(rtk->gsi.ssat[i]==0)continue;*/
				if(rtk->passtime[sta][i][j].st==NULL)continue;
				if(rtk->passtime[sta][i][j].et==NULL)continue;
				for(k=0;k<rtk->passtime[sta][i][j].n;k++){

					/* 衛星種別 */
					if( rtk->ssat[i].sys == SYS_GPS ) {
						fprintf( fp, "GPS");
					}else if( rtk->ssat[i].sys == SYS_GLO ) {
						fprintf( fp, "GLO");
					}else if( rtk->ssat[i].sys == SYS_QZS ) {
						fprintf( fp, "QZS");
					}else if( rtk->ssat[i].sys == SYS_GAL ) {
						fprintf( fp, "GAL");
					}
					satsys(i+1,&prn);
					/* 衛星番号 */
					/*fprintf( fp, ",%02d",satno(rtk->ssat[i].sys,i)+1);*/
					fprintf( fp, ",%02d",prn);

					/* 周波数帯 */
					switch( j ){
						case 0:
 #if 0
							if( popt->ionoopt != IONOOPT_IFLC ) {
								fprintf( fp, ",L1" );
							}
							else {
								fprintf( fp, ",LC" );
							}
#endif
							fprintf( fp, ",L1" );

							break;
						case 1:
							fprintf( fp, ",L2" );
							break;
						case 2:
							fprintf( fp, ",L5" );
							break;
						default:
							break;
					}

/*					if(!(rtk->passtime[i][j].st[k].time)&&!(rtk->passtime[i][j].st[k].sec)&&!(rtk->passtime[i][j].et[k].time)&&!(rtk->passtime[i][j].et[k].sec)){  */
					if((rtk->passtime[sta][i][j].st[k].time)&&(rtk->passtime[sta][i][j].et[k].time)){
						time_inf = gmtime(&(rtk->passtime[sta][i][j].st[k].time));
						fprintf( fp, ",%04d/%02d/%02d %02d:%02d:%02d", time_inf->tm_year + 1900, time_inf->tm_mon+1, time_inf->tm_mday, time_inf->tm_hour, time_inf->tm_min, time_inf->tm_sec );

						durationtime = timediff(rtk->passtime[sta][i][j].et[k], rtk->passtime[sta][i][j].st[k]);
						fprintf( fp, ",%d\n",(int)floor(durationtime));
					}else{
					fprintf( fp, ",NULL\n");
					}
				}
			}
		}

		/* cvsファイルクローズ */
		fclose( fp );
		return 0;

}

static double calc_baseline(rtk_t *rtk){

	double a,fe,L,delta,sigma,sigma_d,delta_d,xi,xi_d,eta,eta_d,u1,u2,x,y,c,epsilon;
	double theta,g,h,sigma_s,J,K,gamma_s,gamma,zeta,zeta_d,D,E,F,G;
	int i;

	double n0,A,B,D12;


	double blh1[3], blh2[3];
	/* 起点 */
	cov2ecef(rtk->gsi.rb,NULL,blh1,NULL,RE_GRS80,FE_GRS80,XYZ_TO_PLH);
	/* 終点 */
	cov2ecef(rtk->gsi.rr,NULL,blh2,NULL,RE_GRS80,FE_GRS80,XYZ_TO_PLH);

	a = RE_GRS80;
	fe = FE_GRS80;


		L = blh2[1] - blh1[1];
		if(L<-PI) L = 2.0*PI+L;
		if(L>PI) L = L-2.0*PI;
		delta = blh2[0] - blh1[0];
		sigma = blh1[0] + blh2[0];
		sigma_d = atan(((1.0-fe)*sin(sigma))/(cos(sigma)+fe*(2.0-fe)*sin(blh1[0])*sin(blh2[0])));
		delta_d = atan(((1.0-fe)*sin(delta))/(cos(delta)-fe*(2.0-fe)*sin(blh1[0])*sin(blh2[0])));
		xi = cos(sigma_d/2.0);
		xi_d = sin(sigma_d/2.0);
		eta = sin(delta_d/2.0);
		eta_d = cos(delta_d/2.0);
		u1 = atan((1.0-fe)*tan(blh1[0]));
		u2 = atan((1.0-fe)*tan(blh2[0]));
		x = sin(u1)*sin(u2);
		y = cos(u1)*cos(u2);
		c = y*cos(L)+x;
		epsilon = fe * (2.0-fe) / ((1.0-fe)*(1.0-fe));

		trace(5,"u1=%f u2=%f \n", u1,u2);

		trace(5,"x=%f y=%f c=%f epsilon=%f \n", x,y,c,epsilon);


		if(c<0){
			//線長計算不可
			return 0.0;
		}else{
			theta = L*(1+fe*y);

			for(i=0;1;i++){
				g = SQRT((eta*eta*cos(theta/2.0)*cos(theta/2.0))+(xi*xi*sin(theta/2.0)*sin(theta/2.0)));
				h = SQRT((eta_d*eta_d*cos(theta/2.0)*cos(theta/2.0))+(xi_d*xi_d*sin(theta/2.0)*sin(theta/2.0)));
				sigma_s = 2.0*atan(g/h);
				J = 2.0*g*h;
				if(J ==0) break;
				K = (h*h)-(g*g);
				gamma_s = y * sin(theta) / J;
				gamma = 1-(gamma_s*gamma_s);
				zeta = (gamma*K)-(2.0*x);
				zeta_d = zeta+x;
				D = (fe*(1.0+fe)/4.0)-(3.0/16.0)*fe*fe*gamma;
				E = (1.0-D*gamma)*fe*gamma_s*(sigma_s+D*J*(zeta+D*K*(2.0*zeta*zeta-gamma*gamma)));
				F = theta-L-E;
				G = fe*gamma_s*gamma_s*(1.0-2.0*D*gamma)+fe*zeta_d*(sigma_s/J)*(1.0-D*gamma+fe*gamma_s*gamma_s/2.0)+fe*fe*zeta*zeta_d/4.0;

				theta = theta-F/(1.0-G);
				if((fabs(F)<1.e-15)||(10000<i)){
					break;
				}
			}
		}
		n0 = (epsilon*gamma)/((SQRT(1.0+(epsilon*gamma))+1.0)*(SQRT(1.0+(epsilon*gamma))+1.0));
		A = (1.0+n0)*(1.0+(5.0/4.0)*n0*n0);
		B = (epsilon*(1.0-(3.0/8.0)*n0*n0))/((SQRT(1.0+(epsilon*gamma))+1.0)*(SQRT(1.0+(epsilon*gamma))+1.0));
		D12 = (1.0-fe)*a*A*(sigma_s-B*J*(zeta-B/4.0*(K*(gamma*gamma-2.0*zeta*zeta)-(B/6.0)*zeta*(1.0-4.0*K*K)*(3.0*gamma*gamma-4.0*zeta*zeta))));
	return D12;
}


static int make_single_report(int sta,char* path_str,char* file_str,const prcopt_t *popt, rtk_t *rtk)
{
	//手簿ファイルの作成

    int us;
	char path[1024]={0};
	FILE *fp;
	double blh[3], enu[3];
	gtime_t gtm;
	struct tm *time_inf;
	int i,j,prn;
    int nf=popt->ionoopt==IONOOPT_IFLC?1:popt->nfreq;
    int f;
    int minsat;
    double bline;
    double satflg[3];  
    
	if( popt->mode == PMODE_STATIC ) {
		/* 単一基線解析結果ファイル：観測手簿・スタティック用のファイルオープン */
		strcpy( path, path_str);

		if(sta==0){
			strcat( path, ".rover_FieldbookStatic.out" );
		} else{
			strcat( path, ".base_FieldbookStatic.out" );

		}

		if(( fp = fopen( path, "w" )) == NULL ) {
/*
			showmsg("error : cannot open solution output file");
			trace(1,"cannot open solution output file\n");
*/
			return( 0 );
		}

        /* 情報出力 */
		fprintf( fp, "CD_FileType=%s\n", "0" ); /* ファイルタイプ */
		fprintf( fp, "FS_Title=%s\n", "GNSS測量手簿" ); /* タイトル */
		fprintf( fp, "FS_ObservationPoint=\n"); /* 地区名 */

		if( popt->rectype[sta] == '\0' ){
			fprintf( fp, "FS_ReceiverName=%s\n", stas[sta].rectype ); /* 受信機名称 */
		}
		else{
			fprintf( fp, "FS_ReceiverName=%s\n", popt->rectype[sta] ); /* 受信機名称 */
		}

		fprintf( fp, "FS_ReceiverSerial=%s\n", stas[sta].recsno ); /* 受信機シリアル */

		if( strchr( popt->anttype[sta], '*' ) != NULL ){
			fprintf( fp, "FS_AntennaName=%s\n", stas[sta].antdes ); /* アンテナ名称 */
		}
		else{
			fprintf( fp, "FS_AntennaName=%s\n", popt->anttype[sta] ); /* アンテナ名称 */
		}

		fprintf( fp, "FS_AntennaSerial=%s\n", stas[sta].antsno ); /* アンテナシリアル */

		switch( popt->rovpos ){
			case 0:
				fprintf( fp, "FS_AntennaHeight=%f\n", popt->antdel[sta][2] ); /* アンテナ底面高 */
				break;
			case 1:
			case 2:
				fprintf( fp, "FS_AntennaHeight=%f\n", 0.0 ); /* アンテナ底面高 */
				break;
			case 3:
				switch( stas[sta].deltype ){
					case 0:
						fprintf( fp, "FS_AntennaHeight=%f\n", stas[sta].del[2] ); /* アンテナ底面高 */
						break;
					case 1:
						ecef2pos(stas[sta].pos,blh);
						ecef2enu(blh,stas[0].del,enu);
						fprintf( fp, "FS_AntennaHeight=%f\n", enu[2] ); /* アンテナ底面高 */
						break;
					default:
						break;
				}
				break;
			default:
				break;
		}

		fprintf( fp, "FS_SessionName=\n"); /* セッション名 */

		gtm = gpst2utc(obss.data[0].time);
		time_inf = gmtime(&(gtm.time));
		fprintf( fp, "FS_ObservationStart=%04d/%02d/%02d %02d:%02d:%02d UTC\n", time_inf->tm_year + 1900, time_inf->tm_mon+1, time_inf->tm_mday, time_inf->tm_hour, time_inf->tm_min, time_inf->tm_sec ); /* 開始 */

		gtm = gpst2utc(obss.data[obss.n-1].time);
		time_inf = gmtime(&(gtm.time));
		fprintf( fp, "FS_ObservationEnd=%04d/%02d/%02d %02d:%02d:%02d UTC\n", time_inf->tm_year + 1900, time_inf->tm_mon+1, time_inf->tm_mday, time_inf->tm_hour, time_inf->tm_min, time_inf->tm_sec ); /* 終了 */

		fprintf( fp, "FS_ObservationInterval=%d\n", rtk->gsi.tint ); /* データ取得間隔 */
		fprintf( fp, "FS_ElevationAngle=%f\n", popt->elmin*R2D ); /* 最低高度角 */
//		fprintf( fp, "FS_MinimumSatellites=%d\n", rtk->gsi.minsat ); /* 最少衛星数 */
        
        bline=calc_baseline(rtk);

		//最小衛星数が設定ファイルに存在しない場合
		if(popt->minsat==0){
			satflg[0]=0;
			satflg[1]=0;
			satflg[2]=0;

			for(i=0;i<MAXSAT;i++){
				if(rtk->gsi.ssat[sta][i]==0)continue;
				satsys(i+1,&prn);
				if( rtk->ssat[i].sys == SYS_GPS ) {
					satflg[0]=1;
				}else if( rtk->ssat[i].sys == SYS_GLO ) {
					satflg[1]=1;
				}else if( rtk->ssat[i].sys == SYS_QZS ) {
					satflg[2]=1;
				}
			}
			//GPS+GLO+QZSまたはGPS+GLOの場合
			if((satflg[0]==1 && satflg[1]==1 && satflg[2]==1) ||
			(satflg[0]==1 && satflg[1]==1 && satflg[2]==0)) {
				if(bline >NASELINE_LENGTH){
					minsat=MIN_SATELLITE_10K_ALL;
				}else{
					minsat=MIN_SATELLITE_ALL;

				}


			}else if( (satflg[0]==1 && satflg[1]==0 && satflg[2]==0) ||
					  (satflg[0]==1 && satflg[1]==0 && satflg[2]==1) ||
					  (satflg[0]==0 && satflg[1]==1 && satflg[2]==0) ) {
				if(bline >NASELINE_LENGTH){
					minsat=MIN_SATELLITE_10K;
				}else{
					minsat=MIN_SATELLITE;

				}
			}else{
				minsat=MIN_SATELLITE;
			}

		}else{
			minsat=popt->minsat;
        }

		fprintf( fp, "FS_MinimumSatellites=%d\n", popt->minsat ); /* 最少衛星数 */


		fprintf( fp, "FS_ReceivingStatus=%s%s\n", file_str, ".FieldbookStatic.out.TrackingStatus.csv" ); /* 電波の受信状況 */

		if(sta==0){
			fprintf( fp, "FS_ReceivingStatus=%s%s\n", file_str, ".rover_FieldbookStatic.out.TrackingStatus.csv" ); /* 電波の受信状況 */
		} else{
			fprintf( fp, "FS_ReceivingStatus=%s%s\n", file_str, ".base_FieldbookStatic.out.TrackingStatus.csv" ); /* 電波の受信状況 */

		}



		for(i=0;i<MAXSAT;i++){
			if(rtk->gsi.ssat[sta][i]==0)continue;
			satsys(i+1,&prn);
			if( rtk->ssat[i].sys == SYS_GPS ) {
				fprintf( fp, "FS_StatusGPS[%d]=OK\n", prn ); /* 衛星の状態(GPS) */
			}else if( rtk->ssat[i].sys == SYS_GLO ) {
				fprintf( fp, "FS_StatusGLONASS[%d]=OK\n", prn ); /* 衛星の状態(GLONASS) */
			}else if( rtk->ssat[i].sys == SYS_QZS ) {
				prn=(prn-MINPRNQZS)+1;
				fprintf( fp, "FS_StatusQZSS[%d]=OK\n", prn ); /* 衛星の状態(QZSS) */
			}else if( rtk->ssat[i].sys == SYS_GAL ) {
				fprintf( fp, "FS_StatusGalileo[%d]=OK\n", prn ); /* 衛星の状態(Galileo) */
			}
		}


 #if 0
		for( i = MINPRNGPS; i < MAXPRNGPS; i++ ) {
			if( rtk->gsi.ssat[satno(SYS_GPS,i)] ){
				fprintf( fp, "FS_StatusGPS[%d]=OK\n", i ); /* 衛星の状態(GPS) */
			}
			else {
				fprintf( fp, "FS_StatusGPS[%d]=NG\n", i ); /* 衛星の状態(GPS) */
			}
		}

		for( i = MINPRNGLO; i < MAXPRNGLO; i++ ) {
			if( rtk->gsi.ssat[satno(SYS_GLO,i)] ){
				fprintf( fp, "FS_StatusGLONASS[%d]=OK\n", i ); /* 衛星の状態(GLONASS) */
			}
			else {
				fprintf( fp, "FS_StatusGLONASS[%d]=NG\n", i ); /* 衛星の状態(GLONASS) */
			}
		}

		for( i = MINPRNQZS; i < MAXPRNQZS; i++ ) {
			if( rtk->gsi.ssat[satno(SYS_QZS,i)] ){
				fprintf( fp, "FS_StatusQZSS[%d]=OK\n", i ); /* 衛星の状態(QZSS) */
			}
			else {
				fprintf( fp, "FS_StatusQZSS[%d]=NG\n", i ); /* 衛星の状態(QZSS) */
			}
		}

		for( i = MINPRNGAL; i < MAXPRNGAL; i++ ) {
			if( rtk->gsi.ssat[satno(SYS_GAL,i)] ){
				fprintf( fp, "FS_StatusGalileo[%d]=OK\n", i ); /* 衛星の状態(Galileo) */
			}
			else {
				fprintf( fp, "FS_StatusGalileo[%d]=NG\n", i ); /* 衛星の状態(Galileo) */
			}
		}
 #endif
		/* 単一基線解析結果ファイル：観測手簿・スタティック用のファイルクローズ */
		fclose( fp );
    }
	else if( popt->mode == PMODE_KINEMA ) {
        /* 単一基線解析結果ファイル：観測手簿・キネマティック用のファイルオープン */
		strcpy( path, path_str);

        strcat( path, ".FieldbookKinema.out" );

		if(( fp = fopen( path, "w" )) == NULL ) {
/*
			showmsg("error : cannot open solution output file");
			trace(1,"cannot open solution output file\n");
*/
			return( 0 );
		}

        /* 情報出力 */
		gtm = gpst2utc(obss.data[0].time);
		time_inf = gmtime(&(gtm.time));
		fprintf( fp, "CD_FileType=%s\n", "2" ); /* ファイルタイプ */
		fprintf( fp, "FK_Title=%s\n", "RTK 観測手簿" ); /* タイトル */
		fprintf( fp, "FK_ObservationDate=%04d年%02d月%02d日UTC\n", time_inf->tm_year + 1900, time_inf->tm_mon+1, time_inf->tm_mday );   /* 観測日 */
		fprintf( fp, "FK_SessionName=\n");       /* セッション名 */
		fprintf( fp, "FK_ObservationMethod=\n");       /* 観測方法 */
        
        us = 0;
		fprintf( fp, "FK_ObservationSatellites=" ); /* 使用した周波数 */
        for(j=SYS_GPS;j<SYS_ALL;j<<=1) {
		    if( popt->navsys & j && rtk->gsi.sys & j ) {
                if(us) fprintf( fp, "  " ); /* 使用した周波数 */
                us=1;
			    fprintf( fp, "%s  ", systemstrs[sysind(j)]); /* 使用した衛星システム */
			    fprintf( fp, "%s", freqstrs[rtk->opt.oprfrq[0]]);       /* 使用した周波数 */
                for(i=1;i<nf;i++) fprintf( fp, ",%s", freqstrs[rtk->opt.oprfrq[i]]); /* 使用した周波数 */
		    }
        }
		fprintf( fp, "\n" ); /* 使用した周波数 */

        fprintf( fp, "FK_RefNumber=%s\n", stas[1].marker); /* 観測点(固定点)No */
        fprintf( fp, "FK_RefName=%s\n", stas[1].name ); /* 観測点名称 */
		if( popt->rectype[1] == '\0' ){
			fprintf( fp, "FK_RefReceiverName=%s\n", stas[1].rectype ); /* 受信機名称 */
		}
		else{
			fprintf( fp, "FK_RefReceiverName=%s\n", popt->rectype[1] ); /* 受信機名称 */
		}

		fprintf( fp, "FK_RefReceiverSerial=%s\n", stas[1].recsno ); /* 受信機シリアル */

		if( strchr( popt->anttype[1], '*' ) != NULL ){
			fprintf( fp, "FK_RefAntennaName=%s\n", stas[1].antdes ); /* アンテナ名称 */
		}
		else{
			fprintf( fp, "FK_RefAntennaName=%s\n", popt->anttype[1] ); /* アンテナ名称 */
		}
		fprintf( fp, "FK_RefAntennaSerial=%s\n", stas[1].antsno ); /* アンテナシリアル */
		fprintf( fp, "FK_RefObservationInterval=%d\n", rtk->gsi.tint ); /* データ取得間隔 */
		fprintf( fp, "FS_RefElevationAngle=%f\n", popt->elmin*R2D ); /* 最低高度角 */
		switch( popt->rovpos ){
			case 0:
				fprintf( fp, "FK_RefAntennaHeight=%f\n", popt->antdel[1][2] ); /* アンテナ底面高 */
				break;
			case 1:
			case 2:
				fprintf( fp, "FK_RefAntennaHeight=%f\n", 0.0 ); /* アンテナ底面高 */
				break;
			case 3:
				switch( stas[1].deltype ){
					case 0:
						fprintf( fp, "FK_RefAntennaHeight=%f\n", stas[1].del[2] ); /* アンテナ底面高 */
						break;
					case 1:
						ecef2pos(stas[1].pos,blh);
						ecef2enu(blh,stas[0].del,enu);
						fprintf( fp, "FK_RefAntennaHeight=%f\n", enu[2] ); /* アンテナ底面高 */
						break;
					default:
						break;
				}
				break;
			default:
				break;
		}

        fprintf( fp, "FK_RovNumber=%s\n", stas[0].marker); /* 観測点(固定点)No */
        fprintf( fp, "FK_RovName=%s\n", stas[0].name ); /* 観測点名称 */
		if( popt->rectype[0] == '\0' ){
			fprintf( fp, "FK_RovReceiverName=%s\n", stas[0].rectype ); /* 受信機名称 */
		}
		else{
			fprintf( fp, "FK_RovReceiverName=%s\n", popt->rectype[0] ); /* 受信機名称 */
		}

		fprintf( fp, "FK_RovReceiverSerial=%s\n", stas[0].recsno ); /* 受信機シリアル */

		if( strchr( popt->anttype[0], '*' ) != NULL ){
			fprintf( fp, "FK_RovAntennaName=%s\n", stas[0].antdes ); /* アンテナ名称 */
		}
		else{
			fprintf( fp, "FK_RovAntennaName=%s\n", popt->anttype[0] ); /* アンテナ名称 */
		}
		fprintf( fp, "FK_RovAntennaSerial=%s\n", stas[0].antsno ); /* アンテナシリアル */
		fprintf( fp, "FK_RovObservationInterval=%d\n", rtk->gsi.tint ); /* データ取得間隔 */
		fprintf( fp, "FS_RovElevationAngle=%f\n", popt->elmin*R2D ); /* 最低高度角 */
		switch( popt->rovpos ){
			case 0:
				fprintf( fp, "FK_RovAntennaHeight=%f\n", popt->antdel[0][2] ); /* アンテナ底面高 */
				break;
			case 1:
			case 2:
				fprintf( fp, "FK_RovAntennaHeight=%f\n", 0.0 ); /* アンテナ底面高 */
				break;
			case 3:
				switch( stas[0].deltype ){
					case 0:
						fprintf( fp, "FK_RovAntennaHeight=%f\n", stas[0].del[2] ); /* アンテナ底面高 */
						break;
					case 1:
						ecef2pos(stas[0].pos,blh);
						ecef2enu(blh,stas[0].del,enu);
						fprintf( fp, "FK_RovAntennaHeight=%f\n", enu[2] ); /* アンテナ底面高 */
						break;
					default:
						break;
				}
				break;
			default:
				break;
		}
		fprintf( fp, "FK_ObservationStart=%02d:%02d:%02d\n", time_inf->tm_hour, time_inf->tm_min, time_inf->tm_sec ); /* 観測開始時間 */

		gtm = gpst2utc(obss.data[obss.n-1].time);
		time_inf = gmtime(&(gtm.time));
		fprintf( fp, "FK_ObservationEnd=%02d:%02d:%02d\n", time_inf->tm_hour, time_inf->tm_min, time_inf->tm_sec ); /* 観測終了時間 */

		fprintf( fp, "FK_NumberSatellites=%d\n", popt->minsat ); /* 共通衛星数 */
		fprintf( fp, "FK_Remark=\n" ); /* 備考 */


		/* 単一基線解析結果ファイル：観測手簿・キネマティック用のファイルクローズ */
		fclose( fp );
    }

	return 1;

}


/*  -----------------------------------------------------
* 単一基線解析結果出力
* args   : prcopt_t *popt     I  処理条件
*          rtk_t *rtk         I  RTK構造体
* return : int     ステータス(0:異常,1:正常)
* note   :
*-------------------------------------------------------*/
static int save_sbres(const prcopt_t *popt, rtk_t *rtk)
{
	FILE *fp;
	char path[1024]={0};

	/* 単一基線解析結果ファイル：観測手簿・スタティック用 */
	int i,j,k,us,f,prn;
	double blh[3], enu[3];
	/* 単一基線解析結果ファイル：観測記簿・スタティック用 */
	gtime_t gtm;
	struct tm *time_inf;
	double pres, temp;
	double blh1[3], blh2[3];
	double dms1[3]={0},dms2[3]={0};
	double d12[3],d21[3];
    double* u=NULL;
	double* Q=NULL;
    double* uQ=NULL;
	double* uQu=NULL;
	double alpha12,s12,v12,alpha21,s21,v21;
	double durationtime;

	double a,fe,L,delta,sigma,sigma_d,delta_d,xi,xi_d,eta,eta_d,u1,u2,x,y,c,epsilon;
	double theta,g,h,sigma_s,J,K,gamma_s,gamma,zeta,zeta_d,D,E,F,G;
	double n0,A,B,D12,cs;
	double alpha;
	double d_alpha2;
	double alpha1;
	double alpha2;


	char *ret;
	char path_str[1024]={0};
	char file_str[1024]={0};
	int nf=popt->ionoopt==IONOOPT_IFLC?1:popt->nfreq;

	/* メイン画面で指定した出力先、及びファイル名を取得 */
	ret = strrchr( fsb, '.' );
	if(ret) {
		strncpy( path_str, fsb, ret - fsb );
		path_str[ret - fsb] = '\0';
	}
	//パス取り除く
	ret = strrchr( path_str, '\\' );
	if(ret)	strcpy(file_str,ret+1);
	else strcpy(file_str,path_str);

	if( popt->mode == PMODE_STATIC ) {

		/*観測衛星の状態を取得*/
		setPassData(rtk);

		for(i=0;i<2;i++){
			/* 単一基線解析結果ファイル：観測手簿・スタティック用のファイル作成 */
			make_single_report(i,path_str, file_str, popt, rtk);
		}

		/*衛星情報CVSファイル作成*/
		for (i=0;i<2;i++) {
			makeTrackingStatus(i, path_str, popt,rtk);

		}

		/*FIX解がでていない場合、処理を終了する。*/
		if(rtk->gsi.rb[0]==0.0||rtk->gsi.rb[1]==0.0|| rtk->gsi.rb[2]==0.0 ||
		   rtk->gsi.rr[0]==0.0 || rtk->gsi.rr[1]==0.0 || rtk->gsi.rr[2]==0.0 ){

			showmsg("no fix data: not make RecordStatic file");
			trace(3,"no fix data: not make RecordStatic file\n");
			return( 1 );
		}

		/* 単一基線解析結果ファイル：観測記簿・スタティック用のファイルオープン */

		strcpy( path, path_str);
		strcat( path, ".RecordStatic.out" );
		if(( fp = fopen( path, "w" )) == NULL ) {
/*
			showmsg("error : cannot open solution output file");
			trace(1,"cannot open solution output file\n");
*/
			return( 0 );
		}

        /* 情報出力 */
		fprintf( fp, "CD_FileType=%s\n", "1" ); /* ファイルタイプ */
		fprintf( fp, "RS_Title=%s\n", "GNSS測量記簿" ); /* タイトル */
		fprintf( fp, "RS_AnalysisSoftware=%s\n", "GSIPOST ver1.0.2" ); /* 解析したソフトウェア */

		switch( popt->sateph ){
			case EPHOPT_BRDC:
				fprintf( fp, "RS_AnalysisNavigation=%s\n", "放送暦" ); /* 使用した軌道情報 */
				break;
			case EPHOPT_PREC:
				fprintf( fp, "RS_AnalysisNavigation=%s\n", "精密暦" ); /* 使用した軌道情報 */
				break;
			case EPHOPT_SBAS:
				fprintf( fp, "RS_AnalysisNavigation=%s\n", "放送暦＋SBAS補正情報" ); /* 使用した軌道情報 */
				break;
			case EPHOPT_SSRAPC:
				fprintf( fp, "RS_AnalysisNavigation=%s\n", "放送暦＋SSR APC" ); /* 使用した軌道情報 */
				break;
			case EPHOPT_SSRCOM:
				fprintf( fp, "RS_AnalysisNavigation=%s\n", "放送暦＋SSR COM" ); /* 使用した軌道情報 */
				break;
			case EPHOPT_LEX:
				fprintf( fp, "RS_AnalysisNavigation=%s\n", "QZSS LEX放送暦" ); /* 使用した軌道情報 */
				break;
			default:
				break;
		}

		fprintf( fp, "RS_AnalysisEllipsoid=%s\n", "GRS-80" ); /* 使用した楕円体 */
        us = 0;
		fprintf( fp, "RS_AnalysisSatellites=" ); /* 使用した周波数 */
        for(j=SYS_GPS;j<SYS_ALL;j<<=1) {
		    if( popt->navsys & j && rtk->gsi.sys & j ) {
                if(us) fprintf( fp, "  " ); /* 使用した周波数 */
                us=1;
			    fprintf( fp, "%s  ", systemstrs[sysind(j)]); /* 使用した衛星システム */
			    fprintf( fp, "%s", freqstrs[rtk->opt.oprfrq[0]]);       /* 使用した周波数 */
                for(i=1;i<nf;i++) fprintf( fp, ",%s", freqstrs[rtk->opt.oprfrq[i]]); /* 使用した周波数 */
		    }
        }
		fprintf( fp, "\n" ); /* 使用した周波数 */

        switch( popt->mode ){
			case PMODE_SINGLE:
				fprintf( fp, "RS_BaselineAnalysis=%s\n", "単独測位" ); /* 基線解析モード */
				break;
			case PMODE_DGPS:
				fprintf( fp, "RS_BaselineAnalysis=%s\n", "DGPS/DGNSS" ); /* 基線解析モード */
				break;
			case PMODE_KINEMA:
				fprintf( fp, "RS_BaselineAnalysis=%s\n", "キネマティック測位" ); /* 基線解析モード */
				break;
			case PMODE_STATIC:
				fprintf( fp, "RS_BaselineAnalysis=%s\n", "スタティック測位" ); /* 基線解析モード */
				break;
			case PMODE_MOVEB:
				fprintf( fp, "RS_BaselineAnalysis=%s\n", "既知局の移動" ); /* 基線解析モード */
				break;
			case PMODE_FIXED:
				fprintf( fp, "RS_BaselineAnalysis=%s\n", "既知局・新局固定" ); /* 基線解析モード */
				break;
			case PMODE_PPP_KINEMA:
				fprintf( fp, "RS_BaselineAnalysis=%s\n", "キネマティックPPP" ); /* 基線解析モード */
				break;
			case PMODE_PPP_STATIC:
				fprintf( fp, "RS_BaselineAnalysis=%s\n", "スタティックPPP" ); /* 基線解析モード */
				break;
			case PMODE_MULTI:
				fprintf( fp, "RS_BaselineAnalysis=%s\n", "複数基線解析" ); /* 基線解析モード */
				break;
			//	kawauchi
//			case PMODE_PPPAR:
//				fprintf( fp, "RS_BaselineAnalysis=%s\n", "PPP-AR" ); /* 基線解析モード */
//				break;
			default:
				break;
		}

		fprintf( fp, "RS_SessionName=\n"); /* セッション名 */

		/* 解析使用データ */
		gtm = gpst2utc(obss.data[0].time);
		time_inf = gmtime(&(gtm.time));
		fprintf( fp, "RS_AnalysisStart=%04d/%02d/%02d %02d:%02d:%02d UTC\n", time_inf->tm_year + 1900, time_inf->tm_mon+1, time_inf->tm_mday, time_inf->tm_hour, time_inf->tm_min, time_inf->tm_sec ); /* 開始 */

		gtm = gpst2utc(obss.data[obss.n-1].time);
		time_inf = gmtime(&(gtm.time));
		fprintf( fp, "RS_AnalysisEnd=%04d/%02d/%02d %02d:%02d:%02d UTC\n", time_inf->tm_year + 1900, time_inf->tm_mon+1, time_inf->tm_mday, time_inf->tm_hour, time_inf->tm_min, time_inf->tm_sec ); /* 終了 */

		fprintf( fp, "RS_ElevationAngle=%f\n", popt->elmin*R2D ); /* 最低高度角 */

		if(stas[0].pos[0]==0 || stas[0].pos[1]==0){
			cov2ecef(rtk->gsi.rr,NULL,blh,NULL,RE_GRS80,FE_GRS80,XYZ_TO_PLH);
		}else{
			ecef2pos(stas[0].pos,blh);
		}
//		pres=1013.25*pow(1.0-2.2557E-5*blh[2],5.2568);
//		fprintf( fp, "RS_Pressure=%f\n", pres); /* 気圧 */
		fprintf( fp, "RS_Pressure=%d\n", RPT_PRESSURE); /* 気圧 */

/*		temp=15.0-6.5E-3*blh[2]+273.16;*/
//		temp=15.0-6.5E-3*blh[2];
//		fprintf( fp, "RS_Temperature=%f\n", temp); /* 温度 */
		fprintf( fp, "RS_Temperature=%d\n", RPT_TEMPERATURE); /* 温度 */

		fprintf( fp, "RS_Humidity=%d\n", RPT_HUMIDITY); /* 湿度 */

        /* 観測点１ */
		fprintf( fp, "RS_Observation1=%s(%s)\n", stas[1].name, stas[1].name3 ); /* 観測点１ */

		if( popt->rectype[1] == '\0' ){
			fprintf( fp, "RS_Obs1ReceiverName=%s\n", stas[1].rectype ); /* 受信機名 */
		}
		else{
			fprintf( fp, "RS_Obs1ReceiverName=%s\n", popt->rectype[1] ); /* 受信機名 */
		}

		fprintf( fp, "RS_Obs1ReceiverSerial=%s\n", stas[1].recsno ); /* 受信機シリアル */

		if( strchr( popt->anttype[1], '*' ) != NULL ){
			fprintf( fp, "RS_Obs1AntennaName=%s\n", stas[1].antdes ); /* アンテナ名 */
		}
		else{
			fprintf( fp, "RS_Obs1AntennaName=%s\n", popt->anttype[1] ); /* アンテナ名 */
		}

		fprintf( fp, "RS_Obs1AntennaSerial=%s\n", stas[1].antsno ); /* アンテナシリアル */


		if (strcmp(stas[1].pcvr.date,"")!=0) {
			fprintf( fp, "RS_Obs1AntennaPCV=%s(%s)\n", "あり", stas[1].pcvr.date ); /* PCV補正 */
		}else{
			fprintf( fp, "RS_Obs1AntennaPCV=%s\n", "なし" ); /* PCV補正 */
		}

        switch( popt->rovpos ){
			case 0:
				fprintf( fp, "RS_Obs1AntennaHeight=%f\n", popt->antdel[1][2] ); /* アンテナ底面高 */
				break;
			case 1:
			case 2:
				fprintf( fp, "RS_Obs1AntennaHeight=%f\n", 0.0 ); /* アンテナ底面高 */
				break;
			case 3:
				switch( stas[1].deltype ){
					case 0:
						fprintf( fp, "RS_Obs1AntennaHeight=%f\n", stas[1].del[2] ); /* アンテナ底面高 */
						break;
					case 1:
						ecef2pos(stas[1].pos,blh);
						ecef2enu(blh,stas[1].del,enu);
						fprintf( fp, "RS_Obs1AntennaHeight=%f\n", enu[2] ); /* アンテナ底面高 */
						break;
					default:
						break;
				}
				break;
			default:
				break;
		}

		/* 観測点２ */
		fprintf( fp, "RS_Observation2=%s(%s)\n", stas[0].name, stas[0].name3 ); /* 観測点１ */

		if( popt->rectype[0] == '\0' ){
			fprintf( fp, "RS_Obs2ReceiverName=%s\n", stas[0].rectype ); /* 受信機名 */
		}
		else{
			fprintf( fp, "RS_Obs2ReceiverName=%s\n", popt->rectype[0] ); /* 受信機名 */
		}

		fprintf( fp, "RS_Obs2ReceiverSerial=%s\n", stas[0].recsno ); /* 受信機シリアル */

		if( strchr( popt->anttype[0], '*' ) != NULL ){
			fprintf( fp, "RS_Obs2AntennaName=%s\n", stas[0].antdes ); /* アンテナ名 */
		}
		else{
			fprintf( fp, "RS_Obs2AntennaName=%s\n", popt->anttype[0] ); /* アンテナ名 */
		}

		fprintf( fp, "RS_Obs2AntennaSerial=%s\n", stas[0].antsno ); /* アンテナシリアル */

		if (strcmp(stas[0].pcvr.date,"")!=0) {
			fprintf( fp, "RS_Obs2AntennaPCV=%s(%s)\n", "あり", stas[0].pcvr.date ); /* PCV補正 */
		}else{
			fprintf( fp, "RS_Obs2AntennaPCV=%s\n", "なし" ); /* PCV補正 */
		}

		switch( popt->rovpos ){
			case 0:
				fprintf( fp, "RS_Obs2AntennaHeight=%f\n", popt->antdel[0][2] ); /* アンテナ底面高 */
				break;
			case 1:
			case 2:
				fprintf( fp, "RS_Obs2AntennaHeight=%f\n", 0.0 ); /* アンテナ底面高 */
				break;
			case 3:
				switch( stas[0].deltype ){
					case 0:
						fprintf( fp, "RS_Obs2AntennaHeight=%f\n", stas[0].del[2] ); /* アンテナ底面高 */
						break;
					case 1:
						ecef2pos(stas[0].pos,blh);
						ecef2enu(blh,stas[0].del,enu);
						fprintf( fp, "RS_Obs2AntennaHeight=%f\n", enu[2] ); /* アンテナ底面高 */
						break;
					default:
						break;
				}
				break;
			default:
				break;
		}

		/* 起点 */
		cov2ecef(rtk->gsi.rb,NULL,blh1,NULL,RE_GRS80,FE_GRS80,XYZ_TO_PLH);
		degtodms(blh1[0]*R2D,dms1);
		fprintf( fp, "RS_StartLatitude=%.0f° %02.0f' %07.4f\"\n",fabs(dms1[0]),dms1[1],dms1[2]); /* 緯度 */

		degtodms(blh1[1]*R2D,dms2);
		fprintf( fp, "RS_StartLongitude=%.0f° %02.0f' %07.4f\"\n",fabs(dms2[0]),dms2[1],dms2[2]); /* 経度 */

		fprintf( fp, "RS_StartEllipsoidHeight=%f\n", blh1[2] ); /* 楕円体高 */
		fprintf( fp, "RS_StartCoordinateX=%f\n", rtk->gsi.rb[0] ); /* 座標値X */
		fprintf( fp, "RS_StartCoordinateY=%f\n", rtk->gsi.rb[1] ); /* 座標値Y */
		fprintf( fp, "RS_StartCoordinateZ=%f\n", rtk->gsi.rb[2] ); /* 座標値Z */

		/* 終点 */
		cov2ecef(rtk->gsi.rr,NULL,blh2,NULL,RE_GRS80,FE_GRS80,XYZ_TO_PLH);
		degtodms(blh2[0]*R2D,dms1);
		fprintf( fp, "RS_EndLatitude=%.0f° %02.0f' %07.4f\"\n",fabs(dms1[0]),dms1[1],dms1[2]); /* 緯度 */

		degtodms(blh2[1]*R2D,dms2);
		fprintf( fp, "RS_EndLongitude=%.0f° %02.0f' %07.4f\"\n",fabs(dms2[0]),dms2[1],dms2[2]); /* 経度 */

		fprintf( fp, "RS_EndEllipsoidHeight=%f\n", blh2[2] ); /* 楕円体高 */
		fprintf( fp, "RS_EndCoordinateX=%f\n", rtk->gsi.rr[0] ); /* 座標値X */
		fprintf( fp, "RS_EndCoordinateY=%f\n", rtk->gsi.rr[1] ); /* 座標値Y */
		fprintf( fp, "RS_EndCoordinateZ=%f\n", rtk->gsi.rr[2] ); /* 座標値Z */

		/* 解析結果 */
		fprintf( fp, "RS_AnalysisQuality=%s\n", "FIX" ); /* 解の種類 */
		fprintf( fp, "RS_AnalysisRatio=%f\n", rtk->gsi.ratio ); /* バイアス決定比 */

		/* 観測点１－２ */
		d12[0] = rtk->gsi.rr[0] - rtk->gsi.rb[0];
		d12[1] = rtk->gsi.rr[1] - rtk->gsi.rb[1];
		d12[2] = rtk->gsi.rr[2] - rtk->gsi.rb[2];
		fprintf( fp, "RS_Obs12DeltaX=%f\n", d12[0] ); /* DX */
		fprintf( fp, "RS_Obs12DeltaY=%f\n", d12[1] ); /* DY */
		fprintf( fp, "RS_Obs12DeltaZ=%f\n", d12[2] ); /* DZ */
		fprintf( fp, "RS_Obs12Distance=%f\n", norm(d12,3) ); /* 斜距離 */

		/* 標準偏差 */
		fprintf( fp, "RS_SDDeltaX=%e\n", SQRT(rtk->gsi.qr[0]) ); /* DX */
		fprintf( fp, "RS_SDDeltaY=%e\n", SQRT(rtk->gsi.qr[1]) ); /* DY */
		fprintf( fp, "RS_SDDeltaZ=%e\n", SQRT(rtk->gsi.qr[2]) ); /* DZ */

		u = mat(3, 1);
		if (!normv3(d12,u)){
		fprintf( fp, "RS_SDDistance=%e\n", 0.0 ); /* 斜距離 */
		}else{

		//u[0] = d12[0];
		//u[1] = d12[1];
		//u[2] = d12[2];

		Q = mat(3, 3);
/*
		Q[0] = rtk->gsi.qr[0];
		Q[4] = rtk->gsi.qr[3];
		Q[8] = rtk->gsi.qr[5];
		Q[1] = Q[3] = rtk->gsi.qr[1];
		Q[5] = Q[7] = rtk->gsi.qr[4];
		Q[2] = Q[6] = rtk->gsi.qr[2];
*/
		Q[0] = rtk->gsi.qr[0];
		Q[4] = rtk->gsi.qr[1];
		Q[8] = rtk->gsi.qr[2];
		Q[1] = Q[3] = rtk->gsi.qr[3];
		Q[5] = Q[7] = rtk->gsi.qr[4];
		Q[2] = Q[6] = rtk->gsi.qr[5];


		uQ = mat(3, 1);
		matmul("NN", 1, 3, 3, 1.0, u, Q, 0.0, uQ);

		uQu = mat(1, 1);
		matmul("NT", 1, 1, 3, 1.0, uQ, u, 0.0, uQu);

		fprintf( fp, "RS_SDDistance=%e\n", SQRT(fabs(*uQu)) ); /* 斜距離 */
		}

		trace(5," blh1[0]=%f blh1[1]=%f blh1[2]=%f \n", blh1[0], blh1[1], blh1[2]);
		trace(5," blh2[0]=%f blh2[1]=%f blh2[2]=%f \n", blh2[0], blh2[1], blh2[2]);

		a = RE_GRS80;
		fe = FE_GRS80;
		L = blh2[1] - blh1[1];
		if(L<-PI) L = 2.0*PI+L;
		if(L>PI) L = L-2.0*PI;
		delta = blh2[0] - blh1[0];
		sigma = blh1[0] + blh2[0];
		sigma_d = atan(((1.0-fe)*sin(sigma))/(cos(sigma)+fe*(2.0-fe)*sin(blh1[0])*sin(blh2[0])));
		delta_d = atan(((1.0-fe)*sin(delta))/(cos(delta)-fe*(2.0-fe)*sin(blh1[0])*sin(blh2[0])));
		xi = cos(sigma_d/2.0);
		xi_d = sin(sigma_d/2.0);
		eta = sin(delta_d/2.0);
		eta_d = cos(delta_d/2.0);
		u1 = atan((1.0-fe)*tan(blh1[0]));
		u2 = atan((1.0-fe)*tan(blh2[0]));
		x = sin(u1)*sin(u2);
		y = cos(u1)*cos(u2);
		c = y*cos(L)+x;
		epsilon = fe * (2.0-fe) / ((1.0-fe)*(1.0-fe));

		trace(5,"u1=%f u2=%f \n", u1,u2);

		trace(5,"x=%f y=%f c=%f epsilon=%f \n", x,y,c,epsilon);


		if(c<0){
			//fprintf( fp, "RS_GeodesicDistance=%f\n", 0.0 ); /* 測地線長計算不可パターン */
		}else{
			theta = L*(1+fe*y);

			for(i=0;1;i++){
				g = SQRT((eta*eta*cos(theta/2.0)*cos(theta/2.0))+(xi*xi*sin(theta/2.0)*sin(theta/2.0)));
				h = SQRT((eta_d*eta_d*cos(theta/2.0)*cos(theta/2.0))+(xi_d*xi_d*sin(theta/2.0)*sin(theta/2.0)));
				sigma_s = 2.0*atan(g/h);
				J = 2.0*g*h;
				if(J ==0) break;
				K = (h*h)-(g*g);
				gamma_s = y * sin(theta) / J;
				gamma = 1-(gamma_s*gamma_s);
				zeta = (gamma*K)-(2.0*x);
				zeta_d = zeta+x;
				D = (fe*(1.0+fe)/4.0)-(3.0/16.0)*fe*fe*gamma;
				E = (1.0-D*gamma)*fe*gamma_s*(sigma_s+D*J*(zeta+D*K*(2.0*zeta*zeta-gamma*gamma)));
				F = theta-L-E;
				G = fe*gamma_s*gamma_s*(1.0-2.0*D*gamma)+fe*zeta_d*(sigma_s/J)*(1.0-D*gamma+fe*gamma_s*gamma_s/2.0)+fe*fe*zeta*zeta_d/4.0;

				theta = theta-F/(1.0-G);
				if((fabs(F)<1.e-15)||(10000<i)){
					break;
				}
			}
		}

/*
		alpha12=atan((-sin(blh1[1])*(rtk->gsi.rr[0]-rtk->gsi.rb[0])+cos(blh1[1])*(rtk->gsi.rr[1]-rtk->gsi.rb[1]))
		/(-sin(blh1[0])*cos(blh1[1])*(rtk->gsi.rr[0]-rtk->gsi.rb[0])-sin(blh1[0])*sin(blh1[1])*(rtk->gsi.rr[1]-rtk->gsi.rb[1])+cos(blh1[0])*(rtk->gsi.rr[2]-rtk->gsi.rb[2])));
		degtodms(alpha12*R2D,dms1);
		fprintf( fp, "RS_Obs12AzimuthAngle=%.0f° %02.0f' %07.4f\"\n",fabs(dms1[0]),dms1[1],dms1[2]); /* 方位角 * /
*/

		alpha = atan2(xi*tan(theta/2.),eta);
		d_alpha2 = atan2(xi_d*tan(theta/2.),eta_d);
		alpha1 = alpha - d_alpha2;
		alpha2 = alpha + d_alpha2;
		trace(5,"α=%f α12=%f α21=%f (deg.)\n", alpha*180./PI, alpha1*180./PI, alpha21*180./PI);

		alpha21 = PI + alpha2;
		if(alpha1 < 0) alpha1+=2*PI;
		if(alpha21 < 0) alpha21+=2*PI;
		if(alpha1 > 2*PI) alpha1-=2*PI;
		if(alpha21 > 2*PI) alpha21-=2*PI;
		degtodms(alpha1*R2D,dms1);
		fprintf( fp, "RS_Obs12AzimuthAngle=%.0f° %02.0f' %07.4f\"\n",fabs(dms1[0]),dms1[1],dms1[2]); /* 方位角 */

		s12=SQRT((d12[0]*d12[0])+(d12[1]*d12[1])+(d12[2]*d12[2]));

		cs=(cos(blh1[0])*cos(blh1[1])*(d12[0])+cos(blh1[0])*sin(blh1[1])*(d12[1])+sin(blh1[0])*(d12[2]))/s12;

		if(-1<cs && cs<1){
			v12=asin(cs);
		}else{
			v12=asin(0);
		}
		degtodms(v12*R2D,dms1);

		trace(5,"d12[0]=%f d12[1]=%f d12[2]=%f \n", d12[0], d12[1], d12[2]);
		trace(5,"blh1[0]=%f blh1[1]=%f blh1[2]=%f \n", blh1[0], blh1[1], blh1[2]);
		trace(5,"s12=%f cs=%f v12=%f v12*R2D=%f dms=%f,%f,%f\n", s12, cs, v12,v12*R2D,dms1[0],dms1[1],dms1[2]);
		fprintf( fp, "RS_Obs12HeightAngle=%.0f° %02.0f' %07.4f\"\n",dms1[0],dms1[1],dms1[2]); /* 高度角 */

		n0 = (epsilon*gamma)/((SQRT(1.0+(epsilon*gamma))+1.0)*(SQRT(1.0+(epsilon*gamma))+1.0));
		A = (1.0+n0)*(1.0+(5.0/4.0)*n0*n0);
		B = (epsilon*(1.0-(3.0/8.0)*n0*n0))/((SQRT(1.0+(epsilon*gamma))+1.0)*(SQRT(1.0+(epsilon*gamma))+1.0));
		D12 = (1.0-fe)*a*A*(sigma_s-B*J*(zeta-B/4.0*(K*(gamma*gamma-2.0*zeta*zeta)-(B/6.0)*zeta*(1.0-4.0*K*K)*(3.0*gamma*gamma-4.0*zeta*zeta))));

		trace(5,"n0=%f A=%f B=%f D12=%f\n", n0,A,B,D12);

		if(c<0){
			fprintf( fp, "RS_GeodesicDistance=%f\n", 0.0 ); /* 測地線長計算不可パターン */
		}else{
			fprintf( fp, "RS_GeodesicDistance=%f\n", D12 ); /* 測地線長 */
		}

		fprintf( fp, "RS_EllipsoidDistance=%f\n", blh2[2]-blh1[2] ); /* 楕円体比高 */

		/* 観測点２－１ */
/*
		alpha21=atan((-sin(blh2[1])*(rtk->gsi.rb[0]-rtk->gsi.rr[0])+cos(blh2[1])*(rtk->gsi.rb[1]-rtk->gsi.rr[1]))
		/(-sin(blh2[0])*cos(blh2[1])*(rtk->gsi.rb[0]-rtk->gsi.rr[0])-sin(blh2[0])*sin(blh2[1])*(rtk->gsi.rb[1]-rtk->gsi.rr[1])+cos(blh2[0])*(rtk->gsi.rb[2]-rtk->gsi.rr[2])));
*/
		degtodms(alpha21*R2D,dms1);

		fprintf( fp, "RS_Obs21AzimuthAngle=%.0f° %02.0f' %07.4f\"\n",fabs(dms1[0]),dms1[1],dms1[2]); /* 方位角 */

		d21[0] = rtk->gsi.rb[0] - rtk->gsi.rr[0];
		d21[1] = rtk->gsi.rb[1] - rtk->gsi.rr[1];
		d21[2] = rtk->gsi.rb[2] - rtk->gsi.rr[2];

		s21=SQRT((d21[0]*d21[0])+(d21[1]*d21[1])+(d21[2]*d21[2]));

		cs = (cos(blh2[0])*cos(blh2[1])*(d21[0])+cos(blh2[0])*sin(blh2[1])*(d21[1])+sin(blh2[0])*(d21[2]))/s21;

		if(-1<cs && cs<1){
			v21=asin(cs);
		}else{
			v21=asin(0);
		}
		degtodms(v21*R2D,dms1);
		trace(5,"d21[0]=%f d21[1]=%f d21[2]=%f \n", d21[0], d21[1], d21[2]);
		trace(5,"blh2[0]=%f blh2[1]=%f blh2[2]=%f \n", blh2[0], blh2[1], blh2[2]);
		trace(5,"s21=%f cs=%f v21=%f v12*R2D=%f dms=%f,%f,%f\n", s21, cs, v21,v21*R2D,dms1[0],dms1[1],dms1[2]);

		fprintf( fp, "RS_Obs21HeightAngle=%.0f° %02.0f' %07.4f\"\n",dms1[0],dms1[1],dms1[2]); /* 高度角 */

		/* 分散・共分散 */
		fprintf( fp, "RS_VarCovDXDX=%e\n", rtk->gsi.qr[0] ); /* DX*DX */
		fprintf( fp, "RS_VarCovDXDY=%e\n", rtk->gsi.qr[3] ); /* DX*DY */
		fprintf( fp, "RS_VarCovDXDZ=%e\n", rtk->gsi.qr[5] ); /* DX*DZ */
		fprintf( fp, "RS_VarCovDYDY=%e\n", rtk->gsi.qr[1] ); /* DY*DY */
		fprintf( fp, "RS_VarCovDYDZ=%e\n", rtk->gsi.qr[4] ); /* DY*DZ */
		fprintf( fp, "RS_VarCovDZDZ=%e\n", rtk->gsi.qr[2] ); /* DZ*DZ */

		fprintf( fp, "RS_NumberEpochUse=%d\n", rtk->gsi.ndata ); /* 使用したデータ数 */

		rtk->gsi.nrej=0;
		for (i=0;i<MAXSAT;i++) for (f=0;f<nf;f++) {
		  rtk->gsi.nrej+=rtk->ssat[i].rejc[f];
			trace(5,"nrej: sat[%d]f[%d] rejc[%d]nrej[%d]\n",i,f,rtk->ssat[i].rejc[f],rtk->gsi.nrej);
		}

		fprintf( fp, "RS_NumberEpochReject=%d\n", rtk->gsi.nrej ); /* 棄却したデータ数 */
		if(rtk->gsi.ndata==0){
			fprintf( fp, "RS_RejectRate=%.1f\n", 0.0); /* 棄却率 */
		}else{
			fprintf( fp, "RS_RejectRate=%.2f\n", (double)(rtk->gsi.nrej)/(double)(rtk->gsi.nrej+rtk->gsi.ndata)*100 ); /* 棄却率 */
		}
		fprintf( fp, "RS_EpochInterval=%d\n", rtk->gsi.tint ); /* 使用したデータ数 */
		fprintf( fp, "RS_RMS=%f\n", SQRT(rtk->gsi.var) ); /* RMS */
		fprintf( fp, "RS_RDOP=%f\n", rtk->gsi.rdop ); /* RDOP */

		/* 単一基線解析結果ファイル：観測記簿・スタティック用 */
		SAFE_FREE(u);
		SAFE_FREE(Q);
		SAFE_FREE(uQ);
		SAFE_FREE(uQu);

		/* 単一基線解析結果ファイル：観測記簿・スタティック用のファイルクローズ */
		fclose( fp );
	}
	else if( popt->mode == PMODE_KINEMA ) {
		/* 単一基線解析結果ファイル：観測手簿・キネマティック用のファイルオープン */
//		strcpy( path, popt->mopt.ofdir );
//		if( path[ strlen( path ) - 1 ] != '\\' ) strcat( path, "\\" );
//		strcat( path, outfbase );
//		strcat( path, ".FieldbookKinematic.out" );
//		if(( fp = fopen( path, "w" )) == NULL ) {
///*
//			showmsg("error : cannot open solution output file");
//			trace(1,"cannot open solution output file\n");
//*/
//			return( 0 );
//		}
//
//		/* 情報出力 */
//
//		/* 単一基線解析結果ファイル：観測手簿・キネマティック用のファイルクローズ */
//		fclose( fp );
//

        /*観測衛星の状態を取得*/
		setPassData(rtk);

		/* 単一基線解析結果ファイル：観測手簿・キネマティック用のファイル作成 */
		make_single_report(0,path_str, file_str, popt, rtk);
		
		/*FIX解がでていない場合、処理を終了する。*/
		if(rtk->gsi.rb[0]==0.0||rtk->gsi.rb[1]==0.0|| rtk->gsi.rb[2]==0.0 ||
		   rtk->gsi.rr[0]==0.0 || rtk->gsi.rr[1]==0.0 || rtk->gsi.rr[2]==0.0 ){

			showmsg("no fix data: not make RecordStatic file");
			trace(3,"no fix data: not make RecordStatic file\n");
			return( 1 );
		}

//		/* 単一基線解析結果ファイル：観測記簿・キネマティック用のファイルオープン */
//		strcpy( path, popt->mopt.ofdir );
//		if( path[ strlen( path ) - 1 ] != '\\' ) strcat( path, "\\" );
//		strcat( path, outfbase );
//		strcat( path, ".RecordKinematic.out" );
//		if(( fp = fopen( path, "w" )) == NULL ) {
///*
//			showmsg("error : cannot open solution output file");
//			trace(1,"cannot open solution output file\n");
//*/
//			return( 0 );
//		}
//
//		/* 情報出力 */
//
//		/* 単一基線解析結果ファイル：観測記簿・キネマティック用のファイルクローズ */
//		fclose( fp );

        /* 単一基線解析結果ファイル：観測記簿・キネマティック用のファイルオープン */

		strcpy( path, path_str);
		strcat( path, ".RecordKinematic.out" );
		if(( fp = fopen( path, "w" )) == NULL ) {
/*
			showmsg("error : cannot open solution output file");
			trace(1,"cannot open solution output file\n");
*/
			return( 0 );
		}

        /* 情報出力 */
		fprintf( fp, "CD_FileType=%s\n", "3" ); /* ファイルタイプ */
		fprintf( fp, "RK_Title=%s\n", "RTK 観測記簿" ); /* タイトル */
		fprintf( fp, "RK_AnalysisSoftware=%s%s\n", "GSIPOST ver", VER_RTKLIB ); /* 解析したソフトウェア */

		fprintf( fp, "RK_AnalysisEllipsoid=%s\n", "GRS-80" ); /* 使用した楕円体 */

        us = 0;
		fprintf( fp, "FK_ObservationSatellites=" ); /* 使用した周波数 */
        for(j=SYS_GPS;j<SYS_ALL;j<<=1) {
		    if( popt->navsys & j && rtk->gsi.sys & j ) {
                if(us) fprintf( fp, "  " ); /* 使用した周波数 */
                us=1;
			    fprintf( fp, "%s  ", systemstrs[sysind(j)]); /* 使用した衛星システム */
			    fprintf( fp, "%s", freqstrs[rtk->opt.oprfrq[0]]);       /* 使用した周波数 */
                for(i=1;i<nf;i++) fprintf( fp, ",%s", freqstrs[rtk->opt.oprfrq[i]]); /* 使用した周波数 */
		    }
        }
		fprintf( fp, "\n" ); /* 使用した周波数 */
        
		fprintf( fp, "RK_AnalysisMethod=\n" ); /* 解析方法 */
        fprintf( fp, "RK_BaselineAnalysis=%s\n", "キネマティック測位" ); /* 基線解析モード */
        fprintf( fp, "RK_SessionName=\n"); /* セッション名 */

		/* 解析使用データ */
		gtm = gpst2utc(obss.data[0].time);
		time_inf = gmtime(&(gtm.time));
		fprintf( fp, "RK_AnalysisStart=%04d/%02d/%02d %02d:%02d:%02d UTC\n", time_inf->tm_year + 1900, time_inf->tm_mon+1, time_inf->tm_mday, time_inf->tm_hour, time_inf->tm_min, time_inf->tm_sec ); /* 開始 */

		gtm = gpst2utc(obss.data[obss.n-1].time);
		time_inf = gmtime(&(gtm.time));
		fprintf( fp, "RK_AnalysisEnd=%04d/%02d/%02d %02d:%02d:%02d UTC\n", time_inf->tm_year + 1900, time_inf->tm_mon+1, time_inf->tm_mday, time_inf->tm_hour, time_inf->tm_min, time_inf->tm_sec ); /* 終了 */

		fprintf( fp, "RK_ElevationAngle=%f\n", popt->elmin*R2D ); /* 最低高度角 */

		if(stas[0].pos[0]==0 || stas[0].pos[1]==0){
			cov2ecef(rtk->gsi.rr,NULL,blh,NULL,RE_GRS80,FE_GRS80,XYZ_TO_PLH);
		}else{
			ecef2pos(stas[0].pos,blh);
		}
		fprintf( fp, "RK_Pressure=%d\n", RPT_PRESSURE); /* 気圧 */
		fprintf( fp, "RK_Temperature=%d\n", RPT_TEMPERATURE); /* 温度 */
		fprintf( fp, "RK_Humidity=%d\n", RPT_HUMIDITY); /* 湿度 */

        /* 観測点１ */
		fprintf( fp, "RK_Observation1=%s(%s)\n", stas[1].name, stas[1].name3 ); /* 観測点１ */

		if( popt->rectype[1] == '\0' ){
			fprintf( fp, "RK_Obs1ReceiverName=%s\n", stas[1].rectype ); /* 受信機名 */
		}
		else{
			fprintf( fp, "RK_Obs1ReceiverName=%s\n", popt->rectype[1] ); /* 受信機名 */
		}

		fprintf( fp, "RK_Obs1ReceiverSerial=%s\n", stas[1].recsno ); /* 受信機シリアル */

        switch( popt->rovpos ){
			case 0:
				fprintf( fp, "RK_Obs1AntennaHeight=%f\n", popt->antdel[1][2] ); /* アンテナ底面高 */
				break;
			case 1:
			case 2:
				fprintf( fp, "RK_Obs1AntennaHeight=%f\n", 0.0 ); /* アンテナ底面高 */
				break;
			case 3:
				switch( stas[1].deltype ){
					case 0:
						fprintf( fp, "RK_Obs1AntennaHeight=%f\n", stas[1].del[2] ); /* アンテナ底面高 */
						break;
					case 1:
						ecef2pos(stas[1].pos,blh);
						ecef2enu(blh,stas[1].del,enu);
						fprintf( fp, "RK_Obs1AntennaHeight=%f\n", enu[2] ); /* アンテナ底面高 */
						break;
					default:
						break;
				}
				break;
			default:
				break;
		}

		/* 観測点２ */
		fprintf( fp, "RK_Observation2=%s(%s)\n", stas[0].name, stas[0].name3 ); /* 観測点１ */

		if( popt->rectype[0] == '\0' ){
			fprintf( fp, "RK_Obs2ReceiverName=%s\n", stas[0].rectype ); /* 受信機名 */
		}
		else{
			fprintf( fp, "RK_Obs2ReceiverName=%s\n", popt->rectype[0] ); /* 受信機名 */
		}

		fprintf( fp, "RK_Obs2ReceiverSerial=%s\n", stas[0].recsno ); /* 受信機シリアル */

		switch( popt->rovpos ){
			case 0:
				fprintf( fp, "RK_Obs2AntennaHeight=%f\n", popt->antdel[0][2] ); /* アンテナ底面高 */
				break;
			case 1:
			case 2:
				fprintf( fp, "RK_Obs2AntennaHeight=%f\n", 0.0 ); /* アンテナ底面高 */
				break;
			case 3:
				switch( stas[0].deltype ){
					case 0:
						fprintf( fp, "RK_Obs2AntennaHeight=%f\n", stas[0].del[2] ); /* アンテナ底面高 */
						break;
					case 1:
						ecef2pos(stas[0].pos,blh);
						ecef2enu(blh,stas[0].del,enu);
						fprintf( fp, "RK_Obs2AntennaHeight=%f\n", enu[2] ); /* アンテナ底面高 */
						break;
					default:
						break;
				}
				break;
			default:
				break;
		}

        /* 観測点１ */
        cov2ecef(popt->rb,NULL,blh1,NULL,RE_GRS80,FE_GRS80,XYZ_TO_PLH);
		degtodms(blh1[0]*R2D,dms1);
		degtodms(blh1[1]*R2D,dms2);
        fprintf( fp, "RK_Obs1Latitude=%.0f° %02.0f' %07.4f\n", dms1[0],dms1[1],dms1[2]); /* 緯度 */
        fprintf( fp, "RK_Obs1Longitude=%.0f° %02.0f' %07.4f\n", dms2[0],dms2[1],dms2[2]); /* 経度 */
        fprintf( fp, "RK_Obs1EllipsoidHeight=%f\n", blh1[2]); /* 楕円体高 */
        fprintf( fp, "RK_Obs1CoordinateX=%f\n",popt->rb[0]); /* 座標値Ｘ */
        fprintf( fp, "RK_Obs1CoordinateY=%f\n",popt->rb[1]); /* 座標値Ｙ */
        fprintf( fp, "RK_Obs1CoordinateZ=%f\n",popt->rb[2]); /* 座標値Ｚ */
        /* 観測点２ */
		cov2ecef(rtk->gsi.rr_fixAve,NULL,blh2,NULL,RE_GRS80,FE_GRS80,XYZ_TO_PLH);
		degtodms(blh2[0]*R2D,dms1);
		degtodms(blh2[1]*R2D,dms2);
        fprintf( fp, "RK_Obs2Latitude=%.0f° %02.0f' %07.4f\n", dms1[0],dms1[1],dms1[2]); /* 緯度 */
        fprintf( fp, "RK_Obs2Longitude=%.0f° %02.0f' %07.4f\n", dms2[0],dms2[1],dms2[2]); /* 経度 */
        fprintf( fp, "RK_Obs2EllipsoidHeight=%f\n", blh2[2]); /* 楕円体高 */
        fprintf( fp, "RK_Obs2CoordinateX=%f\n", rtk->gsi.rr_fixAve[0]); /* 座標値Ｘ */
        fprintf( fp, "RK_Obs2CoordinateY=%f\n", rtk->gsi.rr_fixAve[1]); /* 座標値Ｙ */
        fprintf( fp, "RK_Obs2CoordinateZ=%f\n", rtk->gsi.rr_fixAve[2]); /* 座標値Ｚ */

		/* 解析結果 */
		fprintf( fp, "RK_AnalysisQuality=%s\n", "FIX" ); /* 解の種類 */
		/* 観測点１－２ */
		d12[0] = rtk->gsi.rr[0] - rtk->gsi.rb[0];
		d12[1] = rtk->gsi.rr[1] - rtk->gsi.rb[1];
		d12[2] = rtk->gsi.rr[2] - rtk->gsi.rb[2];
		fprintf( fp, "RK_Obs12DeltaX=%f\n", d12[0] ); /* DX */
		fprintf( fp, "RK_Obs12DeltaY=%f\n", d12[1] ); /* DY */
		fprintf( fp, "RK_Obs12DeltaZ=%f\n", d12[2] ); /* DZ */
		fprintf( fp, "RK_Obs12Distance=%f\n", norm(d12,3) ); /* 斜距離 */
		fprintf( fp, "RK_VarCovDXDX=%e\n", rtk->gsi.qr[0] ); /* DX*DX */
		fprintf( fp, "RK_VarCovDXDY=%e\n", rtk->gsi.qr[3] ); /* DX*DY */
		fprintf( fp, "RK_VarCovDXDZ=%e\n", rtk->gsi.qr[5] ); /* DX*DZ */
		fprintf( fp, "RK_VarCovDYDY=%e\n", rtk->gsi.qr[1] ); /* DY*DY */
		fprintf( fp, "RK_VarCovDYDZ=%e\n", rtk->gsi.qr[4] ); /* DY*DZ */
		fprintf( fp, "RK_VarCovDZDZ=%e\n", rtk->gsi.qr[2] ); /* DZ*DZ */
		fprintf( fp, "RK_RMS=%f\n", SQRT(rtk->gsi.var) ); /* RMS */
		fprintf( fp, "RK_RDOP=%f\n", rtk->gsi.rdop ); /* RDOP */

		/* 単一基線解析結果ファイル：観測記簿・キネマティック用のファイルクローズ */
		fclose( fp );
	}

	return( 1 );
}

extern void readpreceph(char **infile, int n, const prcopt_t *prcopt,
						nav_t *nav, sbs_t *sbs, lex_t *lex)
{
    seph_t seph0={0};
    int i;
    char *ext;
    
    trace(3,"readpreceph: n=%d\n",n);
    
    nav->ne=nav->nemax=0;
    nav->nc=nav->ncmax=0;
    sbs->n =sbs->nmax =0;
    lex->n =lex->nmax =0;
    nav->ns=nav->nsmax=0;

    /* read precise ephemeris files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        readsp3(infile[i],nav,0);
    }
    /* read precise clock files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
		readrnxc(infile[i],nav);
    }
    /* read sbas message files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        sbsreadmsg(infile[i],prcopt->sbassatsel,sbs);
    }
    /* read lex message files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        lexreadmsg(infile[i],0,lex);
    }
    /* allocate sbas ephemeris */
    nav->ns=nav->nsmax=NSATSBS*2;
    if (!(nav->seph=(seph_t *)malloc(sizeof(seph_t)*nav->ns))) {
         showmsg("error : sbas ephem memory allocation");
         trace(1,"error : sbas ephem memory allocation");
         return;
    }
    for (i=0;i<nav->ns;i++) nav->seph[i]=seph0;
    
    /* set rtcm file and initialize rtcm struct */
    rtcm_file[0]=rtcm_path[0]='\0'; fp_rtcm=NULL;
    
    for (i=0;i<n;i++) {
        if ((ext=strrchr(infile[i],'.'))&&
            (!strcmp(ext,".rtcm3")||!strcmp(ext,".RTCM3"))) {
            strcpy(rtcm_file,infile[i]);
        init_rtcm(&rtcm);
        break;
        }
    }
}


/* output reference position -------------------------------------------------*/
static void outrpos(FILE *fp, const double *r, const solopt_t *opt)
{
    double pos[3],dms1[3],dms2[3];
    const char *sep=opt->sep;
    
    trace(3,"outrpos :\n");
    
    if (opt->posf==SOLF_LLH||opt->posf==SOLF_ENU) {
        ecef2pos(r,pos);
        if (opt->degf) {
            deg2dms(pos[0]*R2D,dms1);
            deg2dms(pos[1]*R2D,dms2);
            fprintf(fp,"%3.0f%s%02.0f%s%08.5f%s%4.0f%s%02.0f%s%08.5f%s%10.4f",
                    dms1[0],sep,dms1[1],sep,dms1[2],sep,dms2[0],sep,dms2[1],
                    sep,dms2[2],sep,pos[2]);
        }
        else {
            fprintf(fp,"%13.9f%s%14.9f%s%10.4f",pos[0]*R2D,sep,pos[1]*R2D,
                    sep,pos[2]);
        }
    }
    else if (opt->posf==SOLF_XYZ) {
        fprintf(fp,"%14.4f%s%14.4f%s%14.4f",r[0],sep,r[1],sep,r[2]);
    }
    else if (opt->posf==SOLF_RVA) {
        fprintf(fp,"%14.4f%s%14.4f%s%14.4f",r[0],sep,r[1],sep,r[2]);
    }
}
/* output header -------------------------------------------------------------*/
static void outheader(FILE *fp, char **file, int n, const prcopt_t *popt,
                      const solopt_t *sopt)
{
    const char *s1[]={"GPST","UTC","JST"};
    gtime_t ts,te;
    double t1,t2;
    int i,j,w1,w2;
    char s2[32],s3[32];
    
    trace(3,"outheader: n=%d\n",n);
    
    if (sopt->posf==SOLF_NMEA) return;
    
    if (sopt->outhead) {
        if (!*sopt->prog) {
            fprintf(fp,"%s program   : GSIPOST ver.%s\n",COMMENTH,VER_RTKLIB);
        }
        else {
            fprintf(fp,"%s program   : %s\n",COMMENTH,sopt->prog);
        }
        for (i=0;i<n;i++) {
            fprintf(fp,"%s inp file  : %s\n",COMMENTH,file[i]);
        }
        for (i=0;i<obss.n;i++)    if (obss.data[i].rcv==1) break;
        for (j=obss.n-1;j>=0;j--) if (obss.data[j].rcv==1) break;
        if (j<i) {fprintf(fp,"\n%s no rover obs data\n",COMMENTH); return;}
        ts=obss.data[i].time;
        te=obss.data[j].time;
        t1=time2gpst(ts,&w1);
        t2=time2gpst(te,&w2);
        if (sopt->times>=1) ts=gpst2utc(ts);
        if (sopt->times>=1) te=gpst2utc(te);
        if (sopt->times==2) ts=timeadd(ts,9*3600.0);
        if (sopt->times==2) te=timeadd(te,9*3600.0);
        time2str(ts,s2,1);
        time2str(te,s3,1);
        fprintf(fp,"%s obs start : %s %s (week%04d %8.1fs)\n",COMMENTH,s2,s1[sopt->times],w1,t1);
        fprintf(fp,"%s obs end   : %s %s (week%04d %8.1fs)\n",COMMENTH,s3,s1[sopt->times],w2,t2);
    }
    if (sopt->outopt) {
        outprcopt(fp,popt);
    }
    if (PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED&&popt->mode!=PMODE_MOVEB) {
        fprintf(fp,"%s ref pos   :",COMMENTH);
        outrpos(fp,popt->rb,sopt);
        fprintf(fp,"\n");
    }
    if (sopt->outhead||sopt->outopt) fprintf(fp,"%s\n",COMMENTH);

    outsolhead(fp,popt,sopt);
}

/* process positioning -------------------------------------------------------*/
static void procpos(FILE *fp, const prcopt_t *popt, const solopt_t *sopt,
					int mode)
{
	gtime_t time={0};
	sol_t sol={{0}};
	rtk_t rtk={{0}};
	obsd_t obs[MAXOBS]={{0}};
	gtime_t t0={{0}};
	double rb[3]={0};
	char buff[MAXSOLMSG+1]={0};
	int i,nobs,n,solstatic,pri[]={0,1,2,3,4,5,1,6};
	int ini=0;

	trace(3,"procpos : mode=%d\n",mode);

	solstatic=sopt->solstatic&&
			  (popt->mode==PMODE_STATIC||popt->mode==PMODE_PPP_STATIC);

	rtkinit(&rtk,popt);
	rtcm_path[0]='\0';
	while ((nobs=inputobs(obs,rtk.sol.stat,popt))>=0) {
		t0=obs[0].time;
		if(0==ini) {
			setsnxtstart(&t0);
			ini=1;
		}
		/* exclude satellites */
        for (i=n=0;i<nobs;i++) {
            if ((satsys(obs[i].sat,NULL)&popt->navsys)&&
				popt->exsats[obs[i].sat-1]!=1) obs[n++]=obs[i];
		}
		if (n<=0) {
			continue;
		}
		if (!rtkpos(&rtk,obs,n,&navs)) continue;

//		if((popt->mode!=PMODE_STATIC) && (popt->mode!=PMODE_PPP_STATIC)) {
			outsolsnx(&rtk, obs);
//		}

		if (mode==0) { /* forward/backward */
			if (!solstatic) {
				outsol(fp,&rtk.sol,rtk.rb,popt,sopt);
			}
			else if (time.time==0||pri[rtk.sol.stat]<=pri[sol.stat]) {
                sol=rtk.sol;
                for (i=0;i<3;i++) rb[i]=rtk.rb[i];
                if (time.time==0||timediff(rtk.sol.time,time)<0.0) {
                    time=rtk.sol.time;
				}
			}
        }
        else if (!revs) { /* combined-forward */
            if (isolf>=nepoch) return;
            solf[isolf]=rtk.sol;
            for (i=0;i<3;i++) rbf[i+isolf*3]=rtk.rb[i];
            isolf++;
        }
        else { /* combined-backward */
            if (isolb>=nepoch) return;
            solb[isolb]=rtk.sol;
            for (i=0;i<3;i++) rbb[i+isolb*3]=rtk.rb[i];
			isolb++;
		}
		outsolcrinexrec(&obs[0].time,stas[0].name,rtk.sol.dtr[0]);
	}
	if (mode==0&&solstatic&&time.time!=0.0) {
		sol.time=time;
		outsol(fp,&sol,rb,popt,sopt);
	}
	/*単一基線解析結果ファイルを出力するコードを追加する。*/
	if(sopt->sbresout==SBROUTOPT_NEW) {
		if (!((popt->soltype==2)&&(revs==0))) {
			if (((popt->mode==PMODE_KINEMA)||(popt->mode==PMODE_STATIC))
					&&(popt->antestmode==ANTEST_MODE_NONE)) {
				if (!(save_sbres(popt,&rtk))) {
					showmsg("error : save relative positioning result.");
					trace(1,"error : save relative positioning result.\n");
				}
			}
		}
	}

	if((popt->isb==ISBOPT_EST) || (popt->isb==ISBOPT_EST_P) || (popt->isb==ISBOPT_EST_L) || (popt->isb==ISBOPT_EST_0M)) {
		outisbtable(popt, sopt, &rtk.sol, NULL);
	}

	outgl2table(popt, sopt, &rtk.sol, NULL);
	setsnxtend(&obs[0].time);
//	if(popt->mode == PMODE_PPP_STATIC) {
//		outsolsnxs(&rtk, obs);
//	}
	rtkfree(&rtk);
}


/* validation of combined solutions ------------------------------------------*/
static int valcomb(const sol_t *solf, const sol_t *solb)
{
    double dr[3],var[3];
    int i;
    char tstr[32];
    
    trace(3,"valcomb :\n");
    
    /* compare forward and backward solution */
    for (i=0;i<3;i++) {
        dr[i]=solf->rr[i]-solb->rr[i];
        var[i]=solf->qr[i]+solb->qr[i];
    }
    for (i=0;i<3;i++) {
        if (dr[i]*dr[i]<=16.0*var[i]) continue; /* ok if in 4-sigma */
        
        time2str(solf->time,tstr,2);
        trace(2,"degrade fix to float: %s dr=%.3f %.3f %.3f std=%.3f %.3f %.3f\n",
              tstr+11,dr[0],dr[1],dr[2],SQRT(var[0]),SQRT(var[1]),SQRT(var[2]));
        return 0;
    }
    return 1;
}
/* combine forward/backward solutions and output results ---------------------*/
static void combres(FILE *fp, const prcopt_t *popt, const solopt_t *sopt)
{
    gtime_t time={0};
    sol_t sols={{0}},sol={{0}};
    double tt,Qf[9],Qb[9],Qs[9],rbs[3]={0},rb[3]={0};
    int i,j,k,solstatic,pri[]={0,1,2,3,4,5,1,6};
    
    trace(3,"combres : isolf=%d isolb=%d\n",isolf,isolb);
    
    solstatic=sopt->solstatic&&
              (popt->mode==PMODE_STATIC||popt->mode==PMODE_PPP_STATIC);
    
    for (i=0,j=isolb-1;i<isolf&&j>=0;i++,j--) {
        
        if ((tt=timediff(solf[i].time,solb[j].time))<-DTTOL) {
            sols=solf[i];
            for (k=0;k<3;k++) rbs[k]=rbf[k+i*3];
            j++;
        }
        else if (tt>DTTOL) {
            sols=solb[j];
            for (k=0;k<3;k++) rbs[k]=rbb[k+j*3];
            i--;
        }
        else if (solf[i].stat<solb[j].stat) {
            sols=solf[i];
            for (k=0;k<3;k++) rbs[k]=rbf[k+i*3];
        }
        else if (solf[i].stat>solb[j].stat) {
            sols=solb[j];
            for (k=0;k<3;k++) rbs[k]=rbb[k+j*3];
        }
        else {
            sols=solf[i];
            for (k=0;k<3;k++) rbs[k]=rbf[k+i*3];
            sols.time=timeadd(sols.time,-tt/2.0);
            
            if ((popt->mode==PMODE_KINEMA||popt->mode==PMODE_MOVEB)&&
                sols.stat==SOLQ_FIX) {
                
                /* degrade fix to float if validation failed */
                if (!valcomb(solf+i,solb+j)) sols.stat=SOLQ_FLOAT;
            }
            for (k=0;k<3;k++) {
                Qf[k+k*3]=solf[i].qr[k];
                Qb[k+k*3]=solb[j].qr[k];
            }
            Qf[1]=Qf[3]=solf[i].qr[3];
            Qf[5]=Qf[7]=solf[i].qr[4];
            Qf[2]=Qf[6]=solf[i].qr[5];
            Qb[1]=Qb[3]=solb[j].qr[3];
            Qb[5]=Qb[7]=solb[j].qr[4];
            Qb[2]=Qb[6]=solb[j].qr[5];
            
            if (smoother(solf[i].rr,Qf,solb[j].rr,Qb,3,sols.rr,Qs)) continue;
            
            sols.qr[0]=(float)Qs[0];
            sols.qr[1]=(float)Qs[4];
            sols.qr[2]=(float)Qs[8];
            sols.qr[3]=(float)Qs[1];
            sols.qr[4]=(float)Qs[5];
            sols.qr[5]=(float)Qs[2];
        }
        if (!solstatic) {
            outsol(fp,&sols,rbs,popt,sopt);
        }
        else if (time.time==0||pri[sols.stat]<=pri[sol.stat]) {
            sol=sols;
            for (k=0;k<3;k++) rb[k]=rbs[k];
            if (time.time==0||timediff(sols.time,time)<0.0) {
                time=sols.time;
            }
        }
    }
    if (solstatic&&time.time!=0.0) {
        sol.time=time;
        outsol(fp,&sol,rb,popt,sopt);
    }
}
/* read prec ephemeris, sbas data, lex data, tec grid and open rtcm ----------*/


/* free prec ephemeris and sbas data -----------------------------------------*/
static void freepreceph(nav_t *nav, sbs_t *sbs, lex_t *lex)
{
    int i;
    
    trace(3,"freepreceph:\n");
    
    free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
    free(nav->pclk); nav->pclk=NULL; nav->nc=nav->ncmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
    free(sbs->msgs); sbs->msgs=NULL; sbs->n =sbs->nmax =0;
    free(lex->msgs); lex->msgs=NULL; lex->n =lex->nmax =0;
    for (i=0;i<nav->nt;i++) {
        free(nav->tec[i].data);
        free(nav->tec[i].rms );
    }
    free(nav->tec ); nav->tec =NULL; nav->nt=nav->ntmax=0;
    
#ifdef EXTSTEC
    stec_free(nav);
#endif
    
    if (fp_rtcm) fclose(fp_rtcm);
    free_rtcm(&rtcm);
}

/* station position from file ------------------------------------------------*/
extern int getstapos(const char *file, char *name, double *r, char *name3)
{
	FILE *fp;
	char buff[256],sname[256],*p,*q;
	double pos[3];
	double xyz[3]={0};
	char stf[512]={0};
	int i;

	trace(3,"getstapos: file=%s name=%s\n",file,name);

	if (!(fp=fopen(file,"r"))) {
		trace(1,"station position file open error: %s\n",file);
		return 0;
	}
	if(!fgets(buff,sizeof(buff),fp)) {
		trace(2,"no station position: %s %s\n",name,file);
	}
	if (0==strncmp(buff,"SITE",4)) {
		fgets(buff,sizeof(buff),fp);
		while (fgets(buff,sizeof(buff),fp)) {
			if(strlen(buff) < 250) continue;

			strncpy(sname, buff, 4);
			for (p=sname,q=name;*p&&*q;p++,q++) {
				if (toupper((int)*p)!=toupper((int)*q)) break;
			}
			if ((!*p)&&(!*q)) {
				strncpy(name3,buff+6,16);
				name3[16]='\0';
				trim(name3);
				sprintf(stf, "%13.13s\0", buff+140 );
				xyz[0] = atof(stf);
				sprintf(stf, "%13.13s\0", buff+155 );
				xyz[1] = atof(stf);
				sprintf(stf, "%13.13s\0", buff+170 );
				xyz[2] = atof(stf);
				for(i=0;i<3;++i) r[i]=xyz[i];
				fclose(fp);
				return 1;
			}
		}
	}
	else {
		while (fgets(buff,sizeof(buff),fp)) {
			if ((p=strchr(buff,'%'))) *p='\0';

			if (sscanf(buff,"%lf %lf %lf %s %s",pos,pos+1,pos+2,sname,name3)<5) continue;

			for (p=sname,q=name;*p&&*q;p++,q++) {
				if (toupper((int)*p)!=toupper((int)*q)) break;
			}
			if ((!*p)&&(!*q)) {
				pos[0]*=D2R;
				pos[1]*=D2R;
				pos2ecef(pos,r);
				fclose(fp);
				return 1;
			}
		}
	}
	fclose(fp);
	trace(1,"no station position: %s %s\n",name,file);
	return 0;
}

/* average of single position ------------------------------------------------*/
static int avepos(double *ra, int rcv, const obs_t *obs, const nav_t *nav,
                  const prcopt_t *opt)
{
	obsd_t data[MAXOBS]={{0}};
	gtime_t ts={0};
    sol_t sol={{0}};
    int i,j,n=0,m,iobs;
    char msg[128];
    
    trace(3,"avepos: rcv=%d obs.n=%d\n",rcv,obs->n);
    
    for (i=0;i<3;i++) ra[i]=0.0;
    
    for (iobs=0;(m=nextobsf(obs,&iobs,rcv))>0;iobs+=m) {
        
        for (i=j=0;i<m&&i<MAXOBS;i++) {
            data[j]=obs->data[iobs+i];
            if ((satsys(data[j].sat,NULL)&opt->navsys)&&
                opt->exsats[data[j].sat-1]!=1) j++;
        }
        if (j<=0||!screent(data[0].time,ts,ts,1.0)) continue; /* only 1 hz */
        
        if (!pntpos(data,j,nav,opt,&sol,NULL,NULL,msg)) continue;
        
        for (i=0;i<3;i++) ra[i]+=sol.rr[i];
        n++;
    }
    if (n<=0) {
        trace(1,"no average of base station position\n");
        return 0;
    }
    for (i=0;i<3;i++) ra[i]/=n;
    return 1;
}

/* antenna phase center position ---------------------------------------------*/
static int antpos(prcopt_t *opt, int rcvno, const obs_t *obs, const nav_t *nav,
                  const sta_t *sta, const char *posfile)
{
    double *rr=rcvno==1?opt->ru:opt->rb,del[3],pos[3],dr[3]={0};
    int i,postype=rcvno==1?opt->rovpos:opt->refpos;
    char *name;
	char name3[MAXANT];
    
    trace(3,"antpos  : rcvno=%d\n",rcvno);
    
    if (postype==1) { /* avarage of single position */
        if (!avepos(rr,rcvno,obs,nav,opt)) {
            showmsg("error : station pos computation");
            return 0;
        }
    }
    else if (postype==2) { /* read from position file */
        name=stas[rcvno==1?0:1].name;
		if (!getstapos(posfile,name,rr,name3)) {
            showmsg("error : no position of %s in %s",name,posfile);
            return 0;
        }
		strcpy(stas[rcvno==1?0:1].name3,name3);
    }
    else if (postype==3) { /* get from rinex header */
        if (norm(stas[rcvno==1?0:1].pos,3)<=0.0) {
            showmsg("error : no position in rinex header");
            trace(1,"no position position in rinex header\n");
            return 0;
        }
        /* antenna delta */
        if (stas[rcvno==1?0:1].deltype==0) { /* enu */
            for (i=0;i<3;i++) del[i]=stas[rcvno==1?0:1].del[i];
            del[2]+=stas[rcvno==1?0:1].hgt;
            ecef2pos(stas[rcvno==1?0:1].pos,pos);
            enu2ecef(pos,del,dr);
        }
        else { /* xyz */
            for (i=0;i<3;i++) dr[i]=stas[rcvno==1?0:1].del[i];
        }
        for (i=0;i<3;i++) rr[i]=stas[rcvno==1?0:1].pos[i]+dr[i];
    }
    return 1;
}

/* set antenna parameters ----------------------------------------------------*/
static void setpcv(gtime_t time, prcopt_t *popt, nav_t *nav, const pcvs_t *pcvs,
				   const pcvs_t *pcvr, sta_t *sta)
{
    pcv_t *pcv;
    double pos[3],del[3];
    int i,j,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;
    char id[64];
    
    /* set satellite antenna parameters */
    for (i=0;i<MAXSAT;i++) {
		if (!(satsys(i+1,NULL)&popt->navsys)) continue;
        if (!(pcv=searchpcv(i+1,"",time,pcvs))) {
            satno2id(i+1,id);
            trace(2,"no satellite antenna pcv: %s\n",id);
            continue;
        }
        nav->pcvs[i]=*pcv;
    }
    for (i=0;i<(mode?2:1);i++) {
        if (!strcmp(popt->anttype[i],"*")) { /* set by station parameters */
            strcpy(popt->anttype[i],sta[i].antdes);
            if (sta[i].deltype==1) { /* xyz */
                if (norm(sta[i].pos,3)>0.0) {
                    ecef2pos(sta[i].pos,pos);
                    ecef2enu(pos,sta[i].del,del);
                    for (j=0;j<3;j++) popt->antdel[i][j]=del[j];
                }
            }
            else { /* enu */
                for (j=0;j<3;j++) popt->antdel[i][j]=stas[i].del[j];
            }
        }
        if (!(pcv=searchpcv(0,popt->anttype[i],time,pcvr))) {
            trace(2,"no receiver antenna pcv: %s\n",popt->anttype[i]);
            *popt->anttype[i]='\0';
            continue;
        }
        strcpy(popt->anttype[i],pcv->type);
		popt->pcvr[i]=*pcv;
		sta[i].pcvr=*pcv;
    }
}

/* read ocean tide loading parameters ----------------------------------------*/
static void readotl(prcopt_t *popt, const char *file, const sta_t *sta)
{
    int i,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;
    
    for (i=0;i<(mode?2:1);i++) {
        readblq(file,sta[i].name,popt->odisp[i]);
    }
}

/* write header to output file -----------------------------------------------*/
static int outhead(const char *outfile, char **infile, int n,
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
    outheader(fp,infile,n,popt,sopt);
    
    if (*outfile) fclose(fp);
    
    return 1;
}

/* execute processing session ------------------------------------------------*/
static int execses(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
				   const solopt_t *sopt, const filopt_t *fopt, int flag,
				   char **infile, const int *index, int n, char *outfile)
{
    FILE *fp;
    prcopt_t popt_=*popt;
    char tracefile[1024],statfile[1024],snxfile[1024],clcfile[1024],ionfile[1024];
	fcb_t *fcb;
    int i;
    int mode;
    mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;

    trace(3,"execses : n=%d outfile=%s\n",n,outfile);
    
    /* open debug trace */
    if (flag&&sopt->trace>0) {
        if (*outfile) {
            strcpy(tracefile,outfile);
            strcat(tracefile,".trace");
        }
        else {
            strcpy(tracefile,fopt->trace);
        }
        traceclose();
        traceopen(tracefile);
        tracelevel(sopt->trace);
	}
    /*改修概要*/
	if((popt_.mode==PMODE_PPP_KINEMA) || (popt_.mode==PMODE_PPP_STATIC) || (popt_.mode==PMODE_PPP_FIXED)) {
	    if(popt_.modepppar==PPPAR_FCB){
		    if(!(fcb=loadfcb(((prcopt_t *)popt)))){
			    showmsg("error : read fcb file");
			    trace(1,"fcb read error: %s\n",popt_.mopt.iffcb);
				return 0;
		    }
		    navs.fcb=fcb;
		}
    }
	
    /* read obs and nav data */
	if (!readobsnav(ts,te,ti,infile,index,n,&popt_,&obss,&navs,stas)) return 0;
    
	/* set antenna paramters */
	if (popt_.sateph==EPHOPT_BRDC||popt_.sateph==EPHOPT_PREC||popt_.sateph==EPHOPT_SSRCOM) {
		setpcv(obss.n>0?obss.data[0].time:timeget(),&popt_,&navs,&pcvss,&pcvsr,
               stas);
    }
    /* quarter sysle shift paramters */
    for (i=0;i<(mode?2:1);i++) {
		setL2Csft(popt_.phasshft, popt_.rectype[i], &navs.sfts, stas+i);
    }
    /* dcb */
    for (i=0;i<(mode?2:1);i++) {
	 //	setdcb(&navs, popt_.rectype[i], stas+i);
		setdcb(&navs, stas+i);
    }
    /* ISB */
//	for (i=0;i<(mode?2:1);i++) {
//		setisb(&navs, popt_.rectype[i], stas+i);
////		setisb(&navs, stas[i].rectype, stas+i);
//	//	setisb(&navs, stas+i);
//	}

	/* ISB */
	if     (mode == 0) setisb(&navs, popt_.rectype[0], ""              , stas, NULL  );
	else if(mode == 1) setisb(&navs, popt_.rectype[0], popt_.rectype[1], stas, stas+1);
//	if     (mode == 0) setisb(0, &navs, popt_.rectype[0], ""              , stas, NULL  );
//	else if(mode == 1) setisb(1, &navs, popt_.rectype[0], popt_.rectype[1], stas, stas+1);


    /* read ocean tide loading parameters */
	if (popt_.mode>PMODE_SINGLE&&fopt->blq) {
		readotl(&popt_,fopt->blq,stas);
    }
	/* rover/reference fixed position */
    if (popt_.mode==PMODE_FIXED) {
        if (!antpos(&popt_,1,&obss,&navs,stas,fopt->stapos)) {
            freeobsnav(&obss,&navs);
            return 0;
        }
    }
	else if (PMODE_DGPS<=popt_.mode&&popt_.mode<=PMODE_STATIC) {
		if (!antpos(&popt_,2,&obss,&navs,stas,fopt->stapos)) {
            freeobsnav(&obss,&navs);
			return 0;
		}
    }
    /* open solution statistics */
    if (flag&&sopt->sstat>0) {
        strcpy(statfile,outfile);
        strcat(statfile,".stat");
        rtkclosestat();
		rtkopenstat(statfile,sopt->sstat);
	}
	if((sopt->recclout)
		&& ((popt->mode==PMODE_SINGLE)||(popt->mode==PMODE_PPP_KINEMA)||(popt->mode==PMODE_PPP_STATIC)||(popt->mode==PMODE_PPP_FIXED))) {
		strcpy(clcfile,outfile);
        strcat(clcfile,".clc");
        rtkclosecrinex();
		rtkopencrinex(popt,sopt,NULL,clcfile,stas);
	}
	if((sopt->ionout)
		&& ((popt->mode>=PMODE_DGPS) && (popt->mode<=PMODE_FIXED))) {
		strcpy(ionfile,outfile);
		strcat(ionfile,".ion");
	//	rtkclosecrinex();
//		rtkopencrinex(popt,sopt,NULL,clcfile);

//		fp_rnxion = fopen(ionfile,"w");
//		outrnxobsh_(fp_rnxion,&ropt,nav);
		openrnxobsh_ion(popt,ionfile);
	}
    /* write header to output file */
	if (flag&&!outhead(outfile,infile,n,&popt_,sopt)) {
        freeobsnav(&obss,&navs);
        return 0;
    }
    iobsu=iobsr=isbs=ilex=revs=aborts=0;
    
	if(sopt->possnxout) {
        strcpy(snxfile,outfile);
        strcat(snxfile,".snx");
    //    rtkclosesnx(&te);
		if(!rtkopensnx(&popt_,NULL,stas,snxfile,infile,n)) {
			return 0;
		}
	}
	if (popt_.mode==PMODE_SINGLE||popt_.soltype==0) {
        if ((fp=openfile(outfile))) {
			procpos(fp,&popt_,sopt,0); /* forward */
            fclose(fp);
        }
    }
    else if (popt_.soltype==1) {
        if ((fp=openfile(outfile))) {
            revs=1; iobsu=iobsr=obss.n-1; isbs=sbss.n-1; ilex=lexs.n-1;
			procpos(fp,&popt_,sopt,0); /* backward */
            fclose(fp);
        }
    }
    else { /* combined */
		solf=(sol_t *)malloc(sizeof(sol_t)*nepoch);
        solb=(sol_t *)malloc(sizeof(sol_t)*nepoch);
        rbf=(double *)malloc(sizeof(double)*nepoch*3);
        rbb=(double *)malloc(sizeof(double)*nepoch*3);
        
        if (solf&&solb) {
            isolf=isolb=0;
			procpos(NULL,&popt_,sopt,1); /* forward */
            revs=1; iobsu=iobsr=obss.n-1; isbs=sbss.n-1; ilex=lexs.n-1;
			procpos(NULL,&popt_,sopt,1); /* backward */
            
            /* combine forward/backward solutions */
            if (!aborts&&(fp=openfile(outfile))) {
                combres(fp,&popt_,sopt);
                fclose(fp);
            }
        }
        else showmsg("error : memory allocation");
        free(solf);
        free(solb);
        free(rbf);
        free(rbb);
    }
    /* free obs and nav data */
	freeobsnav(&obss,&navs);

	rtkclosesnx(NULL);
	rtkclosecrinex();
    
	return aborts?1:0;
}

/* execute processing session for each rover ---------------------------------*/
static int execses_r(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                     const solopt_t *sopt, const filopt_t *fopt, int flag,
                     char **infile, const int *index, int n, char *outfile,
                     const char *rov)
{
    gtime_t t0={0};
    int i,stat=0;
    char *ifile[MAXINFILE],ofile[1024],*rov_,*p,*q,s[64]="";
    
    trace(3,"execses_r: n=%d outfile=%s\n",n,outfile);
    
    for (i=0;i<n;i++) if (strstr(infile[i],"%r")) break;
    
    if (i<n) { /* include rover keywords */
        if (!(rov_=(char *)malloc(strlen(rov)+1))) return 0;
		strcpy(rov_,rov);
        
        for (i=0;i<n;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                free(rov_); for (;i>=0;i--) free(ifile[i]);
                return 0;
            }
        }
        for (p=rov_;;p=q+1) { /* for each rover */
            if ((q=strchr(p,' '))) *q='\0';
            
            if (*p) {
                strcpy(proc_rov,p);
                if (ts.time) time2str(ts,s,0); else *s='\0';
                if (checkbrk("reading    : %s",s)) {
                    stat=1;
                    break;
                }
                for (i=0;i<n;i++) reppath(infile[i],ifile[i],t0,p,"");
                reppath(outfile,ofile,t0,p,"");
                
                /* execute processing session */
                stat=execses(ts,te,ti,popt,sopt,fopt,flag,ifile,index,n,ofile);
            }
            if (stat==1||!q) break;
        }
        free(rov_); for (i=0;i<n;i++) free(ifile[i]);
    }
    else {
        /* execute processing session */
        stat=execses(ts,te,ti,popt,sopt,fopt,flag,infile,index,n,outfile);
    }
    return stat;
}
/* execute processing session for each base station --------------------------*/
static int execses_b(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
					 const solopt_t *sopt, const filopt_t *fopt, int flag,
					 char **infile, const int *index, int n, char *outfile,
					 const char *rov, const char *base)
{
    gtime_t t0={0};
    int i,stat=0;
    char *ifile[MAXINFILE],ofile[1024],*base_,*p,*q,s[64];
    
    trace(3,"execses_b: n=%d outfile=%s\n",n,outfile);
    
    /* read prec ephemeris and sbas data */
	readpreceph(infile,n,popt,&navs,&sbss,&lexs);

#ifndef RTKPOS4PCV
	/*複数基線解析？*/
	if(popt->mode==PMODE_MULTI){
		/*EPHOPT_PRECまたはEPHOPT_SSRCOM*/
		if(popt->sateph==EPHOPT_PREC){
			/*複数基線解析*/
			stat=mbstatic(ts,te,ti,(prcopt_t *)popt,(solopt_t *)sopt,(filopt_t *)fopt,
					flag,infile,(int *)index,n,(char *)rov,(char *)base,
					&navs,&pcvsr,&pcvss);
			if(stat== 0) {
				showmsg("done");
			}
			return 0;
		}else{
			/*精密暦,SBASなどメモリ解放*/
			showmsg( "mbs start err :set ephemeris option: precise" );
			trace( 1,"mbs start err :set ephemeris option: precise\n" );

			freepreceph(&navs,&sbss,&lexs);
			return 0;
		}
	}
#endif

    for (i=0;i<n;i++) if (strstr(infile[i],"%b")) break;
    
    if (i<n) { /* include base station keywords */
        if (!(base_=(char *)malloc(strlen(base)+1))) {
            freepreceph(&navs,&sbss,&lexs);
            return 0;
        }
        strcpy(base_,base);
        
        for (i=0;i<n;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                free(base_); for (;i>=0;i--) free(ifile[i]);
                freepreceph(&navs,&sbss,&lexs);
                return 0;
            }
        }
        for (p=base_;;p=q+1) { /* for each base station */
            if ((q=strchr(p,' '))) *q='\0';
            
            if (*p) {
                strcpy(proc_base,p);
                if (ts.time) time2str(ts,s,0); else *s='\0';
                if (checkbrk("reading    : %s",s)) {
                    stat=1;
                    break;
                }
                for (i=0;i<n;i++) reppath(infile[i],ifile[i],t0,"",p);
                reppath(outfile,ofile,t0,"",p);
                
                stat=execses_r(ts,te,ti,popt,sopt,fopt,flag,ifile,index,n,ofile,rov);
            }
            if (stat==1||!q) break;
        }
        free(base_); for (i=0;i<n;i++) free(ifile[i]);
    }
    else {
        stat=execses_r(ts,te,ti,popt,sopt,fopt,flag,infile,index,n,outfile,rov);
    }
    /* free prec ephemeris and sbas data */
    freepreceph(&navs,&sbss,&lexs);
    
    return stat;
}




/* post-processing positioning -------------------------------------------------
* post-processing positioning
* args   : gtime_t ts       I   processing start time (ts.time==0: no limit)
*        : gtime_t te       I   processing end time   (te.time==0: no limit)
*          double ti        I   processing interval  (s) (0:all)
*          double tu        I   processing unit time (s) (0:all)
*          prcopt_t *popt   I   processing options
*          solopt_t *sopt   I   solution options
*          filopt_t *fopt   I   file options
*          char   **infile  I   input files (see below)
*          int    n         I   number of input files
*          char   *outfile  I   output file ("":stdout, see below)
*          char   *rov      I   rover id list        (separated by " ")
*          char   *base     I   base station id list (separated by " ")
* return : status (0:ok,0>:error,1:aborted)
* notes  : input files should contain observation data, navigation data, precise 
*          ephemeris/clock (optional), sbas log file (optional), ssr message
*          log file (optional) and tec grid file (optional). only the first 
*          observation data file in the input files is recognized as the rover
*          data.
*
*          the type of an input file is recognized by the file extention as ]
*          follows:
*              .sp3,.SP3,.eph*,.EPH*: precise ephemeris (sp3c)
*              .sbs,.SBS,.ems,.EMS  : sbas message log files (rtklib or ems)
*              .lex,.LEX            : qzss lex message log files
*              .rtcm3,.RTCM3        : ssr message log files (rtcm3)
*              .*i,.*I              : tec grid files (ionex)
*              others               : rinex obs, nav, gnav, hnav, qnav or clock
*
*          inputs files can include wild-cards (*). if an file includes
*          wild-cards, the wild-card expanded multiple files are used.
*
*          inputs files can include keywords. if an file includes keywords,
*          the keywords are replaced by date, time, rover id and base station
*          id and multiple session analyses run. refer reppath() for the
*          keywords.
*
*          the output file can also include keywords. if the output file does
*          not include keywords. the results of all multiple session analyses
*          are output to a single output file.
*
*          ssr corrections are valid only for forward estimation.
*-----------------------------------------------------------------------------*/
extern int postpos(gtime_t ts, gtime_t te, double ti, double tu,
				   const prcopt_t *popt, const solopt_t *sopt,
				   const filopt_t *fopt, char **infile, int n, char *outfile,
				   const char *rov, const char *base)
{
    gtime_t tts,tte,ttte;
    double tunit,tss;
    int i,j,k,nf,stat=0,week,flag=1,index[MAXINFILE]={0};
    char *ifile[MAXINFILE],ofile[1024],*ext;


	if(!checkopts(popt,sopt,&fopt))
	{
        return -3;
	}

    trace(3,"postpos : ti=%.0f tu=%.0f n=%d outfile=%s\n",ti,tu,n,outfile);
    
    /* 単一基線結果ファイル出力先をメイン画面から取得 */
	strcpy( fsb, outfile );

    /* open processing session */
    if (!openses(popt,sopt,fopt,&navs,&pcvss,&pcvsr)) return -1;
    
    if (ts.time!=0&&te.time!=0&&tu>=0.0) {
        if (timediff(te,ts)<0.0) {
            showmsg("error : no period");
			closeses(&navs,&pcvss,&pcvsr);
            return 0;
        }
        for (i=0;i<MAXINFILE;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                for (;i>=0;i--) free(ifile[i]);
                closeses(&navs,&pcvss,&pcvsr);
                return -1;
            }
        }
        if (tu==0.0||tu>86400.0*MAXPRCDAYS) tu=86400.0*MAXPRCDAYS;
        settspan(ts,te);
        tunit=tu<86400.0?tu:86400.0;
        tss=tunit*(int)floor(time2gpst(ts,&week)/tunit);
        
        for (i=0;;i++) { /* for each periods */
            tts=gpst2time(week,tss+i*tu);
            tte=timeadd(tts,tu-DTTOL);
            if (timediff(tts,te)>0.0) break;
            if (timediff(tts,ts)<0.0) tts=ts;
            if (timediff(tte,te)>0.0) tte=te;
            
            strcpy(proc_rov ,"");
            strcpy(proc_base,"");
            if (checkbrk("reading    : %s",time_str(tts,0))) {
                stat=1;
                break;
            }
            for (j=k=nf=0;j<n;j++) {
                
                ext=strrchr(infile[j],'.');
                
                if (ext&&(!strcmp(ext,".rtcm3")||!strcmp(ext,".RTCM3"))) {
                    strcpy(ifile[nf++],infile[j]);
                }
                else {
                    /* include next day precise ephemeris or rinex brdc nav */
                    ttte=tte;
                    if (ext&&(!strcmp(ext,".sp3")||!strcmp(ext,".SP3")||
                              !strcmp(ext,".eph")||!strcmp(ext,".EPH"))) {
                        ttte=timeadd(ttte,3600.0);
                    }
                    else if (strstr(infile[j],"brdc")) {
                        ttte=timeadd(ttte,7200.0);
                    }
                    nf+=reppaths(infile[j],ifile+nf,MAXINFILE-nf,tts,ttte,"","");
                }
                while (k<nf) index[k++]=j;
                
                if (nf>=MAXINFILE) {
                    trace(2,"too many input files. trancated\n");
                    break;
                }
            }
            if (!reppath(outfile,ofile,tts,"","")&&i>0) flag=0;
            
            /* execute processing session */
			stat=execses_b(tts,tte,ti,popt,sopt,fopt,flag,ifile,index,nf,ofile,rov,base);
            
            if (stat==1) break;
        }
        for (i=0;i<MAXINFILE;i++) free(ifile[i]);
    }
    else if (ts.time!=0) {
        for (i=0;i<n&&i<MAXINFILE;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                for (;i>=0;i--) free(ifile[i]); return -1;
            }
            reppath(infile[i],ifile[i],ts,"","");
            index[i]=i;
        }
        reppath(outfile,ofile,ts,"","");
        
        /* execute processing session */
        stat=execses_b(ts,te,ti,popt,sopt,fopt,1,ifile,index,n,ofile,rov,base);
        
        for (i=0;i<n&&i<MAXINFILE;i++) free(ifile[i]);
    }
	else {
        for (i=0;i<n;i++) index[i]=i;
        
        /* execute processing session */
		stat=execses_b(ts,te,ti,popt,sopt,fopt,1,infile,index,n,outfile,rov,base);
    }
    /* close processing session */
	closeses(&navs,&pcvss,&pcvsr);
    
    return stat;
}
