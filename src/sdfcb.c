/*------------------------------------------------------------------------------
* sdfcb.c : SD-FCB estimation
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
#include <stdlib.h>
#include <math.h>
#include "rtklib.h"

/* 一重差パス組み合わせの作成 */
void pairsdmv(mbs_t *mbs, fcb_t *fcb)
{
	/* 処理フロー No.0 */
	/* 変数定義 */
	int i,j;		/* ループカウンタ */
	double d;       /* Melbourne-Wubbena線型結合一重差 */
	gtime_t ti, tj; /* passデータエポック */
	sdd_t sdd={{0}};/* 一重差データ */
	pass_t pass=mbs->pass;

	/* 処理フロー No.1 */
	/* 衛星ペアループ */
	for(i=0;i<pass.n -1;i++)
	{
		if(satsys(pass.data[i].sat,NULL)!=SYS_GPS) continue;

		for(j=i+1;j<pass.n -1;j++)
		{
			/* GPS以外の衛星は対象外 */
			if(satsys(pass.data[j].sat,NULL)!=SYS_GPS) continue;
			/* 同一衛星の一重差は計算しない（当然） */
			if( pass.data[i].sat == pass.data[j].sat ) continue;
			/* パス（観測が連続している区間）が重複しない衛星ペアは対象外 */
			if(0==check_oltime(&pass.data[i].ts,&pass.data[i].te,&pass.data[j].ts,&pass.data[j].te,0.0,&ti,&tj)) continue;
			/* 局が異なるペアの一重差は計算しない（配列dataは局、衛星、時刻でソート済み） */
			if(pass.data[i].sta != pass.data[j].sta)
			{
				break;
			}


			/* 処理フロー No.2 */
			/* 一重差諸量の計算 */
			sdd.s1 = pass.data[i].sat;
			sdd.s2 = pass.data[j].sat;
			sdd.r = pass.data[i].sta;
				/* Melbourne-Wubbena線型結合(L1-L2) */
			d = pass.data[i].mw -pass.data[j].mw;
			sdd.mwv = pass.data[i].mwv+pass.data[j].mwv;/* MW分散(cycle) */
			sdd.mwf = fmod(d,1.0);						/* MW少数部分[0.0:1.0)(cycle) */
			if(sdd.mwf<0) sdd.mwf+=1.0;
			sdd.mwi = floor(d);			/* MW整数部分(cycle) */
				/* 電離層フリー線形結合(L1-L2) */
			sdd.lc = pass.data[i].lcamb - pass.data[j].lcamb;	/* LCアンビギュイティ(m) */
			sdd.lcv = pass.data[i].lcambv + pass.data[j].lcambv;/* LCアンビギュイティ分散(m) */
			ti = pass.data[i].ts;
			tj = pass.data[j].ts;
			sdd.ts = timediff(ti,tj)>0.0?ti:tj;
			ti = pass.data[i].te;
			tj = pass.data[j].te;
			sdd.te = timediff(ti,tj)>0.0?ti:tj;
			sdd.exc = mbs->stas.sta[pass.data[i].sta].cctype;

			/* 処理フロー No.3 */
			/* 一重差データの追加 */
			addsdpass(fcb,&sdd);
		}
	}
}

///* 一重差パス組み合わせの作成 */
//void pairsdmv(pass_t pass, fcb_t *fcb)
//{
//	/* 処理フロー No.0 */
//	/* 変数定義 */
//	int i,j;		/* ループカウンタ */
//	double d;       /* Melbourne-Wubbena線型結合一重差 */
//	gtime_t ti, tj; /* passデータエポック */
//	sdd_t sdd={{0}};/* 一重差データ */
//
//	/* 処理フロー No.1 */
//	/* 衛星ペアループ */
//	for(i=0;i<pass.n -1;i++)
//	{
//		if(satsys(pass.data[i].sat,NULL)!=SYS_GPS) continue;
//
//		for(j=i+1;j<pass.n -1;j++)
//		{
//			/* GPS以外の衛星は対象外 */
//			if(satsys(pass.data[j].sat,NULL)!=SYS_GPS) continue;
//			/* 同一衛星の一重差は計算しない（当然） */
//			if( pass.data[i].sat == pass.data[j].sat ) continue;
//			/* パス（観測が連続している区間）が重複しない衛星ペアは対象外 */
//			if(0==check_oltime(&pass.data[i].ts,&pass.data[i].te,&pass.data[j].ts,&pass.data[j].te,0.0,&ti,&tj)) continue;
//			/* 局が異なるペアの一重差は計算しない（配列dataは局、衛星、時刻でソート済み） */
//			if(pass.data[i].sta != pass.data[j].sta)
//			{
//				break;
//			}
//
//
//			/* 処理フロー No.2 */
//			/* 一重差諸量の計算 */
//			sdd.s1 = pass.data[i].sat;
//			sdd.s2 = pass.data[j].sat;
//			sdd.r = pass.data[i].sta;
//				/* Melbourne-Wubbena線型結合(L1-L2) */
//			d = pass.data[i].mw -pass.data[j].mw;
//			sdd.mwv = pass.data[i].mwv+pass.data[j].mwv;/* MW分散(cycle) */
//			sdd.mwf = fmod(d,1.0);						/* MW少数部分[0.0:1.0)(cycle) */
//			if(sdd.mwf<0) sdd.mwf+=1.0;
//			sdd.mwi = floor(d);			/* MW整数部分(cycle) */
//				/* 電離層フリー線形結合(L1-L2) */
//			sdd.lc = pass.data[i].lcamb - pass.data[j].lcamb;	/* LCアンビギュイティ(m) */
//			sdd.lcv = pass.data[i].lcambv + pass.data[j].lcambv;/* LCアンビギュイティ分散(m) */
//			ti = pass.data[i].ts;
//			tj = pass.data[j].ts;
//			sdd.ts = timediff(ti,tj)>0.0?ti:tj;
//			ti = pass.data[i].te;
//			tj = pass.data[j].te;
//			sdd.te = timediff(ti,tj)>0.0?ti:tj;
//
//			/* 処理フロー No.3 */
//			/* 一重差データの追加 */
//			addsdpass(fcb,&sdd);
//		}
//	}
//}


/* WL-SD-FCB 平均値及び分散の算出 */
void calwlsdfcb(fcb_t *fcb/*, prcopt *popt*/)
{
	/* 処理フロー No.0 */
	/* 変数初期化 */
	int ni, nj, nij;/* 配列インデックス */
	int i;			/* ループカウンタ */
	for(i=0;i<NSATGPS*(NSATGPS-1)/2;i++){
		fcb->wlfcb[i].n = 0;
		fcb->wlfcb[i].ave = 0.0;
		fcb->wlfcb[i].var = 0.0;
	}

	/* 処理フロー No.1 */
	/* MW結合一重差の小数部の平均と分散の算出(全局で計算した統計値) */
	for(i=0;i<fcb->n;i++){
		if ( fcb->sdd[i].exc == 1 ) continue;
		ni=fcb->sdd[i].s1;
		nj=fcb->sdd[i].s2;
		nij=(ni-1)*(2*NSATGPS-ni)/2+nj-ni-1;
		fcb->wlfcb[nij].ave += fcb->sdd[i].mwf;
		fcb->wlfcb[nij].var += pow(fcb->sdd[i].mwf,2);
		fcb->wlfcb[nij].n ++;
		/*debug*///printf("debug: r=%d s1=%d s2=%d mwf=%8.6f\n",fcb->sdd[i].r,fcb->sdd[i].s1,fcb->sdd[i].s2,fcb->sdd[i].mwf);
	}
	for(i=0;i<NSATGPS*(NSATGPS-1)/2;i++){
		if (fcb->wlfcb[i].n == 0 ) continue;
		fcb->wlfcb[i].ave /= (double)fcb->wlfcb[i].n;
		fcb->wlfcb[i].var /= (double)fcb->wlfcb[i].n;
		fcb->wlfcb[i].var -= pow(fcb->wlfcb[i].ave,2);
	}
}

/* WLアンビギュイティFIX判定 */
int fixwlsdfcb(fcb_t *fcb, prcopt_t *popt)
{
	/* 処理フロー No.0 */
	/* 変数初期化 */
	int fixflg=0;/* FIX判定フラグ FIX率80％? 0:OK 1:NG  */
	int ni, nj, nij;/* 配列インデックス */
	int i, j;			/* ループカウンタ */
	double nw, nwfix, nwv;
	double p;
	int nrfix=0, nr=0, nf=0;

	/* 処理フロー No.1 */
	/* 一重差データ整数化判定 */
	for( i=0;i<fcb->n;i++){
		if ( fcb->sdd[i].exc == 1 ) continue;
		ni=fcb->sdd[i].s1;
		nj=fcb->sdd[i].s2;
		nij=(ni-1)*(2*NSATGPS-ni)/2+nj-ni-1;
		nw = fcb->sdd[i].mwi + fcb->sdd[i].mwf - fcb->wlfcb[nij].ave;  /* ...(10) */
		nwfix = floor(nw+0.5);
		nwv = fcb->sdd[i].mwv + fcb->wlfcb[nij].var;
		if ( fabs(nwv) < 1.0e-30 )
		{   /* TODO ゼロ割りエラー出力
			showmsg("error: hogehoge");
			trace(1,"error: hogehoge");
			*/
			fcb->sdd[i].exc = 1;
			continue;
		}
		fcb->sdd[i].mwfix = check_intdd( popt->mopt.minconfw , nw, nwfix, nwv, &p);
		if(fcb->sdd[i].mwfix) nrfix++;
		else 
		{
			//printf("debug: r=%d s1=%d s2=%d f=%8.6f fave=%8.6f\n",fcb->sdd[i].r,fcb->sdd[i].s1,fcb->sdd[i].s2,fcb->sdd[i].mwf,fcb->wlfcb[nij].ave);
		}
		nr ++;
		/*
		*  局毎にFIX 率を判定。80%より少なければ棄却フラグをたてて、
		*  FIX 判定値をNG にセット
		*/
		if ( i+1==fcb->n || fcb->sdd[i].r != fcb->sdd[i+1].r )
		{
			if ( (double)nrfix / (double)nr < 0.8 )
			{
				for(j=i;j>i-nr;j--) fcb->sdd[j].exc = 1;
				fixflg = 1;
			}
			else
			{
				nf ++; /* 80%以上FIX した局数をカウント */
			}
			nrfix = 0;
			nr = 0;
		}
	}

	/* 処理フロー No.2 */
	/* FIX局数判定 */
	if ( nf == 0 ) fixflg = 9;

	/* 処理フロー No.3 */
	/* 終了処理 */
	return fixflg;
}

/* パス毎のNL-SD-FCBの算出 */
void calnlsdfcb(fcb_t *fcb/*, prcopt_t popt*/)
{
	/* 処理フロー No.0 */
	/* 変数初期化 */
	double b1k;
	double fai1, fai1v;
	int fain;
	int i;
	int j=0,k=0;
	gtime_t ti,tj;
	double nw;
	int ni, nj, nij;

	/* 処理フロー No.1 */
	/* 処理判定 */
	for (i=0;i<fcb->n;i++) {


		if ( fcb->sdd[i].exc == 1 ) continue;
		if ( fcb->sdd[i].mwfix == 0 ) continue;
		
		fcb->sdd[i].b1fix = 1;

		/* 処理フロー No.2 */
		/* b1kの算出 */
		ni=fcb->sdd[i].s1;
		nj=fcb->sdd[i].s2;
		nij=(ni-1)*(2*NSATGPS-ni)/2+nj-ni-1;
		nw = fcb->sdd[j].mwi + fcb->sdd[j].mwf - fcb->wlfcb[nij].ave;
		b1k = (FREQ1+FREQ2)/FREQ1 * fcb->sdd[i].lc - FREQ2/(FREQ1-FREQ2) * nw;

		/* 処理フロー No.3 */
		/* b1kの整数・小数部分離 */
		fcb->sdd[i].b1i = b1k>=0?floor(b1k+0.5):ceil(b1k-0.5);
		fcb->sdd[i].b1f = b1k - fcb->sdd[i].b1i;
	}
	/* 処理フロー No.4 */
	/* 終了 */
	return;
}

/* NL-SD-FCBのタイムテーブル作成 */
void makenlfcbtbl(fcb_t *fcb,double tintv)
{
	/* 処理フロー No.0 */
	/* 変数初期化 */
	int ni, nj;
	int nij; /* 衛星ペア インデックス */
	//double tintv;
	gtime_t ts, te, tsdum, tedum;
	//int nlcnt[NSATGPS*(NSATGPS-1)/2][fcb->ntfcb]={0};
	int **nlcnt;
	int i, j, k; /* ループカウンタ */
	int intnum; /* NL-FCB区間数 */
	int maxpnum;  /* 衛星ペア総数 */
	intnum=fcb->ntfcb;
	maxpnum=NSATGPS*(NSATGPS-1)/2;

	nlcnt = (int**)malloc(sizeof(int*)*maxpnum);
	for (i=0;i<maxpnum;i++)
	{
	   nlcnt[i] = (int*)malloc(sizeof(int)*intnum);
	}

	/* 処理フロー No.1 */
	/* NL-SD-FCB テーブル時刻設定 */
	//tintv = timediff( fcb->te, fcb->ts )/(double)fcb->ntfcb;
	for (i=0;i<intnum;i++)
	{
		ts = timeadd( fcb->ts, (double)i*tintv );
		te = timeadd( fcb->ts, (double)(i+1)*tintv );
		for(nij=0;nij<maxpnum;nij++)
		{
			fcb->nlfcb[nij*intnum+i].ts = ts;
			fcb->nlfcb[nij*intnum+i].te = te;
			fcb->nlfcb[nij*intnum+i].ave = 0.0;
			fcb->nlfcb[nij*intnum+i].var = 0.0;
		}
	}

	/* 処理フロー No.2 */
	/* NL-SD-FCBの平均・分散計算 */
	for (i=0;i<fcb->n;i++)/* 一重差ループ */
	{
		if ( fcb->sdd[i].b1fix == 0 ) continue;
		ni=fcb->sdd[i].s1;
		nj=fcb->sdd[i].s2;
		nij=(ni-1)*(2*NSATGPS-ni)/2+nj-ni-1;
		for (j=0;j<intnum;j++)
		{
			if ( check_oltime(	&fcb->sdd[i].ts,
								&fcb->sdd[i].te,
								&fcb->nlfcb[nij*intnum+j].ts,
								&fcb->nlfcb[nij*intnum+j].te,
								tintv/2.0,
								&tsdum,
								&tedum ) )
			{
				fcb->nlfcb[nij*intnum+j].ave += fcb->sdd[i].b1f;
				fcb->nlfcb[nij*intnum+j].var += pow(fcb->sdd[i].b1f,2 );
				nlcnt[nij][j] ++;
			}
		}
	}
	for (i=0;i<maxpnum;i++)
	{
		for (j=0;j<intnum;j++) {
			if(nlcnt[i][j]==0) continue;
			fcb->nlfcb[i*intnum+j].ave /= nlcnt[i][j];
			fcb->nlfcb[i*intnum+j].var /= nlcnt[i][j];
			fcb->nlfcb[i*intnum+j].var -= pow(fcb->nlfcb[i*intnum+j].ave,2);
		}
	}

	for (i=0;i<maxpnum;i++)
	{
		free(nlcnt[i]);
	}
	free(nlcnt);

	/* 処理フロー No.3 */
	/* 終了 */
	return;
}

/* SD-FCBファイル出力 */
int outputfcb(prcopt_t *popt, filopt_t *fopt, fcb_t *fcb, solopt_t *sopt)
{
	/* 処理フロー No.0 */
	/* 変数定義 */
	FILE *fp;      /* 出力ファイルポインタ */
	char ts[32]={"\0"}, te[32]={"\0"}; /* エポック文字列 */
	int nwl;/* WL-FCBデータ数 */
	int nnl;/* NL-FCN区間数 */
	int itvlnl;/* NL-FCB区間間隔 */
	int i,j;/* ループカウンタ */

	if(sopt->fcbout==FCBOUTOPT_OFF) return 0;

	/* 処理フロー No.1 */
	/* ファイルオープン */
	if(NULL==(fp=fopen(sopt->fcb,"w")))
	{
		showmsg("error: fcb file open");
		trace(1,"error: fcb file open");
		return 1;
	}

	/* コメント書き込み */
	fprintf(fp,"# file-satantfile=%s\n",fopt->satantp);
	fprintf(fp,"# file-rcvantfile=%s\n",fopt->rcvantp);

	/* 処理フロー No.2 */
	/* ヘッダ書き込み */
	time2str(fcb->ts,ts,6);
	time2str(fcb->te,te,6);
	nwl=NSATGPS*(NSATGPS-1)/2;
//	fprintf(fp,"%s %s %d %d %d\n",
//		ts,te,nwl,fcb->ntfcb,popt->nlfcbitvl);
	fprintf(fp,"%s %s %d %d %d\n",
		ts,te,nwl,fcb->ntfcb,popt->mopt.tifcb);

	fprintf(fp,"\n");

	/* 処理フロー No.3 */
	/* WL-FCB書き込み */
	for(i=0;i<nwl;i++)
	{
		 fprintf(fp,"%d %.6f %.6f\n",
			fcb->wlfcb[i].n, fcb->wlfcb[i].ave, fcb->wlfcb[i].var);
	}
	fprintf(fp,"\n");

	/* 処理フロー No.4 */
	/* NL-FCB書き込み */
	nnl=fcb->ntfcb;
	for(i=0;i<nnl;i++)
	{
		for(j=0;j<nwl;j++)
		{
			time2str(fcb->nlfcb[nnl*j+i].ts,ts,6);
			fprintf(fp,"%s %.6f %.6f\n",
			ts, fcb->nlfcb[nnl*j+i].ave, fcb->nlfcb[nnl*j+i].var);
		}
		fprintf(fp,"\n");
	}

	/* 処理フロー No.5 */
	/* ファイルクローズ */
	if(NULL!=fp)
	{
		fclose(fp);
	}

	/* 処理フロー No.6 */
	/* 終了処理 */
	return 0;
}

/* Cross-Correlation Receiver テーブルファイル読み込み */
int readCC(const char *file, char rcvname[MAXRECTYP][256]){

	FILE *fp;
	char buff[1024],ori[34]="Cross-Correlation Receivers Table";
	int line=1, nod=0;
	int len = 0;
	const int nType = 20;
	const int nBias = 22;

	trace(3,"readL2C: file=%s\n",file);

	if((fp=fopen(file,"r"))==NULL){
		trace(1,"cross-correlation receiver table file open error: %s",file);
		return -1;
	}

	while(fgets(buff,sizeof(buff),fp)){
		/*タイトルチェック*/
		if(line==1){
//			/*先頭から比較*/
//			if(((strncmp(buff,ori,32)))!=0){
//				showmsg("cross-correlation receiver table file title error: %s",file);
//				trace(1,"cross-correlation receiver table file title error: %s\n",file);
//				fclose(fp);
//				return -1;
//			}
		}
		if(line>4){
			len = strlen(buff);
			if(buff[0] != '\n') {
				strncpy(rcvname[nod], buff, len-1);
				trim(rcvname[nod]);
				/* RINEX Headerと同じフォーマット(A20)かチェック */
				if(strlen(rcvname[nod])>20)
				{
					showmsg("invalid receiver type name: file=%s line=%d\n",file,line);
					trace(2,"invalid receiver type name: file=%s line=%d\n",file,line);
				}
				nod++;
			}

			/*データ数が最大値以上かどうか*/
			if(nod>=MAXRECTYP){
				showmsg("cross-correlation receiver table file size error: %s",file);
				trace(1,"cross-correlation receiver table file size error: %s\n",file);
				fclose(fp);
				return -1;
			}
		}
		line++;
	}
	fclose(fp);
	trace(4,"number of readCC data: %d\n",nod);
	if(nod <1){
		showmsg("no cross-correlation receiver data error: %s\n",file);
		trace(1,"no cross-correlation receiver data error: %s\n",file);
	}
	return nod;
}

/* cross-correlation receiver チェック */
int checkccrcv(char *file, mbs_t *mbs)
{
	char ccrcv[MAXRECTYP][256];
	int nod=0;
	int i,j;

	if(NULL==file) return 1;

	nod=readCC(file, ccrcv);

	if(nod>0)
	{
		for(i=0;i<mbs->stas.nsta;i++)
		{
			for(j=0;j<nod;j++)
			{
				if(0==strcmp(mbs->stas.sta[i].rectype, ccrcv[j]))
				{
					mbs->stas.sta[i].cctype=1;
					break;
                }
			}
		}
	}

	return 0;
}

/* SD-FCB算出 */
int calsdfcb(mbs_t *mbs, prcopt_t *popt, filopt_t *fopt, solopt_t *sopt)
{
	/* 処理フロー No.0 */
	/* 変数定義 */
	fcb_t *fcb; /* FCB構造体 */
	int status;	/* 処理ステータス 0:正常 1:異常 */
	int nt;		/* NL-SD-FCBテーブル区間数 */
	double elipse;/* 測位処理経過時間（秒） */
	double dt;/* NL-SD-FCB推定区間（秒） */

	/* アンテナ補正ファイル入力チェック */
	if( NULL==fopt->satantp)
	{
		showmsg("error : no satantfile file");
		trace(1,"no satantfile file\n");
		return 1;
	}
	if(NULL==fopt->rcvantp)
	{
		showmsg("error : no rcvantfile file");
		trace(1,"no rcvantfile file\n");
		return 1;
	}

	/* Cross-Correlation タイプ受信機フラグ設定 */
	checkccrcv(fopt->cc,mbs);


	/* 処理フロー No.1 */
	/* FCB格納領域確保 */
	elipse = timediff(mbs->te, mbs->ts);
//	dt=(double)(popt->nlfcbitvl);
	dt=(double)(popt->mopt.tifcb);
	nt=(int)(elipse/dt)+1;
	fcb=fcbnew(nt);
	fcb->ts=mbs->ts;
	fcb->te=mbs->te;

	/* 処理フロー No.2 */
	/* 一重差パスの作成 */
	//pairsdmv(mbs->pass, fcb);
	pairsdmv(mbs, fcb);

	/* 処理フロー No.3 */
	/* WL-SD-FCB 平均及び分散の算出 */
	while(1)
	{
		calwlsdfcb(fcb);

		/* 処理フロー No.4 */
		/* WLアンビギュイティFIX判定 */
		status = fixwlsdfcb(fcb, popt);
		if(0==status) break;
		else if(1==status) continue;
		else if(9==status) return 1;
		else return 1;
	}

	/* 処理フロー No.5 */
	/* パス毎のNL-SD-FCBの算出 */
	calnlsdfcb(fcb/* , pass*/);

	/* 処理フロー No.6 */
	/* NL-SD-FCBのタイムテーブル作成 */
	makenlfcbtbl(fcb,dt);

	/* 処理フロー No.7 */
	/* SD-FCBファイル出力 */
	outputfcb(popt, fopt, fcb, sopt);

	/* 処理フロー No.8 */
	/* 終了処理 */
	return status;
}


