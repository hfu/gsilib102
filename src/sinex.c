/*------------------------------------------------------------------------------
* sinex.c : sinex functions
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

static FILE *fp_snx=NULL;        /* sinex file pointer */
static FILE *fp_snx_tm=NULL;        /* sinex file pointer */
static FILE *fp_snx_pos=NULL;        /* sinex file pointer */
static gtime_t tstart;
static gtime_t tend;
static int tstartset=0;
static int tendset=0;
static unsigned int estPrmInd=0;
static unsigned int solInd=0;
static const char* frmVer = "2.10";
static const char* fileAgency = "GSI";
static const char* dataAgency = "GSI";
static const int cnstrCode = 0;
static const char* solCont = " S          ";
static prcopt_t *popt=NULL;
static solopt_t *sopt=NULL;
static mbs_t *mbs=NULL;
static sta_t *sta;
static int solSta;
static sol_t solLast;

static void outsolsnxs();
static int outsolsnxm();

extern int setsnxtstart(gtime_t* t) {
	int ret = 0;
	if((NULL!=t) && ((t->time!=0) || (t->sec!=0.0)) && (tstartset == 0)) {
		tstart = *t;
		tstartset = 1;
		ret = 1;
	}
	return ret;
}

extern int setsnxtend(gtime_t* t) {
	int ret = 1;
	tend = *t;
	tendset = 1;
	return ret;
}

extern int rtkopensnx(const prcopt_t *opt, const mbs_t *mbs0, const sta_t *sta0, const char *file, char **infile, int n)
{
	int i;

	trace(3,"rtkopensnx: file=%s\n",file);
	estPrmInd=0;
	solInd=0;

	tstartset=0;
	tendset=0;

	if (!(fp_snx=fopen(file,"w"))) {
		trace(1,"rtkopensnx: file open error file=%s\n",file);
		fp_snx=NULL;
		return 0;
	}
	if (!(fp_snx_tm=tmpfile())) {
		trace(1,"rtkopensnx: file open error file=%s\n",file);
		fclose(fp_snx);
		fp_snx=NULL;
		fp_snx_tm=NULL;
		return 0;
	}
	if (!(fp_snx_pos=tmpfile())) {
		trace(1,"rtkopensnx: file open error file=%s\n",file);
		fclose(fp_snx);
		fclose(fp_snx_tm);
		fp_snx=NULL;
		fp_snx_tm=NULL;
		fp_snx_pos=NULL;
		return 0;
	}

	popt = opt;
	mbs = mbs0;
	sta = sta0;


	if(
	(popt->mode==PMODE_MULTI) && (mbs)) {
		setsnxtstart(&(mbs->ts));
		setsnxtend(&(mbs->te));
	}

	fprintf(fp_snx,"%%=SNX %-4s %-3s --:---:----- %-3s --:---:----- --:---:----- P ----- %1d%-12s\n"
		, frmVer
		, fileAgency
		, dataAgency
		, cnstrCode
		, solCont);

	fprintf(fp_snx,"+FILE/REFERENCE\n");
	fprintf(fp_snx,"DESCRIPTION        %-60s\n",ANAL_CENTER_FULL);
	fprintf(fp_snx,"SOFTWARE           %-60s\n",PROGRAM_NAME_VER);
	for(i=0;i<n;i++) {
		fprintf(fp_snx,"INPUT              %-60s\n",infile[i]);
	}
	fprintf(fp_snx,"-FILE/REFERENCE\n");


	fprintf(fp_snx_tm,"+SOLUTION/EPOCHS\n");
	fprintf(fp_snx_tm,"*CODE PT SOLN T _DATA_START_ __DATA_END__ _MEAN_EPOCH_\n");


	fprintf(fp_snx_pos,"+SOLUTION/ESTIMATE\n");
	fprintf(fp_snx_pos,"*INDEX TYPE__ CODE PT SOLN _REF_EPOCH__ UNIT S __ESTIMATED VALUE____ _STD_DEV___\n");


	fflush(fp_snx);
	fflush(fp_snx_pos);
    
	return 1;
}

int timedoy(const gtime_t* t) {
	double ep0[6]={0.0 ,1.0 ,1.0 ,0.0 ,0.0, 0.0};
	double ep[6];
	double diff;
	int doy;
	time2epoch(*t, ep);
	ep0[0] = ep[0];
	doy = (int)floor(timediff(*t, epoch2time(ep0))/86400.0) + 1;
	return doy;
}

int snxtimefrm(const gtime_t* t, char *frm) {

	int week;
	double ep[6];
	int doy;
	double sod;
//	double tow
	gtime_t t0;

	t0.time=t->time;
	t0.sec=t->sec;

	if(t0.sec > 0.5) {
	  t0.time += 1;
	}
	t0.sec = 0.0;
//	tow=time2gpst(t0,&week);
	doy=timedoy(&t0);
	time2epoch(t0, ep);
	sod = ep[3]*3600.0 + ep[4]*60.0 + ep[5];

	sprintf(frm, "%02d:%03d:%05.0f",(int)ep[0]%100, doy, sod);
}


static int fortranformat(const unsigned int w,const unsigned int p, double v, char* str) {
    char* buf = (char*)calloc(w*2, sizeof(char));
    int x;
    int w0;
    char t;
    if(w-p>5) {
        w0 = sprintf(buf, "%-+*.*E", w, p-1, v*10);
        t = buf[1];
        buf[1] = buf[2];
        buf[2] = t;
        x = (buf+w0-1) - strchr(buf, 'E');
        while(x>3) {
            if     (strstr(buf,"E+0")) repstr(buf,"E+0", "E+");
            else if(strstr(buf,"E-0")) repstr(buf,"E-0", "E-");
            else break;
            x = (buf+strlen(buf)-1) - strchr(buf, 'E');
        }
        if(v>=0) buf[0] = '0';
        strcpy(str, buf);
    }
	else if((w-p==5) && (v>=0)) {
		w0 = sprintf(buf, "%-*.*E" , w, p-1, v*10);
        t = buf[0];
        buf[0] = buf[1];
        buf[1] = t;
        x = (buf+w0-1) - strchr(buf, 'E');
        while(x>3) {
			if     (strstr(buf,"E+0")) repstr(buf,"E+0", "E+");
            else if(strstr(buf,"E-0")) repstr(buf,"E-0", "E-");
            else break;
            x = (buf+strlen(buf)-1) - strchr(buf, 'E');
        }
		strcpy(str, buf);
    }
    else {
        w0 = sprintf(buf, "%-*.*E" , w, p-1, v*10);
        t = buf[1];
        buf[1] = buf[2];
        buf[2] = t;
		x = (buf+w0-1) - strchr(buf, 'E');
        while(x>2) {
            if     (strstr(buf,"E+0")) repstr(buf,"E+0", "E+");
            else if(strstr(buf,"E-0")) repstr(buf,"E-0", "E-");
            else break;
            x = (buf+strlen(buf)-1) - strchr(buf, 'E');
        }
        strcpy(str, buf);
    }
    free(buf);
}

extern void rtkclosesnx(const gtime_t* te)
{
	gtime_t gtn;
    double eps[6],epe[6],epn[6];
    int doys,doye,doyn;
	double sods,sode,sodn;
	char buf[256]={0};
	int i;
	double antpos[3],dms1[3],dms2[3];
	char frms[13] = {0},frme[13] = {0},frmn[13] = {0};

	trace(3,"rtkclosesnx:\n");
	if (fp_snx) {

		fprintf(fp_snx,"+SITE/ID\n");
		fprintf(fp_snx,"*CODE PT __DOMES__ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_\n");
		if(popt->mode==PMODE_MULTI) {
			for(i=0;i<mbs->stas.nsta;i++) {
				ecef2pos(mbs->stas.sta[i].pos,antpos);
				antpos[2]-=geoidh(antpos);
				antpos[0]*=R2D;
				antpos[1]*=R2D;
				if(antpos[1]<0) antpos[1]+=360.0;
				deg2dms(antpos[0],dms1);
				deg2dms(antpos[1],dms2);
				fprintf(fp_snx," %-4s  A           P %-22s %3.0f %02.0f %04.1f %3.0f %02.0f %04.1f %7.1f\n"
					, mbs->stas.sta[i].name2, "", dms2[0],dms2[1],dms2[2],dms1[0],dms1[1],dms1[2],antpos[2]);
			}
		}
		else if((popt->mode==PMODE_FIXED) || (popt->mode==PMODE_PPP_FIXED)) {
			ecef2pos(popt->ru,antpos);
			antpos[2]-=geoidh(antpos);
			antpos[0]*=R2D;
			antpos[1]*=R2D;
			if(antpos[1]<0) antpos[1]+=360.0;
			deg2dms(antpos[0],dms1);
			deg2dms(antpos[1],dms2);
			fprintf(fp_snx," %-4s  A           P %-22s %3.0f %02.0f %04.1f %3.0f %02.0f %04.1f %7.1f\n"
				, sta[0].name, "", dms2[0],dms2[1],dms2[2],dms1[0],dms1[1],dms1[2],antpos[2]);
		}
		else {
			ecef2pos(sta[0].pos,antpos);
			antpos[2]-=geoidh(antpos);
			antpos[0]*=R2D;
			antpos[1]*=R2D;
			if(antpos[1]<0) antpos[1]+=360.0;
			deg2dms(antpos[0],dms1);
			deg2dms(antpos[1],dms2);
			fprintf(fp_snx," %-4s  A           P %-22s %3.0f %02.0f %04.1f %3.0f %02.0f %04.1f %7.1f\n"
				, sta[0].name, "", dms2[0],dms2[1],dms2[2],dms1[0],dms1[1],dms1[2],antpos[2]);
		}
		fprintf(fp_snx,"-SITE/ID\n");


		if((popt->mode==PMODE_STATIC) || (popt->mode == PMODE_PPP_STATIC)) {
			outsolsnxs();
		}
		else if(popt->mode==PMODE_MULTI) {
			outsolsnxm();
		}

		if (fp_snx_tm) {
			fprintf(fp_snx_tm,"-SOLUTION/EPOCHS\n");
			rewind(fp_snx_tm);
			while(NULL!=fgets(buf,sizeof(buf), fp_snx_tm)) {
				fprintf(fp_snx, "%s", buf);
			}
			fclose(fp_snx_tm);
			fp_snx_tm=NULL;
		}

		if (fp_snx_pos) {
			fprintf(fp_snx_pos,"-SOLUTION/ESTIMATE\n");
			rewind(fp_snx_pos);
			while(NULL!=fgets(buf,sizeof(buf), fp_snx_pos)) {
				fprintf(fp_snx, "%s", buf);
			}
			fclose(fp_snx_pos);
			fp_snx_pos=NULL;
		}


		fprintf(fp_snx,"%%ENDSNX\n");
		rewind(fp_snx);
		gtn = pctime();

		if(NULL!=te) tend = *te;
		snxtimefrm(&tstart, frms);
		snxtimefrm(&tend, frme);
		snxtimefrm(&gtn, frmn);

		fprintf(fp_snx,"%%=SNX %-4s %-3s %s %-3s %s %s P %05d %1d%-12s\n"
			, frmVer
			, fileAgency
			, frmn
			, dataAgency
			, frms
			, frme
			, estPrmInd
			, cnstrCode
			, solCont);
		fclose(fp_snx);
	}
	fp_snx=NULL;
	popt=NULL;
	mbs=NULL;
	sta=NULL;
}

extern void outsolsnx(rtk_t *rtk, const obsd_t *obs)
{
	double vela[3]={0},acca[3]={0},x[3],p[3];
	int i;
	char pos[100] = {0}, sgm[100] = {0};
	char frm[13] = {0};

	double P[9];

	setsnxtstart(&rtk->sol.time);
	setsnxtend(&rtk->sol.time);


	if (!fp_snx) return;

	if((popt->mode==PMODE_STATIC) || (popt->mode==PMODE_PPP_STATIC)) {
		if((rtk->sol.stat == 1)) {
			solLast = rtk->sol;
			solSta = 1;
		}
		else if(solSta == 0) {
			solLast = rtk->sol;
		}
		return;
	}

	trace(3,"outsolsnx:\n");

	soltocov(&rtk->sol,P);

	snxtimefrm(&rtk->sol.time, frm);

	++solInd;

 //	for (i=0;i<3;i++) x[i]=rtk->xa[i];
	for (i=0;i<3;i++) x[i]=rtk->sol.rr[i];
	p[0]=P[0];
	p[1]=P[4];
	p[2]=P[8];

	fprintf(fp_snx_tm," %-4s  A %4d P %s %s %s\n", sta->name, solInd, frm, frm, frm);

	fortranformat(21, 15, x[0], pos);
	fortranformat(11,  6, p[0], sgm);
	fprintf(fp_snx_pos," %5d STAX   %-4s  A %4d %s m    %1d %-21s %-11s\n"
		, ++estPrmInd, sta->name, solInd, frm, cnstrCode, pos, sgm);
	fortranformat(21, 15, x[1], pos);
	fortranformat(11,  6, p[1], sgm);
	fprintf(fp_snx_pos," %5d STAY   %-4s  A %4d %s m    %1d %-21s %-11s\n"
		, ++estPrmInd, sta->name, solInd, frm, cnstrCode, pos, sgm);
	fortranformat(21, 15, x[2], pos);
	fortranformat(11,  6, p[2], sgm);
	fprintf(fp_snx_pos," %5d STAZ   %-4s  A %4d %s m    %1d %-21s %-11s\n"
		, ++estPrmInd, sta->name, solInd, frm, cnstrCode, pos, sgm);
}

//extern void outsolsnxs(rtk_t *rtk, const obsd_t *obs)
static void outsolsnxs()
{
	double vela[3]={0},acca[3]={0},x[3],p[3];
	int i;
	gtime_t tmean;
	char pos[100] = {0}, sgm[100] = {0};
	char frms[13] = {0};
	char frme[13] = {0};
	char frmm[13] = {0};
	double P[9];

	if (!fp_snx) return;

	trace(3,"outsolsnx:\n");

	tmean = timeadd(tstart, timediff(tend, tstart)/2.0);
	snxtimefrm(&tstart, frms);
	snxtimefrm(&tend, frme);
	snxtimefrm(&tmean, frmm);

	fprintf(fp_snx_tm," %-4s  A    1 P %s %s %s\n", sta->name, frms, frme, frmm);

	for (i=0;i<3;i++) x[i]=solLast.rr[i];
	soltocov(&solLast,P);
	p[0]=P[0];
	p[1]=P[4];
	p[2]=P[8];

	fortranformat(21, 15, x[0], pos);
	fortranformat(11,  6, p[0], sgm);
	fprintf(fp_snx_pos," %5d STAX   %-4s  A    1 %s m    %1d %-21s %-11s\n"
		, ++estPrmInd, sta->name, frmm, cnstrCode, pos, sgm);
	fortranformat(21, 15, x[1], pos);
	fortranformat(11,  6, p[1], sgm);
	fprintf(fp_snx_pos," %5d STAY   %-4s  A    1 %s m    %1d %-21s %-11s\n"
		, ++estPrmInd, sta->name, frmm, cnstrCode, pos, sgm);
	fortranformat(21, 15, x[2], pos);
	fortranformat(11,  6, p[2], sgm);
	fprintf(fp_snx_pos," %5d STAZ   %-4s  A    1 %s m    %1d %-21s %-11s\n"
		, ++estPrmInd, sta->name, frmm, cnstrCode, pos, sgm);
}

static int outsolsnxm()
{
    int i;
	gtime_t tmean;
	char    pos[100] = {0}, sgm[100] = {0};
	char frms[13] = {0};
	char frme[13] = {0};
	char frmm[13] = {0};

	if (!fp_snx) return 0;

	tmean = timeadd(mbs->ts, timediff(mbs->te, mbs->ts)/2.0);
	snxtimefrm(&tstart, frms);
	snxtimefrm(&tend, frme);
	snxtimefrm(&tmean, frmm);

	trace(3,"outsolsnx:\n");

	for(i=0;i<mbs->stas.nsta;i++) {
		fprintf(fp_snx_tm," %-4s  A    1 P %s %s %s\n", mbs->stas.sta[i].name2
			, frms, frme, frmm);
	}

	for(i=0;i<mbs->stas.nsta;i++) {
		fortranformat(21, 15, mbs->param.par[i*3], pos);
		fortranformat(11,  6, mbs->param.pg[i*3], sgm);
		fprintf(fp_snx_pos," %5d STAX   %-4s  A    1 %s m    %1d %-21s %-11s\n"
			, ++estPrmInd, mbs->stas.sta[i].name2, frmm, cnstrCode, pos, sgm);
		fortranformat(21, 15, mbs->param.par[i*3+1], pos);
		fortranformat(11,  6, mbs->param.pg[i*3+1], sgm);
		fprintf(fp_snx_pos," %5d STAY   %-4s  A    1 %s m    %1d %-21s %-11s\n"
			, ++estPrmInd, mbs->stas.sta[i].name2, frmm, cnstrCode, pos, sgm);
		fortranformat(21, 15, mbs->param.par[i*3+2], pos);
		fortranformat(11,  6, mbs->param.pg[i*3+2], sgm);
		fprintf(fp_snx_pos," %5d STAZ   %-4s  A    1 %s m    %1d %-21s %-11s\n"
			, ++estPrmInd, mbs->stas.sta[i].name2, frmm, cnstrCode, pos, sgm);
	}

    return 1;
}