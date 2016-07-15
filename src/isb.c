/*------------------------------------------------------------------------------
* isb.c : isb (inter system bias) functions
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

#define NINCISB     262144

static int addisbdata(nav_t *nav, const isb_t *data)
{
	isb_t *isb_data;

	if (nav->nimax<=nav->ni) {
		if (nav->nimax<=0) nav->nimax=NINCISB; else nav->nimax*=2;
		if (!(isb_data=(isb_t *)realloc(nav->isb,sizeof(isb_t)*nav->nimax))) {
			trace(1,"addisbdata: memalloc error n=%dx%d\n",sizeof(isb_t),nav->nimax);
			free(nav->isb); nav->isb=NULL; nav->ni=nav->nimax=0;
            return -1;
        }
		nav->isb=isb_data;
    }
	nav->isb[nav->ni++]=*data;
	return 1;
}

/* read isb parameters file --------------------------------------------------*/
static int readisbf(const char *file, nav_t *nav)
{
	FILE *fp;
	double bias;
	char buff[256];
	int sys,i,s_code,res;
	char stationname[21]={0};
	char stationname0[21]={0};
	char stationname_base[21]={0};
	char stationname0_base[21]={0};
	isb_t data={{0}};
    int nbuf = 0;
    int freq = 0;
    char type = '\0';
    int n = 0;
    int ifreq = 0;
	int itype = 0; /* =0:phase,1:code*/
	int ver;

	trace(3,"readisbf: file=%s\n",file);

	if (!(fp=fopen(file,"r"))) {
		trace(2,"isb parameters file open error: %s\n",file);
		return 0;
	}
	while (fgets(buff,sizeof(buff),fp)) {
		++n;
		if(4>n) continue;

		chop(buff);
		if(4==n) {
			if(strcmp(buff,"******************** * * * **********************") == 0){
				ver = 0;
			}
			else if(strcmp(buff,"******************** ******************** * * * **********************") == 0){
				ver = 1;
			}
			continue;
		}

		nbuf = strlen(buff);
		if(28 > nbuf) {
			trace(2,"isb parameter error: %s (%s:%d)\n",buff,file,n);
			fprintf(stderr,"isb parameter error: %s (%s:%d)\n",buff,file,n);
			continue;
		}

		if(ver==0) {
			setstr(stationname, buff, 20);
			bias = str2num(buff, 27, nbuf - 27);
			freq = str2num(buff, 23, 1);
			type = buff[25];
			sys = code2sys(buff[21]);

			switch(sys){
				case SYS_GPS:
					s_code = NSYSGPS-1;
					break;
				case SYS_GLO:
					s_code = NSYSGPS+NSYSGLO-1;
					break;
				case SYS_GAL:
					s_code = NSYSGPS+NSYSGLO+NSYSGAL-1;
					break;
				case SYS_QZS:
					s_code = NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS-1;
					break;
				case SYS_CMP:
					s_code = NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS+NSYSCMP-1;
					break;
				default:
					s_code = -1;
					break;
			}

			if(    (' ' != buff[20]) || (' ' != buff[22]) || (' ' != buff[24]) || (' ' != buff[26])
				|| (('P' != type) && ('L' != type))
				|| ((1 != freq) && (2 != freq) && (5 != freq))
				|| ( -1 == s_code)){
				trace(2,"isb parameter error: %s (%s:%d)\n",buff,file,n);
				fprintf(stderr,"isb parameter error: %s (%s:%d)\n",buff,file,n);
				continue;
			}

			ifreq = 5==freq?2:freq -1;
			itype = 'L'==type?0:1;

			if(stationname[0]!='\0') {
				strcpy(stationname0, stationname);
			}
			if(stationname0[0]!='\0') {
				res=0;
				for(i=0;i<nav->ni;i++)
				{
					if(strcmp( nav->isb[i].sta_name,stationname0)==0){
						/*àÍívÉfÅ[É^Ç†ÇË ÉoÉCÉAÉXèÓïÒÇÃÇ›í«â¡ */
						nav->isb[i].gsb[s_code][ifreq][itype]=bias*1E-9*CLIGHT;
						res=1;
						break;
					}
				}
				/*ã«ñ¢ìoò^ êVãKí«â¡*/
				if(res==0){
					strcpy(data.sta_name, stationname0);
					data.gsb[s_code][ifreq][itype]=bias*1E-9*CLIGHT;
					addisbdata(nav,&data);
					data.gsb[s_code][ifreq][itype]=0.0;
				}
			}
		}

		else if(ver==1) {
			setstr(stationname, buff, 20);
			setstr(stationname_base, buff+21, 20);
			bias = str2num(buff, 48, nbuf - 48);
			freq = str2num(buff, 44, 1);
			type = buff[46];
			sys = code2sys(buff[42]);

			switch(sys){
				case SYS_GPS:
					s_code = NSYSGPS-1;
					break;
				case SYS_GLO:
					s_code = NSYSGPS+NSYSGLO-1;
					break;
				case SYS_GAL:
					s_code = NSYSGPS+NSYSGLO+NSYSGAL-1;
					break;
				case SYS_QZS:
					s_code = NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS-1;
					break;
				case SYS_CMP:
					s_code = NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS+NSYSCMP-1;
					break;
				default:
					s_code = -1;
					break;
			}

			if(    (' ' != buff[20]) || (' ' != buff[41]) || (' ' != buff[43]) || (' ' != buff[45]) || (' ' != buff[47])
				|| (('P' != type) && ('L' != type))
				|| ((1 != freq) && (2 != freq) && (5 != freq))
				|| ( -1 == s_code)){
				trace(2,"isb parameter error: %s (%s:%d)\n",buff,file,n);
				fprintf(stderr,"isb parameter error: %s (%s:%d)\n",buff,file,n);
				continue;
			}

			ifreq = 5==freq?2:freq -1;
			itype = 'L'==type?0:1;

			if(stationname[0]!='\0') {
				strcpy(stationname0, stationname);
				stationname0_base[0] = '\0';
				if(stationname_base[0]!='\0') {
					strcpy(stationname0_base, stationname_base);
				}
			}
			else if(stationname0[0]!='\0') {
				if(stationname_base[0]!='\0') {
					strcpy(stationname0_base, stationname_base);
				}
			}

			if(stationname0[0]!='\0') {
				res=0;
				for(i=0;i<nav->ni;i++)
				{
					if(    (strcmp( nav->isb[i].sta_name     ,stationname0     )==0)
						&& (strcmp( nav->isb[i].sta_name_base,stationname0_base)==0)){
						/*àÍívÉfÅ[É^Ç†ÇË ÉoÉCÉAÉXèÓïÒÇÃÇ›í«â¡ */
						nav->isb[i].gsb[s_code][ifreq][itype]=bias*1E-9*CLIGHT;
						res=1;
						break;
					}
				}
				/*ã«ñ¢ìoò^ êVãKí«â¡*/
				if(res==0){
					strcpy(data.sta_name     , stationname0     );
					strcpy(data.sta_name_base, stationname0_base);
					data.gsb[s_code][ifreq][itype]=bias*1E-9*CLIGHT;
					addisbdata(nav,&data);
					data.gsb[s_code][ifreq][itype]=0.0;
					if(stationname0_base[0]!='\0') {
						strcpy(data.sta_name_base, stationname0     );
						strcpy(data.sta_name     , stationname0_base);
						data.gsb[s_code][ifreq][itype]=-bias*1E-9*CLIGHT;
						addisbdata(nav,&data);
						data.gsb[s_code][ifreq][itype]=0.0;
					}
				}
			}
		}
	}
	fclose(fp);
    
	return 1;
}
/* read isb parameters ---------------------------------------------------------
* read inter-system bias (isb) parameters
* args   : char   *file       I   isb parameters file (wild-card * expanded)
*          nav_t  *nav        IO  navigation data
* return : status (1:ok,0:error)
* notes  : 
*-----------------------------------------------------------------------------*/
extern int readisb(const char *file, nav_t *nav)
{
	int i,j,n;
	char *efiles[MAXEXFILE]={0};
    
	trace(3,"readisb : file=%s\n",file);
    
	for (i=0;i<MAXSAT;i++) for (j=0;j<3;j++) {
		nav->cbias[i][j]=0.0;
	}
	for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
			for (i--;i>=0;i--) free(efiles[i]);
            return 0;
        }
    }
    nav->ni=nav->nimax=0;

    n=expath(file,efiles,MAXEXFILE);
    
    for (i=0;i<n;i++) {
		readisbf(efiles[i],nav);
    }
	for (i=0;i<MAXEXFILE;i++) free(efiles[i]);
    
    return 1;
}


void setisb(const nav_t *nav, const char* rectype0, const char* rectype1, sta_t *sta0, sta_t *sta1)
{
	int i,j,k,m;
	int r=0;
	isb_t* pisb[2]={0};

	if(sta0!=NULL) {
		for(i=0;i<NSYSISB;++i) {
			for(j=0;j<NFREQ;++j) {
				for(k=0;k<2;++k) sta0->isb[i][j][k] = 0.0;
			}
		}
	}
	if(sta1!=NULL) {
		for(i=0;i<NSYSISB;++i) {
			for(j=0;j<NFREQ;++j) {
				for(k=0;k<2;++k) sta1->isb[i][j][k] = 0.0;
			}
		}
	}

	if(rectype0!=NULL) {
		for(m=0;m<nav->ni;++m) {
			if ( !strcmp(rectype0, nav->isb[m].sta_name)) {
				if ( !strcmp(rectype1, nav->isb[m].sta_name_base)) {
					for(i=0;i<NSYSISB;++i) {
						for(j=0;j<NFREQ;++j) {
							for(k=0;k<2;++k) sta0->isb[i][j][k] = nav->isb[m].gsb[i][j][k];
						}
					}
					r=1;
					break;
				}
			}
			if(rectype1!=NULL) {
				if ( !strcmp(rectype1, nav->isb[m].sta_name)) {
					if ( !strcmp(rectype0, nav->isb[m].sta_name_base)) {
						for(i=0;i<NSYSISB;++i) {
							for(j=0;j<NFREQ;++j) {
								for(k=0;k<2;++k) sta0->isb[i][j][k] = -nav->isb[m].gsb[i][j][k];
							}
						}
						r=1;
						break;
					}
				}
			}
			if(r==0) {
				if ( !strcmp(rectype0, nav->isb[m].sta_name)) {
					if ( !strcmp(nav->isb[m].sta_name_base, "")) {
						pisb[0] = &nav->isb[m];
					}
				}
				else if ( !strcmp(rectype1, nav->isb[m].sta_name)) {
					if ( !strcmp(nav->isb[m].sta_name_base, "")) {
						pisb[1] = &nav->isb[m];
					}
				}
			}
		}
	}
	if(r==0) {
		if(pisb[0] != NULL) {
			for(i=0;i<NSYSISB;++i) {
				for(j=0;j<NFREQ;++j) {
					for(k=0;k<2;++k) sta0->isb[i][j][k] = pisb[0]->gsb[i][j][k];
				}
			}
		}

		if(pisb[1] != NULL) {
			for(i=0;i<NSYSISB;++i) {
				for(j=0;j<NFREQ;++j) {
					for(k=0;k<2;++k) sta1->isb[i][j][k] = pisb[1]->gsb[i][j][k];
				}
			}
		}
	}
}

void chk_isb(char sysno, const prcopt_t *opt, const sta_t *sta, double y[NFREQ][2])
{
	int i,j;
//	int isys = sysind(sysno) - 2;
	int isys = sysind(sysno) - 1;

	for(i=0;i<3;++i) {
		for(j=0;j<2;++j) y[i][j] = 0.0;
    }

	if(isys < 0) return;
    if(opt->isb!=ISBOPT_TABLE) return;

    for(i=0;i<3;++i) {
        for(j=0;j<2;++j) {
            y[i][j] = sta->isb[isys][i][j];
        }
    }
	return;
}

