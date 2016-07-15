/*------------------------------------------------------------------------------
* options.c : options functions
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

/* system options buffer -----------------------------------------------------*/
static prcopt_t prcopt_;
static solopt_t solopt_;
static filopt_t filopt_;
static int antpostype_[2];
static double elmask_,elmaskar_,elmaskhold_;
static double antpos_[2][3];
static char exsats_[1024];
static char snrmask_[NFREQ][1024];

/* system options table ------------------------------------------------------*/
#define SWTOPT  "0:off,1:on"
#define MODOPT  "0:single,1:dgps,2:kinematic,3:static,4:movingbase,5:fixed,6:ppp-kine,7:ppp-static,8:ppp-fixed,9:mb-static"

#define FRQOPT  "1:l1,2:l1+l2,3:l1+l2+l5,4:l1+l2+l5+l6,5:l1+l2+l5+l6+l7,2:l1+l5"
#define FRQOPT2 "1:l1,2:l2,3:l5,12:l1+l2,13:l1+l5,23:l2+l5,123:l1+l2+l5,1:L1,2:L2,3:L5,12:L1+L2,13:L1+L5,23:L2+L5,123:L1+L2+L5"
#define TYPOPT  "0:forward,1:backward,2:combined"
#define IONOPT  "0:off,1:brdc,2:sbas,3:dual-freq,4:est-stec,5:ionex-tec,6:qzs-brdc,7:qzs-lex,8:vtec_sf,9:vtec_ef,10:gtec"
#define TRPOPT  "0:off,1:saas,2:sbas,3:est-ztd,4:est-ztdgrad"
#define EPHOPT  "0:brdc,1:precise,2:brdc+sbas,3:brdc+ssrapc,4:brdc+ssrcom"
#define NAVOPT  "1:gps+2:sbas+4:glo+8:gal+16:qzs+32:comp"
#define GAROPT  "0:off,1:on,2:autocal,3:use IFB table"

#define SOLOPT  "0:llh,1:xyz,2:enu,3:nmea,4:rva"
#define TSYOPT  "0:gpst,1:utc,2:jst"
#define TFTOPT  "0:tow,1:hms"
#define DFTOPT  "0:deg,1:dms"
#define HGTOPT  "0:ellipsoidal,1:geodetic"
#define GEOOPT  "0:internal,1:egm96,2:egm08_2.5,3:egm08_1,4:gsi2000"
#define STAOPT  "0:all,1:single"
#define STSOPT  "0:off,1:state,2:residual"
#define OWAOPT  "0:off,1:new,2:append"
#define OUTOPT  "0:off,1:out"
#define ARMOPT  "1:continuous,2:instantaneous,3:fix-and-hold"
#define RAROPT  "0:off,1:LAMBDA,2:LC"
#define PAROPT  "0:off,1:cnes,2:cnes-ils,3:fcb,4:est"
#define POSOPT  "0:llh,1:xyz,2:single,3:posfile,4:rinexhead,5:rtcm"
#define ERROPT	"0:user settings,1:table"
#define TSCOPT	"0:off,1:on"
#define PHAOPT	"0:off,1:table,2:rinexrtcm"
#define GL2OPT	"0:off,1:table,2:est"       //  table:DCB file
#define L2COPT	"0:L2P,1:L2C"
#define ISBOPT	"0:off,1:table,2:est,3:est-P,4:est-L,5:est(L:estimate only)"       //  table:ISB file
#define DIFOPT	"0:inall,1:exc-glo"         //  default:1
//#define GLTOPT	"0:off,1:est"               //  default:1)
#define IDTOPT	"0:set,1:finlname,2:rinex-header"               //  default:1)

#define CHIOPT	"0:off,1:on"                //  default:1



opt_t sysopts[]={
    {"pos1-posmode",    3,  (void *)&prcopt_.mode,       MODOPT, 1},
	{"pos1-frequency",  3,  (void *)&prcopt_.nfreq,      FRQOPT, 0 },
	{"pos1-freqs",      4,  (void *)&prcopt_.oprfrq,     FRQOPT2, 1 },
	{"pos1-l2cprior",  	3,  (void *)&prcopt_.l2cprior,   L2COPT, 1 },
	{"pos1-gpscodepriL1",2, (void *)&prcopt_.codepri[ISYSGPS][0],"", 0},
	{"pos1-gpscodepriL2",2, (void *)&prcopt_.codepri[ISYSGPS][1],"", 0},
	{"pos1-gpscodepriL5",2, (void *)&prcopt_.codepri[ISYSGPS][2],"", 0},
	{"pos1-glocodepriL1",2, (void *)&prcopt_.codepri[ISYSGLO][0],"", 1},
	{"pos1-glocodepriL2",2, (void *)&prcopt_.codepri[ISYSGLO][1],"", 1},
	{"pos1-qzscodepriL1",2, (void *)&prcopt_.codepri[ISYSQZS][0],"", 0},
	{"pos1-qzscodepriL2",2, (void *)&prcopt_.codepri[ISYSQZS][1],"", 0},
	{"pos1-qzscodepriL5",2, (void *)&prcopt_.codepri[ISYSQZS][2],"", 0},
	{"pos1-galcodepriL1",2, (void *)&prcopt_.codepri[ISYSGAL][0],"", 0},
	{"pos1-galcodepriL2",2, (void *)&prcopt_.codepri[ISYSGAL][1],"", 0},
    {"pos1-soltype",    3,  (void *)&prcopt_.soltype,    TYPOPT, 1 },
    {"pos1-elmask",     1,  (void *)&elmask_,            "deg", 1  },
    {"pos1-snrmask_r",  3,  (void *)&prcopt_.snrmask.ena[0],SWTOPT, 1},
    {"pos1-snrmask_b",  3,  (void *)&prcopt_.snrmask.ena[1],SWTOPT, 1},
    {"pos1-snrmask_L1", 2,  (void *)snrmask_[0],         ""     , 1},
    {"pos1-snrmask_L2", 2,  (void *)snrmask_[1],         ""     , 1},
    {"pos1-snrmask_L5", 2,  (void *)snrmask_[2],         ""     , 1},
    {"pos1-dynamics",   3,  (void *)&prcopt_.dynamics,   SWTOPT, 1 },
    {"pos1-tidecorr",   3,  (void *)&prcopt_.tidecorr,   SWTOPT, 1 },
    {"pos1-ionoopt",    3,  (void *)&prcopt_.ionoopt,    IONOPT, 1 },
    {"pos1-tropopt",    3,  (void *)&prcopt_.tropopt,    TRPOPT, 1 },
	{"pos1-tsyscorr",   3,  (void *)&prcopt_.tsyscorr,   TSCOPT, 1 },
    {"pos1-sateph",     3,  (void *)&prcopt_.sateph,     EPHOPT, 1 },
    {"pos1-posopt1",    3,  (void *)&prcopt_.posopt[0],  SWTOPT, 1 },
    {"pos1-posopt2",    3,  (void *)&prcopt_.posopt[1],  SWTOPT, 1 },
    {"pos1-posopt3",    3,  (void *)&prcopt_.posopt[2],  SWTOPT, 1 },
    {"pos1-posopt4",    3,  (void *)&prcopt_.posopt[3],  SWTOPT, 1 },
    {"pos1-posopt5",    3,  (void *)&prcopt_.posopt[4],  SWTOPT, 1 },
	{"pos1-exclsats",   2,  (void *)exsats_,             "prn ...", 1},
    {"pos1-navsys",     0,  (void *)&prcopt_.navsys,     NAVOPT, 1 },
    
	{"pos2-armode",     3,  (void *)&prcopt_.modear,     ARMOPT, 1 },
	{"pos2-rtkarmode",  3,  (void *)&prcopt_.modertkar,  RAROPT, 1 },
    {"pos2-ppparmode",  3,  (void *)&prcopt_.modepppar,  PAROPT, 1 },
    {"pos2-gloarmode",  3,  (void *)&prcopt_.glomodear,  GAROPT, 1 },
    {"pos2-arthres",    1,  (void *)&prcopt_.thresar[0], "", 1     },
    {"pos2-arlockcnt",  0,  (void *)&prcopt_.minlock,    "", 1     },
    {"pos2-arelmask",   1,  (void *)&elmaskar_,          "deg", 1  },
    {"pos2-arminfix",   0,  (void *)&prcopt_.minfix,     "", 1     },
    {"pos2-elmaskhold", 1,  (void *)&elmaskhold_,        "deg", 1  },
    {"pos2-aroutcnt",   0,  (void *)&prcopt_.maxout,     ""   , 1  },
	{"pos2-phasshft",   3,  (void *)&prcopt_.phasshft,   PHAOPT, 1 },
    {"pos2-maxage",     1,  (void *)&prcopt_.maxtdiff,   "s"   , 1 },
    {"pos2-syncsol",    3,  (void *)&prcopt_.syncsol,    SWTOPT, 1 },
    {"pos2-slipthres",  1,  (void *)&prcopt_.thresslip,  "m"   , 1 },
    {"pos2-rejionno",   1,  (void *)&prcopt_.maxinno,    "m"   , 1 },
    {"pos2-rejgdop",    1,  (void *)&prcopt_.maxgdop,    ""    , 1 },
    {"pos2-niter",      0,  (void *)&prcopt_.niter,      ""    , 1 },
    {"pos2-baselen",    1,  (void *)&prcopt_.baseline[0],"m"   , 1 },
    {"pos2-basesig",    1,  (void *)&prcopt_.baseline[1],"m"   , 1 },

	{"pos2-gpsl2bias",  3,  (void *)&prcopt_.gl2bias,    GL2OPT, 1 },
	{"pos2-isb",        3,  (void *)&prcopt_.isb,        ISBOPT, 1 },
    {"pos2-diff",       3,  (void *)&prcopt_.diff,       DIFOPT, 1 },
//    {"pos2-glot",       3,  (void *)&prcopt_.glt,        GLTOPT, 1 },

	{"pos3-phacycfile", 2,  (void *)&prcopt_.mopt.ifpcs,   ""  , 1 },
	{"pos3-gloifbfile", 2,  (void *)&prcopt_.mopt.ififb,   ""  , 1 },
	{"pos3-errmodfile", 2,  (void *)&prcopt_.mopt.iferr,   ""  , 1 },

	{"pos3-estsatclock",3,  (void *)&prcopt_.mopt.estsatclk,SWTOPT, 1},
	{"pos3-estsatfcb",  3,  (void *)&prcopt_.mopt.estsatfcb,SWTOPT, 1},
	{"pos3-semidcpara", 2,  (void *)&prcopt_.mopt.ifsdp, ""     , 1},
	{"pos3-soldirfile", 2,  (void *)&prcopt_.mopt.ofdir, ""     , 1},

	{"pos3-estinttroze",1,  (void *)&prcopt_.mopt.tiztd,  "sec" , 1},
	{"pos3-estinttroew",1,  (void *)&prcopt_.mopt.tigra,  "sec" , 1},
	{"pos3-rwsigzen",   1,  (void *)&prcopt_.mopt.sigtrop[0],"m", 1},
	{"pos3-rwsigew",    1,  (void *)&prcopt_.mopt.sigtrop[1],"m", 1},
	{"pos3-rwsigns",    1,  (void *)&prcopt_.mopt.sigtrop[2],"m", 1},
	{"pos3-thresholdo", 1,  (void *)&prcopt_.mopt.cpomcth,   "m", 1},
	{"pos3-thresholdc", 1,  (void *)&prcopt_.mopt.promcth,   "m", 1},
	{"pos3-maxbaseline",1,  (void *)&prcopt_.mopt.maxdddist, "m", 1},
	{"pos3-judgevalwl", 1,  (void *)&prcopt_.mopt.wlddfix,   "%", 1},
	{"pos3-judgevall1", 1,  (void *)&prcopt_.mopt.l1ddfix,   "%", 1},
	{"pos3-weightddlc", 1,  (void *)&prcopt_.mopt.wdd,  "m"     , 1},
	{"pos3-concriite",  1,  (void *)&prcopt_.mopt.itrconv,   "" , 1},
	{"pos3-maxitera",   0,  (void *)&prcopt_.mopt.itrmax,    "" , 1},
	{"pos3-fcbfile",    2,  (void *)&prcopt_.mopt.iffcb, ""     , 1},

	{"pos3-temstofile", 2,  (void *)&prcopt_.mopt.eptmp,     "" , 1},
	{"pos3-nlfcb",      1,  (void *)&prcopt_.mopt.tifcb,  "sec" , 1},
	{"pos3-minsd",      1,  (void *)&prcopt_.mopt.minpass,"sec" , 1},
	{"pos3-mindd",      1,  (void *)&prcopt_.mopt.minddpass,"sec", 1},
	{"pos3-maxsd",      1,  (void *)&prcopt_.mopt.maxsigw,"cycle", 1},
	{"pos3-fixwl",      1,  (void *)&prcopt_.mopt.minconfw,   "", 1},
	{"pos3-fixnl",      1,  (void *)&prcopt_.mopt.minconf1,   "", 1},
	{"pos3-mobstadn",   1,  (void *)&prcopt_.mopt.sigr[0],    "", 1},
	{"pos3-mobstade",   1,  (void *)&prcopt_.mopt.sigr[1],    "", 1},
	{"pos3-mobstadu",   1,  (void *)&prcopt_.mopt.sigr[2],    "", 1},
	{"pos3-basestadn",  1,  (void *)&prcopt_.mopt.sigb[0],    "", 1},
	{"pos3-basestade",  1,  (void *)&prcopt_.mopt.sigb[1],    "", 1},
	{"pos3-basestadu",  1,  (void *)&prcopt_.mopt.sigb[2],    "", 1},

	{"out-solformat",   3,  (void *)&solopt_.posf,       SOLOPT , 1},
	{"out-outhead",     3,  (void *)&solopt_.outhead,    SWTOPT , 1},// ok
	{"out-outopt",      3,  (void *)&solopt_.outopt,     SWTOPT , 1},// okヘッダー出力時のオプション出力ONOFF
	{"out-timesys",     3,  (void *)&solopt_.times,      TSYOPT , 1},// ok時系
	{"out-timeform",    3,  (void *)&solopt_.timef,      TFTOPT , 1},// 時刻フォーマット
	{"out-timendec",    0,  (void *)&solopt_.timeu,      ""     , 1},//   ok
	{"out-degform",     3,  (void *)&solopt_.degf,       DFTOPT , 1},//okdeg
	{"out-fieldsep",    2,  (void *) solopt_.sep,        ""     , 1},//okセパレータ
	{"out-height",      3,  (void *)&solopt_.height,     HGTOPT , 1},// okgeodetic height
	{"out-geoid",       3,  (void *)&solopt_.geoid,      GEOOPT , 1}, // okジオイド
	{"out-solstatic",   3,  (void *)&solopt_.solstatic,  STAOPT , 1}, //ok
	{"out-nmeaintv1",   1,  (void *)&solopt_.nmeaintv[0],"s"    , 1}, // ok
	{"out-nmeaintv2",   1,  (void *)&solopt_.nmeaintv[1],"s"    , 1}, // ok
	{"out-outstat",     3,  (void *)&solopt_.sstat,      STSOPT , 1},

	{"out-isbout",      3,  (void *)&solopt_.isbout,     OWAOPT , 1},
	{"out-isbfile",     2,  (void *)&solopt_.isbfile,     ""    , 1},
	{"out-gl2out",      3,  (void *)&solopt_.gl2out,     OWAOPT , 1},
	{"out-gl2file",     2,  (void *)&solopt_.gl2file,     ""    , 1},
	{"out-sbres",       3,  (void *)&solopt_.sbresout,   OUTOPT , 1},
	{"out-possinex",    3,  (void *)&solopt_.possnxout,  OUTOPT , 1},
	{"out-possinexfile",2,  (void *)&solopt_.possinexfile,""    , 1},
    {"out-ion",         3,  (void *)&solopt_.ionout,     OUTOPT , 1},
	{"out-ionfile",     2,  (void *)&solopt_.ionfile,     ""    , 1},
	{"out-trop",        3,  (void *)&solopt_.tropout,    OUTOPT , 1},
	{"out-recclock",    3,  (void *)&solopt_.recclout,   OUTOPT , 1},
	{"out-clcfile",     2,  (void *)&solopt_.clcfile,     ""    , 1},
    {"out-satclock",    3,  (void *)&solopt_.satclout,   OUTOPT , 1},
    
    {"stats-errmodel",  3,  (void *)&prcopt_.errmodel,   ERROPT , 1},
    {"stats-eratio1",   1,  (void *)&prcopt_.eratio[0],  ""     , 1},
    {"stats-eratio2",   1,  (void *)&prcopt_.eratio[1],  ""     , 1},
    {"stats-eratio3",   1,  (void *)&prcopt_.eratio[2],  ""     , 1},
    {"stats-errphase",  1,  (void *)&prcopt_.err[1],     "m"    , 1},
    {"stats-errphaseel",1,  (void *)&prcopt_.err[2],     "m"    , 1},
    {"stats-dcbratio",	1,  (void *)&prcopt_.dcbratio,   "", 1},
    {"stats-errphasebl",1,  (void *)&prcopt_.err[3],     "m/10km", 1},
    {"stats-errdoppler",1,  (void *)&prcopt_.err[4],     "Hz"   , 1},
    {"stats-stdbias",   1,  (void *)&prcopt_.std[0],     "m"    , 1},
    {"stats-stdiono",   1,  (void *)&prcopt_.std[1],     "m"    , 1},
    {"stats-stdtrop",   1,  (void *)&prcopt_.std[2],     "m"    , 1},
    {"stats-stdgrad",   1,  (void *)&prcopt_.std[3],     "m"    , 1},
    {"stats-stdclkr",   1,  (void *)&prcopt_.std[4],     "m"    , 1},
    {"stats-stdclks",   1,  (void *)&prcopt_.std[5],     "m"    , 1},
    {"stats-stdissb",   1,  (void *)&prcopt_.std[6],     "m"    , 1},
    {"stats-prnaccelh", 1,  (void *)&prcopt_.prn[3],     "m/s^2", 1},
    {"stats-prnaccelv", 1,  (void *)&prcopt_.prn[4],     "m/s^2", 1},
    {"stats-prnbias",   1,  (void *)&prcopt_.prn[0],     "m"    , 1},
    {"stats-prniono",   1,  (void *)&prcopt_.prn[1],     "m"    , 1},
    {"stats-prntrop",   1,  (void *)&prcopt_.prn[2],     "m"    , 1},
    {"stats-prnisbl",   1,  (void *)&prcopt_.prn[5],     "m"    , 1},
    {"stats-prnisbp",   1,  (void *)&prcopt_.prn[6],     "m"    , 1},
    {"stats-prngl2" ,   1,  (void *)&prcopt_.prn[7],     "m"    , 1},
    {"stats-clkstab",   1,  (void *)&prcopt_.sclkstab,   "s/s"  , 1},
    {"stats-chisqr",    3,  (void *)&prcopt_.chisqr,     CHIOPT , 0},
    
    {"ant1-postype",    3,  (void *)&antpostype_[0],     POSOPT , 1},
    {"ant1-pos1",       1,  (void *)&antpos_[0][0],      "deg|m", 1},
    {"ant1-pos2",       1,  (void *)&antpos_[0][1],      "deg|m", 1},
    {"ant1-pos3",       1,  (void *)&antpos_[0][2],      "m|m"  , 1},
    {"ant1-anttype",    2,  (void *)prcopt_.anttype[0],  ""     , 1},
    {"ant1-antdele",    1,  (void *)&prcopt_.antdel[0][0],"m"   , 1},
    {"ant1-antdeln",    1,  (void *)&prcopt_.antdel[0][1],"m"   , 1},
    {"ant1-antdelu",    1,  (void *)&prcopt_.antdel[0][2],"m"   , 1},
    {"ant1-rectype",    2,  (void *)prcopt_.rectype[0],  ""     , 1},
    {"ant1-idtype",     3,  (void *)&prcopt_.idtype[0],  IDTOPT , 1},
    {"ant1-id4char",    2,  (void *)prcopt_.id4char[0],  ""     , 1},
    
    {"ant2-postype",    3,  (void *)&antpostype_[1],     POSOPT , 1},
    {"ant2-pos1",       1,  (void *)&antpos_[1][0],      "deg|m", 1},
    {"ant2-pos2",       1,  (void *)&antpos_[1][1],      "deg|m", 1},
    {"ant2-pos3",       1,  (void *)&antpos_[1][2],      "m|m"  , 1},
    {"ant2-anttype",    2,  (void *)prcopt_.anttype[1],  ""     , 1},
    {"ant2-antdele",    1,  (void *)&prcopt_.antdel[1][0],"m"   , 1},
    {"ant2-antdeln",    1,  (void *)&prcopt_.antdel[1][1],"m"   , 1},
    {"ant2-antdelu",    1,  (void *)&prcopt_.antdel[1][2],"m"   , 1},
    {"ant2-rectype",    2,  (void *)prcopt_.rectype[1],  ""     , 1},
    {"ant2-idtype",     3,  (void *)&prcopt_.idtype[1],  IDTOPT , 1},
    {"ant2-id4char",    2,  (void *)prcopt_.id4char[1],  ""     , 1},
    
    {"misc-timeinterp", 3,  (void *)&prcopt_.intpref,    SWTOPT , 1},
    {"misc-sbasatsel",  0,  (void *)&prcopt_.sbassatsel, "0:all", 1},
    {"misc-rnxopt1",    2,  (void *)prcopt_.rnxopt[0],   ""     , 1},
    {"misc-rnxopt2",    2,  (void *)prcopt_.rnxopt[1],   ""     , 1},
    
    {"file-satantfile", 2,  (void *)&filopt_.satantp,    ""     , 1},
    {"file-rcvantfile", 2,  (void *)&filopt_.rcvantp,    ""     , 1},
    {"file-staposfile", 2,  (void *)&filopt_.stapos,     ""     , 1},
    {"file-geoidfile",  2,  (void *)&filopt_.geoid,      ""     , 1},
    {"file-ionofile",   2,  (void *)&filopt_.iono,       ""     , 1},
    {"file-dcbfile",    2,  (void *)&filopt_.dcb,        ""     , 1},
    {"file-cirtfile",   2,  (void *)&filopt_.cirtfile,   ""     , 1},
	{"file-eopfile",    2,  (void *)&filopt_.eop,        ""     , 1},
    {"file-blqfile",    2,  (void *)&filopt_.blq,        ""     , 1},
    {"file-tempdir",    2,  (void *)&filopt_.tempdir,    ""     , 1},
    {"file-geexefile",  2,  (void *)&filopt_.geexe,      ""     , 1},
	{"file-solstatfile",2,  (void *)&filopt_.solstat,    ""     , 1},
	{"file-tracefile",  2,  (void *)&filopt_.trace,      ""     , 1},
	{"file-isbfile",    2,  (void *)&filopt_.isb,        ""     , 1},

	/* 以下、GSILIB v1.0から追加されたパラメタ */
	{"out-fcbout",      3,  (void *)&solopt_.fcbout,     SWTOPT , 1},
	{"out-fcb",			2,  (void *)&solopt_.fcb,        ""     , 1},
	{"out-nlfcbintv",   0,  (void *)&prcopt_.nlfcbitvl,  "s"    , 1},
	{"misc-baselist",   2,  (void *)&solopt_.baselist,   ""     , 1},
	{"misc-roverlist",  2,  (void *)&solopt_.roverlist,  ""     , 1},
	{"file-ccfile",     2,  (void *)&filopt_.cc,         ""     , 1},

    {"",0,NULL,"", 0} /* terminator */
};

/* enum to string ------------------------------------------------------------*/
static int enum2str(char *s, const char *comment, int val)
{
    char str[32],*p,*q;
    int n;
    
    n=sprintf(str,"%d:",val);
    if (!(p=strstr(comment,str))) {
        return sprintf(s,"%d",val);
    }
    if (!(q=strchr(p+n,','))&&!(q=strchr(p+n,')'))) {
        strcpy(s,p+n);
        return strlen(p+n);
    }
    strncpy(s,p+n,q-p-n); s[q-p-n]='\0';
    return q-p-n;
}
/* flag to string ------------------------------------------------------------*/
static int flag2str(char *s, const char *comment, int* val)
{
    int num=0;
    int i=0;
	while(-1<val[i] && val[i]<10){
		num*=10;
		num+=val[i]+1;
		i++;
	}
    return enum2str(s, comment, num);
}
/* string to enum ------------------------------------------------------------*/
static int str2enum(const char *str, const char *comment, int *val)
{
    const char *p;
    char s[32];
    
    for (p=comment;;p++) {
       if (!(p=strstr(p,str))) break;
       if (*(p-1)!=':') continue;
       for (p-=2;'0'<=*p&&*p<='9';p--) ;
       return sscanf(p+1,"%d",val)==1;
    }
    sprintf(s,"%30.30s:",str);
	trim(s);
    if ((p=strstr(comment,s))) { /* number */
        return sscanf(p,"%d",val)==1;
    }
    return 0;
}
/* string to enum ------------------------------------------------------------*/
extern int str2flag(const char *str, const char *comment, int *val)
{
    const char *p;
    char s[32],vs[32]={0};
    int v=-1;
    int i=0,j=0,f=0,n=0,ns=0;
     freq_ind ch[32];
    char* c;
    
    for (p=comment;;p++) {
       if (!(p=strstr(p,str))) break;
       if (*(p-1)!=':') continue;
       for (p-=2;'0'<=*p&&*p<='9';p--) ;
       if(sscanf(p+1,"%d",&v)!=1) return 0;
       else break;
    }
    sprintf(s,"%30.30s:",str);
    trim(s);
    if ((p=strstr(comment,s))) { /* number */
        if(sscanf(p,"%d",&v)!=1) {
           return 0;
        }
    }
    if(v==-1) return 0;

    i=0;
    sprintf(vs,"%d",v);
    ns = strlen(vs);
    for(i=0;i<ns;i++) {
        f=vs[i]-0x30;
        if(f>0) val[j++]=f-1;
    }
    val[j]=NOSELECT;

    return 1;
}
/* search option ---------------------------------------------------------------
* search option record
* args   : char   *name     I  option name
*          opt_t  *opts     I  options table
*                              (terminated with table[i].name="")
* return : option record (NULL: not found)
*-----------------------------------------------------------------------------*/
extern opt_t *searchopt(const char *name, const opt_t *opts)
{
    int i;
    
    trace(3,"searchopt: name=%s\n",name);
    
    for (i=0;*opts[i].name;i++) {
     //   if (strstr(opts[i].name,name)) return (opt_t *)(opts+i);
        if (!strcmp(opts[i].name,name)) return (opt_t *)(opts+i);
    }
    return NULL;
}
/* string to option value ------------------------------------------------------
* convert string to option value
* args   : opt_t  *opt      O  option
*          char   *str      I  option value string
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int str2opt(opt_t *opt, const char *str)
{
    switch (opt->format) {
        case 0: *(int    *)opt->var=atoi(str); break;
        case 1: *(double *)opt->var=atof(str); break;
        case 2: strcpy((char *)opt->var,str);  break;
		case 3: return str2enum(str,opt->comment,(int *)opt->var);
		case 4: return str2flag(str,opt->comment,(int *)opt->var);
        default: return 0;
	}
    return 1;
}
/* option value to string ------------------------------------------------------
* convert option value to string
* args   : opt_t  *opt      I  option
*          char   *str      O  option value string
* return : length of output string
*-----------------------------------------------------------------------------*/
extern int opt2str(const opt_t *opt, char *str)
{
	char *p=str;

	trace(3,"opt2str : name=%s\n",opt->name);

	switch (opt->format) {
		case 0: p+=sprintf(p,"%d"   ,*(int   *)opt->var); break;
		case 1: p+=sprintf(p,"%.15g",*(double*)opt->var); break;
		case 2: p+=sprintf(p,"%s"   , (char  *)opt->var); break;
		case 3: p+=enum2str(p,opt->comment,*(int *)opt->var); break;
		case 4: p+=flag2str(p,opt->comment,(int *)opt->var); break;
    }
    return (int)(p-str);
}
/* option to string -------------------------------------------------------------
* convert option to string (keyword=value # comment)
* args   : opt_t  *opt      I  option
*          char   *buff     O  option string
* return : length of output string
*-----------------------------------------------------------------------------*/
extern int opt2buf(const opt_t *opt, char *buff)
{
    char *p=buff;
    int n;
    
    trace(3,"opt2buf : name=%s\n",opt->name);
    
    p+=sprintf(p,"%-18s =",opt->name);
    p+=opt2str(opt,p);
    if (*opt->comment) {
        if ((n=(int)(buff+30-p))>0) p+=sprintf(p,"%*s",n,"");
        p+=sprintf(p," # (%s)",opt->comment);
    }
    return (int)(p-buff);
}
/* load options ----------------------------------------------------------------
* load options from file
* args   : char   *file     I  options file
*          opt_t  *opts     IO options table
*                              (terminated with table[i].name="")
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int loadopts(const char *file, opt_t *opts)
{
    FILE *fp;
    opt_t *opt;
    char buff[2048],*p;
    int n=0;
    
    trace(3,"loadopts: file=%s\n",file);
    
    if (!(fp=fopen(file,"r"))) {
        trace(1,"loadopts: options file open error (%s)\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        n++;
        chop(buff);
        
        if (buff[0]=='\0') continue;
        
        if (!(p=strstr(buff,"="))) {
            fprintf(stderr,"invalid option %s (%s:%d)\n",buff,file,n);
            continue;
        }
        *p++='\0';
        chop(buff);
        if (!(opt=searchopt(buff,opts))) {
            fprintf(stderr,"invalid option %s (%s:%d)\n",buff,file,n);
            continue;
        }
        if (!str2opt(opt,p)) {
            fprintf(stderr,"invalid option value %s (%s:%d)\n",buff,file,n);
            continue;
        }
    }
    fclose(fp);
    
    return 1;
}
/* load options ----------------------------------------------------------------
* load options from file
* args   : char   *file     I  options file
*          opt_t  *rcvpts     IO receiver options table
*          opt_t  *syspts     IO system options table
*                              (terminated with table[i].name="")
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int loadopts4rcv(const char *file, opt_t *rcvpts, opt_t *sysopts)
{
    FILE *fp;
    opt_t *opt;
    char buff[2048],*p;
    int n=0;
    
    trace(3,"loadopts: file=%s\n",file);
    
    if (!(fp=fopen(file,"r"))) {
        trace(1,"loadopts: options file open error (%s)\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        n++;
        chop(buff);
        
        if (buff[0]=='\0') continue;
        
        if (!(p=strstr(buff,"="))) {
            fprintf(stderr,"invalid option %s (%s:%d)\n",buff,file,n);
            continue;
        }
        *p++='\0';
        chop(buff);
        if (!(opt=searchopt(buff,rcvpts)) && !(opt=searchopt(buff,sysopts))) {
            fprintf(stderr,"invalid option %s (%s:%d)\n",buff,file,n);
            continue;
        }
        if (!str2opt(opt,p)) {
            fprintf(stderr,"invalid option value %s (%s:%d)\n",buff,file,n);
            continue;
        }
    }
    fclose(fp);
    
    return 1;
}
/* save options to file --------------------------------------------------------
* save options to file
* args   : char   *file     I  options file
*          char   *mode     I  write mode ("w":overwrite,"a":append);
*          char   *comment  I  header comment (NULL: no comment)
*          opt_t  *opts     I  options table
*                              (terminated with table[i].name="")
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int saveopts(const char *file, const char *mode, const char *comment,
                    const opt_t *opts)
{
    FILE *fp;
    char buff[2048];
    int i;
    
    trace(3,"saveopts: file=%s mode=%s\n",file,mode);
    
    if (!(fp=fopen(file,mode))) {
        trace(1,"saveopts: options file open error (%s)\n",file);
        return 0;
    }
    if (comment) fprintf(fp,"# %s\n\n",comment);
    
    for (i=0;*opts[i].name;i++) {
        if(opts[i].update==0) continue;
        opt2buf(opts+i,buff);
        fprintf(fp,"%s\n",buff);
    }
    fclose(fp);
    return 1;
}
/* system options buffer to options ------------------------------------------*/
static void buff2sysopts(void)
{
    double pos[3],*rr;
    char buff[1024],*p,*id;
    int i,j,sat,*ps;
    
    prcopt_.elmin     =elmask_    *D2R;
    prcopt_.elmaskar  =elmaskar_  *D2R;
    prcopt_.elmaskhold=elmaskhold_*D2R;
    
    for (i=0;i<2;i++) {
        ps=i==0?&prcopt_.rovpos:&prcopt_.refpos;
        rr=i==0?prcopt_.ru:prcopt_.rb;
        
        if (antpostype_[i]==0) { /* lat/lon/hgt */
            *ps=0;
            pos[0]=antpos_[i][0]*D2R;
            pos[1]=antpos_[i][1]*D2R;
            pos[2]=antpos_[i][2];
            pos2ecef(pos,rr);
        }
        else if (antpostype_[i]==1) { /* xyz-ecef */
            *ps=0;
            rr[0]=antpos_[i][0];
            rr[1]=antpos_[i][1];
            rr[2]=antpos_[i][2];
        }
        else *ps=antpostype_[i]-1;
    }
    /* excluded satellites */
    for (i=0;i<MAXSAT;i++) prcopt_.exsats[i]=0;
    if (exsats_[0]!='\0') {
        strcpy(buff,exsats_);
        for (p=strtok(buff," ");p;p=strtok(NULL," ")) {
            if (*p=='+') id=p+1; else id=p;
            if (!(sat=satid2no(id))) continue;
            prcopt_.exsats[sat-1]=*p=='+'?2:1;
        }
    }
    /* snrmask */
    for (i=0;i<NFREQ;i++) {
        for (j=0;j<9;j++) prcopt_.snrmask.mask[i][j]=0.0;
        strcpy(buff,snrmask_[i]);
        for (p=strtok(buff,","),j=0;p&&j<9;p=strtok(NULL,",")) {
            prcopt_.snrmask.mask[i][j++]=atof(p);
        }
    }
}
/* options to system options buffer ------------------------------------------*/
static void sysopts2buff(void)
{
    double pos[3],*rr;
    char id[32],*p;
    int i,j,sat,*ps;
    
    elmask_    =prcopt_.elmin     *R2D;
    elmaskar_  =prcopt_.elmaskar  *R2D;
    elmaskhold_=prcopt_.elmaskhold*R2D;
    
    for (i=0;i<2;i++) {
        ps=i==0?&prcopt_.rovpos:&prcopt_.refpos;
        rr=i==0?prcopt_.ru:prcopt_.rb;
        
        if (*ps==0) {
            antpostype_[i]=0;
            ecef2pos(rr,pos);
            antpos_[i][0]=pos[0]*R2D;
            antpos_[i][1]=pos[1]*R2D;
            antpos_[i][2]=pos[2];
        }
        else antpostype_[i]=*ps+1;
    }
    /* excluded satellites */
    exsats_[0]='\0';
    for (sat=1,p=exsats_;sat<=MAXSAT&&p-exsats_<(int)sizeof(exsats_)-32;sat++) {
        if (prcopt_.exsats[sat-1]) {
            satno2id(sat,id);
            p+=sprintf(p,"%s%s%s",p==exsats_?"":" ",
                       prcopt_.exsats[sat-1]==2?"+":"",id);
        }
    }
    /* snrmask */
    for (i=0;i<NFREQ;i++) {
        snrmask_[i][0]='\0';
        p=snrmask_[i];
        for (j=0;j<9;j++) {
            p+=sprintf(p,"%s%.0f",j>0?",":"",prcopt_.snrmask.mask[i][j]);
        }
    }
}
/* reset system options to default ---------------------------------------------
* reset system options to default
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
extern void resetsysopts(void)
{
    int i,j;
    
    trace(3,"resetsysopts:\n");
    
	prcopt_=prcopt_default;
    solopt_=solopt_default;
    filopt_.satantp[0]='\0';
    filopt_.rcvantp[0]='\0';
    filopt_.stapos [0]='\0';
    filopt_.geoid  [0]='\0';
    filopt_.dcb    [0]='\0';
    filopt_.blq    [0]='\0';
    filopt_.solstat[0]='\0';
	filopt_.trace  [0]='\0';
	filopt_.cc     [0]='\0';
	solopt_.fcb    [0]='\0';
    for (i=0;i<2;i++) antpostype_[i]=0;
    elmask_=15.0;
    elmaskar_=0.0;
    elmaskhold_=0.0;
    for (i=0;i<2;i++) for (j=0;j<3;j++) {
        antpos_[i][j]=0.0;
    }
    exsats_[0] ='\0';
	prcopt_.mopt = mbsopt_default;
}

/* get system options ----------------------------------------------------------
* get system options
* args   : prcopt_t *popt   IO processing options (NULL: no output)
*          solopt_t *sopt   IO solution options   (NULL: no output)
*          folopt_t *fopt   IO file options       (NULL: no output)
* return : none
* notes  : to load system options, use loadopts() before calling the function
*-----------------------------------------------------------------------------*/
extern void getsysopts(prcopt_t *popt, solopt_t *sopt, filopt_t *fopt)
{
    trace(3,"getsysopts:\n");
    
    buff2sysopts();
    if (popt) *popt=prcopt_;
    if (sopt) *sopt=solopt_;
    if (fopt) *fopt=filopt_;
}
/* set system options ----------------------------------------------------------
* set system options
* args   : prcopt_t *prcopt I  processing options (NULL: default)
*          solopt_t *solopt I  solution options   (NULL: default)
*          filopt_t *filopt I  file options       (NULL: default)
* return : none
* notes  : to save system options, use saveopts() after calling the function
*-----------------------------------------------------------------------------*/
extern void setsysopts(const prcopt_t *prcopt, const solopt_t *solopt,
                       const filopt_t *filopt)
{
	trace(3,"setsysopts:\n");
    
    resetsysopts();
    if (prcopt) prcopt_=*prcopt;
    if (solopt) solopt_=*solopt;
    if (filopt) filopt_=*filopt;
    sysopts2buff();
}

extern int checkopts(prcopt_t *popt, solopt_t *sopt, filopt_t *fopt) {
    int i;
    int nf=0;
    int ret=1;


	if((popt->modertkar==RTKAR_LC) || (popt->modertkar==RTKAR_TCAR)) {
		if(popt->ionoopt!=IONOOPT_IFLC) {
			fprintf(stderr,"apply Ionosphere-free Combination \n");
			popt->ionoopt=IONOOPT_IFLC;
        }
        else {
        }
    }


 
    for(i=0;i<6;++i) {
        if(NOSELECT==popt->oprfrq[i]) break;
        ++nf;
    }
    if((popt->nfreq==0) && (nf==0)) {        // not set
        popt->nfreq=2;
        popt->oprfrq[0]=ind_L1;
        popt->oprfrq[1]=ind_L2;
	}
    else if((popt->nfreq!=0) && (nf==0)) {    // rtklib/prototype
        if(popt->nfreq>=1) popt->oprfrq[0]=ind_L1;
        if(popt->nfreq>=2) popt->oprfrq[1]=ind_L2;
        if(popt->nfreq>=3) popt->oprfrq[2]=ind_L5;
		if(popt->nfreq>=4) popt->oprfrq[3]=ind_L6;
        if(popt->nfreq>=5) popt->oprfrq[4]=ind_L7;
        if(popt->nfreq>=6) popt->oprfrq[5]=ind_L8;
    }
    else if((popt->nfreq==0) && (nf!=0)) {
        popt->nfreq=nf;
    }
    else if(popt->nfreq!=nf) {
        popt->nfreq=nf;
        fprintf(stderr,"mismatch between ""pos1-frequency"" and ""pos1-freqs""\n");
    }

    if((popt->ionoopt==IONOOPT_IFLC) && (popt->nfreq==1)) {
        fprintf(stderr,"mismatch between pos1-ionoopt and pos1-freqs\n");
        popt->nfreq=2;
        popt->oprfrq[0]=ind_L1;
        popt->oprfrq[1]=ind_L2;
        popt->oprfrq[2]=NOSELECT;
        popt->oprfrq[3]=NOSELECT;
		popt->oprfrq[4]=NOSELECT;
		popt->oprfrq[5]=NOSELECT;
	}
	if(((popt->mode==PMODE_PPP_KINEMA) || (popt->mode==PMODE_PPP_STATIC) || (popt->mode==PMODE_PPP_STATIC))
		&& (popt->modepppar!=PPPAR_OFF)){
		if(popt->navsys!=SYS_GPS) {
			popt->modepppar=PPPAR_OFF;
			showmsg("%s", "PPPAR only suppots GPS. it cannot be run.\n");
		}
	}

	return ret;
}
