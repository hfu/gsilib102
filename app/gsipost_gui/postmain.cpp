/*------------------------------------------------------------------------------
* postmain.cpp
*
*    Copyright (C) 2014 by Geospatial Information Authority of Japan,
*    All rights reserved.
*
*
*  Original software: RTKLIB ver.2.4.2 p4
*
*    Copyright (C) 2007-2013 by T.Takasu, All rights reserved.
*
*
*  Options : gsipost_gui [-t title][-i file][-r file][-b file][-n file ...]
*                        [-d dir][-o file]
*                        [-ts y/m/d h:m:s][-te y/m/d h:m:s][-ti tint][-tu tunit]
*
*            -t title   window title
*            -i file    ini file path
*            -r file    rinex obs rover file
*            -b file    rinex obs base station file
*            -n file    rinex nav/clk, sp3, ionex or sp3 file
*            -d dir     output directory
*            -o file    output file
*            -ts y/m/d h:m:s time start
*            -te y/m/d h:m:s time end
*            -ti tint   time interval (s)
*            -tu tunit  time unit (hr)
*
* references :
*
* history : 2015/01/08  1.0  new
*-----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <vcl.h>

#include <locale.h>
#include <wchar.h>

#pragma hdrstop

#include "rtklib.h"
#include "postpos.h"
#include "postmain.h"
#include "postopt.h"
#include "kmzconv.h"
#include "refdlg.h"
#include "timedlg.h"
#include "confdlg.h"
#include "keydlg.h"
#include "aboutdlg.h"
#include "viewer.h"

#pragma package(smart_init)
#pragma resource "*.dfm"

TMainForm *MainForm;

#define PRGNAME     L"GSIPOST"
#define MAXHIST     20
//#define GOOGLE_EARTH "C:\\Program Files\\Google\\Google Earth\\googleearth.exe"
#define GOOGLE_EARTH L"C:\\Program Files (x86)\\Google\\Google Earth\\client\\googleearth.exe"

static const wchar_t version[]=L"$Revision: 1.1 $ $Date: 2008/07/17 22:14:45 $";

// global variables ---------------------------------------------------------
static wchar_t rov_ [256]=L"";          // rover name
static wchar_t base_[256]=L"";          // base-station name
extern "C" {

// show message in message area ---------------------------------------------
extern int showmsg(char *format, ...)
{
    va_list arg;
	char buff[1024];
	wchar_t *wbuff;
    if (*format) {
        va_start(arg,format);
		vsprintf(buff,format,arg);
                    va_end(arg);
		if (!(buff==NULL)){
			wchar_t *wch = new wchar_t[1024];
			if (!(wch==NULL)){
				mbstowcs(wch,buff,1024);
				MainForm->ShowMsg(wch);
				delete [] wch;
            }
        }
    }
    else Application->ProcessMessages();
    return !MainForm->BtnExec->Enabled;
}
// set time span of progress bar --------------------------------------------
extern void settspan(gtime_t ts, gtime_t te)
{
	obsts=ts;
	obste=te;
}
// set current time to show progress ----------------------------------------
extern void settime(gtime_t time)
{
    static int i=0;
    double tt;
	if (obste.time!=0&&obsts.time!=0&&(tt=timediff(obste,obsts))>0.0) {
		MainForm->Progress->Position=(int)(timediff(time,obsts)/tt*100.0+0.5);
    }
    if (i++%23==0) Application->ProcessMessages();
}

} // extern "C"

// convert string to double -------------------------------------------------
static double str2dbl(UnicodeString str)
{
    double val=0.0;
    swscanf(str.c_str(),L"%lf",&val);
    return val;
}
// constructor --------------------------------------------------------------
__fastcall TMainForm::TMainForm(TComponent* Owner)
    : TForm(Owner)
{
	int i;
	wchar_t *dotpoint,*file;
	char getfile[1024]={0};

	::GetModuleFileName(NULL,getfile,sizeof(getfile));
	file=char2wchar(getfile,-1);
	dotpoint=wcsrchr(file,L'.');
	if (dotpoint==NULL){
		wcscat(file,L".ini");
	}else{
		wcscpy(dotpoint,L".ini");
	}
	IniFile=file;

	if(file!=NULL){
		free(file);
		file=NULL;
	}

    setlocale( LC_ALL, "" );
	_wsetlocale(LC_ALL, L"" );
	DynamicModel=IonoOpt=TropOpt=RovAntPcv=RefAntPcv=AmbResMethod=0;
	AmbRes=1;
	RovPosType=RefPosType=0;
    OutCntResetAmb=5; LockCntFixAmb=5; FixCntHoldAmb=10;
    MaxAgeDiff=30.0; RejectThres=30.0; RejectGdop=30.0;
	MeasErrR1=MeasErrR2=MeasErrR3=100.0; MeasErr2=0.004; MeasErr3=0.003; MeasErr4=1.0;
	SatClkStab=1E-11; ValidThresAR=3.0;
	RovAntE=RovAntN=RovAntU=RefAntE=RefAntN=RefAntU=0.0;
	for (i=0;i<3;i++) RovPos[i]=0.0;
	for (i=0;i<3;i++) RefPos[i]=0.0;
	EstIntZ=7200.0;EstIntES=43200.0;
	RanWalSigZen=1.0;RanWalSigEW=1.0;RanWalSigNS=1.0;
	MaxBasDD=1.0;WeiDD=1.0;	ConCri=1.0;
	MaxIte=1;
    
    DoubleBuffered=true;
}
// callback on form create --------------------------------------------------
void __fastcall TMainForm::FormCreate(TObject *Sender)
{
    UnicodeString s;
    
    wchar_t *wch = char2wchar(VER_RTKLIB,-1);
    Caption=s.sprintf(L"%s ver.%s",PRGNAME,wch);
	if(wch!=NULL){
		free(wch);
		wch=NULL;
	}
	::DragAcceptFiles(Handle,true);
}
// callback on form show ----------------------------------------------------
void __fastcall TMainForm::FormShow(TObject *Sender)
{
    TComboBox *ifile[]={InputFile3,InputFile4,InputFile5};
    wchar_t *p,*argv[32],buff[1024];
    int argc=0,n=0,inputflag=0;;
    
    wcscpy(buff,GetCommandLineW());
    
    for (p=buff;*p&&argc<32;p++) {
        if (*p==' ') continue;
        if (*p=='"') {
            argv[argc++]=p+1;
            p=wcschr(p+1,'"');
			if (p==NULL) break;
        }
        else {
            argv[argc++]=p;
            p=wcschr(p+1,' ');
			if (p==NULL) break;
        }
        *p='\0';
    }
    for (int i=1;i<argc;i++) { // get ini file option
        if (!wcscmp(argv[i],L"-i")&&i+1<argc) IniFile=argv[++i];
    }
    LoadOpt();
    
    for (int i=1;i<argc;i++) {
        if      (!wcscmp(argv[i],L"-t")&&i+1<argc) Caption=argv[++i];
        else if (!wcscmp(argv[i],L"-r")&&i+1<argc) {
            InputFile1->Text=argv[++i];
            inputflag=1;
        }
        else if (!wcscmp(argv[i],L"-b")&&i+1<argc) InputFile2->Text=argv[++i];
        else if (!wcscmp(argv[i],L"-d")&&i+1<argc) {
            OutDirEna->Checked=true;
            OutDir->Text=argv[++i];
        }
        else if (!wcscmp(argv[i],L"-o")&&i+1<argc) OutputFile->Text=argv[++i];
        else if (!wcscmp(argv[i],L"-n")&&i+1<argc) {
            if (n<3) ifile[n++]->Text=argv[++i];
        }
        else if (!wcscmp(argv[i],L"-ts")&&i+2<argc) {
            TimeStart->Checked=true;
            TimeY1->Text=argv[++i]; TimeH1->Text=argv[++i];
        }
        else if (!wcscmp(argv[i],L"-te")&&i+2<argc) {
            TimeEnd->Checked=true;
            TimeY2->Text=argv[++i]; TimeH2->Text=argv[++i];
        }
        else if (!wcscmp(argv[i],L"-ti")&&i+1<argc) {
            TimeIntF->Checked=true;
            TimeInt->Text=argv[++i];
        }
        else if (!wcscmp(argv[i],L"-tu")&&i+1<argc) {
            TimeUnitF->Checked=true;
            TimeUnit->Text=argv[++i];
        }
    }
    if (inputflag) SetOutFile();
    
    UpdateEnable();
}
// callback on form close ---------------------------------------------------
void __fastcall TMainForm::FormClose(TObject *Sender, TCloseAction &Action)
{
    SaveOpt();
}
// callback on drop files ---------------------------------------------------
void __fastcall TMainForm::DropFiles(TWMDropFiles msg)
{
    POINT point={0};
    int top;
	wchar_t *p,file[1024];
    
	if (DragQueryFile((HDROP)msg.Drop,0xFFFFFFFF,NULL,0)<=0) return;
    DragQueryFileW((HDROP)msg.Drop,0,file,sizeof(file));
    if (!DragQueryPoint((HDROP)msg.Drop,&point)) return;
    
    top=Panel1->Top+Panel4->Top;
    if (point.y<=top+InputFile1->Top+InputFile1->Height) {
		InputFile1->Text=file;
        SetOutFile();
    }
    else if (point.y<=top+InputFile2->Top+InputFile2->Height) {
        InputFile2->Text=file;
    }
    else if (point.y<=top+InputFile3->Top+InputFile3->Height) {
        InputFile3->Text=file;
    }
    else if (point.y<=top+InputFile4->Top+InputFile4->Height) {
        InputFile4->Text=file;
    }
    else if (point.y<=top+InputFile5->Top+InputFile5->Height) {
        InputFile5->Text=file;
    }
}
// callback on button-plot --------------------------------------------------
void __fastcall TMainForm::BtnPlotClick(TObject *Sender)
{
	wchar_t outfile[1024]={0};

	if (OutDirEna->Checked == true) {
		wcscpy(outfile,OutDir->Text.c_str());
		wcscat(outfile, L"\\");
		wcscat(outfile, OutputFile->Text.c_str());
	}
	else {
		wcscpy(outfile,OutputFile->Text.c_str());
	}

	UnicodeString file=FilePath(outfile);
	UnicodeString cmd=L"gsiplot \""+file+"\"";
	if (!ExecCmd(cmd,1)) ShowMsg(L"error : gsiplot execution");
}
// callback on button-view --------------------------------------------------
void __fastcall TMainForm::BtnViewClick(TObject *Sender)
{
	wchar_t outfile[1024]={0};

	if (OutDirEna->Checked == true) {
		wcscpy(outfile,OutDir->Text.c_str());
		wcscat(outfile, L"\\");
		wcscat(outfile, OutputFile->Text.c_str());
	}
	else {
		wcscpy(outfile,OutputFile->Text.c_str());
	}
	ViewFile(FilePath(outfile));
}
// callback on button-to-kml ------------------------------------------------
void __fastcall TMainForm::BtnToKMLClick(TObject *Sender)
{
	wchar_t outfile[1024]={0};

	if (OutDirEna->Checked == true) {
		wcscpy(outfile,OutDir->Text.c_str());
		wcscat(outfile, L"\\");
		wcscat(outfile, OutputFile->Text.c_str());
	}
	else {
		wcscpy(outfile,OutputFile->Text.c_str());
	}

    ConvDialog->Show();
    ConvDialog->SetInput(FilePath(outfile));
}
// callback on button-options -----------------------------------------------
void __fastcall TMainForm::BtnOptionClick(TObject *Sender)
{
    int format=SolFormat;
    if (OptDialog->ShowModal()!=mrOk) return;
    if ((format==SOLF_NMEA)!=(SolFormat==SOLF_NMEA)) {
        SetOutFile();
    }
    UpdateEnable();
}
// callback on button-execute -----------------------------------------------
void __fastcall TMainForm::BtnExecClick(TObject *Sender)
{
	char *p,*wp;
    
	if (BtnExec->Caption=="Abort") {
        BtnExec->Enabled=false;
        return;
    }
    if (InputFile1->Text=="") {
        showmsg("error : no rinex obs file (rover)");
        return;
    }
	if (InputFile2->Text==""&&PMODE_DGPS<=PosMode&&PosMode<=PMODE_FIXED) {
		showmsg("error : no rinex obs file (base station)");
		return;
    }
    if (OutputFile->Text=="") {
        showmsg("error : no output file");
        return;
    }

	wp = wchar2char(OutputFile->Text.c_str(),-1);
	p=strrchr(wp,'.');
    if (p!=NULL) {
		if (!strcmp(p,".obs")||!strcmp(p,".OBS")||!strcmp(p,".nav")||
			!strcmp(p,".NAV")||!strcmp(p,".gnav")||!strcmp(p,".GNAV")||
			!strcmp(p,".gz")||!strcmp(p,".Z")||
			!strcmp(p+3,"o")||!strcmp(p+3,"O")||!strcmp(p+3,"d")||
			!strcmp(p+3,"D")||!strcmp(p+3,"n")||!strcmp(p+3,"N")||
			!strcmp(p+3,"g")||!strcmp(p+3,"G")) {
            showmsg("error : invalid extension of output file (%s)",p);
			if(wp!=NULL){
				free(wp);
				wp=NULL;
			}
            return;
        }
    }
	if(wp!=NULL){
		free(wp);
		wp=NULL;
	}
    showmsg("");
    BtnExec  ->Caption=L"Abort";
    BtnExit  ->Enabled=false;
    BtnView  ->Enabled=false;
    BtnToKML ->Enabled=false;
    BtnPlot  ->Enabled=false;
    BtnOption->Enabled=false;
    Panel1   ->Enabled=false;
    
	if (ExecProc()>=0) {
        AddHist(InputFile1);
        AddHist(InputFile2);
        AddHist(InputFile3);
        AddHist(InputFile4);
        AddHist(InputFile5);
        AddHist(OutputFile);
    }
    if (wcsstr(Message->Caption.c_str(),L"processing")) {
        showmsg("done");
    }
    BtnExec  ->Caption="E&xecute";
    BtnExec  ->Enabled=true;
    BtnExit  ->Enabled=true;
    BtnView  ->Enabled=true;
    BtnToKML ->Enabled=true;
    BtnPlot  ->Enabled=true;
    BtnOption->Enabled=true;
    Panel1   ->Enabled=true;
}
// callback on button-abort -------------------------------------------------
void __fastcall TMainForm::BtnStopClick(TObject *Sender)
{
    showmsg("abort");
}
// callback on button-exit --------------------------------------------------
void __fastcall TMainForm::BtnExitClick(TObject *Sender)
{
    Close();
}
// callback on button-about -------------------------------------------------
void __fastcall TMainForm::BtnAboutClick(TObject *Sender)
{
    UnicodeString prog=PRGNAME;
#ifdef MKL
    prog+="_MKL";
#endif
    AboutDialog->About=prog;
    AboutDialog->IconIndex=1;
    AboutDialog->ShowModal();
}
// callback on button-time-1 ------------------------------------------------
void __fastcall TMainForm::BtnTime1Click(TObject *Sender)
{
    TimeDialog->Time=GetTime1();
    TimeDialog->ShowModal();
}
// callback on button-time-2 ------------------------------------------------
void __fastcall TMainForm::BtnTime2Click(TObject *Sender)
{
    TimeDialog->Time=GetTime2();
    TimeDialog->ShowModal();
}
// callback on button-inputfile-1 -------------------------------------------
void __fastcall TMainForm::BtnInputFile1Click(TObject *Sender)
{
    wchar_t file[1024],*p;
    
    OpenDialog->Title=L"RINEX OBS (Rover) File";
    OpenDialog->FileName=L"";
    OpenDialog->FilterIndex=2;
    if (!OpenDialog->Execute()) return;
    InputFile1->Text=OpenDialog->FileName;
    SetOutFile();
}
// callback on button-inputfile-2 -------------------------------------------
void __fastcall TMainForm::BtnInputFile2Click(TObject *Sender)
{
    OpenDialog->Title=L"RINEX OBS (Base Station) File";
    OpenDialog->FileName=L"";
    OpenDialog->FilterIndex=2;
    if (!OpenDialog->Execute()) return;
    InputFile2->Text=OpenDialog->FileName;
}
// callback on button-inputfile-3 -------------------------------------------
void __fastcall TMainForm::BtnInputFile3Click(TObject *Sender)
{
    OpenDialog->Title=L"RINEX NAV/CLK, SP3, IONEX or SBAS/EMS File";
    OpenDialog->FileName=L"";
    OpenDialog->FilterIndex=3;
    if (!OpenDialog->Execute()) return;
    InputFile3->Text=OpenDialog->FileName;
}
// callback on button-inputfile-4 -------------------------------------------
void __fastcall TMainForm::BtnInputFile4Click(TObject *Sender)
{
    OpenDialog->Title=L"RINEX NAV/CLK, SP3, IONEX or SBAS/EMS File";
    OpenDialog->FileName=L"";
    OpenDialog->FilterIndex=4;
    if (!OpenDialog->Execute()) return;
    InputFile4->Text=OpenDialog->FileName;
}
// callback on button-inputfile-5 -------------------------------------------
void __fastcall TMainForm::BtnInputFile5Click(TObject *Sender)
{
    OpenDialog->Title=L"RINEX NAV/CLK, SP3, IONEX or SBAS/EMS File";
    OpenDialog->FileName=L"";
    OpenDialog->FilterIndex=5;
    if (!OpenDialog->Execute()) return;
    InputFile5->Text=OpenDialog->FileName;
}
// callback on button-outputfile --------------------------------------------
void __fastcall TMainForm::BtnOutputFileClick(TObject *Sender)
{
	wchar_t *ret;

    SaveDialog->Title=L"Output File";
    OpenDialog->FileName=L"";
	if (!SaveDialog->Execute()) return;

	if (OutDirEna->Checked == true) {
		ret = wcsrchr( SaveDialog->FileName.c_str(), '\\' );
		OutputFile->Text=++ret;
	}
	else {
		OutputFile->Text=SaveDialog->FileName;
	}
}
// callback on button-inputview-1 -------------------------------------------
void __fastcall TMainForm::BtnInputView1Click(TObject *Sender)
{
    ViewFile(FilePath(InputFile1->Text));
}
// callback on button-inputview-2 -------------------------------------------
void __fastcall TMainForm::BtnInputView2Click(TObject *Sender)
{
    ViewFile(FilePath(InputFile2->Text));
}
// callback on button-inputview-3 -------------------------------------------
void __fastcall TMainForm::BtnInputView3Click(TObject *Sender)
{
    UnicodeString file=FilePath(InputFile3->Text);
    wchar_t f[1024];
    
    if (file==L"") {
        file=FilePath(InputFile1->Text);
        if (!ObsToNav(file.c_str(),f)) return;
        file=f;
    }
    ViewFile(file);
}
// callback on button-inputview-4 -------------------------------------------
void __fastcall TMainForm::BtnInputView4Click(TObject *Sender)
{
    ViewFile(FilePath(InputFile4->Text));
}
// callback on button-inputview-5 -------------------------------------------
void __fastcall TMainForm::BtnInputView5Click(TObject *Sender)
{
    ViewFile(FilePath(InputFile5->Text));
}
// callback on button-outputview-1 ------------------------------------------
void __fastcall TMainForm::BtnOutputView1Click(TObject *Sender)
{
	wchar_t outfile[1024]={0};

	if (OutDirEna->Checked == true) {
		wcscpy(outfile,OutDir->Text.c_str());
		wcscat(outfile, L"\\");
		wcscat(outfile, OutputFile->Text.c_str());
	}
	else {
		wcscpy(outfile,OutputFile->Text.c_str());
	}

	UnicodeString file=FilePath(outfile)+L".stat";
    FILE *fp=_wfopen(file.c_str(),L"r");
    if (fp) fclose(fp); else return;
    ViewFile(file);
}
// callback on button-outputview-2 ------------------------------------------
void __fastcall TMainForm::BtnOutputView2Click(TObject *Sender)
{
	wchar_t outfile[1024]={0};

	if (OutDirEna->Checked == true) {
		wcscpy(outfile,OutDir->Text.c_str());
		wcscat(outfile, L"\\");
		wcscat(outfile, OutputFile->Text.c_str());
	}
	else {
		wcscpy(outfile,OutputFile->Text.c_str());
	}

    UnicodeString file=FilePath(outfile)+L".trace";
    FILE *fp=_wfopen(file.c_str(),L"r");
    if (fp) fclose(fp); else return;
    ViewFile(file);
}
// callback on button-inputplot-1 -------------------------------------------
void __fastcall TMainForm::BtnInputPlot1Click(TObject *Sender)
{
    UnicodeString files[5],cmd;
    wchar_t navfile[1024];
    
    files[0]=FilePath(InputFile1->Text); /* obs rover */
    files[1]=FilePath(InputFile2->Text); /* obs base */
    files[2]=FilePath(InputFile3->Text);
    files[3]=FilePath(InputFile4->Text);
    files[4]=FilePath(InputFile5->Text);
    
    if (files[2]=="") {
        if (ObsToNav(files[0].c_str(),navfile)) files[2]=navfile;
    }
	cmd=L"gsiplot -r \""+files[0]+L"\" \""+files[2]+L"\" \""+files[3]+L"\" \""+
        files[4]+L"\"";
    
	if (!ExecCmd(cmd,1)) ShowMsg(L"error : gsiplot execution");
}
// callback on button-inputplot-2 -------------------------------------------
void __fastcall TMainForm::BtnInputPlot2Click(TObject *Sender)
{
    UnicodeString files[5],cmd;
    wchar_t navfile[1024],gnavfile[1024];
    
    files[0]=FilePath(InputFile1->Text); /* obs rover */
    files[1]=FilePath(InputFile2->Text); /* obs base */
    files[2]=FilePath(InputFile3->Text);
    files[3]=FilePath(InputFile4->Text);
    files[4]=FilePath(InputFile5->Text);
    
    if (files[2]==L"") {
        if (ObsToNav(files[0].c_str(),navfile)) files[2]=navfile;
    }
	cmd=L"gsiplot -r \""+files[1]+L"\" \""+files[2]+L"\" \""+files[3]+L"\" \""+
        files[4]+L"\"";
    
    if (!ExecCmd(cmd,1)) ShowMsg(L"error : gsiplot execution");
}
// callback on button-output-directory --------------------------------------
void __fastcall TMainForm::BtnOutDirClick(TObject *Sender)
{
    UnicodeString dir=OutDir->Text;
    if (!SelectDirectory(L"Output Directory",L"",dir)) return;
    OutDir->Text=dir;
}
// callback on button keyword -----------------------------------------------
void __fastcall TMainForm::BtnKeywordClick(TObject *Sender)
{
    KeyDialog->Flag=2;
    KeyDialog->Show();
}
// callback on time-start/end check -----------------------------------------
void __fastcall TMainForm::TimeStartClick(TObject *Sender)
{
    UpdateEnable();
}
// callback on time-interval check ------------------------------------------
void __fastcall TMainForm::TimeIntFClick(TObject *Sender)
{
    UpdateEnable();
}
// callback on time-unit check ----------------------------------------------
void __fastcall TMainForm::TimeUnitFClick(TObject *Sender)
{
    UpdateEnable();
}
// callback on time-ymd-1 updown --------------------------------------------
void __fastcall TMainForm::TimeY1UDChangingEx(TObject *Sender,
      bool &AllowChange, short NewValue, TUpDownDirection Direction)
{
    UnicodeString s;
    double ep[]={2000,1,1,0,0,0};
    int p=TimeY1->SelStart,ud=Direction==updUp?1:-1;
    
    swscanf(TimeY1->Text.c_str(),L"%lf/%lf/%lf",ep,ep+1,ep+2);
    if (4<p&&p<8) {
        ep[1]+=ud;
        if (ep[1]<=0) {ep[0]--; ep[1]+=12;}
        else if (ep[1]>12) {ep[0]++; ep[1]-=12;}
    }
    else if (p>7||p==0) ep[2]+=ud; else ep[0]+=ud;
    time2epoch(epoch2time(ep),ep);
    TimeY1->Text=s.sprintf(L"%04.0f/%02.0f/%02.0f",ep[0],ep[1],ep[2]);
    TimeY1->SelStart=p>7||p==0?10:(p>4?7:4);
}
// callback on time-hms-1 updown --------------------------------------------
void __fastcall TMainForm::TimeH1UDChangingEx(TObject *Sender,
      bool &AllowChange, short NewValue, TUpDownDirection Direction)
{
    UnicodeString s;
    int hms[3]={0},sec,p=TimeH1->SelStart,ud=Direction==updUp?1:-1;
    
    swscanf(TimeH1->Text.c_str(),L"%d:%d:%d",hms,hms+1,hms+2);
    if (p>5||p==0) hms[2]+=ud; else if (p>2) hms[1]+=ud; else hms[0]+=ud;
    sec=hms[0]*3600+hms[1]*60+hms[2];
    if (sec<0) sec+=86400; else if (sec>=86400) sec-=86400;
    TimeH1->Text=s.sprintf(L"%02d:%02d:%02d",sec/3600,(sec%3600)/60,sec%60);
    TimeH1->SelStart=p>5||p==0?8:(p>2?5:2);
}
// callback on time-ymd-2 updown --------------------------------------------
void __fastcall TMainForm::TimeY2UDChangingEx(TObject *Sender,
      bool &AllowChange, short NewValue, TUpDownDirection Direction)
{
    UnicodeString s;
    double ep[]={2000,1,1,0,0,0};
    int p=TimeY2->SelStart,ud=Direction==updUp?1:-1;
    
    swscanf(TimeY2->Text.c_str(),L"%lf/%lf/%lf",ep,ep+1,ep+2);
    if (4<p&&p<8) {
        ep[1]+=ud;
        if (ep[1]<=0) {ep[0]--; ep[1]+=12;}
        else if (ep[1]>12) {ep[0]++; ep[1]-=12;}
    }
    else if (p>7||p==0) ep[2]+=ud; else ep[0]+=ud;
    time2epoch(epoch2time(ep),ep);
    TimeY2->Text=s.sprintf(L"%04.0f/%02.0f/%02.0f",ep[0],ep[1],ep[2]);
    TimeY2->SelStart=p>7||p==0?10:(p>4?7:4);
}
// callback on time-hms-2 updown --------------------------------------------
void __fastcall TMainForm::TimeH2UDChangingEx(TObject *Sender,
      bool &AllowChange, short NewValue, TUpDownDirection Direction)
{
    UnicodeString s;
    int hms[3]={0},sec,p=TimeH2->SelStart,ud=Direction==updUp?1:-1;
    
    swscanf(TimeH2->Text.c_str(),L"%d:%d:%d",hms,hms+1,hms+2);
    if (p>5||p==0) hms[2]+=ud; else if (p>2) hms[1]+=ud; else hms[0]+=ud;
    sec=hms[0]*3600+hms[1]*60+hms[2];
    if (sec<0) sec+=86400; else if (sec>=86400) sec-=86400;
    TimeH2->Text=s.sprintf(L"%02d:%02d:%02d",sec/3600,(sec%3600)/60,sec%60);
    TimeH2->SelStart=p>5||p==0?8:(p>2?5:2);
}
// callback on inputfile-1 change -------------------------------------------
void __fastcall TMainForm::InputFile1Change(TObject *Sender)
{
    SetOutFile();
}
// callback on output-directory checked -------------------------------------
void __fastcall TMainForm::OutDirEnaClick(TObject *Sender)
{
	UpdateEnable();
	SetOutFile();
}
// callback on output-directory change --------------------------------------
void __fastcall TMainForm::OutDirChange(TObject *Sender)
{
    SetOutFile();
}
// set output file path -----------------------------------------------------
void __fastcall TMainForm::SetOutFile(void)
{
    wchar_t *p,ofile[1024],ifile[1024];
    FILE *fp;
    if (InputFile1->Text=="") return;
    
	wcscpy(ifile,InputFile1->Text.c_str());
    
    if (OutDirEna->Checked) {
        p=wcsrchr(ifile,'\\');
        if (p!=NULL) p++; else p=ifile;
		swprintf(ofile,L"%s",p);
    }
    else {
        wcscpy(ofile,ifile);
	}
    p=wcsrchr(ofile,'.');
    if (p==NULL) p=ofile+wcslen(ofile);
	wcscpy(p,SolFormat==SOLF_NMEA?L".nmea":L".pos");
    for (p=ofile;*p;p++) if (*p=='*') *p='0';
	OutputFile->Text=ofile;
}
// execute post-processing --------------------------------------------------
int __fastcall TMainForm::ExecProc(void)
{
    FILE *fp;
    prcopt_t prcopt=prcopt_default;
    solopt_t solopt=solopt_default;
    filopt_t filopt={""};
    gtime_t ts={0},te={0};
    double ti=0.0,tu=0.0;
	int i,n=0,stat;
	wchar_t infile_[5][1024]={L""},*infile[5],outfile[1024]={0};
    wchar_t *rov,*base,*p,*q,*r;
    
	prcopt.mopt=mbsopt_default;

	// get processing options
    if (TimeStart->Checked) ts=GetTime1();
    if (TimeEnd  ->Checked) te=GetTime2();
    if (TimeIntF ->Checked) ti=str2dbl(TimeInt ->Text);
    if (TimeUnitF->Checked) tu=str2dbl(TimeUnit->Text)*3600.0;
    
	if (!GetOption(prcopt,solopt,filopt)) return 0;
    
    // set input/output files
    for (i=0;i<5;i++) infile[i]=infile_[i];
    
	wcscpy(infile[n++],InputFile1->Text.c_str());
    
	if ((PMODE_DGPS<=prcopt.mode&&prcopt.mode<=PMODE_FIXED) ||prcopt.mode==PMODE_MULTI) {
		wcscpy(infile[n++],InputFile2->Text.c_str());
	}
    if (InputFile3->Text!=L"") {
        wcscpy(infile[n++],InputFile3->Text.c_str());
    }
	else if ((prcopt.navsys&SYS_GPS)&&!ObsToNav(InputFile1->Text.c_str(),infile[n++])) {
		showmsg("error: no gps navigation data");
		return 0;
	}
    if (InputFile4->Text!="") {
        wcscpy(infile[n++],InputFile4->Text.c_str());
    }
    if (InputFile5->Text!="") {
        wcscpy(infile[n++],InputFile5->Text.c_str());
    }
	if (OutDirEna->Checked == true) {
		wcscpy(outfile,OutDir->Text.c_str());
		wcscat(outfile, L"\\");
		wcscat(outfile, OutputFile->Text.c_str());
	}
	else {
    wcscpy(outfile,OutputFile->Text.c_str());
	}

    // confirm overwrite
    if (!TimeStart->Checked||!TimeEnd->Checked) {
        fp=_wfopen(outfile,L"r");
		if (fp!=NULL) {
            fclose(fp);
            ConfDialog->Label2->Caption=outfile;
            if (ConfDialog->ShowModal()!=mrOk) return 0;
        }
    }
    // set rover and base station list
	rov =new wchar_t [wcslen(RovList .c_str())];
	base=new wchar_t [wcslen(BaseList.c_str())];

	for (p=RovList.c_str(),r=rov;*p;p=q+2) {

		q=wcsstr(p,L"\r\n");
        if (q==NULL) {
            if (*p!='#') wcscpy(r,p); break;
        }
        else if (*p!='#') {
            wcsncpy(r,p,q-p); r+=q-p;
			wcscpy(r++,L" ");
        }
    }
    for (p=BaseList.c_str(),r=base;*p;p=q+2) {
        
        q=wcsstr(p,L"\r\n");
        if (q==NULL) {
			if (*p!='#') wcscpy(r,p); break;
        }
        else if (*p!='#') {
            wcsncpy(r,p,q-p); r+=q-p;
			wcscpy(r++,L" ");
        }
    }
	Progress->Position=0;
    showmsg("reading...");
    
    // post processing positioning
    char *cinfile[5]={""};
    for(int i=0;i<5;i++){
        cinfile[i] = wchar2char(infile[i],1024);
    }

	char *coutfile = wchar2char(outfile,1024);
	char *ch3 = wchar2char(rov,-1);
	char *ch4 = wchar2char(base,-1);

	if ((stat=postpos(ts,te,ti,tu,&prcopt,&solopt,&filopt,cinfile,n,coutfile,
					  ch3,ch4))==1) {
		showmsg("aborted");
	}
	delete [] rov ;
	delete [] base;

	for(int i=0;i<5;i++){
		if(cinfile[i]!=NULL){
			free(cinfile[i]);
			cinfile[i]=NULL;
		}
	}
	if(coutfile!=NULL){
		free(coutfile);
		coutfile=NULL;
	}
	if(ch3!=NULL){
		free(ch3);
		ch3=NULL;
	}
	if(ch4!=NULL){
		free(ch4);
		ch4=NULL;
	}
	return stat;
}
// get processing and solution options --------------------------------------
int __fastcall TMainForm::GetOption(prcopt_t &prcopt, solopt_t &solopt,
                                    filopt_t &filopt)
{
//	char buff[1024],id[32],*p;
	char *buff,id[32],*p;
	int sat,ex;

    // processing options
	prcopt.mode     =PosMode;
	prcopt.soltype  =Solution;
	if(Freqs==1) {
		prcopt.nfreq=1;
		prcopt.oprfrq[0]=0;
		prcopt.oprfrq[1]=-1;
	}
	else if(Freqs==5) {
		prcopt.nfreq=1;
		prcopt.oprfrq[0]=2;
		prcopt.oprfrq[1]=-1;
	}
	else if(Freqs==12) {
		prcopt.nfreq=2;
		prcopt.oprfrq[0]=0;
		prcopt.oprfrq[1]=1;
		prcopt.oprfrq[2]=-1;
	}
	else if(Freqs==15) {
		prcopt.nfreq=2;
		prcopt.oprfrq[0]=0;
		prcopt.oprfrq[1]=2;
		prcopt.oprfrq[2]=-1;
	}
	else if(Freqs==125) {
		prcopt.nfreq=3;
		prcopt.oprfrq[0]=0;
		prcopt.oprfrq[1]=1;
		prcopt.oprfrq[2]=2;
		prcopt.oprfrq[3]=-1;
	}
	prcopt.navsys   =NavSys;

	char *cGloCodePri1 = wchar2char(GloCodePri1.c_str(),MAXFCODE+2);
	strcpy(prcopt.codepri[ISYSGLO][0],cGloCodePri1);
	if(cGloCodePri1!=NULL){
		free(cGloCodePri1);
		cGloCodePri1=NULL;
	}
	char *cGloCodePri2 = wchar2char(GloCodePri2.c_str(),MAXFCODE+2);
	strcpy(prcopt.codepri[ISYSGLO][1],cGloCodePri2);
	if(cGloCodePri2!=NULL){
		free(cGloCodePri2);
		cGloCodePri2=NULL;
	}

    prcopt.elmin    =ElMask*D2R;
    prcopt.snrmask  =SnrMask;
    prcopt.sateph   =SatEphem;
    prcopt.modertkar=AmbResMethod;
    prcopt.modear   =AmbRes;
    prcopt.glomodear=GloAmbRes;
	prcopt.modepppar=PppAmbRes;
	prcopt.maxout   =OutCntResetAmb;
	prcopt.minfix   =FixCntHoldAmb;
    prcopt.minlock  =LockCntFixAmb;
    prcopt.ionoopt  =IonoOpt;
	prcopt.tropopt  =TropOpt;
	prcopt.posopt[0]=PosOpt[0];
	prcopt.posopt[1]=PosOpt[1];
	prcopt.posopt[2]=PosOpt[2];
	prcopt.posopt[3]=PosOpt[3];
    prcopt.posopt[4]=PosOpt[4];
    prcopt.dynamics =DynamicModel;
    prcopt.tidecorr =TideCorr;
    prcopt.niter    =NumIter;
    prcopt.intpref  =IntpRefObs;
    prcopt.sbassatsel=SbasSat;
	prcopt.eratio[0]=MeasErrR1;
	prcopt.eratio[1]=MeasErrR2;
	prcopt.eratio[2]=MeasErrR3;
	prcopt.err[1]   =MeasErr2;
	prcopt.err[2]   =MeasErr3;
    prcopt.err[3]   =MeasErr4;
    prcopt.err[4]   =MeasErr5;
    prcopt.prn[0]   =PrNoise1;
    prcopt.prn[1]   =PrNoise2;
    prcopt.prn[2]   =PrNoise3;
	prcopt.prn[3]   =PrNoise4;
	prcopt.prn[4]   =PrNoise5;
	prcopt.prn[5]   =PrNoise6;
	prcopt.prn[6]   =PrNoise7;
	prcopt.prn[7]   =PrNoise8;
	prcopt.sclkstab =SatClkStab;
	prcopt.thresar[0]=ValidThresAR;
	prcopt.thresar[1]=ThresAR2;
	prcopt.thresar[2]=ThresAR3;
    prcopt.elmaskar =ElMaskAR*D2R;
    prcopt.elmaskhold=ElMaskHold*D2R;
    prcopt.thresslip=SlipThres;
    prcopt.maxtdiff =MaxAgeDiff;
    prcopt.maxgdop  =RejectGdop;
    prcopt.maxinno  =RejectThres;
    if (BaseLineConst) {
        prcopt.baseline[0]=BaseLine[0];
        prcopt.baseline[1]=BaseLine[1];
    }
    else {
        prcopt.baseline[0]=0.0;
        prcopt.baseline[1]=0.0;
    }
	if (PosMode!=PMODE_FIXED&&PosMode!=PMODE_PPP_FIXED) {
        for (int i=0;i<3;i++) prcopt.ru[i]=0.0;
    }
    else if (RovPosType<=2) {
		for (int i=0;i<3;i++) prcopt.ru[i]=RovPos[i];
    }
    else prcopt.rovpos=RovPosType-2; /* 1:single,2:posfile,3:rinex */
    
    if (PosMode==PMODE_SINGLE||PosMode==PMODE_MOVEB) {
        for (int i=0;i<3;i++) prcopt.rb[i]=0.0;
    }
    else if (RefPosType<=2) {
        for (int i=0;i<3;i++) prcopt.rb[i]=RefPos[i];
    }
    else prcopt.refpos=RefPosType-2;
    
    if (RovAntPcv) {
        char *cRovAnt = wchar2char(RovAnt.c_str(),MAXANT);
        strcpy(prcopt.anttype[0],cRovAnt);
		if(cRovAnt!=NULL){
			free(cRovAnt);
			cRovAnt=NULL;
		}
		prcopt.antdel[0][0]=RovAntE;
        prcopt.antdel[0][1]=RovAntN;
		prcopt.antdel[0][2]=RovAntU;
    }
    if (RefAntPcv) {
        char *cRefAnt = wchar2char(RefAnt.c_str(),MAXANT);
        strcpy(prcopt.anttype[1],cRefAnt);
		if(cRefAnt!=NULL){
			free(cRefAnt);
			cRefAnt=NULL;
		}
		prcopt.antdel[1][0]=RefAntE;
        prcopt.antdel[1][1]=RefAntN;
        prcopt.antdel[1][2]=RefAntU;
    }
	if (ExSats!="") { // excluded satellites

		buff = wchar2char(ExSats.c_str(),1024);

		for (p=strtok(buff," ");p;p=strtok(NULL," ")) {
			if (*p=='+') {ex=2; p++;} else ex=1;
            if (!(sat=satid2no(p))) continue;
            prcopt.exsats[sat-1]=ex;
        }

		if(buff!=NULL){
			free(buff);
			buff=NULL;
		}
	}

	prcopt.l2cprior=L2Cod;/* L2 code priority (0:L2P,1:L2C) */
	prcopt.tsyscorr=TimSys;/* time system correction ( 0:off,1:on) */
	prcopt.phasshft=PhaCyc;/* phase cycle shift (0:off,1:table) */
	prcopt.gl2bias=L2CPBias;
	prcopt.errmodel=ErrMod;/* error model (0:user settings,1:table) */

	prcopt.dcbratio=MeasErr6;/* code error ratio(no DCB) */
	prcopt.isb=Isb;/* phase cycle shift (0:off,1:table) */
	prcopt.diff=Diff;/* (0:inall,1:exc-glo) */

	/*BIPMCircularTFile*/
	char *wcirtfile = wchar2char(BIPMCircularTFile.c_str(),MAXSTRPATH);
	strcpy(filopt.cirtfile,wcirtfile);
	if(wcirtfile!=NULL){
		free(wcirtfile);
		wcirtfile=NULL;
	}
	/* receiver types {rover,base} */
	char *wrectype0 = wchar2char(RovRecTyp.c_str(),MAXANT);
	strcpy(prcopt.rectype[0],wrectype0);
	if(wrectype0!=NULL){
		free(wrectype0);
		wrectype0=NULL;
	}
	char *wrectype1 = wchar2char(RefRecTyp.c_str(),MAXANT);
	strcpy(prcopt.rectype[1],wrectype1);
	if(wrectype1!=NULL){
		free(wrectype1);
		wrectype1=NULL;
	}
	/* 1/4cycle phase correction file */
	char *wL2Ccorr = wchar2char(PhaCycFile.c_str(),MAXSTRPATH);
	strcpy(prcopt.mopt.ifpcs,wL2Ccorr);
	if(wL2Ccorr!=NULL){
		free(wL2Ccorr);
		wL2Ccorr=NULL;
	}
	/* GLONASS IFB table file*/
	char *wgloifb = wchar2char(GloIfbFile.c_str(),MAXSTRPATH);
	strcpy(prcopt.mopt.ififb,wgloifb);
	if(wgloifb!=NULL){
		free(wgloifb);
		wgloifb=NULL;
	}
	/* error model file */
	char *werrmodel = wchar2char(ErrModFile.c_str(),MAXSTRPATH);
	strcpy(prcopt.mopt.iferr,werrmodel);
	if(werrmodel!=NULL){
		free(werrmodel);
		werrmodel=NULL;
	}

	/*Estimate Satellite Clock/FCB*/
	prcopt.mopt.estsatclk = EstSatClo;
	prcopt.mopt.estsatfcb = EstSatFCB;

	/*Semi-Dynamic Correction Parameter*/
	char *wSemiDCPara = wchar2char(SemiDCPara.c_str(),MAXSTRPATH);
	strcpy(prcopt.mopt.ifsdp,wSemiDCPara);
	if(wSemiDCPara!=NULL){
		free(wSemiDCPara);
		wSemiDCPara=NULL;
	}
	/*Solution Dir*/
	char *wSolDirFile = wchar2char(SolDirFile.c_str(),MAXSTRPATH);
	strcpy(prcopt.mopt.ofdir,wSolDirFile);
	if(wSolDirFile!=NULL){
		free(wSolDirFile);
		wSolDirFile=NULL;
	}
	/*Satellite FCB*/
	char *wFCBFile = wchar2char(FCBFile.c_str(),MAXSTRPATH);
	strcpy(prcopt.mopt.iffcb,wFCBFile);
	if(wFCBFile!=NULL){
		free(wFCBFile);
		wFCBFile=NULL;
	}
	prcopt.mopt.tiztd      = EstIntZ;
	prcopt.mopt.tigra      = EstIntES;
	prcopt.mopt.sigtrop[0] = RanWalSigZen;
	prcopt.mopt.sigtrop[1] = RanWalSigEW;
	prcopt.mopt.sigtrop[2] = RanWalSigNS;
	prcopt.mopt.cpomcth    = ThrO;
	prcopt.mopt.promcth    = ThrC;
	prcopt.mopt.maxdddist  = MaxBasDD;
	prcopt.mopt.wlddfix    = JudValWL;
	prcopt.mopt.l1ddfix    = JudValL1;
	prcopt.mopt.wdd        = WeiDD;
	prcopt.mopt.itrconv    = ConCri;
	prcopt.mopt.itrmax     = MaxIte;

	/*Temporary storage epoch parameters*/
	char *wTemStoFile = wchar2char(TemStoFile.c_str(),MAXSTRPATH);
	strcpy(prcopt.mopt.eptmp,wTemStoFile);
	if(wTemStoFile!=NULL){
		free(wTemStoFile);
		wTemStoFile=NULL;
	}
	prcopt.mopt.tifcb         = NLFcb;
	prcopt.mopt.minpass       = Minsd;
	prcopt.mopt.minddpass     = Mindd;
	prcopt.mopt.maxsigw       = Maxsd;
	prcopt.mopt.minconfw      = Fixwl;
	prcopt.mopt.minconf1      = Fixnl;
	prcopt.mopt.sigr[0]       = Mobstadn;
	prcopt.mopt.sigr[1]       = Mobstade;
	prcopt.mopt.sigr[2]       = Mobstadu;
	prcopt.mopt.sigb[0]       = Basestadn;
	prcopt.mopt.sigb[1]       = Basestade;
	prcopt.mopt.sigb[2]       = Basestadu;

	for (int i=0;i<7;i++) prcopt.std[i]=std[i];
	prcopt.minsat = minsat;

	// solution options
	solopt.posf     =SolFormat;
	solopt.times    =TimeFormat==0?0:TimeFormat-1;
	solopt.timef    =TimeFormat==0?0:1;
	solopt.timeu    =TimeDecimal<=0?0:TimeDecimal;
	solopt.degf     =LatLonFormat;
	solopt.outhead  =OutputHead;
	solopt.outopt   =OutputOpt;
	solopt.datum    =OutputDatum;
    solopt.height   =OutputHeight;
    solopt.geoid    =OutputGeoid;
	solopt.solstatic=SolStatic;
	solopt.sstat    =DebugStatus;
	solopt.trace    =DebugTrace;
	solopt.isbout   =IsbOut;
	solopt.gl2out   =Gl2Out;
	solopt.fcbout   =FCBOut;
	solopt.possnxout=PossnxOut;
	solopt.ionout   =IonOut;
	solopt.tropout  =TropOut;
	solopt.recclout =RecClOut;
	solopt.satclout =SatClOut;
	solopt.sbresout =StaticOut;

	wchar_t wFieldSep[64];
    memset(wFieldSep,0,64);
    wcscpy(wFieldSep,FieldSep!=L""?FieldSep.c_str():L" ");
    
    char *cFieldSep = wchar2char(wFieldSep,64);
    
    strcpy(solopt.sep,cFieldSep);
	if(cFieldSep!=NULL){
		free(cFieldSep);
		cFieldSep=NULL;
	}
    char *cch = wchar2char(PRGNAME,-1);
    sprintf(solopt.prog,"%s ver.%s",cch,VER_RTKLIB);
	if(cch!=NULL){
		free(cch);
		cch=NULL;
	}
    // file options
	char *wsatantp = wchar2char(SatPcvFile.c_str(),MAXSTRPATH);
	strcpy(filopt.satantp,wsatantp);
	if(wsatantp!=NULL){
		free(wsatantp);
		wsatantp=NULL;
	}
	char *wrcvantp = wchar2char(AntPcvFile.c_str(),MAXSTRPATH);
	strcpy(filopt.rcvantp,wrcvantp);
	if(wrcvantp!=NULL){
		free(wrcvantp);
		wrcvantp=NULL;
	}
	char *wstapos = wchar2char(StaPosFile.c_str(),MAXSTRPATH);
	strcpy(filopt.stapos,wstapos);
	if(wstapos!=NULL){
		free(wstapos);
		wstapos=NULL;
	}
	char *wgeoid = wchar2char(GeoidDataFile.c_str(),MAXSTRPATH);
	strcpy(filopt.geoid,wgeoid);
	if(wgeoid!=NULL){
		free(wgeoid);
		wgeoid=NULL;
	}
    char *wiono = wchar2char(IonoFile.c_str(),MAXSTRPATH);
    strcpy(filopt.iono,wiono);
	if(wiono!=NULL){
		free(wiono);
		wiono=NULL;
	}
	char *wdcb = wchar2char(DCBFile.c_str(),MAXSTRPATH);
    strcpy(filopt.dcb,wdcb);
	if(wdcb!=NULL){
		free(wdcb);
		wdcb=NULL;
	}
	char *wisb = wchar2char(ISBFile.c_str(),MAXSTRPATH);
	strcpy(filopt.isb,wisb);
	if(wisb!=NULL){
		free(wisb);
		wisb=NULL;
	}
	char *wisbout = wchar2char(ISBOutFile.c_str(),MAXSTRPATH);
	strcpy(solopt.isbfile,wisbout);
	if(wisbout!=NULL){
		free(wisbout);
		wisbout=NULL;
	}
	char *wgl2out = wchar2char(GL2OutFile.c_str(),MAXSTRPATH);
	strcpy(solopt.gl2file,wgl2out);
	if(wgl2out!=NULL){
		free(wgl2out);
		wgl2out=NULL;
	}
	char *wfcbout = wchar2char(FCBOutFile.c_str(),MAXSTRPATH);
	strcpy(solopt.fcb,wfcbout);
	if(wfcbout!=NULL){
		free(wfcbout);
		wfcbout=NULL;
	}
	char *weop = wchar2char(EOPFile.c_str(),MAXSTRPATH);
	strcpy(filopt.eop,weop);
	if(weop!=NULL){
		free(weop);
		weop=NULL;
	}
	char *wblq = wchar2char(BLQFile.c_str(),MAXSTRPATH);
	strcpy(filopt.blq,wblq);
	if(wblq!=NULL){
		free(wblq);
		wblq=NULL;
	}
    return 1;
}
// observation file to nav file ---------------------------------------------
int __fastcall TMainForm::ObsToNav(const wchar_t *obsfile, wchar_t *navfile)
{
    wchar_t *p;
    wcscpy(navfile,obsfile);
    p=wcsrchr(navfile,'.');
    if (p==NULL) return 0;
    if      (wcslen(p)==4&&*(p+3)=='o') *(p+3)='*';
    else if (wcslen(p)==4&&*(p+3)=='d') *(p+3)='*';
    else if (wcslen(p)==4&&*(p+3)=='O') *(p+3)='*';
    else if (!wcscmp(p,L".obs")) wcscpy(p,L".*nav");
    else if (!wcscmp(p,L".OBS")) wcscpy(p,L".*NAV");
    else if (!wcscmp(p,L".gz")||!wcscmp(p,L".Z")) {
        if      (*(p-1)=='o') *(p-1)='*';
        else if (*(p-1)=='d') *(p-1)='*';
        else if (*(p-1)=='O') *(p-1)='*';
        else return 0;
    }
    else return 0;
    return 1;
}
// replace file path with keywords ------------------------------------------
UnicodeString __fastcall TMainForm::FilePath(UnicodeString file)
{
    UnicodeString s;
    gtime_t ts={0};
    wchar_t wrov[256]=L"",wbase[256]=L"",*wpath,*p,*q;
    char path[1024]="";
    char *rov, *base;
    if (TimeStart->Checked) ts=GetTime1();
    
    p=RovList.c_str();
    q=wcsstr(p,L"\r\n");
    while(q!=NULL){
        if (*p&&*p!='#') break;
        p=q+2;
        q=wcsstr(p,L"\r\n");
    }
    if (!q) {
        rov = wchar2char(p,256);
    }else{
        wcsncpy(wrov,p,q-p);
        rov = wchar2char(wrov,256);
    }
    
    
    p=BaseList.c_str();
    q=wcsstr(p,L"\r\n");
	while(q!=NULL){
		if (*p&&p[0]!='#') break;
		p=q+2;
		q=wcsstr(p,L"\r\n");
	}

	if (!q) {
		base = wchar2char(p,256);
	}else{
		wcsncpy(wbase,p,q-p);
		base = wchar2char(wbase,256);
	}

	char *cfile = wchar2char(file.c_str(),-1);
	reppath(cfile,path,ts,rov,base);

	if(rov!=NULL){
		free(rov);
		rov=NULL;
	}
	if(base!=NULL){
		free(base);
		base=NULL;
	}
	if(cfile!=NULL){
		free(cfile);
		cfile=NULL;
	}
	wpath = char2wchar(path,1024);
	s=wpath;
	if(wpath!=NULL){
		free(wpath);
		wpath=NULL;
	}
	return s;
}
// read history -------------------------------------------------------------
TStringList * __fastcall TMainForm::ReadList(TIniFile *ini, UnicodeString cat,
    UnicodeString key)
{
    TStringList *list=new TStringList;
    UnicodeString s,item;
    int i;
    
    for (i=0;i<100;i++) {
        item=ini->ReadString(cat,s.sprintf(L"%s_%03d",key.c_str(),i),"");
        if (item!="") list->Add(item); else break;
    }
    return list;
}
// write history ------------------------------------------------------------
void __fastcall TMainForm::WriteList(TIniFile *ini, UnicodeString cat,
    UnicodeString key, TStrings *list)
{
    UnicodeString s;
    int i;
    
    for (i=0;i<list->Count;i++) {
        ini->WriteString(cat,s.sprintf(L"%s_%03d",key.c_str(),i),list->Strings[i]);
    }
}
// add history --------------------------------------------------------------
void __fastcall TMainForm::AddHist(TComboBox *combo)
{
    UnicodeString hist=combo->Text;
    if (hist=="") return;
    TStrings *list=combo->Items;
    int i=list->IndexOf(hist);
    if (i>=0) list->Delete(i);
    list->Insert(0,hist);
    for (int i=list->Count-1;i>=MAXHIST;i--) list->Delete(i);
    combo->ItemIndex=0;
}
// execute command ----------------------------------------------------------
int __fastcall TMainForm::ExecCmd(UnicodeString cmd, int show)
{
    PROCESS_INFORMATION info;
    STARTUPINFOW si={0};
    si.cb=sizeof(si);
    wchar_t *p=cmd.c_str();
    
    if (!CreateProcessW(NULL,p,NULL,NULL,false,show?0:CREATE_NO_WINDOW,NULL,
                       NULL,&si,&info)) return 0;
    CloseHandle(info.hProcess);
    CloseHandle(info.hThread);
    return 1;
}
// view file ----------------------------------------------------------------
void __fastcall TMainForm::ViewFile(UnicodeString file)
{
    TTextViewer *viewer;
    UnicodeString f;
    wchar_t *tmpfile;
    int cstat;
    
    if (file.IsEmpty()) return;
    
    char *ch;
    ch = wchar2char(file.c_str(), -1);
    
    char *ch2= (char *)malloc(1024);
    if (ch2==NULL) return;
    
    cstat=uncompress(ch,ch2);
	if(ch!=NULL){
		free(ch);
		ch=NULL;
	}
    
    tmpfile = char2wchar(ch2, -1);
    
    f=!cstat?file.c_str():tmpfile;
    
    viewer=new TTextViewer(Application);
    viewer->Caption=file;
    viewer->Show();
    viewer->Read(f);
    if (cstat==1) _wremove(tmpfile);
	if(tmpfile!=NULL){
		free(tmpfile);
		tmpfile=NULL;
	}
}
// show message in message area ---------------------------------------------
void __fastcall TMainForm::ShowMsg(wchar_t *msg)
{
	Message->Caption=msg;
	Message->Update();
}
// get time from time-1 -----------------------------------------------------
gtime_t _fastcall TMainForm::GetTime1(void)
{
    double ep[]={2000,1,1,0,0,0};
    
    swscanf(TimeY1->Text.c_str(),L"%lf/%lf/%lf",ep,ep+1,ep+2);
    swscanf(TimeH1->Text.c_str(),L"%lf:%lf:%lf",ep+3,ep+4,ep+5);
    return epoch2time(ep);
}
// get time from time-2 -----------------------------------------------------
gtime_t _fastcall TMainForm::GetTime2(void)
{
    double ep[]={2000,1,1,0,0,0};
    
    swscanf(TimeY2->Text.c_str(),L"%lf/%lf/%lf",ep,ep+1,ep+2);
    swscanf(TimeH2->Text.c_str(),L"%lf:%lf:%lf",ep+3,ep+4,ep+5);
    return epoch2time(ep);
}
// set time to time-1 -------------------------------------------------------
void _fastcall TMainForm::SetTime1(gtime_t time)
{
    UnicodeString s;
    double ep[6];
    
    time2epoch(time,ep);
    TimeY1->Text=s.sprintf(L"%04.0f/%02.0f/%02.0f",ep[0],ep[1],ep[2]);
    TimeH1->Text=s.sprintf(L"%02.0f:%02.0f:%02.0f",ep[3],ep[4],ep[5]);
    TimeY1->SelStart=10; TimeH1->SelStart=10;
}
// set time to time-2 -------------------------------------------------------
void _fastcall TMainForm::SetTime2(gtime_t time)
{
    UnicodeString s;
    double ep[6];
    
    time2epoch(time,ep);
    TimeY2->Text=s.sprintf(L"%04.0f/%02.0f/%02.0f",ep[0],ep[1],ep[2]);
    TimeH2->Text=s.sprintf(L"%02.0f:%02.0f:%02.0f",ep[3],ep[4],ep[5]);
    TimeY2->SelStart=10; TimeH2->SelStart=10;
}
// update enable/disable of widgets -----------------------------------------
void __fastcall TMainForm::UpdateEnable(void)
{
	int moder=0;
	if((PMODE_DGPS<=PosMode&&PosMode<=PMODE_FIXED) || PosMode==PMODE_MULTI){
		moder=PosMode;
	}
    
    LabelInputFile1->Caption=moder?L"RINEX OBS: Rover":L"RINEX OBS";
    InputFile2     ->Enabled=moder;
    BtnInputFile2  ->Enabled=moder;
    BtnInputPlot2  ->Enabled=moder;
    BtnInputView2  ->Enabled=moder;
	BtnOutputView1 ->Enabled=DebugStatus>0;
	BtnOutputView2 ->Enabled=DebugTrace >0;
    LabelInputFile3->Enabled=moder;
	TimeY1         ->Enabled=TimeStart->Checked;
    TimeH1         ->Enabled=TimeStart->Checked;
    TimeY1UD       ->Enabled=TimeStart->Checked;
    TimeH1UD       ->Enabled=TimeStart->Checked;
    BtnTime1       ->Enabled=TimeStart->Checked;
    TimeY2         ->Enabled=TimeEnd  ->Checked;
    TimeH2         ->Enabled=TimeEnd  ->Checked;
    TimeY2UD       ->Enabled=TimeEnd  ->Checked;
    TimeH2UD       ->Enabled=TimeEnd  ->Checked;
	BtnTime2       ->Enabled=TimeEnd  ->Checked;
    TimeInt        ->Enabled=TimeIntF ->Checked;
    LabelTimeInt   ->Enabled=TimeIntF ->Checked;
    TimeUnitF      ->Enabled=TimeStart->Checked&&TimeEnd  ->Checked;
    TimeUnit       ->Enabled=TimeUnitF->Enabled&&TimeUnitF->Checked;
    LabelTimeUnit  ->Enabled=TimeUnitF->Enabled&&TimeUnitF->Checked;
    OutDir         ->Enabled=OutDirEna->Checked;
    BtnOutDir      ->Enabled=OutDirEna->Checked;
    LabelOutDir    ->Enabled=OutDirEna->Checked;
}
// load options from ini file -----------------------------------------------
void __fastcall TMainForm::LoadOpt(void)
{
    TIniFile *ini=new TIniFile(IniFile);
	AnsiString s;
//	char *p;
	wchar_t *p;

    TimeStart->Checked =ini->ReadInteger("set","timestart",   0);
    TimeEnd->Checked   =ini->ReadInteger("set","timeend",     0);
    TimeY1->Text       =ini->ReadString ("set","timey1",      "2000/01/01");
    TimeY1->Text       =ini->ReadString ("set","timey1",      "2000/01/01");
    TimeH1->Text       =ini->ReadString ("set","timeh1",      "00:00:00");
    TimeY2->Text       =ini->ReadString ("set","timey2",      "2000/01/01");
    TimeH2->Text       =ini->ReadString ("set","timeh2",      "00:00:00");
    TimeIntF ->Checked =ini->ReadInteger("set","timeintf",    0);
    TimeInt->Text      =ini->ReadString ("set","timeint",     "0");
    TimeUnitF->Checked =ini->ReadInteger("set","timeunitf",   0);
    TimeUnit->Text     =ini->ReadString ("set","timeunit",    "24");
    InputFile1->Text   =ini->ReadString ("set","inputfile1",  "");
    InputFile2->Text   =ini->ReadString ("set","inputfile2",  "");
    InputFile3->Text   =ini->ReadString ("set","inputfile3",  "");
    InputFile4->Text   =ini->ReadString ("set","inputfile4",  "");
    InputFile5->Text   =ini->ReadString ("set","inputfile5",  "");
    OutDirEna->Checked =ini->ReadInteger("set","outputdirena", 0);
    OutDir->Text       =ini->ReadString ("set","outputdir",   "");
    OutputFile->Text   =ini->ReadString ("set","outputfile",  "");
    
    InputFile1->Items  =ReadList(ini,"hist","inputfile1");
    InputFile2->Items  =ReadList(ini,"hist","inputfile2");
    InputFile3->Items  =ReadList(ini,"hist","inputfile3");
    InputFile4->Items  =ReadList(ini,"hist","inputfile4");
    InputFile5->Items  =ReadList(ini,"hist","inputfile5");
	OutputFile->Items  =ReadList(ini,"hist","outputfile");
    
	PosMode            =ini->ReadInteger("opt","posmode",        0);
	Freqs              =ini->ReadInteger("opt","freqs",         1);
    L2Cod		       =ini->ReadInteger("opt", "l2cprior",      0);
    Solution           =ini->ReadInteger("opt","solution",       0);
	ElMask             =ini->ReadFloat  ("opt","elmask",      15.0);
	SnrMask.ena[0]     =ini->ReadInteger("opt","snrmask_ena1",   0);
	SnrMask.ena[1]     =ini->ReadInteger("opt","snrmask_ena2",   0);
    for (int i=0;i<3;i++) for (int j=0;j<9;j++) {
		SnrMask.mask[i][j]=
			ini->ReadFloat("opt",s.sprintf("snrmask_%d_%d",i+1,j+1),0.0);
	}
    IonoOpt            =ini->ReadInteger("opt","ionoopt",     IONOOPT_BRDC);
	TropOpt            =ini->ReadInteger("opt","tropopt",     TROPOPT_SAAS);
    TimSys			   =ini->ReadInteger("opt","tsyscorr",		 0);
    RcvBiasEst         =ini->ReadInteger("opt","rcvbiasest",     0);
    DynamicModel       =ini->ReadInteger("opt","dynamicmodel",   0);
    TideCorr           =ini->ReadInteger("opt","tidecorr",       0);
    SatEphem           =ini->ReadInteger("opt","satephem",       0);
	ExSats             =ini->ReadString ("opt","exsats",        "");
	NavSys             =ini->ReadInteger("opt","navsys",   SYS_GPS);
	GloCodePri1        =ini->ReadString ("opt","codepriglo1",   "");
	GloCodePri2        =ini->ReadString ("opt","codepriglo2",   "");

	PosOpt[0]          =ini->ReadInteger("opt","posopt1",        0);
	PosOpt[1]          =ini->ReadInteger("opt","posopt2",        0);
    PosOpt[2]          =ini->ReadInteger("opt","posopt3",        0);
    PosOpt[3]          =ini->ReadInteger("opt","posopt4",        0);
    PosOpt[4]          =ini->ReadInteger("opt","posopt5",        0);
    MapFunc            =ini->ReadInteger("opt","mapfunc",        0);

	AmbResMethod       =ini->ReadInteger("opt","ambresmethod",   0);
	AmbRes             =ini->ReadInteger("opt","ambres",         1);
	GloAmbRes          =ini->ReadInteger("opt","gloambres",      1);
	PppAmbRes          =ini->ReadInteger("opt","pppambres",      0);
	ValidThresAR       =ini->ReadFloat  ("opt","validthresar", 3.0);
	ThresAR2           =ini->ReadFloat  ("opt","thresar2",  0.9999);
	ThresAR3           =ini->ReadFloat  ("opt","thresar3",    0.25);
    LockCntFixAmb      =ini->ReadInteger("opt","lockcntfixamb",  0);
    FixCntHoldAmb      =ini->ReadInteger("opt","fixcntholdamb", 10);
    ElMaskAR           =ini->ReadFloat  ("opt","elmaskar",     0.0);
    ElMaskHold         =ini->ReadFloat  ("opt","elmaskhold",   0.0);
    OutCntResetAmb     =ini->ReadInteger("opt","outcntresetbias",5);
    SlipThres          =ini->ReadFloat  ("opt","slipthres",   0.05);
	PhaCyc		   	   =ini->ReadFloat  ("opt","phasshft",   	 0);
	L2CPBias		   =ini->ReadInteger("opt","l2cpbias",   	 0);
	MaxAgeDiff         =ini->ReadFloat  ("opt","maxagediff",  30.0);
	RejectThres        =ini->ReadFloat  ("opt","rejectthres", 30.0);
	RejectGdop         =ini->ReadFloat  ("opt","rejectgdop",  30.0);
	NumIter            =ini->ReadInteger("opt","numiter",        1);
	CodeSmooth         =ini->ReadInteger("opt","codesmooth",     0);
	BaseLine[0]        =ini->ReadFloat  ("opt","baselinelen",  0.0);
	BaseLine[1]        =ini->ReadFloat  ("opt","baselinesig",  0.0);
	BaseLineConst      =ini->ReadInteger("opt","baselineconst",  0);
	Isb	     	   	   =ini->ReadInteger("opt","isb",   	 0);
	Diff     	   	   =ini->ReadInteger("opt","diff",   	 0);

    SolFormat          =ini->ReadInteger("opt","solformat",      0);
    TimeFormat         =ini->ReadInteger("opt","timeformat",     1);
    TimeDecimal        =ini->ReadInteger("opt","timedecimal",    3);
    LatLonFormat       =ini->ReadInteger("opt","latlonformat",   0);
    FieldSep           =ini->ReadString ("opt","fieldsep",      "");
    OutputHead         =ini->ReadInteger("opt","outputhead",     1);
    OutputOpt          =ini->ReadInteger("opt","outputopt",      1);
    OutputDatum        =ini->ReadInteger("opt","outputdatum",    0);
    OutputHeight       =ini->ReadInteger("opt","outputheight",   0);
    OutputGeoid        =ini->ReadInteger("opt","outputgeoid",    0);
	SolStatic          =ini->ReadInteger("opt","solstatic",      0);
    DebugTrace         =ini->ReadInteger("opt","debugtrace",     0);
	DebugStatus        =ini->ReadInteger("opt","debugstatus",    0);
	IsbOut             =ini->ReadInteger("opt","isbout",         0);
	Gl2Out             =ini->ReadInteger("opt","gl2out",         0);
	FCBOut             =ini->ReadInteger("opt","fcbout",         0);

	PossnxOut          =ini->ReadInteger("opt", "possnxout",     0);
	IonOut             =ini->ReadInteger("opt", "ionout",        0);
	TropOut            =ini->ReadInteger("opt", "tropout",       0);
	RecClOut           =ini->ReadInteger("opt", "recclockout",   0);
	SatClOut           =ini->ReadInteger("opt", "satclockout",   0);
	StaticOut          =ini->ReadInteger("opt", "staticout",     0);

    PhaCycFile         =ini->ReadString ("opt","phacycfile",    "");
    GloIfbFile         =ini->ReadString ("opt","gloifbfile",    "");
    ErrModFile	       =ini->ReadString ("opt","errmodfile",	"");
     
    ErrMod			   =ini->ReadInteger("opt","errmodel",       0);
	MeasErrR1          =ini->ReadFloat  ("opt","measeratio1",100.0);
	MeasErrR2          =ini->ReadFloat  ("opt","measeratio2",100.0);
	MeasErrR3          =ini->ReadFloat  ("opt","measeratio3",100.0);
	MeasErr2           =ini->ReadFloat  ("opt","measerr2",   0.003);
	MeasErr3           =ini->ReadFloat  ("opt","measerr3",   0.003);
    MeasErr6		   =ini->ReadFloat  ("opt","dcbratio",    10.0);
    MeasErr4           =ini->ReadFloat  ("opt","measerr4",   0.000);
    MeasErr5           =ini->ReadFloat  ("opt","measerr5",  10.000);
    SatClkStab         =ini->ReadFloat  ("opt","satclkstab", 5E-12);
    PrNoise1           =ini->ReadFloat  ("opt","prnoise1",    1E-4);
    PrNoise2           =ini->ReadFloat  ("opt","prnoise2",    1E-3);
    PrNoise3           =ini->ReadFloat  ("opt","prnoise3",    1E-4);
	PrNoise4           =ini->ReadFloat  ("opt","prnoise4",    1E-1);
	PrNoise5           =ini->ReadFloat  ("opt","prnoise5",    1E-2);
	PrNoise6           =ini->ReadFloat  ("opt","prnoise6",     0.0);
	PrNoise7           =ini->ReadFloat  ("opt","prnoise7",     0.0);
	PrNoise8           =ini->ReadFloat  ("opt","prnoise8",     0.0);
    
    RovPosType         =ini->ReadInteger("opt","rovpostype",     0);
    RefPosType         =ini->ReadInteger("opt","refpostype",     0);
    RovPos[0]          =ini->ReadFloat  ("opt","rovpos1",      0.0);
    RovPos[1]          =ini->ReadFloat  ("opt","rovpos2",      0.0);
    RovPos[2]          =ini->ReadFloat  ("opt","rovpos3",      0.0);
    RefPos[0]          =ini->ReadFloat  ("opt","refpos1",      0.0);
    RefPos[1]          =ini->ReadFloat  ("opt","refpos2",      0.0);
    RefPos[2]          =ini->ReadFloat  ("opt","refpos3",      0.0);
    RovAntPcv          =ini->ReadInteger("opt","rovantpcv",      0);
    RefAntPcv          =ini->ReadInteger("opt","refantpcv",      0);
    RovAnt             =ini->ReadString ("opt","rovant",        "");
    RefAnt             =ini->ReadString ("opt","refant",        "");
    RovRecTyp	       =ini->ReadString ("opt","rovrectyp",		"");
    RefRecTyp    	   =ini->ReadString ("opt","refrectyp",		"");
    RovAntE            =ini->ReadFloat  ("opt","rovante",      0.0);
    RovAntN            =ini->ReadFloat  ("opt","rovantn",      0.0);
    RovAntU            =ini->ReadFloat  ("opt","rovantu",      0.0);
    RefAntE            =ini->ReadFloat  ("opt","refante",      0.0);
    RefAntN            =ini->ReadFloat  ("opt","refantn",      0.0);
    RefAntU            =ini->ReadFloat  ("opt","refantu",      0.0);
    
	RnxOpts1           =ini->ReadString ("opt","rnxopts1",      "");
	RnxOpts2           =ini->ReadString ("opt","rnxopts2",      "");
    
    AntPcvFile         =ini->ReadString ("opt","antpcvfile",    "");
    IntpRefObs         =ini->ReadInteger("opt","intprefobs",     0);
	SbasSat            =ini->ReadInteger("opt","sbassat",        0);
    NetRSCorr          =ini->ReadInteger("opt","netrscorr",      0);
    SatClkCorr         =ini->ReadInteger("opt","satclkcorr",     0);
    SbasCorr           =ini->ReadInteger("opt","sbascorr",       0);
    SbasCorr1          =ini->ReadInteger("opt","sbascorr1",      0);
    SbasCorr2          =ini->ReadInteger("opt","sbascorr2",      0);
    SbasCorr3          =ini->ReadInteger("opt","sbascorr3",      0);
    SbasCorr4          =ini->ReadInteger("opt","sbascorr4",      0);
    SbasCorrFile       =ini->ReadString ("opt","sbascorrfile",  "");
    PrecEphFile        =ini->ReadString ("opt","precephfile",   "");
    SatPcvFile         =ini->ReadString ("opt","satpcvfile",    "");
    StaPosFile         =ini->ReadString ("opt","staposfile",    "");
    GeoidDataFile      =ini->ReadString ("opt","geoiddatafile", "");
	IonoFile           =ini->ReadString ("opt","ionofile",      "");
	DCBFile            =ini->ReadString ("opt","dcbfile",       "");
	ISBFile            =ini->ReadString ("opt","isbfile",       "");
	ISBOutFile         =ini->ReadString ("opt","isboutfile",    "");
	GL2OutFile         =ini->ReadString ("opt","gl2outfile",    "");
	FCBOutFile         =ini->ReadString ("opt","fcboutfile",    "");
	EOPFile            =ini->ReadString ("opt","eopfile",       "");
	BLQFile            =ini->ReadString ("opt","blqfile",       "");
    GoogleEarthFile    =ini->ReadString ("opt","googleearthfile",GOOGLE_EARTH);
    BIPMCircularTFile  =ini->ReadString ("opt","bipmcirculartfile","");
    
	/*Estimate Satellite Clock/FCB*/
	EstSatClo		   =ini->ReadInteger("opt","estsatclo",    0.0);
	EstSatFCB	       =ini->ReadInteger("opt","estsatfcb",    0.0);
	SemiDCPara	       =ini->ReadString ("opt","semidcpara",    "");/*Semi-Dynamic Correction Parameter*/
	SolDirFile	       =ini->ReadString ("opt","soldirfile",    "");/*Solution Dir*/

	EstIntZ		       =ini->ReadFloat  ("opt","estintz",   7200.0);
	EstIntES           =ini->ReadFloat  ("opt","estintes", 43200.0);
	RanWalSigZen       =ini->ReadFloat  ("opt","ranwalsigzen",1E-4);
	RanWalSigEW        =ini->ReadFloat  ("opt","ranwalsigew", 1E-4);
	RanWalSigNS        =ini->ReadFloat  ("opt","ranwalsigns", 1E-7);
	ThrO		       =ini->ReadFloat  ("opt","thro",         5.0);
	ThrC	           =ini->ReadFloat  ("opt","thrc",         5.0);
	MaxBasDD	       =ini->ReadFloat  ("opt","maxbased",5000000.0);
	JudValWL	       =ini->ReadFloat  ("opt","judvalwl", 0.99990);
	JudValL1           =ini->ReadFloat  ("opt","judvall1", 0.99990);
	WeiDD		       =ini->ReadFloat  ("opt","weidd",       1E+6);
	ConCri		       =ini->ReadFloat  ("opt","concri",     0.001);
	MaxIte             =ini->ReadInteger("opt","maxite",         3);
	FCBFile            =ini->ReadString ("opt","fcbfile",       "");/*Satellite FCB*/


	TemStoFile	       =ini->ReadString ("opt","temstofile",    "");
//	NLFcb			   =ini->ReadFloat  ("opt","nlfcb",        900);
	NLFcb			   =ini->ReadInteger("opt","nlfcb",        900);
	Minsd			   =ini->ReadFloat  ("opt","minsd",        600);
	Mindd			   =ini->ReadFloat  ("opt","maxdd",       1200);
	Maxsd			   =ini->ReadFloat  ("opt","maxsd",        0.2);
	Fixwl			   =ini->ReadFloat  ("opt","fixwl",    0.99990);
	Fixnl			   =ini->ReadFloat  ("opt","fixnl",    0.99990);
	Mobstadn		   =ini->ReadFloat  ("opt","mobstadn",     100);
	Mobstade		   =ini->ReadFloat  ("opt","mobstade",     100);
	Mobstadu		   =ini->ReadFloat  ("opt","mobstadu",     100);
	Basestadn		   =ini->ReadFloat  ("opt","basestadn",    100);
	Basestade		   =ini->ReadFloat  ("opt","basestade",    100);
	Basestadu		   =ini->ReadFloat  ("opt","basestadu",    100);

	std[0]		   	   =ini->ReadFloat  ("opt","stdbias",     30.0);
	std[1]		   	   =ini->ReadFloat  ("opt","stdiono",     0.03);
	std[2]		   	   =ini->ReadFloat  ("opt","stdtrop",      0.3);
	std[3]		   	   =ini->ReadFloat  ("opt","stdgrad",     0.01);
	std[4]		   	   =ini->ReadFloat  ("opt","stdclkr",    100.0);
	std[5]		   	   =ini->ReadFloat  ("opt","stdclks",    100.0);
	std[6]		   	   =ini->ReadFloat  ("opt","stdissb",    100.0);

	minsat			   =ini->ReadInteger("opt","minsat",         5);

    RovList="";
	for (int i=0;i<10;i++) {
        RovList +=ini->ReadString("opt",s.sprintf("rovlist%d",i+1),"");
	}
    BaseList="";
    for (int i=0;i<10;i++) {
        BaseList+=ini->ReadString("opt",s.sprintf("baselist%d",i+1),"");
    }
//	for (p=RovList.c_str();*p;p++) {
//		if ((p=strstr(p,"@@"))) strncpy(p,"\r\n",2); else break;
//	}
//	for (p=BaseList.c_str();*p;p++) {
//		if ((p=strstr(p,"@@"))) strncpy(p,"\r\n",2); else break;
//	}
	for (p=RovList.c_str();*p;p++) {
		p=wcsstr(p,L"@@");
		if (p!=NULL) wcsncpy(p,L"\r\n",2); else break;
	}
	for (p=BaseList.c_str();*p;p++) {
		p=wcsstr(p,L"@@");
		if (p!=NULL) wcsncpy(p,L"\r\n",2); else break;
	}
    ExtErr.ena[0]      =ini->ReadInteger("opt","exterr_ena0",    0);
    ExtErr.ena[1]      =ini->ReadInteger("opt","exterr_ena1",    0);
    ExtErr.ena[2]      =ini->ReadInteger("opt","exterr_ena2",    0);
    ExtErr.ena[3]      =ini->ReadInteger("opt","exterr_ena3",    0);
    for (int i=0;i<3;i++) for (int j=0;j<6;j++) {
        ExtErr.cerr[i][j]=ini->ReadFloat("opt",s.sprintf("exterr_cerr%d%d",i,j),0.3);
    }
    for (int i=0;i<3;i++) for (int j=0;j<6;j++) {
        ExtErr.perr[i][j]=ini->ReadFloat("opt",s.sprintf("exterr_perr%d%d",i,j),0.003);
    }
    ExtErr.gloicb[0]   =ini->ReadFloat  ("opt","exterr_gloicb0",0.0);
    ExtErr.gloicb[1]   =ini->ReadFloat  ("opt","exterr_gloicb1",0.0);
    ExtErr.gpsglob[0]  =ini->ReadFloat  ("opt","exterr_gpsglob0",0.0);
    ExtErr.gpsglob[1]  =ini->ReadFloat  ("opt","exterr_gpsglob1",0.0);
    
    ConvDialog->TimeSpan  ->Checked  =ini->ReadInteger("conv","timespan",  0);
    ConvDialog->TimeIntF  ->Checked  =ini->ReadInteger("conv","timeintf",  0);
    ConvDialog->TimeY1    ->Text     =ini->ReadString ("conv","timey1","2000/01/01");
    ConvDialog->TimeH1    ->Text     =ini->ReadString ("conv","timeh1","00:00:00"  );
    ConvDialog->TimeY2    ->Text     =ini->ReadString ("conv","timey2","2000/01/01");
    ConvDialog->TimeH2    ->Text     =ini->ReadString ("conv","timeh2","00:00:00"  );
    ConvDialog->TimeInt   ->Text     =ini->ReadString ("conv","timeint", "0");
    ConvDialog->TrackColor->ItemIndex=ini->ReadInteger("conv","trackcolor",5);
    ConvDialog->PointColor->ItemIndex=ini->ReadInteger("conv","pointcolor",5);
    ConvDialog->OutputAlt ->ItemIndex=ini->ReadInteger("conv","outputalt", 0);
    ConvDialog->OutputTime->ItemIndex=ini->ReadInteger("conv","outputtime",0);
    ConvDialog->AddOffset ->Checked  =ini->ReadInteger("conv","addoffset", 0);
    ConvDialog->Offset1   ->Text     =ini->ReadString ("conv","offset1", "0");
    ConvDialog->Offset2   ->Text     =ini->ReadString ("conv","offset2", "0");
    ConvDialog->Offset3   ->Text     =ini->ReadString ("conv","offset3", "0");
    ConvDialog->Compress  ->Checked  =ini->ReadInteger("conv","compress",  0);
    
    TTextViewer::Color1=(TColor)ini->ReadInteger("viewer","color1",(int)clBlack);
    TTextViewer::Color2=(TColor)ini->ReadInteger("viewer","color2",(int)clWhite);
    TTextViewer::FontD=new TFont;
    TTextViewer::FontD->Name=ini->ReadString ("viewer","fontname","Courier New");
    TTextViewer::FontD->Size=ini->ReadInteger("viewer","fontsize",9);
    delete ini;
}
// save options to ini file -------------------------------------------------
void __fastcall TMainForm::SaveOpt(void)
{
    TIniFile *ini=new TIniFile(IniFile);
    AnsiString s;
//    char *p;
	wchar_t *p;

    ini->WriteInteger("set","timestart",   TimeStart ->Checked?1:0);
    ini->WriteInteger("set","timeend",     TimeEnd   ->Checked?1:0);
    ini->WriteString ("set","timey1",      TimeY1    ->Text);
    ini->WriteString ("set","timeh1",      TimeH1    ->Text);
    ini->WriteString ("set","timey2",      TimeY2    ->Text);
	ini->WriteString ("set","timeh2",      TimeH2    ->Text);
    ini->WriteInteger("set","timeintf",    TimeIntF  ->Checked?1:0);
    ini->WriteString ("set","timeint",     TimeInt   ->Text);
    ini->WriteInteger("set","timeunitf",   TimeUnitF ->Checked?1:0);
    ini->WriteString ("set","timeunit",    TimeUnit  ->Text);
    ini->WriteString ("set","inputfile1",  InputFile1->Text);
    ini->WriteString ("set","inputfile2",  InputFile2->Text);
    ini->WriteString ("set","inputfile3",  InputFile3->Text);
    ini->WriteString ("set","inputfile4",  InputFile4->Text);
    ini->WriteString ("set","inputfile5",  InputFile5->Text);
    ini->WriteInteger("set","outputdirena",OutDirEna ->Checked);
    ini->WriteString ("set","outputdir",   OutDir    ->Text);
	ini->WriteString ("set","outputfile",  OutputFile->Text);
    
    WriteList(ini,"hist","inputfile1",     InputFile1->Items);
    WriteList(ini,"hist","inputfile2",     InputFile2->Items);
	WriteList(ini,"hist","inputfile3",     InputFile3->Items);
    WriteList(ini,"hist","inputfile4",     InputFile4->Items);
    WriteList(ini,"hist","inputfile5",     InputFile5->Items);
    WriteList(ini,"hist","outputfile",     OutputFile->Items);
    
	ini->WriteInteger("opt","posmode",     PosMode     );
	ini->WriteInteger("opt","freqs",       Freqs       );
    ini->WriteInteger("opt","l2cprior",    L2Cod	   );
    ini->WriteInteger("opt","solution",    Solution    );
	ini->WriteFloat  ("opt","elmask",      ElMask      );
	ini->WriteInteger("opt","snrmask_ena1",SnrMask.ena[0]);
    ini->WriteInteger("opt","snrmask_ena2",SnrMask.ena[1]);
    for (int i=0;i<3;i++) for (int j=0;j<9;j++) {
		ini->WriteFloat("opt",s.sprintf("snrmask_%d_%d",i+1,j+1),
						SnrMask.mask[i][j]);
	}
	ini->WriteInteger("opt","ionoopt",     IonoOpt     );
    ini->WriteInteger("opt","tropopt",     TropOpt     );
    ini->WriteInteger("opt","tsyscorr",	   TimSys	   );
    ini->WriteInteger("opt","rcvbiasest",  RcvBiasEst  );
    ini->WriteInteger("opt","dynamicmodel",DynamicModel);
    ini->WriteInteger("opt","tidecorr",    TideCorr    );
    ini->WriteInteger("opt","satephem",    SatEphem    );
	ini->WriteString ("opt","exsats",      ExSats      );
    ini->WriteInteger("opt","navsys",      NavSys      );
	ini->WriteString ("opt","codepriglo1", GloCodePri1 );
	ini->WriteString ("opt","codepriglo2", GloCodePri2 );
    ini->WriteInteger("opt","posopt1",     PosOpt[0]   );
    ini->WriteInteger("opt","posopt2",     PosOpt[1]   );
    ini->WriteInteger("opt","posopt3",     PosOpt[2]   );
    ini->WriteInteger("opt","posopt4",     PosOpt[3]   );
    ini->WriteInteger("opt","posopt5",     PosOpt[4]   );
    ini->WriteInteger("opt","mapfunc",     MapFunc     );

	ini->WriteInteger("opt","ambresmethod",AmbResMethod);
	ini->WriteInteger("opt","ambres",      AmbRes      );
	ini->WriteInteger("opt","gloambres",   GloAmbRes   );
	ini->WriteInteger("opt","pppambres",   PppAmbRes   );
	ini->WriteFloat  ("opt","validthresar",ValidThresAR);
    ini->WriteFloat  ("opt","thresar2",    ThresAR2    );
    ini->WriteFloat  ("opt","thresar3",    ThresAR3    );
    ini->WriteInteger("opt","lockcntfixamb",LockCntFixAmb);
    ini->WriteInteger("opt","fixcntholdamb",FixCntHoldAmb);
    ini->WriteFloat  ("opt","elmaskar",    ElMaskAR    );
    ini->WriteFloat  ("opt","elmaskhold",  ElMaskHold  );
    ini->WriteInteger("opt","outcntresetbias",OutCntResetAmb);
	ini->WriteFloat  ("opt","slipthres",   SlipThres   );
	ini->WriteFloat  ("opt","phasshft",    PhaCyc	   );
	ini->WriteInteger("opt","l2cpbias",    L2CPBias	   );
    ini->WriteFloat  ("opt","maxagediff",  MaxAgeDiff  );
    ini->WriteFloat  ("opt","rejectgdop",  RejectGdop  );
    ini->WriteFloat  ("opt","rejectthres", RejectThres );
    ini->WriteInteger("opt","numiter",     NumIter     );
    ini->WriteInteger("opt","codesmooth",  CodeSmooth  );
    ini->WriteFloat  ("opt","baselinelen", BaseLine[0] );
    ini->WriteFloat  ("opt","baselinesig", BaseLine[1] );
	ini->WriteInteger("opt","baselineconst",BaseLineConst);
	ini->WriteInteger("opt","isb"          ,Isb    	   );
	ini->WriteInteger("opt","diff"         ,Diff   	   );
    
    ini->WriteInteger("opt","solformat",   SolFormat   );
    ini->WriteInteger("opt","timeformat",  TimeFormat  );
    ini->WriteInteger("opt","timedecimal", TimeDecimal );
    ini->WriteInteger("opt","latlonformat",LatLonFormat);
    ini->WriteString ("opt","fieldsep",    FieldSep    );
    ini->WriteInteger("opt","outputhead",  OutputHead  );
    ini->WriteInteger("opt","outputopt",   OutputOpt   );
    ini->WriteInteger("opt","outputdatum", OutputDatum );
    ini->WriteInteger("opt","outputheight",OutputHeight);
    ini->WriteInteger("opt","outputgeoid", OutputGeoid );
    ini->WriteInteger("opt","solstatic",   SolStatic   );
    ini->WriteInteger("opt","debugtrace",  DebugTrace  );
	ini->WriteInteger("opt","debugstatus", DebugStatus );
	ini->WriteInteger("opt","isbout",      IsbOut      );
	ini->WriteInteger("opt","gl2out",      Gl2Out      );
	ini->WriteInteger("opt","fcbout",      FCBOut      );

	ini->WriteInteger("opt","possnxout",   PossnxOut   );
	ini->WriteInteger("opt","ionout",      IonOut      );
	ini->WriteInteger("opt","tropout",     TropOut     );
	ini->WriteInteger("opt","recclockout", RecClOut    );
	ini->WriteInteger("opt","satclockout", SatClOut    );
	ini->WriteInteger("opt","staticout",  StaticOut   );
    
    ini->WriteInteger("opt","errmodel",    ErrMod	   );
    ini->WriteFloat  ("opt","measeratio1", MeasErrR1   );
    ini->WriteFloat  ("opt","measeratio2", MeasErrR2   );
	ini->WriteFloat  ("opt","measeratio3", MeasErrR3   );
	ini->WriteFloat  ("opt","measerr2",    MeasErr2    );
	ini->WriteFloat  ("opt","measerr3",    MeasErr3    );
    ini->WriteFloat  ("opt","dcbratio",    MeasErr6	   );
    ini->WriteFloat  ("opt","measerr4",    MeasErr4    );
    ini->WriteFloat  ("opt","measerr5",    MeasErr5    );
    ini->WriteFloat  ("opt","satclkstab",  SatClkStab  );
    ini->WriteFloat  ("opt","prnoise1",    PrNoise1    );
    ini->WriteFloat  ("opt","prnoise2",    PrNoise2    );
    ini->WriteFloat  ("opt","prnoise3",    PrNoise3    );
    ini->WriteFloat  ("opt","prnoise4",    PrNoise4    );
    ini->WriteFloat  ("opt","prnoise5",    PrNoise5    );
	ini->WriteFloat  ("opt","prnoise6",    PrNoise6    );
	ini->WriteFloat  ("opt","prnoise7",    PrNoise7    );
	ini->WriteFloat  ("opt","prnoise8",    PrNoise8    );

    ini->WriteString ("opt","phacycfile",  PhaCycFile  );
    ini->WriteString ("opt","gloifbfile",  GloIfbFile  );
    ini->WriteString ("opt","errmodfile",  ErrModFile  );
    ini->WriteInteger("opt","rovpostype",  RovPosType  );
    ini->WriteInteger("opt","refpostype",  RefPosType  );
    ini->WriteFloat  ("opt","rovpos1",     RovPos[0]   );
    ini->WriteFloat  ("opt","rovpos2",     RovPos[1]   );
    ini->WriteFloat  ("opt","rovpos3",     RovPos[2]   );
    ini->WriteFloat  ("opt","refpos1",     RefPos[0]   );
    ini->WriteFloat  ("opt","refpos2",     RefPos[1]   );
    ini->WriteFloat  ("opt","refpos3",     RefPos[2]   );
    ini->WriteInteger("opt","rovantpcv",   RovAntPcv   );
    ini->WriteInteger("opt","refantpcv",   RefAntPcv   );
    ini->WriteString ("opt","rovant",      RovAnt      );
    ini->WriteString ("opt","refant",      RefAnt      );
    ini->WriteString ("opt","rovrectyp",   RovRecTyp   );
    ini->WriteString ("opt","refrectyp",   RefRecTyp   );
    ini->WriteFloat  ("opt","rovante",     RovAntE     );
    ini->WriteFloat  ("opt","rovantn",     RovAntN     );
    ini->WriteFloat  ("opt","rovantu",     RovAntU     );
    ini->WriteFloat  ("opt","refante",     RefAntE     );
    ini->WriteFloat  ("opt","refantn",     RefAntN     );
    ini->WriteFloat  ("opt","refantu",     RefAntU     );
    
	ini->WriteString ("opt","rnxopts1",    RnxOpts1    );
	ini->WriteString ("opt","rnxopts2",    RnxOpts2    );
    
    ini->WriteString ("opt","antpcvfile",  AntPcvFile  );
    ini->WriteInteger("opt","intprefobs",  IntpRefObs  );
    ini->WriteInteger("opt","sbassat",     SbasSat     );
    ini->WriteInteger("opt","netrscorr",   NetRSCorr   );
    ini->WriteInteger("opt","satclkcorr",  SatClkCorr  );
    ini->WriteInteger("opt","sbascorr",    SbasCorr    );
    ini->WriteInteger("opt","sbascorr1",   SbasCorr1   );
    ini->WriteInteger("opt","sbascorr2",   SbasCorr2   );
    ini->WriteInteger("opt","sbascorr3",   SbasCorr3   );
    ini->WriteInteger("opt","sbascorr4",   SbasCorr4   );
    ini->WriteString ("opt","sbascorrfile",SbasCorrFile);
    ini->WriteString ("opt","precephfile", PrecEphFile );
    ini->WriteString ("opt","satpcvfile",  SatPcvFile  );
    ini->WriteString ("opt","staposfile",  StaPosFile  );
    ini->WriteString ("opt","geoiddatafile",GeoidDataFile);
	ini->WriteString ("opt","ionofile",    IonoFile     );
	ini->WriteString ("opt","dcbfile",     DCBFile     );
	ini->WriteString ("opt","isbfile",     ISBFile     );
	ini->WriteString ("opt","isboutfile",  ISBOutFile  );
	ini->WriteString ("opt","gl2outfile",  GL2OutFile  );
	ini->WriteString ("opt","fcboutfile",  FCBOutFile  );
	ini->WriteString ("opt","eopfile",     EOPFile     );
	ini->WriteString ("opt","blqfile",     BLQFile     );
    ini->WriteString ("opt","googleearthfile",GoogleEarthFile);
    ini->WriteString ("opt","bipmcirculartfile",BIPMCircularTFile);
    
	/*Estimate Satellite Clock/FCB*/
	ini->WriteInteger("opt","estsatclo",    EstSatClo);
	ini->WriteInteger("opt","estsatfcb",    EstSatFCB);
	ini->WriteString ("opt","semidcpara",   SemiDCPara);/*Semi-Dynamic Correction Parameter*/
	ini->WriteString ("opt","soldirfile",   SolDirFile);/*Solution Dir*/

	ini->WriteFloat  ("opt","estintz",      EstIntZ);
	ini->WriteFloat  ("opt","estintes",     EstIntES);
	ini->WriteFloat  ("opt","ranwalsigzen", RanWalSigZen);
	ini->WriteFloat  ("opt","ranwalsigew",  RanWalSigEW);
	ini->WriteFloat  ("opt","ranwalsigns",  RanWalSigNS);
	ini->WriteFloat  ("opt","thro",         ThrO);
	ini->WriteFloat  ("opt","thrc",         ThrC);
	ini->WriteFloat  ("opt","maxbased",     MaxBasDD);
	ini->WriteFloat  ("opt","judvalwl",     JudValWL);
	ini->WriteFloat  ("opt","judvall1",     JudValL1);
	ini->WriteFloat  ("opt","weidd",        WeiDD);
	ini->WriteFloat  ("opt","concri",       ConCri);
	ini->WriteInteger("opt","maxite",       MaxIte);
	ini->WriteString ("opt","fcbfile",      FCBFile);/*Satellite FCB*/

	ini->WriteString ("opt","temstofile",   TemStoFile);
	ini->WriteFloat  ("opt","nlfcb",        NLFcb);
	ini->WriteFloat  ("opt","minsd",        Minsd);
	ini->WriteFloat  ("opt","maxdd",        Mindd);
	ini->WriteFloat  ("opt","maxsd",        Maxsd);
	ini->WriteFloat  ("opt","fixwl",        Fixwl);
	ini->WriteFloat  ("opt","fixnl",        Fixnl);
	ini->WriteFloat  ("opt","mobstadn",     Mobstadn);
	ini->WriteFloat  ("opt","mobstade",     Mobstade);
	ini->WriteFloat  ("opt","mobstadu",     Mobstadu);
	ini->WriteFloat  ("opt","basestadn",    Basestadn);
	ini->WriteFloat  ("opt","basestade",    Basestade);
	ini->WriteFloat  ("opt","basestadu",    Basestadu);

	ini->WriteFloat  ("opt","stdbias",      std[0]);
	ini->WriteFloat  ("opt","stdiono",      std[1]);
	ini->WriteFloat  ("opt","stdtrop",      std[2]);
	ini->WriteFloat  ("opt","stdgrad",      std[3]);
	ini->WriteFloat  ("opt","stdclkr",      std[4]);
	ini->WriteFloat  ("opt","stdclks",      std[5]);
	ini->WriteFloat  ("opt","stdissb",      std[6]);
	ini->WriteInteger("opt","minsat",       minsat);

//    for (p=RovList.c_str();*p;p++) {
//        if ((p=strstr(p,"\r\n"))) strncpy(p,"@@",2); else break;
//    }
    for (p=RovList.c_str();*p;p++) {
        p=wcsstr(p,L"\r\n");
        if (p!=NULL) wcsncpy(p,L"@@",2); else break;
    }
    for (int i=0;i<10;i++) {
        ini->WriteString("opt",s.sprintf("rovlist%d",i+1),RovList.SubString(i*2000,2000));
	}
//	for (p=BaseList.c_str();*p;p++) {
//		if ((p=strstr(p,"\r\n"))) strncpy(p,"@@",2); else break;
//	}
    for (p=BaseList.c_str();*p;p++) {
        p=wcsstr(p,L"\r\n");
        if (p!=NULL) wcsncpy(p,L"@@",2); else break;
	}
	for (int i=0;i<10;i++) {
		ini->WriteString("opt",s.sprintf("baselist%d",i+1),BaseList.SubString(i*2000,2000));
	}
    ini->WriteInteger("opt","exterr_ena0", ExtErr.ena[0]);
    ini->WriteInteger("opt","exterr_ena1", ExtErr.ena[1]);
    ini->WriteInteger("opt","exterr_ena2", ExtErr.ena[2]);
    ini->WriteInteger("opt","exterr_ena3", ExtErr.ena[3]);
    
    for (int i=0;i<3;i++) for (int j=0;j<6;j++) {
        ini->WriteFloat("opt",s.sprintf("exterr_cerr%d%d",i,j),ExtErr.cerr[i][j]);
    }
    for (int i=0;i<3;i++) for (int j=0;j<6;j++) {
        ini->WriteFloat("opt",s.sprintf("exterr_perr%d%d",i,j),ExtErr.perr[i][j]);
    }
    ini->WriteFloat  ("opt","exterr_gloicb0",ExtErr.gloicb[0]);
    ini->WriteFloat  ("opt","exterr_gloicb1",ExtErr.gloicb[1]);
    ini->WriteFloat  ("opt","exterr_gpsglob0",ExtErr.gpsglob[0]);
    ini->WriteFloat  ("opt","exterr_gpsglob1",ExtErr.gpsglob[1]);
    
    ini->WriteInteger("conv","timespan",   ConvDialog->TimeSpan  ->Checked  );
    ini->WriteString ("conv","timey1",     ConvDialog->TimeY1    ->Text     );
    ini->WriteString ("conv","timeh1",     ConvDialog->TimeH1    ->Text     );
    ini->WriteString ("conv","timey2",     ConvDialog->TimeY2    ->Text     );
    ini->WriteString ("conv","timeh2",     ConvDialog->TimeH2    ->Text     );
    ini->WriteInteger("conv","timeintf",   ConvDialog->TimeIntF  ->Checked  );
    ini->WriteString ("conv","timeint",    ConvDialog->TimeInt   ->Text     );
    ini->WriteInteger("conv","trackcolor", ConvDialog->TrackColor->ItemIndex);
    ini->WriteInteger("conv","pointcolor", ConvDialog->PointColor->ItemIndex);
    ini->WriteInteger("conv","outputalt",  ConvDialog->OutputAlt ->ItemIndex);
    ini->WriteInteger("conv","outputtime", ConvDialog->OutputTime->ItemIndex);
    ini->WriteInteger("conv","addoffset",  ConvDialog->AddOffset ->Checked  );
    ini->WriteString ("conv","offset1",    ConvDialog->Offset1   ->Text     );
    ini->WriteString ("conv","offset2",    ConvDialog->Offset2   ->Text     );
    ini->WriteString ("conv","offset3",    ConvDialog->Offset3   ->Text     );
    ini->WriteInteger("conv","compress",   ConvDialog->Compress  ->Checked  );
    
    ini->WriteInteger("viewer","color1",(int)TTextViewer::Color1  );
    ini->WriteInteger("viewer","color2",(int)TTextViewer::Color2  );
    ini->WriteString ("viewer","fontname",TTextViewer::FontD->Name);
    ini->WriteInteger("viewer","fontsize",TTextViewer::FontD->Size);
    delete ini;
}

/* ------------------------------------------------
* convert char to wchar
* args   : char *str
*          int length
* return : wchar_t* 
* notes  : 
*--------------------------------------------------*/
wchar_t* __fastcall TMainForm::char2wchar(char *str, int length)
{
	int len;
	if(length<0){
		len = ::MultiByteToWideChar(CP_ACP,0,str,-1,NULL,0);
	}else{
		len = length;
	}
	wchar_t *buff = (wchar_t *)malloc(sizeof(wchar_t)*(len+1));
	wmemset(buff,0,len+1);
	::MultiByteToWideChar(CP_ACP,0,str,-1,buff,len+1);
	return buff;
}

/* ------------------------------------------------
* convert wchar to char
* args   : wchar_t *wstr
*          int length
* return : char* 
* notes  : 
*--------------------------------------------------*/
char* __fastcall TMainForm::wchar2char(wchar_t *wstr, int length)
{
	int len;
	if(length<0){
		len = ::WideCharToMultiByte(CP_ACP,0,wstr,-1,NULL,0,NULL,NULL);
	}else{
		len = length;
	}
	char *buff = (char *)malloc(sizeof(char)*(len+1));
	memset(buff,0,len+1);
	::WideCharToMultiByte(CP_ACP,0,wstr,-1,buff,len+1,NULL,NULL);
	return buff;
}


