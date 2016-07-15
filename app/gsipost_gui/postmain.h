/*------------------------------------------------------------------------------
* postmain.h
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
* references :
*
* history : 2015/01/08  1.0  new
*-----------------------------------------------------------------------------*/

#ifndef postmainH
#define postmainH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <Dialogs.hpp>
#include <ComCtrls.hpp>
#include <Buttons.hpp>
#include <FileCtrl.hpp>
#include <inifiles.hpp>

#include "rtklib.h"
//---------------------------------------------------------------------------
class TMainForm : public TForm
{
__published:
	TPanel *Panel1;
	TPanel *Panel2;
	TPanel *Panel3;
	TPanel *Panel4;
	TPanel *Panel5;
	TPanel *Message;
	
	TButton *BtnPlot;
	TButton *BtnView;
	TButton *BtnToKML;
	TButton *BtnOption;
	TButton *BtnExec;
	TButton *BtnExit;
	TButton *BtnInputFile1;
	TButton *BtnInputFile3;
	TButton *BtnInputFile2;
	TButton *BtnInputFile4;
	TButton *BtnInputFile5;
	TButton *BtnOutputFile;
	
	TSpeedButton *BtnTime1;
	TSpeedButton *BtnTime2;
	TSpeedButton *BtnInputPlot1;
	TSpeedButton *BtnInputView1;
	TSpeedButton *BtnInputView3;
	TSpeedButton *BtnInputPlot2;
	TSpeedButton *BtnInputView2;
	TSpeedButton *BtnInputView4;
	TSpeedButton *BtnInputView5;
	TSpeedButton *BtnAbout;
	TSpeedButton *BtnKeyword;
	
	TCheckBox *TimeStart;
	TCheckBox *TimeEnd;
	TCheckBox *TimeIntF;
	TCheckBox *TimeUnitF;
	TEdit *TimeY1;
	TEdit *TimeH1;
	TEdit *TimeY2;
	TEdit *TimeH2;
	TEdit *TimeUnit;
	TUpDown *TimeY1UD;
	TUpDown *TimeH1UD;
	TUpDown *TimeY2UD;
	TUpDown *TimeH2UD;
	TComboBox *TimeInt;
	TComboBox *InputFile1;
	TComboBox *InputFile3;
	TComboBox *InputFile2;
	TComboBox *InputFile4;
	TComboBox *InputFile5;
	TComboBox *OutputFile;
	
	TOpenDialog *OpenDialog;
	TSaveDialog *SaveDialog;
	
	TProgressBar *Progress;
	
	TLabel *Label1;
	TLabel *LabelInputFile1;
	TLabel *LabelInputFile2;
	TLabel *LabelInputFile3;
	TLabel *LabelTimeInt;
	TLabel *LabelTimeUnit;
	TEdit *OutDir;
	TCheckBox *OutDirEna;
	TButton *BtnOutDir;
	TSpeedButton *BtnOutputView2;
	TSpeedButton *BtnOutputView1;
	TLabel *LabelOutDir;

	void __fastcall FormCreate         (TObject *Sender);
	void __fastcall FormShow           (TObject *Sender);
	void __fastcall FormClose          (TObject *Sender, TCloseAction &Action);
	
	void __fastcall BtnPlotClick       (TObject *Sender);
	void __fastcall BtnViewClick       (TObject *Sender);
	void __fastcall BtnToKMLClick      (TObject *Sender);
	void __fastcall BtnOptionClick     (TObject *Sender);
	void __fastcall BtnExecClick       (TObject *Sender);
	void __fastcall BtnStopClick       (TObject *Sender);
	void __fastcall BtnExitClick       (TObject *Sender);
	void __fastcall BtnAboutClick      (TObject *Sender);
	
	void __fastcall BtnTime1Click      (TObject *Sender);
	void __fastcall BtnTime2Click      (TObject *Sender);
	void __fastcall BtnInputFile1Click (TObject *Sender);
	void __fastcall BtnInputFile3Click (TObject *Sender);
	void __fastcall BtnInputFile2Click (TObject *Sender);
	void __fastcall BtnInputFile4Click (TObject *Sender);
	void __fastcall BtnInputFile5Click (TObject *Sender);
	void __fastcall BtnOutputFileClick (TObject *Sender);
	void __fastcall BtnInputView1Click (TObject *Sender);
	void __fastcall BtnInputView3Click (TObject *Sender);
	void __fastcall BtnInputView2Click (TObject *Sender);
	void __fastcall BtnInputView4Click (TObject *Sender);
	void __fastcall BtnInputView5Click (TObject *Sender);
	void __fastcall BtnOutputView1Click(TObject *Sender);
	void __fastcall BtnOutputView2Click(TObject *Sender);
	void __fastcall BtnInputPlot1Click (TObject *Sender);
	void __fastcall BtnInputPlot2Click (TObject *Sender);
	void __fastcall BtnKeywordClick    (TObject *Sender);
	
	void __fastcall TimeStartClick     (TObject *Sender);
	void __fastcall TimeIntFClick      (TObject *Sender);
	void __fastcall TimeUnitFClick     (TObject *Sender);
	void __fastcall TimeH1UDChangingEx (TObject *Sender, bool &AllowChange,
          short NewValue, TUpDownDirection Direction);
	void __fastcall TimeY1UDChangingEx (TObject *Sender, bool &AllowChange,
          short NewValue, TUpDownDirection Direction);
	void __fastcall TimeY2UDChangingEx (TObject *Sender, bool &AllowChange,
          short NewValue, TUpDownDirection Direction);
	void __fastcall TimeH2UDChangingEx (TObject *Sender, bool &AllowChange,
          short NewValue, TUpDownDirection Direction);
	
	void __fastcall InputFile1Change   (TObject *Sender);
	void __fastcall OutDirEnaClick(TObject *Sender);
	void __fastcall BtnOutDirClick(TObject *Sender);
	void __fastcall OutDirChange(TObject *Sender);

private:
	void __fastcall DropFiles          (TWMDropFiles msg); // for files drop
	
	int  __fastcall ExecProc           (void);
		
	int  __fastcall GetOption(prcopt_t &prcopt, solopt_t &solopt, filopt_t &filopt);
	
	int  __fastcall ObsToNav (const wchar_t *obsfile, wchar_t *navfile);
	
	UnicodeString __fastcall FilePath(UnicodeString file);
	TStringList * __fastcall ReadList(TIniFile *ini, UnicodeString cat,
		UnicodeString key);
	void __fastcall WriteList(TIniFile *ini, UnicodeString cat,
		UnicodeString key, TStrings *list);
	void __fastcall AddHist  (TComboBox *combo);
	int __fastcall ExecCmd(UnicodeString cmd, int show);
	
	gtime_t __fastcall GetTime1(void);
	gtime_t __fastcall GetTime2(void);
	void __fastcall SetOutFile(void);
	void __fastcall SetTime1(gtime_t time);
	void __fastcall SetTime2(gtime_t time);
	void __fastcall UpdateEnable(void);
	void __fastcall LoadOpt(void);
	void __fastcall SaveOpt(void);
	
	BEGIN_MESSAGE_MAP
	MESSAGE_HANDLER(WM_DROPFILES,TWMDropFiles,DropFiles);
	END_MESSAGE_MAP(TForm);
public:
	UnicodeString IniFile;
	
	// options
	int PosMode,Freqs,Solution,DynamicModel,IonoOpt,TropOpt,RcvBiasEst;
	int PosOpt[6],MapFunc;
	int L2Cod,TimSys,PhaCyc,L2CPBias,ErrMod,Isb,Diff;
	int NumIter,CodeSmooth,TideCorr;
	int OutCntResetAmb,FixCntHoldAmb,LockCntFixAmb,RovPosType,RefPosType;
	int SatEphem,NavSys;
	int RovAntPcv,RefAntPcv,AmbResMethod,AmbRes,GloAmbRes,PppAmbRes,OutputHead,OutputOpt,OutputDatum;
	int OutputHeight,OutputGeoid,DebugTrace,DebugStatus,BaseLineConst,IsbOut,Gl2Out,FCBOut;
	int PossnxOut,IonOut,TropOut,RecClOut,SatClOut,StaticOut;
	int SolFormat,TimeFormat,LatLonFormat,IntpRefObs,NetRSCorr,SatClkCorr;
	int SbasCorr,SbasCorr1,SbasCorr2,SbasCorr3,SbasCorr4,TimeDecimal;
	int SolStatic,SbasSat;
	double ElMask,MaxAgeDiff,RejectThres,RejectGdop;
	double MeasErrR1,MeasErrR2,MeasErrR3,MeasErr2,MeasErr3,MeasErr4,MeasErr5;
	double SatClkStab,RovAntE,RovAntN,RovAntU,RefAntE,RefAntN,RefAntU;
	double PrNoise1,PrNoise2,PrNoise3,PrNoise4,PrNoise5,PrNoise6,PrNoise7,PrNoise8;
	double ValidThresAR,ElMaskAR,ElMaskHold,SlipThres;
	double ThresAR2,ThresAR3;
	double RovPos[3],RefPos[3],BaseLine[2];
	double MeasErr6;
	double EstIntZ,EstIntES,RanWalSigZen,RanWalSigEW,RanWalSigNS;
	double ThrO,ThrC,MaxBasDD,JudValWL,JudValL1,WeiDD,ConCri;
	int MaxIte;
	int EstSatClo,EstSatFCB;
    int NLFcb;
	double Minsd,Mindd,Maxsd,Fixwl,Fixnl;
	double Mobstadn,Mobstade,Mobstadu,Basestadn,Basestade,Basestadu;
	double std[7];
	snrmask_t SnrMask;

	UnicodeString FieldSep,RovAnt,RefAnt,AntPcvFile,StaPosFile,PrecEphFile;
	UnicodeString RnxOpts1,RnxOpts2;
	UnicodeString NetRSCorrFile1,NetRSCorrFile2,SatClkCorrFile,GoogleEarthFile,BIPMCircularTFile;
	UnicodeString GeoidDataFile,IonoFile,DCBFile,ISBFile,ISBOutFile,GL2OutFile,FCBOutFile;
	UnicodeString EOPFile,BLQFile;
	UnicodeString SbasCorrFile,SatPcvFile,ExSats;
	UnicodeString GloCodePri1,GloCodePri2;
	UnicodeString RovList,BaseList;
	UnicodeString RovRecTyp, RefRecTyp;
	UnicodeString PhaCycFile,GloIfbFile,ErrModFile;
	UnicodeString SemiDCPara,SolDirFile,FCBFile;
	UnicodeString TemStoFile;

	int minsat;

	void __fastcall ViewFile(UnicodeString file);
	void __fastcall ShowMsg(wchar_t *msg);
	__fastcall TMainForm(TComponent* Owner);
	wchar_t* __fastcall char2wchar(char *str, int length);
	char* __fastcall wchar2char(wchar_t *wstr, int length);
public:
	exterr_t ExtErr;
};
//---------------------------------------------------------------------------
extern PACKAGE TMainForm *MainForm;
//---------------------------------------------------------------------------
#endif
