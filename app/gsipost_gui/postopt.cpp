/*------------------------------------------------------------------------------
* postopt.cpp
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

#include <vcl.h>
#pragma hdrstop

#include "postmain.h"
#include "postopt.h"
#include "keydlg.h"
#include "viewer.h"
#include "refdlg.h"

#include "rtklib.h"
#include "maskoptdlg.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TOptDialog *OptDialog;
//---------------------------------------------------------------------------
static int str2int(UnicodeString str)
{
	int val=0.0;
	swscanf(str.c_str(),L"%d",&val);
	return val;
}
//---------------------------------------------------------------------------
static double str2dbl(UnicodeString str)
{
	double val=0.0;
	swscanf(str.c_str(),L"%lf",&val);
	return val;
}
//---------------------------------------------------------------------------
__fastcall TOptDialog::TOptDialog(TComponent* Owner)
    : TForm(Owner)
{
    UnicodeString label,s;
    int nglo=MAXPRNGLO,ngal=MAXPRNGAL,nqzs=MAXPRNQZS,ncmp=MAXPRNCMP;
    
    setlocale( LC_ALL, "" );
	_wsetlocale(LC_ALL, L"" );	
	if (nglo<=0) NavSys2->Enabled=false;
    if (ngal<=0) NavSys3->Enabled=false;
    if (nqzs<=0) NavSys4->Enabled=false;
    if (ncmp<=0) NavSys6->Enabled=false;
#ifdef EXTLEX
    IonoOpt ->Items->Add("QZSS LEX");
    SatEphem->Items->Add("QZSS LEX");
#endif
    UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::FormShow(TObject *Sender)
{
	GetOpt();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnOkClick(TObject *Sender)
{
	SetOpt();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnLoadClick(TObject *Sender)
{
	OpenDialog->Title="Load Options";
	OpenDialog->FilterIndex=4;
	if (!OpenDialog->Execute()) return;
	LoadOpt(OpenDialog->FileName);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnSaveClick(TObject *Sender)
{
	SaveDialog->Title="Save Options";
	SaveDialog->FilterIndex=2;
	if (!SaveDialog->Execute()) return;
	SaveOpt(SaveDialog->FileName);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnStaPosViewClick(TObject *Sender)
{
	if (StaPosFile->Text=="") return;
	TTextViewer *viewer=new TTextViewer(Application);
	viewer->Show();
	viewer->Read(StaPosFile->Text);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnStaPosFileClick(TObject *Sender)
{
	OpenDialog->Title="Station Postion File";
	OpenDialog->FilterIndex=3;
	if (!OpenDialog->Execute()) return;
	StaPosFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::RovPosTypeChange(TObject *Sender)
{
	TEdit *edit[]={RovPos1,RovPos2,RovPos3};
	double pos[3];
	GetPos(RovPosTypeP,edit,pos);
	SetPos(RovPosType->ItemIndex,edit,pos);
	RovPosTypeP=RovPosType->ItemIndex;
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::RefPosTypeChange(TObject *Sender)
{
	TEdit *edit[]={RefPos1,RefPos2,RefPos3};
	double pos[3];
	GetPos(RefPosTypeP,edit,pos);
	SetPos(RefPosType->ItemIndex,edit,pos);
	RefPosTypeP=RefPosType->ItemIndex;
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnRovPosClick(TObject *Sender)
{
	TEdit *edit[]={RovPos1,RovPos2,RovPos3};
	double p[3],pos[3];
	GetPos(RovPosType->ItemIndex,edit,p);
	ecef2pos(p,pos);
	RefDialog->RovPos[0]=pos[0]*R2D;
	RefDialog->RovPos[1]=pos[1]*R2D;
	RefDialog->Pos[2]=pos[2];
	RefDialog->StaPosFile=StaPosFile->Text;
	RefDialog->Left=Left+Width/2-RefDialog->Width/2;
	RefDialog->Top=Top+Height/2-RefDialog->Height/2;
	if (RefDialog->ShowModal()!=mrOk) return;
	pos[0]=RefDialog->Pos[0]*D2R;
	pos[1]=RefDialog->Pos[1]*D2R;
	pos[2]=RefDialog->Pos[2];
	pos2ecef(pos,p);
	SetPos(RovPosType->ItemIndex,edit,p);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnRefPosClick(TObject *Sender)
{
	TEdit *edit[]={RefPos1,RefPos2,RefPos3};
	double p[3],pos[3];
	GetPos(RefPosType->ItemIndex,edit,p);
	ecef2pos(p,pos);
	RefDialog->RovPos[0]=pos[0]*R2D;
	RefDialog->RovPos[1]=pos[1]*R2D;
	RefDialog->RovPos[2]=pos[2];
	RefDialog->StaPosFile=StaPosFile->Text;
	RefDialog->Left=Left+Width/2-RefDialog->Width/2;
	RefDialog->Top=Top+Height/2-RefDialog->Height/2;
	if (RefDialog->ShowModal()!=mrOk) return;
	pos[0]=RefDialog->Pos[0]*D2R;
	pos[1]=RefDialog->Pos[1]*D2R;
	pos[2]=RefDialog->Pos[2];
	pos2ecef(pos,p);
	SetPos(RefPosType->ItemIndex,edit,p);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnSatPcvViewClick(TObject *Sender)
{
	if (SatPcvFile->Text=="") return;
	TTextViewer *viewer=new TTextViewer(Application);
	viewer->Show();
	viewer->Read(SatPcvFile->Text);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnSatPcvFileClick(TObject *Sender)
{
	OpenDialog->Title="Satellite Antenna PCV File";
	OpenDialog->FilterIndex=2;
	if (!OpenDialog->Execute()) return;
	SatPcvFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnAntPcvViewClick(TObject *Sender)
{
	if (AntPcvFile->Text=="") return;
	TTextViewer *viewer=new TTextViewer(Application);
	viewer->Show();
	viewer->Read(AntPcvFile->Text);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnAntPcvFileClick(TObject *Sender)
{
	OpenDialog->Title="Receiver Antenna PCV File";
	OpenDialog->FilterIndex=2;
	if (!OpenDialog->Execute()) return;
	AntPcvFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnGeoidDataFileClick(TObject *Sender)
{
	OpenDialog->Title="Geoid Data File";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	GeoidDataFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnIonoFileClick(TObject *Sender)
{
	OpenDialog->Title="Iono Data File";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	IonoFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnIonoViewClick(TObject *Sender)
{
	if (IonoFile->Text=="") return;
	TTextViewer *viewer=new TTextViewer(Application);
	viewer->Show();
	viewer->Read(IonoFile->Text);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnDCBFileClick(TObject *Sender)
{
	OpenDialog->Title="DCB Data File";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	DCBFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnDCBViewClick(TObject *Sender)
{
	if (DCBFile->Text=="") return;
	TTextViewer *viewer=new TTextViewer(Application);
	viewer->Show();
	viewer->Read(DCBFile->Text);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnISBFileClick(TObject *Sender)
{
	OpenDialog->Title="ISB Data File";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	ISBFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnISBViewClick(TObject *Sender)
{
	if (ISBFile->Text=="") return;
	TTextViewer *viewer=new TTextViewer(Application);
	viewer->Show();
	viewer->Read(ISBFile->Text);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnISBOutFileClick(TObject *Sender)
{
	OpenDialog->Title="ISB Data File";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	ISBOutFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnISBOutViewClick(TObject *Sender)
{
	if (ISBOutFile->Text=="") return;
	TTextViewer *viewer=new TTextViewer(Application);
	viewer->Show();
	viewer->Read(ISBOutFile->Text);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnGL2OutFileClick(TObject *Sender)
{
	OpenDialog->Title="L2P-L2C Data File";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	GL2OutFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnGL2OutViewClick(TObject *Sender)
{
	if (GL2OutFile->Text=="") return;
	TTextViewer *viewer=new TTextViewer(Application);
	viewer->Show();
	viewer->Read(GL2OutFile->Text);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnFCBOutFileClick(TObject *Sender)
{
	OpenDialog->Title="FCB Data File";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	FCBOutFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnFCBOutViewClick(TObject *Sender)
{
	if (FCBOutFile->Text=="") return;
	TTextViewer *viewer=new TTextViewer(Application);
	viewer->Show();
	viewer->Read(FCBOutFile->Text);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnEOPFileClick(TObject *Sender)
{
	OpenDialog->Title="EOP Data File";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	EOPFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnEOPViewClick(TObject *Sender)
{
	if (EOPFile->Text=="") return;
	TTextViewer *viewer=new TTextViewer(Application);
	viewer->Show();
	viewer->Read(EOPFile->Text);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnBLQFileClick(TObject *Sender)
{
	OpenDialog->Title="BLQ Data File";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	BLQFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnBLQViewClick(TObject *Sender)
{
	if (BLQFile->Text=="") return;
	TTextViewer *viewer=new TTextViewer(Application);
	viewer->Show();
	viewer->Read(BLQFile->Text);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnGoogleEarthFileClick(TObject *Sender)
{
	OpenDialog->Title="Google Earth Exe File";
	OpenDialog->FilterIndex=5;
	if (!OpenDialog->Execute()) return;
	GoogleEarthFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::FreqChange(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::IonoOptChange(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::TropOptChange(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::DynamicModelChange(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::SatEphemChange(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::SolFormatChange(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::IsbOutChange(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::PosModeChange(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::SatEphemClick(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::NavSys2Click(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::AmbResChange(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::RovAntPcvClick(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::NetRSCorrClick(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::SatClkCorrClick(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::RovPosClick(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::RefPosClick(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::SbasCorrClick(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::OutputHeightClick(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BaselineConstClick(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::RovAntClick(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------

void __fastcall TOptDialog::RefAntClick(TObject *Sender)
{
	UpdateEnable();
}

//---------------------------------------------------------------------------

void __fastcall TOptDialog::L2CPBiasChange(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------

void __fastcall TOptDialog::IsbChange(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------

void __fastcall TOptDialog::Gl2OutChange(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------

void __fastcall TOptDialog::FCBOutChange(TObject *Sender)
{
	UpdateEnable();
}
//---------------------------------------------------------------------------

void __fastcall TOptDialog::GetOpt(void)
{
	TEdit *editu[]={RovPos1,RovPos2,RovPos3};
	TEdit *editr[]={RefPos1,RefPos2,RefPos3};
	UnicodeString s;
	PosMode		 ->ItemIndex	=MainForm->PosMode;
	if     (MainForm->Freqs==  1) Freqs->ItemIndex=0;
	else if(MainForm->Freqs==  5) Freqs->ItemIndex=1;
	else if(MainForm->Freqs== 12) Freqs->ItemIndex=2;
	else if(MainForm->Freqs== 15) Freqs->ItemIndex=3;
	else if(MainForm->Freqs==125) Freqs->ItemIndex=4;
	Solution	 ->ItemIndex	=MainForm->Solution;
	ElMask		 ->Text			=s.sprintf(L"%.0f",MainForm->ElMask);
	SnrMask						=MainForm->SnrMask;
	DynamicModel ->ItemIndex	=MainForm->DynamicModel;
	TideCorr	 ->ItemIndex	=MainForm->TideCorr;
	IonoOpt		 ->ItemIndex	=MainForm->IonoOpt;
	TropOpt		 ->ItemIndex	=MainForm->TropOpt;
	SatEphem	 ->ItemIndex	=MainForm->SatEphem;
	ExSats	     ->Text			=MainForm->ExSats;
	NavSys1	     ->Checked		=MainForm->NavSys&SYS_GPS;
	NavSys2	     ->Checked		=MainForm->NavSys&SYS_GLO;
	NavSys3	     ->Checked		=MainForm->NavSys&SYS_GAL;
	NavSys4	     ->Checked		=MainForm->NavSys&SYS_QZS;
	NavSys5	     ->Checked		=MainForm->NavSys&SYS_SBS;
	NavSys6	     ->Checked		=MainForm->NavSys&SYS_CMP;
	GloCodePri1->Text           =MainForm->GloCodePri1;
	GloCodePri2->Text           =MainForm->GloCodePri2;
	PosOpt1	     ->Checked		=MainForm->PosOpt[0];
	PosOpt2	     ->Checked		=MainForm->PosOpt[1];
	PosOpt3	     ->Checked		=MainForm->PosOpt[2];
	PosOpt4	     ->Checked		=MainForm->PosOpt[3];
	PosOpt5	     ->Checked		=MainForm->PosOpt[4];

	if(MainForm->AmbRes==0) {
		AmbResMethod->ItemIndex = 0;
		AmbRes->ItemIndex = 1;
	}
	else {
		AmbResMethod->ItemIndex = MainForm->AmbResMethod;
		AmbRes->ItemIndex = MainForm->AmbRes-1;
	}
	GloAmbRes	 ->ItemIndex	=MainForm->GloAmbRes;
	PppAmbRes	 ->ItemIndex	=MainForm->PppAmbRes;
	ValidThresAR ->Text			=s.sprintf(L"%.3g",MainForm->ValidThresAR);
	ThresAR2     ->Text			=s.sprintf(L"%.9g",MainForm->ThresAR2);
	ThresAR3     ->Text			=s.sprintf(L"%.3g",MainForm->ThresAR3);
	OutCntResetAmb->Text		=s.sprintf(L"%d",MainForm->OutCntResetAmb);
	FixCntHoldAmb->Text			=s.sprintf(L"%d",MainForm->FixCntHoldAmb);
	LockCntFixAmb->Text			=s.sprintf(L"%d",MainForm->LockCntFixAmb);
	ElMaskAR	 ->Text			=s.sprintf(L"%.0f",MainForm->ElMaskAR);
	ElMaskHold	 ->Text			=s.sprintf(L"%.0f",MainForm->ElMaskHold);
	MaxAgeDiff	 ->Text			=s.sprintf(L"%.1f",MainForm->MaxAgeDiff);
	RejectGdop   ->Text			=s.sprintf(L"%.1f",MainForm->RejectGdop);
	RejectThres  ->Text			=s.sprintf(L"%.1f",MainForm->RejectThres);
	SlipThres	 ->Text			=s.sprintf(L"%.3f",MainForm->SlipThres);
	NumIter		 ->Text			=s.sprintf(L"%d",  MainForm->NumIter);
	BaselineLen	 ->Text			=s.sprintf(L"%.3f",MainForm->BaseLine[0]);
	BaselineSig	 ->Text			=s.sprintf(L"%.3f",MainForm->BaseLine[1]);
	BaselineConst->Checked		=MainForm->BaseLineConst;
	
	SolFormat	 ->ItemIndex	=MainForm->SolFormat;
	TimeFormat	 ->ItemIndex	=MainForm->TimeFormat;
	TimeDecimal	 ->Text			=s.sprintf(L"%d",MainForm->TimeDecimal);
	LatLonFormat ->ItemIndex	=MainForm->LatLonFormat;
	FieldSep	 ->Text			=MainForm->FieldSep;
	OutputHead	 ->ItemIndex	=MainForm->OutputHead;
	OutputOpt	 ->ItemIndex	=MainForm->OutputOpt;
	OutputDatum  ->ItemIndex	=MainForm->OutputDatum;
	OutputHeight ->ItemIndex	=MainForm->OutputHeight;
	OutputGeoid  ->ItemIndex	=MainForm->OutputGeoid;
	SolStatic    ->ItemIndex	=MainForm->SolStatic;
	DebugTrace	 ->ItemIndex	=MainForm->DebugTrace;
	DebugStatus	 ->ItemIndex	=MainForm->DebugStatus;
	IsbOut	     ->ItemIndex	=MainForm->IsbOut;
	Gl2Out	     ->ItemIndex	=MainForm->Gl2Out;
	FCBOut	     ->ItemIndex	=MainForm->FCBOut;

	OutPosSinex	 ->ItemIndex	=MainForm->PossnxOut;
	OutIon	     ->ItemIndex	=MainForm->IonOut;
	OutTrop	     ->ItemIndex	=MainForm->TropOut;
	OutReceiverClock->ItemIndex	=MainForm->RecClOut;
	OutSatelliteClock->ItemIndex=MainForm->SatClOut;
	OutStatic->ItemIndex	    =MainForm->StaticOut;

	MeasErrR1	 ->Text			=s.sprintf(L"%.1f",MainForm->MeasErrR1);
	MeasErrR2	 ->Text			=s.sprintf(L"%.1f",MainForm->MeasErrR2);
	MeasErrR3	 ->Text			=s.sprintf(L"%.1f",MainForm->MeasErrR3);
	MeasErr2	 ->Text			=s.sprintf(L"%.3f",MainForm->MeasErr2);
	MeasErr3	 ->Text			=s.sprintf(L"%.3f",MainForm->MeasErr3);
	MeasErr4	 ->Text			=s.sprintf(L"%.3f",MainForm->MeasErr4);
	MeasErr5	 ->Text			=s.sprintf(L"%.3f",MainForm->MeasErr5);
	SatClkStab	 ->Text			=s.sprintf(L"%.2E",MainForm->SatClkStab);
	PrNoise1	 ->Text			=s.sprintf(L"%.2E",MainForm->PrNoise1);
	PrNoise2	 ->Text			=s.sprintf(L"%.2E",MainForm->PrNoise2);
	PrNoise3	 ->Text			=s.sprintf(L"%.2E",MainForm->PrNoise3);
	PrNoise4	 ->Text			=s.sprintf(L"%.2E",MainForm->PrNoise4);
	PrNoise5	 ->Text			=s.sprintf(L"%.2E",MainForm->PrNoise5);
	PrNoise6	 ->Text			=s.sprintf(L"%.2E",MainForm->PrNoise6);
	PrNoise7	 ->Text			=s.sprintf(L"%.2E",MainForm->PrNoise7);
	PrNoise8	 ->Text			=s.sprintf(L"%.2E",MainForm->PrNoise8);
	
	RovAntPcv	 ->Checked		=MainForm->RovAntPcv;
	RefAntPcv	 ->Checked		=MainForm->RefAntPcv;
	RovAnt		 ->Text			=MainForm->RovAnt;
	RefAnt		 ->Text			=MainForm->RefAnt;
	RovAntE		 ->Text			=s.sprintf(L"%.4f",MainForm->RovAntE);
	RovAntN		 ->Text			=s.sprintf(L"%.4f",MainForm->RovAntN);
	RovAntU		 ->Text			=s.sprintf(L"%.4f",MainForm->RovAntU);
	RefAntE		 ->Text			=s.sprintf(L"%.4f",MainForm->RefAntE);
	RefAntN		 ->Text			=s.sprintf(L"%.4f",MainForm->RefAntN);
	RefAntU		 ->Text			=s.sprintf(L"%.4f",MainForm->RefAntU);
	AntPcvFile	 ->Text			=MainForm->AntPcvFile;

	RnxOpts1	 ->Text			=MainForm->RnxOpts1;
	RnxOpts2	 ->Text			=MainForm->RnxOpts2;
	
	IntpRefObs	 ->ItemIndex	=MainForm->IntpRefObs;
	SbasSat		 ->Text			=s.sprintf(L"%d",MainForm->SbasSat);
	SatPcvFile   ->Text			=MainForm->SatPcvFile;
	StaPosFile	 ->Text			=MainForm->StaPosFile;
	GeoidDataFile->Text			=MainForm->GeoidDataFile;
	IonoFile	 ->Text			=MainForm->IonoFile;
	DCBFile		 ->Text			=MainForm->DCBFile;
	ISBFile		 ->Text			=MainForm->ISBFile;
	ISBOutFile	->Text			=MainForm->ISBOutFile;
	GL2OutFile	->Text			=MainForm->GL2OutFile;
	FCBOutFile	->Text			=MainForm->FCBOutFile;
	EOPFile		 ->Text			=MainForm->EOPFile;
	BLQFile		 ->Text			=MainForm->BLQFile;
	GoogleEarthFile->Text		=MainForm->GoogleEarthFile;
	RovPosType	 ->ItemIndex	=MainForm->RovPosType;
	RefPosType	 ->ItemIndex	=MainForm->RefPosType;
	RovPosTypeP					=RovPosType->ItemIndex;
	RefPosTypeP					=RefPosType->ItemIndex;
	SetPos(RovPosType->ItemIndex,editu,MainForm->RovPos);
	SetPos(RefPosType->ItemIndex,editr,MainForm->RefPos);
	ReadAntList();
	
	RovList		 ->Text			=MainForm->RovList;
	BaseList	 ->Text			=MainForm->BaseList;

	L2Cod        ->ItemIndex=MainForm->L2Cod;/* L2 code priority (0:L2P,1:L2C) */
	TimSys       ->ItemIndex=MainForm->TimSys;/* time system correction ( 0:off,1:on) */
	PhaCyc       ->ItemIndex=MainForm->PhaCyc;/* phase cycle shift (0:off,1:table) */
	L2CPBias     ->ItemIndex=MainForm->L2CPBias;
	ErrMod       ->ItemIndex=MainForm->ErrMod;/* error model (0:user settings,1:table) */
	MeasErr6	 ->Text     =s.sprintf(L"%.1f",MainForm->MeasErr6);/* code error ratio(no DCB) */
	Isb          ->ItemIndex=MainForm->Isb;/* ISB (0:off,1:table,2:est,3:est(code),4:est(phase)) */
	Diff         ->ItemIndex=MainForm->Diff;

	/* receiver types {rover,base} */
	RovRecTyp		->Text=MainForm->RovRecTyp;
	RefRecTyp		->Text=MainForm->RefRecTyp;

	PhaCycFile   ->Text			=MainForm->PhaCycFile;
	GloIfbFile	 ->Text			=MainForm->GloIfbFile;
	ErrModFile   ->Text			=MainForm->ErrModFile;
    
    BIPMCircularTFile->Text		=MainForm->BIPMCircularTFile;/*BIPMCircularTFile*/

	/*Estimate Satellite Clock/FCB*/
	EstSatClo	 ->ItemIndex	=MainForm->EstSatClo;
	EstSatFCB	 ->ItemIndex	=MainForm->EstSatFCB;

	SemiDCPara	 ->Text			=MainForm->SemiDCPara;/*Semi-Dynamic Correction Parameter*/
	SolDirFile	 ->Text			=MainForm->SolDirFile;/*Solution Dir*/
	FCBFile		 ->Text			=MainForm->FCBFile;/*Satellite FCB*/

	if(MainForm->EstIntZ>0){
		EstIntZ		->Text     =s.sprintf(L"%.0f",MainForm->EstIntZ);
	}else{
		EstIntZ		->Text     =s.sprintf(L"7200",MainForm->EstIntZ);
	}
	if(MainForm->EstIntES>0){
		EstIntES	->Text     =s.sprintf(L"%.0f",MainForm->EstIntES);
	}else{
		EstIntES	->Text     =s.sprintf(L"43200",MainForm->EstIntES);
	}
	if(MainForm->RanWalSigZen>0){
		RanWalSigZen->Text     =s.sprintf(L"%.2E",MainForm->RanWalSigZen);
	}else{
		RanWalSigZen->Text     =s.sprintf(L"%.2E",1.0E-4);
	}
	if(MainForm->RanWalSigEW>0){
		RanWalSigEW	->Text     =s.sprintf(L"%.2E",MainForm->RanWalSigEW);
	}else{
		RanWalSigEW	->Text     =s.sprintf(L"%.2E",1.0E-4);
	}
	if(MainForm->RanWalSigNS>0){
		RanWalSigNS	->Text     =s.sprintf(L"%.2E",MainForm->RanWalSigNS);
	}else{
		RanWalSigNS	->Text     =s.sprintf(L"%.2E",1.0E-7);
	}
	if(MainForm->ThrO>0.0){
		ThrO		->Text     =s.sprintf(L"%.1f",MainForm->ThrO);
	}else{

		ThrO		->Text     =s.sprintf(L"%.1f",5.0);
	}

	if(MainForm->ThrC>0.0){
		ThrC		->Text     =s.sprintf(L"%.1f",MainForm->ThrC);
	}else{

		ThrC		->Text     =s.sprintf(L"%.1f",5.0);
	}

	if(MainForm->JudValWL<=1.0){
	JudValWL	->Text     =s.sprintf(L"%.5f",MainForm->JudValWL);
	}else{
		JudValWL	->Text     =s.sprintf(L"%.5f",0.99990);
	}
	if(MainForm->JudValL1<=1.0){
	JudValL1	->Text     =s.sprintf(L"%.5f",MainForm->JudValL1);
	}else{
		JudValL1	->Text     =s.sprintf(L"%.5f",0.99990);
	}
	if(MainForm->ConCri>0){
		ConCri		->Text     =s.sprintf(L"%.2E",MainForm->ConCri);
	}else{
		ConCri		->Text     =s.sprintf(L"%.2E",0.001);
    }
	if(MainForm->MaxIte>0){
		MaxIte		->Text     =s.sprintf(L"%d",MainForm->MaxIte);
	}else{
		MaxIte		->Text     =s.sprintf(L"%d",1);
	}

	TemStoFile	 ->Text			=MainForm->TemStoFile;/*Temporary storage epoch parameters*/

//	NLFcb		 ->Text     		=s.sprintf(L"%.0f",MainForm->NLFcb);
	NLFcb		 ->Text     		=s.sprintf(L"%d",MainForm->NLFcb);
	Minsd		 ->Text     		=s.sprintf(L"%.0f",MainForm->Minsd);
	Mindd		 ->Text     		=s.sprintf(L"%.0f",MainForm->Mindd);

	if(MainForm->Fixwl<=1.0){
	Fixwl		 ->Text     		=s.sprintf(L"%.5f",MainForm->Fixwl);
	}else{
		Fixwl	->Text     =s.sprintf(L"%.5f",0.99990);
	}
	if(MainForm->Fixnl<=1.0){
	Fixnl		 ->Text     		=s.sprintf(L"%.5f",MainForm->Fixnl);
	}else{
		Fixnl	->Text     =s.sprintf(L"%.5f",0.99990);
	}
	Mobstadn		 ->Text     	=s.sprintf(L"%.4f",MainForm->Mobstadn);
	Mobstade		 ->Text     	=s.sprintf(L"%.4f",MainForm->Mobstade);
	Mobstadu		 ->Text     	=s.sprintf(L"%.4f",MainForm->Mobstadu);
	Basestadn		 ->Text     	=s.sprintf(L"%.4f",MainForm->Basestadn);
	Basestade		 ->Text     	=s.sprintf(L"%.4f",MainForm->Basestade);
	Basestadu		 ->Text     	=s.sprintf(L"%.4f",MainForm->Basestadu);

	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::SetOpt(void)
{
	TEdit *editu[]={RovPos1,RovPos2,RovPos3};
	TEdit *editr[]={RefPos1,RefPos2,RefPos3};
	
	MainForm->PosMode		=PosMode	->ItemIndex;
	if     (Freqs->ItemIndex==0) MainForm->Freqs=1;
	else if(Freqs->ItemIndex==1) MainForm->Freqs=5;
	else if(Freqs->ItemIndex==2) MainForm->Freqs=12;
	else if(Freqs->ItemIndex==3) MainForm->Freqs=15;
	else if(Freqs->ItemIndex==4) MainForm->Freqs=125;
	MainForm->Solution		=Solution   ->ItemIndex;
	MainForm->ElMask		=str2dbl(ElMask	->Text);
//	MainForm->SnrMask		=str2dbl(SnrMask->Text);
	MainForm->SnrMask		=SnrMask;
	MainForm->DynamicModel	=DynamicModel->ItemIndex;
	MainForm->TideCorr		=TideCorr	->ItemIndex;
	MainForm->IonoOpt	  	=IonoOpt	->ItemIndex;
	MainForm->TropOpt	  	=TropOpt	->ItemIndex;
	MainForm->SatEphem	  	=SatEphem	->ItemIndex;
	MainForm->ExSats	  	=ExSats		->Text;
	MainForm->NavSys	  	=0;
	if (NavSys1->Checked) MainForm->NavSys|=SYS_GPS;
	if (NavSys2->Checked) MainForm->NavSys|=SYS_GLO;
	if (NavSys3->Checked) MainForm->NavSys|=SYS_GAL;
	if (NavSys4->Checked) MainForm->NavSys|=SYS_QZS;
	if (NavSys5->Checked) MainForm->NavSys|=SYS_SBS;
	if (NavSys6->Checked) MainForm->NavSys|=SYS_CMP;
	MainForm->GloCodePri1   =GloCodePri1->Text;
	MainForm->GloCodePri2   =GloCodePri2->Text;
	MainForm->PosOpt[0]	  	=PosOpt1	->Checked;
	MainForm->PosOpt[1]	  	=PosOpt2	->Checked;
	MainForm->PosOpt[2]	  	=PosOpt3	->Checked;
	MainForm->PosOpt[3]	  	=PosOpt4	->Checked;
	MainForm->PosOpt[4]	  	=PosOpt5	->Checked;


	MainForm->AmbResMethod	=AmbResMethod->ItemIndex;
	MainForm->AmbRes	  	=AmbRes		->ItemIndex+1;
	MainForm->GloAmbRes	  	=GloAmbRes	->ItemIndex;
	MainForm->PppAmbRes	  	=PppAmbRes	->ItemIndex;
	MainForm->ValidThresAR	=str2dbl(ValidThresAR->Text);
	MainForm->ThresAR2		=str2dbl(ThresAR2->Text);
	MainForm->ThresAR3		=str2dbl(ThresAR3->Text);
	MainForm->OutCntResetAmb=OutCntResetAmb->Text.ToInt();
	MainForm->FixCntHoldAmb =FixCntHoldAmb->Text.ToInt();
	MainForm->OutCntResetAmb=OutCntResetAmb->Text.ToInt();
	MainForm->LockCntFixAmb	=LockCntFixAmb->Text.ToInt();
	MainForm->ElMaskAR	  	=ElMaskAR   ->Text.ToInt();
	MainForm->ElMaskHold  	=ElMaskHold ->Text.ToInt();
	MainForm->MaxAgeDiff  	=str2dbl(MaxAgeDiff ->Text);
	MainForm->RejectGdop 	=str2dbl(RejectGdop ->Text);
	MainForm->RejectThres 	=str2dbl(RejectThres->Text);
	MainForm->SlipThres   	=str2dbl(SlipThres  ->Text);
	MainForm->NumIter	  	=NumIter	  ->Text.ToInt();
	MainForm->BaseLine[0]  	=str2dbl(BaselineLen->Text);
	MainForm->BaseLine[1]  	=str2dbl(BaselineSig->Text);
	MainForm->BaseLineConst	=BaselineConst->Checked;
	
	MainForm->SolFormat   	=SolFormat  ->ItemIndex;
	MainForm->TimeFormat  	=TimeFormat ->ItemIndex;
	MainForm->TimeDecimal  	=str2dbl(TimeDecimal->Text);
	MainForm->LatLonFormat	=LatLonFormat->ItemIndex;
	MainForm->FieldSep	  	=FieldSep   ->Text;
	MainForm->OutputHead  	=OutputHead ->ItemIndex;
	MainForm->OutputOpt   	=OutputOpt  ->ItemIndex;
	MainForm->OutputDatum 	=OutputDatum->ItemIndex;
	MainForm->OutputHeight	=OutputHeight->ItemIndex;
	MainForm->OutputGeoid 	=OutputGeoid->ItemIndex;
	MainForm->SolStatic	 	=SolStatic  ->ItemIndex;
	MainForm->DebugTrace  	=DebugTrace ->ItemIndex;
	MainForm->DebugStatus  	=DebugStatus->ItemIndex;
	MainForm->IsbOut       	=IsbOut->ItemIndex;
	MainForm->Gl2Out       	=Gl2Out->ItemIndex;
	MainForm->FCBOut       	=FCBOut->ItemIndex;

	MainForm->PossnxOut     =OutPosSinex->ItemIndex;
	MainForm->IonOut       	=OutIon->ItemIndex;
	MainForm->TropOut       =OutTrop->ItemIndex;
	MainForm->RecClOut      =OutReceiverClock->ItemIndex;
	MainForm->SatClOut      =OutSatelliteClock->ItemIndex;
	MainForm->StaticOut     =OutStatic->ItemIndex;

	MainForm->MeasErrR1	  =str2dbl(MeasErrR1  ->Text);
	MainForm->MeasErrR2	  =str2dbl(MeasErrR2  ->Text);
	MainForm->MeasErrR3	  =str2dbl(MeasErrR3  ->Text);
	MainForm->MeasErr2	  =str2dbl(MeasErr2   ->Text);
	MainForm->MeasErr3	  =str2dbl(MeasErr3   ->Text);
	MainForm->MeasErr4	  =str2dbl(MeasErr4   ->Text);
	MainForm->MeasErr5	  =str2dbl(MeasErr5   ->Text);
	MainForm->SatClkStab  =str2dbl(SatClkStab ->Text);
	MainForm->PrNoise1	  =str2dbl(PrNoise1   ->Text);
	MainForm->PrNoise2	  =str2dbl(PrNoise2   ->Text);
	MainForm->PrNoise3	  =str2dbl(PrNoise3   ->Text);
	MainForm->PrNoise4	  =str2dbl(PrNoise4   ->Text);
	MainForm->PrNoise5	  =str2dbl(PrNoise5   ->Text);
	MainForm->PrNoise6	  =str2dbl(PrNoise6   ->Text);
	MainForm->PrNoise7	  =str2dbl(PrNoise7   ->Text);
	MainForm->PrNoise8	  =str2dbl(PrNoise8   ->Text);
	
	MainForm->RovAntPcv   =RovAntPcv	->Checked;
	MainForm->RefAntPcv   =RefAntPcv	->Checked;
	MainForm->RovAnt	  =RovAnt		->Text;
	MainForm->RefAnt	  =RefAnt		->Text;
	MainForm->RovAntE	  =str2dbl(RovAntE	->Text);
	MainForm->RovAntN	  =str2dbl(RovAntN	->Text);
	MainForm->RovAntU	  =str2dbl(RovAntU	->Text);
	MainForm->RefAntE	  =str2dbl(RefAntE	->Text);
	MainForm->RefAntN	  =str2dbl(RefAntN	->Text);
	MainForm->RefAntU	  =str2dbl(RefAntU	->Text);

	MainForm->RnxOpts1	  =RnxOpts1		->Text;
	MainForm->RnxOpts2	  =RnxOpts2		->Text;
	
	MainForm->IntpRefObs  =IntpRefObs	->ItemIndex;
	MainForm->SbasSat     =SbasSat		->Text.ToInt();
	MainForm->AntPcvFile  =AntPcvFile	->Text;
	MainForm->SatPcvFile  =SatPcvFile	->Text;
	MainForm->StaPosFile  =StaPosFile	->Text;
	MainForm->GeoidDataFile=GeoidDataFile->Text;
	MainForm->IonoFile    =IonoFile		->Text;
	MainForm->DCBFile     =DCBFile		->Text;
	MainForm->ISBFile     =ISBFile		->Text;
	MainForm->ISBOutFile  =ISBOutFile	->Text;
	MainForm->GL2OutFile  =GL2OutFile	->Text;
	MainForm->FCBOutFile  =FCBOutFile	->Text;
	MainForm->BLQFile     =EOPFile		->Text;
	MainForm->EOPFile     =BLQFile		->Text;
	MainForm->GoogleEarthFile=GoogleEarthFile->Text;
	MainForm->RovPosType  =RovPosType	->ItemIndex;
	MainForm->RefPosType  =RefPosType	->ItemIndex;
	GetPos(RovPosType->ItemIndex,editu,MainForm->RovPos);
	GetPos(RefPosType->ItemIndex,editr,MainForm->RefPos);
	
	MainForm->RovList	  =RovList		->Text;
	MainForm->BaseList	  =BaseList		->Text;
	
	MainForm->L2Cod		=	L2Cod        ->ItemIndex;/* L2 code priority (0:L2P,1:L2C) */
	MainForm->TimSys	=	TimSys       ->ItemIndex;/* time system correction ( 0:off,1:on) */
	MainForm->PhaCyc	=	PhaCyc       ->ItemIndex;/* phase cycle shift (0:off,1:table) */
    MainForm->L2CPBias	=	L2CPBias     ->ItemIndex;
	MainForm->ErrMod	=	ErrMod       ->ItemIndex;/* error model (0:user settings,1:table) */
	MainForm->MeasErr6  =str2dbl(MeasErr6  ->Text);/* code error ratio(no DCB) */
	MainForm->Isb	    =	Isb          ->ItemIndex;/* Isb (0:off,1:table,2:est,3:est(code),4:est(phase)) */
	MainForm->Diff	    =	Diff          ->ItemIndex;

	/* receiver types {rover,base} */
	MainForm->RovRecTyp	  =RovRecTyp		->Text;
	MainForm->RefRecTyp	  =RefRecTyp		->Text;

	MainForm->PhaCycFile  =PhaCycFile		->Text;/* 1/4cycle phase correction file */
	MainForm->GloIfbFile  =GloIfbFile		->Text;/* GLONASS IFB table file */
	MainForm->ErrModFile  =ErrModFile		->Text;/* error model file */

    MainForm->BIPMCircularTFile=BIPMCircularTFile->Text;/*BIPMCircularTFile*/

	/*Estimate Satellite Clock/FCB*/
	MainForm->EstSatClo	  =EstSatClo->ItemIndex;
	MainForm->EstSatFCB	  =EstSatFCB->ItemIndex;

	MainForm->SemiDCPara	 =SemiDCPara->Text;/*Semi-Dynamic Correction Parameter*/
	MainForm->SolDirFile	 =SolDirFile->Text;/*Solution Dir*/
	MainForm->FCBFile		 =FCBFile->Text;/*Satellite FCB*/

	MainForm->EstIntZ      = str2dbl(EstIntZ		->Text);
	MainForm->EstIntES     = str2dbl(EstIntES		->Text);
	MainForm->RanWalSigZen = str2dbl(RanWalSigZen	->Text);
	MainForm->RanWalSigEW  = str2dbl(RanWalSigEW	->Text);
	MainForm->RanWalSigNS  = str2dbl(RanWalSigNS	->Text);
	MainForm->ThrO         = str2dbl(ThrO			->Text);
	MainForm->ThrC         = str2dbl(ThrC			->Text);
	MainForm->JudValWL     = str2dbl(JudValWL		->Text);
	MainForm->JudValL1     = str2dbl(JudValL1		->Text);
	MainForm->ConCri       = str2dbl(ConCri			->Text);
	MainForm->MaxIte       = str2dbl(MaxIte			->Text);

	MainForm->TemStoFile   = TemStoFile	 ->Text;/*Temporary storage epoch parameters*/

	MainForm->NLFcb        = str2dbl(NLFcb		->Text);
	MainForm->Minsd        = str2dbl(Minsd		->Text);
	MainForm->Mindd        = str2dbl(Mindd		->Text);
	MainForm->Fixwl        = str2dbl(Fixwl		->Text);
	MainForm->Fixnl        = str2dbl(Fixnl		->Text);
	MainForm->Mobstadn     = str2dbl(Mobstadn	->Text);
	MainForm->Mobstade     = str2dbl(Mobstade	->Text);
	MainForm->Mobstadu     = str2dbl(Mobstadu	->Text);
	MainForm->Basestadn    = str2dbl(Basestadn	->Text);
	MainForm->Basestade    = str2dbl(Basestade	->Text);
	MainForm->Basestadu    = str2dbl(Basestadu	->Text);

	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::LoadOpt(UnicodeString file)
{
	TEdit *editu[]={RovPos1,RovPos2,RovPos3};
	TEdit *editr[]={RefPos1,RefPos2,RefPos3};
	UnicodeString s;
	char buff[1024]="",*p,id[32];
	int sat;
	prcopt_t prcopt=prcopt_default;
	solopt_t solopt=solopt_default;
	filopt_t filopt={""};
	prcopt.mopt=mbsopt_default;


	resetsysopts();
	char *ch = wchar2char(file.c_str(),-1);
	if (!loadopts(ch,sysopts)){
		if(ch!=NULL){
			free(ch);
			ch=NULL;
		}
		return;
	}
	if(ch!=NULL){
		free(ch);
		ch=NULL;
	}

	getsysopts(&prcopt,&solopt,&filopt);

	PosMode		 ->ItemIndex	=prcopt.mode;
	if(prcopt.oprfrq[0]==0) {
		Freqs->ItemIndex = 0;
		if(prcopt.oprfrq[1]==1) {
			Freqs->ItemIndex = 2;
			if(prcopt.oprfrq[2]==2) {
				Freqs->ItemIndex = 4;
			}
		}
		else if(prcopt.oprfrq[1]==2) {
			Freqs->ItemIndex = 3;
		}
	}
	else if(prcopt.oprfrq[0]==2) {
		Freqs->ItemIndex = 1;
	}
	Solution	 ->ItemIndex	=prcopt.soltype;
	ElMask		 ->Text			=s.sprintf(L"%.0f",prcopt.elmin*R2D);
	SnrMask						=prcopt.snrmask;
	DynamicModel ->ItemIndex	=prcopt.dynamics;
	TideCorr	 ->ItemIndex	=prcopt.tidecorr;
	IonoOpt		 ->ItemIndex	=prcopt.ionoopt;
	TropOpt		 ->ItemIndex	=prcopt.tropopt;
	SatEphem	 ->ItemIndex	=prcopt.sateph;
	ExSats	     ->Text			="";
	for (sat=1,p=buff;sat<=MAXSAT;sat++) {
		if (!prcopt.exsats[sat-1]) continue;
		satno2id(sat,id);
		p+=sprintf(p,"%s%s%s",p==buff?"":" ",prcopt.exsats[sat-1]==2?"+":"",id);
	}
	wchar_t *wbuff = char2wchar(buff,1024);
	ExSats		 ->Text			=wbuff;
	if(wbuff!=NULL){
		free(wbuff);
		wbuff=NULL;
	}
	NavSys1	     ->Checked		=prcopt.navsys&SYS_GPS;
	NavSys2	     ->Checked		=prcopt.navsys&SYS_GLO;
	NavSys3	     ->Checked		=prcopt.navsys&SYS_GAL;
	NavSys4	     ->Checked		=prcopt.navsys&SYS_QZS;
	NavSys5	     ->Checked		=prcopt.navsys&SYS_SBS;
	NavSys6	     ->Checked		=prcopt.navsys&SYS_CMP;

	GloCodePri1->Text           =prcopt.codepri[ISYSGLO][0];
	GloCodePri2->Text           =prcopt.codepri[ISYSGLO][1];

	PosOpt1	     ->Checked		=prcopt.posopt[0];
	PosOpt2	     ->Checked		=prcopt.posopt[1];
	PosOpt3	     ->Checked		=prcopt.posopt[2];
	PosOpt4	     ->Checked		=prcopt.posopt[3];
	PosOpt5	     ->Checked		=prcopt.posopt[4];

	if(prcopt.modear==0) {
		AmbResMethod->ItemIndex	=0;
		AmbRes	 ->ItemIndex	=1;
	}
	else {
		AmbResMethod->ItemIndex	=prcopt.modertkar;
		AmbRes	 ->ItemIndex	=prcopt.modear-1;
	}
	GloAmbRes	 ->ItemIndex	=prcopt.glomodear;
	PppAmbRes	 ->ItemIndex	=prcopt.modepppar;
	ValidThresAR ->Text			=s.sprintf(L"%.3g",prcopt.thresar[0]);
	ThresAR2	 ->Text			=s.sprintf(L"%.9g",prcopt.thresar[1]);
	ThresAR3	 ->Text			=s.sprintf(L"%.3g",prcopt.thresar[2]);
	OutCntResetAmb->Text		=s.sprintf(L"%d"  ,prcopt.maxout   );
	FixCntHoldAmb->Text			=s.sprintf(L"%d"  ,prcopt.minfix   );
	LockCntFixAmb  ->Text		=s.sprintf(L"%d"  ,prcopt.minlock  );
	ElMaskAR	 ->Text			=s.sprintf(L"%.0f",prcopt.elmaskar*R2D);
	ElMaskHold	 ->Text			=s.sprintf(L"%.0f",prcopt.elmaskhold*R2D);
	MaxAgeDiff	 ->Text			=s.sprintf(L"%.1f",prcopt.maxtdiff );
	RejectGdop   ->Text			=s.sprintf(L"%.1f",prcopt.maxgdop  );
	RejectThres  ->Text			=s.sprintf(L"%.1f",prcopt.maxinno  );
	SlipThres	 ->Text			=s.sprintf(L"%.3f",prcopt.thresslip);
	NumIter		 ->Text			=s.sprintf(L"%d",  prcopt.niter    );
	BaselineLen	 ->Text			=s.sprintf(L"%.3f",prcopt.baseline[0]);
	BaselineSig	 ->Text			=s.sprintf(L"%.3f",prcopt.baseline[1]);
	BaselineConst->Checked		=prcopt.baseline[0]>0.0;

	SolFormat	 ->ItemIndex	=solopt.posf;
	TimeFormat	 ->ItemIndex	=solopt.timef==0?0:solopt.times+1;
	TimeDecimal	 ->Text			=s.sprintf(L"%d",solopt.timeu);
	LatLonFormat ->ItemIndex	=solopt.degf;
	FieldSep	 ->Text			=solopt.sep;
	OutputHead	 ->ItemIndex	=solopt.outhead;
	OutputOpt	 ->ItemIndex	=solopt.outopt;
	OutputDatum  ->ItemIndex	=solopt.datum;
	OutputHeight ->ItemIndex	=solopt.height;
	OutputGeoid  ->ItemIndex	=solopt.geoid;
	SolStatic    ->ItemIndex	=solopt.solstatic;
	NmeaIntv1	 ->Text			=s.sprintf(L"%.2g",solopt.nmeaintv[0]);
	NmeaIntv2	 ->Text			=s.sprintf(L"%.2g",solopt.nmeaintv[1]);
	DebugTrace	 ->ItemIndex	=solopt.trace;
	DebugStatus	 ->ItemIndex	=solopt.sstat;
	IsbOut	     ->ItemIndex	=solopt.isbout;
	Gl2Out	     ->ItemIndex	=solopt.gl2out;
	FCBOut	     ->ItemIndex	=solopt.fcbout;

	OutPosSinex  ->ItemIndex    =solopt.possnxout;
	OutIon       ->ItemIndex    =solopt.ionout;
	OutTrop      ->ItemIndex    =solopt.tropout;
	OutReceiverClock->ItemIndex =solopt.recclout;
	OutSatelliteClock->ItemIndex=solopt.satclout;
	OutStatic->ItemIndex        =solopt.sbresout;
	
	MeasErrR1	 ->Text			=s.sprintf(L"%.1f",prcopt.eratio[0]);
	MeasErrR2	 ->Text			=s.sprintf(L"%.1f",prcopt.eratio[1]);
	MeasErrR3	 ->Text			=s.sprintf(L"%.1f",prcopt.eratio[2]);
	MeasErr2	 ->Text			=s.sprintf(L"%.3f",prcopt.err[1]);
	MeasErr3	 ->Text			=s.sprintf(L"%.3f",prcopt.err[2]);
	MeasErr4	 ->Text			=s.sprintf(L"%.3f",prcopt.err[3]);
	MeasErr5	 ->Text			=s.sprintf(L"%.3f",prcopt.err[4]);
	SatClkStab	 ->Text			=s.sprintf(L"%.2E",prcopt.sclkstab);
	PrNoise1	 ->Text			=s.sprintf(L"%.2E",prcopt.prn[0]);
	PrNoise2	 ->Text			=s.sprintf(L"%.2E",prcopt.prn[1]);
	PrNoise3	 ->Text			=s.sprintf(L"%.2E",prcopt.prn[2]);
	PrNoise4	 ->Text			=s.sprintf(L"%.2E",prcopt.prn[3]);
	PrNoise5	 ->Text			=s.sprintf(L"%.2E",prcopt.prn[4]);
	PrNoise6	 ->Text			=s.sprintf(L"%.2E",prcopt.prn[5]);
	PrNoise7	 ->Text			=s.sprintf(L"%.2E",prcopt.prn[6]);
	PrNoise8	 ->Text			=s.sprintf(L"%.2E",prcopt.prn[7]);
	
	RovAntPcv	 ->Checked		=*prcopt.anttype[0];
	RefAntPcv	 ->Checked		=*prcopt.anttype[1];
	RovAnt		 ->Text			=prcopt.anttype[0];
	RefAnt		 ->Text			=prcopt.anttype[1];
	RovAntE		 ->Text			=s.sprintf(L"%.4f",prcopt.antdel[0][0]);
	RovAntN		 ->Text			=s.sprintf(L"%.4f",prcopt.antdel[0][1]);
	RovAntU		 ->Text			=s.sprintf(L"%.4f",prcopt.antdel[0][2]);
	RefAntE		 ->Text			=s.sprintf(L"%.4f",prcopt.antdel[1][0]);
	RefAntN		 ->Text			=s.sprintf(L"%.4f",prcopt.antdel[1][1]);
	RefAntU		 ->Text			=s.sprintf(L"%.4f",prcopt.antdel[1][2]);

	RnxOpts1	 ->Text			=prcopt.rnxopt[0];
	RnxOpts2	 ->Text			=prcopt.rnxopt[1];
	
	IntpRefObs	 ->ItemIndex	=prcopt.intpref;
	SbasSat		 ->Text			=s.sprintf(L"%d",prcopt.sbassatsel);
	RovPosType	 ->ItemIndex	=prcopt.rovpos==0?0:prcopt.rovpos+2;
	RefPosType	 ->ItemIndex	=prcopt.refpos==0?0:prcopt.refpos+2;
	RovPosTypeP					=RovPosType->ItemIndex;
	RefPosTypeP					=RefPosType->ItemIndex;
	SetPos(RovPosType->ItemIndex,editu,prcopt.ru);
	SetPos(RefPosType->ItemIndex,editr,prcopt.rb);
	
	SatPcvFile ->Text			=filopt.satantp;
	AntPcvFile ->Text			=filopt.rcvantp;
	StaPosFile ->Text			=filopt.stapos;
	GeoidDataFile->Text			=filopt.geoid;
	IonoFile   ->Text			=filopt.iono;
	DCBFile	   ->Text			=filopt.dcb;
	ISBFile	   ->Text			=filopt.isb;
	ISBOutFile ->Text			=solopt.isbfile;
	GL2OutFile ->Text			=solopt.gl2file;
	FCBOutFile ->Text			=solopt.fcb;
	EOPFile	   ->Text			=filopt.eop;
	BLQFile	   ->Text			=filopt.blq;
	GoogleEarthFile->Text		=filopt.geexe;
	
	L2Cod        ->ItemIndex=prcopt.l2cprior;/* L2 code priority (0:L2P,1:L2C) */
	TimSys       ->ItemIndex=prcopt.tsyscorr;/* time system correction ( 0:off,1:on) */
	PhaCyc       ->ItemIndex=prcopt.phasshft;/* phase cycle shift (0:off,1:table) */
    L2CPBias     ->ItemIndex=prcopt.gl2bias;
	ErrMod       ->ItemIndex=prcopt.errmodel;/* error model (0:user settings,1:table) */
	MeasErr6	 ->Text     =s.sprintf(L"%.1f",prcopt.dcbratio);/* code error ratio(no DCB) */
	Isb          ->ItemIndex=prcopt.isb;/* ISB (0:off,1:table,2:est,3:est(code),4:est(phase)) */
	Diff         ->ItemIndex=prcopt.diff;

	RovRecTyp  ->Text		=prcopt.rectype[0];/* receiver types {rover,base} */
	RefRecTyp  ->Text		=prcopt.rectype[1];/* receiver types {rover,base} */

	PhaCycFile	 ->Text=prcopt.mopt.ifpcs;/* 1/4cycle phase correction file */
	GloIfbFile	 ->Text=prcopt.mopt.ififb;/* GLONASS IFB table file */
	ErrModFile	 ->Text=prcopt.mopt.iferr;/* error model file */

	BIPMCircularTFile->Text		=filopt.cirtfile;/*BIPMCircularTFile*/

	/*Estimate Satellite Clock/FCB*/
	EstSatClo	 ->ItemIndex	=prcopt.mopt.estsatclk;
	EstSatFCB	 ->ItemIndex	=prcopt.mopt.estsatfcb;

	SemiDCPara	 ->Text			=prcopt.mopt.ifsdp;/*Semi-Dynamic Correction Parameter*/
	SolDirFile	 ->Text			=prcopt.mopt.ofdir;/*Solution Dir*/
	FCBFile		 ->Text			=prcopt.mopt.iffcb;/*Satellite FCB*/

	EstIntZ		->Text     =s.sprintf(L"%.0f",prcopt.mopt.tiztd);
	EstIntES	->Text     =s.sprintf(L"%.0f",prcopt.mopt.tigra);
	RanWalSigZen->Text     =s.sprintf(L"%.2E",prcopt.mopt.sigtrop[0]);
	RanWalSigEW	->Text     =s.sprintf(L"%.2E",prcopt.mopt.sigtrop[1]);
	RanWalSigNS	->Text     =s.sprintf(L"%.2E",prcopt.mopt.sigtrop[2]);
	ThrO		->Text     =s.sprintf(L"%.1f",prcopt.mopt.cpomcth);
	ThrC		->Text     =s.sprintf(L"%.1f",prcopt.mopt.promcth);
	JudValWL	->Text     =s.sprintf(L"%.5f",prcopt.mopt.wlddfix);
	JudValL1	->Text     =s.sprintf(L"%.5f",prcopt.mopt.l1ddfix);
	ConCri		->Text     =s.sprintf(L"%.4f",prcopt.mopt.itrconv);
	MaxIte		->Text     =s.sprintf(L"%d",prcopt.mopt.itrmax);

	TemStoFile	 ->Text			=prcopt.mopt.eptmp;/*Temporary storage epoch parameters*/

 //	NLFcb		 ->Text     		=s.sprintf(L"%.0f",prcopt.mopt.tifcb);
	NLFcb		 ->Text     		=s.sprintf(L"%d",prcopt.mopt.tifcb);
	Minsd		 ->Text     		=s.sprintf(L"%.0f",prcopt.mopt.minpass);
	Mindd		 ->Text     		=s.sprintf(L"%.0f",prcopt.mopt.minddpass);
	Fixwl		 ->Text     		=s.sprintf(L"%.5f",prcopt.mopt.minconfw);
	Fixnl		 ->Text     		=s.sprintf(L"%.5f",prcopt.mopt.minconf1);
	Mobstadn		 ->Text     	=s.sprintf(L"%.4f",prcopt.mopt.sigr[0]);
	Mobstade		 ->Text     	=s.sprintf(L"%.4f",prcopt.mopt.sigr[1]);
	Mobstadu		 ->Text     	=s.sprintf(L"%.4f",prcopt.mopt.sigr[2]);
	Basestadn		 ->Text     	=s.sprintf(L"%.4f",prcopt.mopt.sigb[0]);
	Basestade		 ->Text     	=s.sprintf(L"%.4f",prcopt.mopt.sigb[1]);
	Basestadu		 ->Text     	=s.sprintf(L"%.4f",prcopt.mopt.sigb[2]);

	ReadAntList();
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::SaveOpt(UnicodeString file)
{
	TEdit *editu[]={RovPos1,RovPos2,RovPos3};
	TEdit *editr[]={RefPos1,RefPos2,RefPos3};
	wchar_t wbuff[1024],id[32];
	char *buff,*p,comment[256],s[64];
	int sat,ex;
	prcopt_t prcopt=prcopt_default;
	solopt_t solopt=solopt_default;
	filopt_t filopt={""};
	prcopt.mopt=mbsopt_default;
	
	prcopt.mode		=PosMode	 ->ItemIndex;
	if(Freqs->ItemIndex==0) {
		prcopt.nfreq=1;
		prcopt.oprfrq[0]=0;
		prcopt.oprfrq[1]=-1;
	}
	else if(Freqs->ItemIndex==1) {
		prcopt.nfreq=1;
		prcopt.oprfrq[0]=2;
		prcopt.oprfrq[1]=-1;
	}
	else if(Freqs->ItemIndex==2) {
		prcopt.nfreq=2;
		prcopt.oprfrq[0]=0;
		prcopt.oprfrq[1]=1;
		prcopt.oprfrq[2]=-1;
	}
	else if(Freqs->ItemIndex==3) {
		prcopt.nfreq=2;
		prcopt.oprfrq[0]=0;
		prcopt.oprfrq[1]=2;
		prcopt.oprfrq[2]=-1;
	}
	else if(Freqs->ItemIndex==4) {
		prcopt.nfreq=3;
		prcopt.oprfrq[0]=0;
		prcopt.oprfrq[1]=1;
		prcopt.oprfrq[2]=2;
		prcopt.oprfrq[3]=-1;
	}
	prcopt.soltype	=Solution	 ->ItemIndex;
	prcopt.elmin	=str2dbl(ElMask	->Text)*D2R;
//	prcopt.snrmin	=str2dbl(SnrMask->Text);
	prcopt.snrmask	=SnrMask;
	prcopt.dynamics	=DynamicModel->ItemIndex;
	prcopt.tidecorr	=TideCorr	 ->ItemIndex;
	prcopt.ionoopt	=IonoOpt	 ->ItemIndex;
	prcopt.tropopt	=TropOpt	 ->ItemIndex;
	prcopt.sateph	=SatEphem	 ->ItemIndex;
	if (ExSats->Text!=L"") {
	buff = wchar2char(ExSats->Text.c_str(),1024);

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
	prcopt.navsys	= (NavSys1->Checked?SYS_GPS:0)|
					  (NavSys2->Checked?SYS_GLO:0)|
					  (NavSys3->Checked?SYS_GAL:0)|
					  (NavSys4->Checked?SYS_QZS:0)|
					  (NavSys5->Checked?SYS_SBS:0)|
					  (NavSys6->Checked?SYS_CMP:0);

	char *cGloCodePri1 = wchar2char(GloCodePri1->Text.c_str(),MAXFCODE+2);
	strcpy(prcopt.codepri[ISYSGLO][0],cGloCodePri1);
	if(cGloCodePri1!=NULL){
		free(cGloCodePri1);
		cGloCodePri1=NULL;
	}
	char *cGloCodePri2 = wchar2char(GloCodePri2->Text.c_str(),MAXFCODE+2);
	strcpy(prcopt.codepri[ISYSGLO][1],cGloCodePri2);
	if(cGloCodePri2!=NULL){
		free(cGloCodePri2);
		cGloCodePri2=NULL;
	}

	prcopt.posopt[0]=PosOpt1	->Checked;
	prcopt.posopt[1]=PosOpt2	->Checked;
	prcopt.posopt[2]=PosOpt3	->Checked;
	prcopt.posopt[3]=PosOpt4	->Checked;
	prcopt.posopt[4]=PosOpt5	->Checked;
	prcopt.modertkar=AmbResMethod->ItemIndex;
	prcopt.modear	=AmbRes		->ItemIndex+1;
	prcopt.glomodear=GloAmbRes->ItemIndex;
	prcopt.modepppar    =PppAmbRes->ItemIndex;
	prcopt.thresar[0]=str2dbl(ValidThresAR->Text);
	prcopt.thresar[1]=str2dbl(ThresAR2->Text);
	prcopt.thresar[2]=str2dbl(ThresAR3->Text);
	prcopt.maxout	=str2dbl(OutCntResetAmb->Text);
	prcopt.minfix	=str2dbl(FixCntHoldAmb->Text);
	prcopt.minlock	=str2dbl(LockCntFixAmb->Text);
	prcopt.elmaskar	=str2dbl(ElMaskAR	->Text)*D2R;
	prcopt.elmaskhold=str2dbl(ElMaskHold->Text)*D2R;
	prcopt.maxtdiff	=str2dbl(MaxAgeDiff	->Text);
	prcopt.maxgdop	=str2dbl(RejectGdop ->Text);
	prcopt.maxinno	=str2dbl(RejectThres->Text);
	prcopt.thresslip=str2dbl(SlipThres	->Text);
	prcopt.niter	=str2dbl(NumIter	->Text);
	if (prcopt.mode==PMODE_MOVEB&&BaselineConst->Checked) {
		prcopt.baseline[0]=str2dbl(BaselineLen->Text);
		prcopt.baseline[1]=str2dbl(BaselineSig->Text);
	}
	solopt.posf		=SolFormat	->ItemIndex;
	solopt.timef	=TimeFormat	->ItemIndex==0?0:1;
	solopt.times	=TimeFormat	->ItemIndex==0?0:TimeFormat->ItemIndex-1;
	solopt.timeu	=str2dbl(TimeDecimal ->Text);
	solopt.degf		=LatLonFormat->ItemIndex;
	
	char *cFieldSep = wchar2char(FieldSep->Text.c_str(),64);
	strcpy(solopt.sep,cFieldSep);
	if(cFieldSep!=NULL){
		free(cFieldSep);
		cFieldSep=NULL;
	}
	solopt.outhead	=OutputHead	 ->ItemIndex;
	solopt.outopt	=OutputOpt	 ->ItemIndex;
	solopt.datum	=OutputDatum ->ItemIndex;
	solopt.height	=OutputHeight->ItemIndex;
	solopt.geoid	=OutputGeoid ->ItemIndex;
	solopt.solstatic=SolStatic   ->ItemIndex;
	solopt.nmeaintv[0]=str2dbl(NmeaIntv1->Text);
	solopt.nmeaintv[1]=str2dbl(NmeaIntv2->Text);
	solopt.trace	=DebugTrace	 ->ItemIndex;
	solopt.sstat	=DebugStatus ->ItemIndex;
	solopt.isbout	=IsbOut      ->ItemIndex;
	solopt.gl2out	=Gl2Out      ->ItemIndex;
	solopt.fcbout	=FCBOut      ->ItemIndex;

	solopt.possnxout=OutPosSinex ->ItemIndex;
	solopt.ionout   =OutIon      ->ItemIndex;
	solopt.tropout  =OutTrop     ->ItemIndex;
	solopt.recclout =OutReceiverClock->ItemIndex;
	solopt.satclout =OutSatelliteClock->ItemIndex;
	solopt.sbresout =OutStatic->ItemIndex;
	
	prcopt.eratio[0]=str2dbl(MeasErrR1->Text);
	prcopt.eratio[1]=str2dbl(MeasErrR2->Text);
	prcopt.eratio[2]=str2dbl(MeasErrR3->Text);
	prcopt.err[1]	=str2dbl(MeasErr2->Text);
	prcopt.err[2]	=str2dbl(MeasErr3->Text);
	prcopt.err[3]	=str2dbl(MeasErr4->Text);
	prcopt.err[4]	=str2dbl(MeasErr5->Text);
	prcopt.sclkstab	=str2dbl(SatClkStab->Text);
	prcopt.prn[0]	=str2dbl(PrNoise1->Text);
	prcopt.prn[1]	=str2dbl(PrNoise2->Text);
	prcopt.prn[2]	=str2dbl(PrNoise3->Text);
	prcopt.prn[3]	=str2dbl(PrNoise4->Text);
	prcopt.prn[4]	=str2dbl(PrNoise5->Text);
	prcopt.prn[5]	=str2dbl(PrNoise6->Text);
	prcopt.prn[6]	=str2dbl(PrNoise7->Text);
	prcopt.prn[7]	=str2dbl(PrNoise8->Text);
	
	char *cRovAnt = wchar2char(RovAnt->Text.c_str(),MAXANT);
	strcpy(prcopt.anttype[0],cRovAnt);
	if(cRovAnt!=NULL){
		free(cRovAnt);
		cRovAnt=NULL;
	}
	char *cRefAnt = wchar2char(RefAnt->Text.c_str(),MAXANT);
	strcpy(prcopt.anttype[1],cRefAnt);
	if(cRefAnt!=NULL){
		free(cRefAnt);
		cRefAnt=NULL;
	}
	prcopt.antdel[0][0]=str2dbl(RovAntE->Text);
	prcopt.antdel[0][1]=str2dbl(RovAntN->Text);
	prcopt.antdel[0][2]=str2dbl(RovAntU->Text);
	prcopt.antdel[1][0]=str2dbl(RefAntE->Text);
	prcopt.antdel[1][1]=str2dbl(RefAntN->Text);
	prcopt.antdel[1][2]=str2dbl(RefAntU->Text);

	prcopt.intpref	=IntpRefObs->ItemIndex;
	prcopt.sbassatsel=SbasSat->Text.ToInt();
	prcopt.rovpos=RovPosType->ItemIndex<3?0:RovPosType->ItemIndex-2;
	prcopt.refpos=RefPosType->ItemIndex<3?0:RefPosType->ItemIndex-2;
	if (prcopt.rovpos==0) GetPos(RovPosType->ItemIndex,editu,prcopt.ru);
	if (prcopt.refpos==0) GetPos(RefPosType->ItemIndex,editr,prcopt.rb);

	char *cRnxOpt1 = wchar2char(RnxOpts1->Text.c_str(),MAXSTRPATH);
	strcpy(prcopt.rnxopt[0],cRnxOpt1);
	if(cRnxOpt1!=NULL){
		free(cRnxOpt1);
		cRnxOpt1=NULL;
	}

	char *cRnxOpt2 = wchar2char(RnxOpts2->Text.c_str(),MAXSTRPATH);
	strcpy(prcopt.rnxopt[0],cRnxOpt2);
	if(cRnxOpt2!=NULL){
		free(cRnxOpt2);
		cRnxOpt2=NULL;
	}

	char *cSatPcvFile = wchar2char(SatPcvFile->Text.c_str(),MAXSTRPATH);
	strcpy(filopt.satantp,cSatPcvFile);
	if(cSatPcvFile!=NULL){
		free(cSatPcvFile);
		cSatPcvFile=NULL;
	}
	char *cAntPcvFile = wchar2char(AntPcvFile->Text.c_str(),MAXSTRPATH);
	strcpy(filopt.rcvantp,cAntPcvFile);
	if(cAntPcvFile!=NULL){
		free(cAntPcvFile);
		cAntPcvFile=NULL;
	}
	char *cStaPosFile = wchar2char(StaPosFile->Text.c_str(),MAXSTRPATH);
	strcpy(filopt.stapos,cStaPosFile);
	if(cStaPosFile!=NULL){
		free(cStaPosFile);
		cStaPosFile=NULL;
	}
	char *cGeoidDataFile = wchar2char(GeoidDataFile->Text.c_str(),MAXSTRPATH);
	strcpy(filopt.geoid,cGeoidDataFile);
	if(cGeoidDataFile!=NULL){
		free(cGeoidDataFile);
		cGeoidDataFile=NULL;
	}
	char *cDCBFile = wchar2char(DCBFile->Text.c_str(),MAXSTRPATH);
	strcpy(filopt.dcb,cDCBFile);
	if(cDCBFile!=NULL){
		free(cDCBFile);
		cDCBFile=NULL;
	}
	char *cISBFile = wchar2char(ISBFile->Text.c_str(),MAXSTRPATH);
	strcpy(filopt.isb,cISBFile);
	if(cISBFile!=NULL){
		free(cISBFile);
		cISBFile=NULL;
	}
	char *cISBOutFile = wchar2char(ISBOutFile->Text.c_str(),MAXSTRPATH);
	strcpy(solopt.isbfile,cISBOutFile);
	if(cISBOutFile!=NULL){
		free(cISBOutFile);
		cISBOutFile=NULL;
	}
	char *cGL2OutFile = wchar2char(GL2OutFile->Text.c_str(),MAXSTRPATH);
	strcpy(solopt.gl2file,cGL2OutFile);
	if(cGL2OutFile!=NULL){
		free(cGL2OutFile);
		cGL2OutFile=NULL;
	}
	char *cFCBOutFile = wchar2char(FCBOutFile->Text.c_str(),MAXSTRPATH);
	strcpy(solopt.fcb,cFCBOutFile);
	if(cFCBOutFile!=NULL){
		free(cFCBOutFile);
		cFCBOutFile=NULL;
	}
	char *cEOPFile = wchar2char(EOPFile->Text.c_str(),MAXSTRPATH);
	strcpy(filopt.eop,cEOPFile);
	if(cEOPFile!=NULL){
		free(cEOPFile);
		cEOPFile=NULL;
	}
	char *cBLQFile = wchar2char(BLQFile->Text.c_str(),MAXSTRPATH);
	strcpy(filopt.blq,cBLQFile);
	if(cBLQFile!=NULL){
		free(cBLQFile);
		cBLQFile=NULL;
	}
	char *cIonoFile = wchar2char(IonoFile->Text.c_str(),MAXSTRPATH);
	strcpy(filopt.iono,cIonoFile);
	if(cIonoFile!=NULL){
		free(cIonoFile);
		cIonoFile=NULL;
	}
	char *cGoogleEarthFile = wchar2char(GoogleEarthFile->Text.c_str(),MAXSTRPATH);
	strcpy(filopt.geexe,cGoogleEarthFile);
	if(cGoogleEarthFile!=NULL){
		free(cGoogleEarthFile);
		cGoogleEarthFile=NULL;
	}
	prcopt.l2cprior	=	L2Cod        ->ItemIndex;/* L2 code priority (0:L2P,1:L2C) */
	prcopt.tsyscorr	=	TimSys       ->ItemIndex;/* time system correction ( 0:off,1:on) */
	prcopt.phasshft	=	PhaCyc       ->ItemIndex;/* phase cycle shift (0:off,1:table) */
    prcopt.gl2bias	=	L2CPBias     ->ItemIndex;
	prcopt.errmodel	=	ErrMod       ->ItemIndex;/* error model (0:user settings,1:table) */
	prcopt.dcbratio =str2dbl(MeasErr6  ->Text);/* code error ratio(no DCB) */
	prcopt.isb	    =	Isb          ->ItemIndex;/* ISB (0:off,1:table,2:est,3:est(code),4:est(phase)) */
	prcopt.diff	    =	Diff         ->ItemIndex;

	/* receiver types {rover,base} */
	char *cRovRecTyp = wchar2char(RovRecTyp->Text.c_str(),-1);
	strcpy(prcopt.rectype[0],cRovRecTyp);
	if(cRovRecTyp!=NULL){
		free(cRovRecTyp);
		cRovRecTyp=NULL;
	}
	char *cRefRecTyp = wchar2char(RefRecTyp->Text.c_str(),-1);
	strcpy(prcopt.rectype[1],cRefRecTyp);
	if(cRefRecTyp!=NULL){
		free(cRefRecTyp);
		cRefRecTyp=NULL;
	}
	/* 1/4cycle phase correction file */
	char *cPhaCycFile = wchar2char(PhaCycFile->Text.c_str(),-1);
	strcpy(prcopt.mopt.ifpcs,cPhaCycFile);
	if(cPhaCycFile!=NULL){
		free(cPhaCycFile);
		cPhaCycFile=NULL;
	}
	/* GLONASS IFB table file */
	char *cGloIfbFile = wchar2char(GloIfbFile->Text.c_str(),-1);
	strcpy(prcopt.mopt.ififb,cGloIfbFile);
	if(cGloIfbFile!=NULL){
		free(cGloIfbFile);
		cGloIfbFile=NULL;
	}
	/* error model file */
	char *cErrModFile = wchar2char(ErrModFile->Text.c_str(),-1);
	strcpy(prcopt.mopt.iferr,cErrModFile);
	if(cErrModFile!=NULL){
		free(cErrModFile);
		cErrModFile=NULL;
	}
	/*BIPMCircularTFile*/
	char *cBIPMCircularTFile = wchar2char(BIPMCircularTFile->Text.c_str(),-1);
	strcpy(filopt.cirtfile,cBIPMCircularTFile);
	if(cBIPMCircularTFile!=NULL){
		free(cBIPMCircularTFile);
		cBIPMCircularTFile=NULL;
	}
	/*Estimate Satellite Clock/FCB*/
	prcopt.mopt.estsatclk = EstSatClo	 ->ItemIndex;
	prcopt.mopt.estsatfcb = EstSatFCB	 ->ItemIndex;

	/*Semi-Dynamic Correction Parameter*/
	char *cSemiDCPara = wchar2char(SemiDCPara->Text.c_str(),-1);
	strcpy(prcopt.mopt.ifsdp,cSemiDCPara);
	if(cSemiDCPara!=NULL){
		free(cSemiDCPara);
		cSemiDCPara=NULL;
	}
	/*Solution Dir*/
	char *cSolDirFile = wchar2char(SolDirFile->Text.c_str(),-1);
	strcpy(prcopt.mopt.ofdir,cSolDirFile);
	if(cSolDirFile!=NULL){
		free(cSolDirFile);
		cSolDirFile=NULL;
	}
	/*Satellite FCB*/
	char *cFCBFile = wchar2char(FCBFile->Text.c_str(),-1);
	strcpy(prcopt.mopt.iffcb,cFCBFile);
	if(cFCBFile!=NULL){
		free(cFCBFile);
		cFCBFile=NULL;
	}
	prcopt.mopt.tiztd      = str2dbl(EstIntZ	 ->Text);
	prcopt.mopt.tigra      = str2dbl(EstIntES	 ->Text);
	prcopt.mopt.sigtrop[0] = str2dbl(RanWalSigZen->Text);
	prcopt.mopt.sigtrop[1] = str2dbl(RanWalSigEW ->Text);
	prcopt.mopt.sigtrop[2] = str2dbl(RanWalSigNS ->Text);
	prcopt.mopt.cpomcth    = str2dbl(ThrO		 ->Text);
	prcopt.mopt.promcth    = str2dbl(ThrC		 ->Text);
	prcopt.mopt.wlddfix    = str2dbl(JudValWL	 ->Text);
	prcopt.mopt.l1ddfix    = str2dbl(JudValL1	 ->Text);
	prcopt.mopt.itrconv    = str2dbl(ConCri		 ->Text);
	prcopt.mopt.itrmax     = str2int(MaxIte		 ->Text);

	/*Temporary storage epoch parameters*/
	char *cTemStoFile = wchar2char(TemStoFile->Text.c_str(),-1);
	strcpy(prcopt.mopt.eptmp,cTemStoFile);
	if(cTemStoFile!=NULL){
		free(cTemStoFile);
		cTemStoFile=NULL;
	}
	prcopt.mopt.tifcb         = str2int(NLFcb		->Text);
	prcopt.mopt.minpass       = str2dbl(Minsd		->Text);
	prcopt.mopt.minddpass     = str2dbl(Mindd		->Text);
	prcopt.mopt.minconfw      = str2dbl(Fixwl		->Text);
	prcopt.mopt.minconf1      = str2dbl(Fixnl		->Text);
	prcopt.mopt.sigr[0]       = str2dbl(Mobstadn	->Text);
	prcopt.mopt.sigr[1]       = str2dbl(Mobstade	->Text);
	prcopt.mopt.sigr[2]       = str2dbl(Mobstadu	->Text);
	prcopt.mopt.sigb[0]       = str2dbl(Basestadn	->Text);
	prcopt.mopt.sigb[1]       = str2dbl(Basestade	->Text);
	prcopt.mopt.sigb[2]       = str2dbl(Basestadu	->Text);

	time2str(utc2gpst(timeget()),s,0);
	sprintf(comment,"gsipost options (%s, v.%s)",s,VER_RTKLIB);
	setsysopts(&prcopt,&solopt,&filopt);

	char *cfile = wchar2char(file.c_str(),-1);
	if (!saveopts(cfile,"w",comment,sysopts)){
		if(cfile!=NULL){
			free(cfile);
			cfile=NULL;
		}
		return;
	}
	if(cfile!=NULL){
		free(cfile);
		cfile=NULL;
	}
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::UpdateEnable(void)
{
	int rel=PMODE_DGPS<=PosMode->ItemIndex&&PosMode->ItemIndex<=PMODE_FIXED;
	int rtk=PMODE_KINEMA<=PosMode->ItemIndex&&PosMode->ItemIndex<=PMODE_FIXED;
	int ppp=PosMode->ItemIndex>=PMODE_PPP_KINEMA;
	int mbl=PosMode->ItemIndex==PMODE_MULTI;
	int ar=rtk||ppp;
	
//	Freqs          ->Enabled=rel;
	Freqs          ->Enabled=(rel||ppp)&&(IonoOpt->ItemIndex!=IONOOPT_IFLC);
	Solution       ->Enabled=rel||ppp;
	DynamicModel   ->Enabled=rel;
	TideCorr       ->Enabled=rel||ppp;
	//IonoOpt        ->Enabled=!ppp;
	PosOpt1        ->Enabled=ppp;
	PosOpt2        ->Enabled=ppp;
	PosOpt3        ->Enabled=ppp;
	PosOpt4        ->Enabled=ppp;

	AmbResMethod   ->Enabled=rtk;
	AmbRes         ->Enabled=rtk;
	GloAmbRes      ->Enabled=ar&&AmbResMethod->ItemIndex>0&&NavSys2->Checked;
	PppAmbRes      ->Enabled=ppp;
//	ValidThresAR   ->Enabled=ar&&AmbRes->ItemIndex>=1&&AmbRes->ItemIndex<3;
	ValidThresAR   ->Enabled=(rtk&&AmbResMethod->ItemIndex!=RTKAR_OFF)||(ppp&&PppAmbRes->ItemIndex!=PPPAR_OFF);
	ThresAR2	   ->Enabled=ar&&AmbRes->ItemIndex>=3;
	ThresAR3	   ->Enabled=ar&&AmbRes->ItemIndex>=3;
//	LockCntFixAmb  ->Enabled=ar&&AmbRes->ItemIndex>=0;
	LockCntFixAmb  ->Enabled=ar&&AmbResMethod->ItemIndex!=RTKAR_OFF;
//	ElMaskAR       ->Enabled=ar&&AmbRes->ItemIndex>=0;
	ElMaskAR       ->Enabled=ar&&AmbResMethod->ItemIndex!=RTKAR_OFF;
	OutCntResetAmb ->Enabled=ar;
	FixCntHoldAmb  ->Enabled=ar&&AmbRes->ItemIndex==2;
	ElMaskHold     ->Enabled=ar&&AmbRes->ItemIndex==2;
	LabelHold      ->Enabled=ar&&AmbRes->ItemIndex==2;
	SlipThres      ->Enabled=ar;
	MaxAgeDiff     ->Enabled=rel;
	RejectThres    ->Enabled=rel||ppp;
	NumIter        ->Enabled=rel||ppp;
	BaselineConst  ->Enabled=PosMode->ItemIndex==PMODE_MOVEB;
	BaselineLen    ->Enabled=BaselineConst->Checked&&PosMode->ItemIndex==PMODE_MOVEB;
	BaselineSig    ->Enabled=BaselineConst->Checked&&PosMode->ItemIndex==PMODE_MOVEB;
	
	OutputHead     ->Enabled=SolFormat->ItemIndex!=3;
	OutputOpt      ->Enabled=SolFormat->ItemIndex!=3;
	TimeFormat     ->Enabled=SolFormat->ItemIndex!=3;
	TimeDecimal    ->Enabled=SolFormat->ItemIndex!=3;
	LatLonFormat   ->Enabled=SolFormat->ItemIndex==0;
	FieldSep       ->Enabled=SolFormat->ItemIndex!=3;
	OutputDatum    ->Enabled=SolFormat->ItemIndex==0;
	OutputHeight   ->Enabled=SolFormat->ItemIndex==0;
	OutputGeoid    ->Enabled=SolFormat->ItemIndex==0&&OutputHeight->ItemIndex==1;
	SolStatic      ->Enabled=PosMode->ItemIndex==PMODE_STATIC||
							 PosMode->ItemIndex==PMODE_PPP_STATIC;

	BtnISBOutView  ->Enabled=IsbOut->ItemIndex!=0;
	BtnISBOutFile  ->Enabled=IsbOut->ItemIndex!=0;
	ISBOutFile     ->Enabled=IsbOut->ItemIndex!=0;

	BtnGL2OutView  ->Enabled=Gl2Out->ItemIndex!=0;
	BtnGL2OutFile  ->Enabled=Gl2Out->ItemIndex!=0;
	GL2OutFile     ->Enabled=Gl2Out->ItemIndex!=0;

	BtnFCBOutView  ->Enabled=FCBOut->ItemIndex!=0;
	BtnFCBOutFile  ->Enabled=FCBOut->ItemIndex!=0;
	FCBOutFile     ->Enabled=FCBOut->ItemIndex!=0;

	OutSatelliteClock->Enabled=mbl;

	RovAntPcv      ->Enabled=rel||ppp;
	RovAnt         ->Enabled=(rel||ppp)&&RovAntPcv->Checked;
	RovAntE        ->Enabled=(rel||ppp)&&RovAntPcv->Checked;
	RovAntN        ->Enabled=(rel||ppp)&&RovAntPcv->Checked;
	RovAntU        ->Enabled=(rel||ppp)&&RovAntPcv->Checked;
	LabelRovAntD   ->Enabled=(rel||ppp)&&RovAntPcv->Checked;
	RefAntPcv      ->Enabled=rel;
	RefAnt         ->Enabled=rel&&RefAntPcv->Checked;
	RefAntE        ->Enabled=rel&&RefAntPcv->Checked;
	RefAntN        ->Enabled=rel&&RefAntPcv->Checked;
	RefAntU        ->Enabled=rel&&RefAntPcv->Checked;
	LabelRefAntD   ->Enabled=rel&&RefAntPcv->Checked;
	
	RovPosType     ->Enabled=PosMode->ItemIndex==PMODE_FIXED||PosMode->ItemIndex==PMODE_PPP_FIXED;
	RovPos1        ->Enabled=RovPosType->Enabled&&RovPosType->ItemIndex<=2;
	RovPos2        ->Enabled=RovPosType->Enabled&&RovPosType->ItemIndex<=2;
	RovPos3        ->Enabled=RovPosType->Enabled&&RovPosType->ItemIndex<=2;
	BtnRovPos      ->Enabled=RovPosType->Enabled&&RovPosType->ItemIndex<=2;
	
	RefPosType     ->Enabled=rel&&PosMode->ItemIndex!=PMODE_MOVEB;
	RefPos1        ->Enabled=RefPosType->Enabled&&RefPosType->ItemIndex<=2;
	RefPos2        ->Enabled=RefPosType->Enabled&&RefPosType->ItemIndex<=2;
	RefPos3        ->Enabled=RefPosType->Enabled&&RefPosType->ItemIndex<=2;
	BtnRefPos      ->Enabled=RefPosType->Enabled&&RefPosType->ItemIndex<=2;

	MeasErrR1	   ->Enabled=ErrMod->ItemIndex==0;/* error model==(1:table) disenabled */
	MeasErrR2	   ->Enabled=ErrMod->ItemIndex==0;/* error model==(1:table) disenabled */
	MeasErrR3	   ->Enabled=ErrMod->ItemIndex==0;/* error model==(1:table) disenabled */
	MeasErr2	   ->Enabled=ErrMod->ItemIndex==0;/* error model==(1:table) disenabled */
	MeasErr3	   ->Enabled=ErrMod->ItemIndex==0;/* error model==(1:table) disenabled */
	Label6		   ->Enabled=ErrMod->ItemIndex==0;/* error model==(1:table) disenabled */
	Label7		   ->Enabled=ErrMod->ItemIndex==0;/* error model==(1:table) disenabled */
#if 0
	RovRecTyp	   ->Enabled=GloAmbRes->ItemIndex==GLO_ARMODE_IFB;/* GLONASS AR mode==(3:use IFB table) */
	RefRecTyp	   ->Enabled=GloAmbRes->ItemIndex==GLO_ARMODE_IFB;/* GLONASS AR mode==(3:use IFB table) */
	Label58		   ->Enabled=GloAmbRes->ItemIndex==GLO_ARMODE_IFB;/* GLONASS AR mode==(3:use IFB table) */
	Label59		   ->Enabled=GloAmbRes->ItemIndex==GLO_ARMODE_IFB;/* GLONASS AR mode==(3:use IFB table) */
	PhaCycFile	   ->Enabled=PhaCyc->ItemIndex==1;/* phase cycle shift==(1:table) */
	GloIfbFile	   ->Enabled=GloAmbRes->ItemIndex==GLO_ARMODE_IFB;/* GLONASS AR mode==(3:use IFB table) */
	ErrModFile	   ->Enabled=ErrMod->ItemIndex==1;/* error model==(1:table) */
	IsbFile	       ->Enabled=Isb->ItemIndex==1;/* isb==(1:table) */

	PhaCycBtn	   ->Enabled=PhaCyc->ItemIndex==1;/* phase cycle shift==(1:table) */
	GloIfbBtn	   ->Enabled=GloAmbRes->ItemIndex==GLO_ARMODE_IFB;/* GLONASS AR mode==(3:use IFB table) */
	ErrModBtn	   ->Enabled=ErrMod->ItemIndex==1;/* error model==(1:table) */
	ISbBtn	       ->Enabled=Isb->ItemIndex==1;/* Isb==(1:table) */

	PhaCycSBtn	   ->Enabled=PhaCyc->ItemIndex==1;/* phase cycle shift==(1:table) */
	GloIfbSBtn	   ->Enabled=GloAmbRes->ItemIndex==GLO_ARMODE_IFB;/* GLONASS AR mode==(3:use IFB table) */
	ErrModSBtn	   ->Enabled=ErrMod->ItemIndex==1;/* error model==(1:table) */
	ISbSBtn   	   ->Enabled=ISb->ItemIndex==1;/* Isb==(1:table) */
#endif

//	FCBBtn	   	   ->Enabled=PosMode->ItemIndex==PMODE_PPPAR;/* positioning mode: PPP-AR */
//	FCBFile	   	   ->Enabled=PosMode->ItemIndex==PMODE_PPPAR;/* positioning mode: PPP-AR */
//	FCBLabel	   ->Enabled=PosMode->ItemIndex==PMODE_PPPAR;/* positioning mode: PPP-AR */



}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::GetPos(int type, TEdit **edit, double *pos)
{
	double p[3]={0},dms1[3]={0},dms2[3]={0};
	
	if (type==1) { /* lat/lon/height dms/m */
		swscanf(edit[0]->Text.c_str(),L"%lf %lf %lf",dms1,dms1+1,dms1+2);
		swscanf(edit[1]->Text.c_str(),L"%lf %lf %lf",dms2,dms2+1,dms2+2);
		p[0]=(dms1[0]<0?-1:1)*(fabs(dms1[0])+dms1[1]/60+dms1[2]/3600)*D2R;
		p[1]=(dms1[0]<0?-1:1)*(fabs(dms2[0])+dms2[1]/60+dms2[2]/3600)*D2R;
		p[2]=str2dbl(edit[2]->Text);
		pos2ecef(p,pos);
	}
	else if (type==2) { /* x/y/z-ecef */
		pos[0]=str2dbl(edit[0]->Text);
		pos[1]=str2dbl(edit[1]->Text);
		pos[2]=str2dbl(edit[2]->Text);
	}
	else {
		p[0]=str2dbl(edit[0]->Text)*D2R;
		p[1]=str2dbl(edit[1]->Text)*D2R;
		p[2]=str2dbl(edit[2]->Text);
		pos2ecef(p,pos);
	}
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::SetPos(int type, TEdit **edit, double *pos)
{
	UnicodeString s;
	double p[3],dms1[3],dms2[3],s1,s2;
	
	if (type==1) { /* lat/lon/height dms/m */
		ecef2pos(pos,p); s1=p[0]<0?-1:1; s2=p[1]<0?-1:1;
		p[0]=fabs(p[0])*R2D+1E-12; p[1]=fabs(p[1])*R2D+1E-12;
		dms1[0]=floor(p[0]); p[0]=(p[0]-dms1[0])*60.0;
		dms1[1]=floor(p[0]); dms1[2]=(p[0]-dms1[1])*60.0;
		dms2[0]=floor(p[1]); p[1]=(p[1]-dms2[0])*60.0;
		dms2[1]=floor(p[1]); dms2[2]=(p[1]-dms2[1])*60.0;
		edit[0]->Text=s.sprintf(L"%.0f %02.0f %09.6f",s1*dms1[0],dms1[1],dms1[2]);
		edit[1]->Text=s.sprintf(L"%.0f %02.0f %09.6f",s2*dms2[0],dms2[1],dms2[2]);
		edit[2]->Text=s.sprintf(L"%.4f",p[2]);
	}
	else if (type==2) { /* x/y/z-ecef */
		edit[0]->Text=s.sprintf(L"%.4f",pos[0]);
		edit[1]->Text=s.sprintf(L"%.4f",pos[1]);
		edit[2]->Text=s.sprintf(L"%.4f",pos[2]);
	}
	else {
		ecef2pos(pos,p);
		edit[0]->Text=s.sprintf(L"%.9f",p[0]*R2D);
		edit[1]->Text=s.sprintf(L"%.9f",p[1]*R2D);
		edit[2]->Text=s.sprintf(L"%.4f",p[2]);
	}
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::ReadAntList(void)
{
	TStringList *list;
	pcvs_t pcvs={0};
	char *p;
	
	char *ch = wchar2char(AntPcvFile->Text.c_str(),-1);
	if (!readpcv(ch,&pcvs)){
		if(ch!=NULL){
			free(ch);
			ch=NULL;
		}
		return;
	}
	if(ch!=NULL){
		free(ch);
		ch=NULL;
	}
	
	list=new TStringList;
	list->Add("");
	list->Add("*");
	
	for (int i=0;i<pcvs.n;i++) {
		if (pcvs.pcv[i].sat) continue;
		if ((p=strchr(pcvs.pcv[i].type,' '))) *p='\0';
		if (i>0&&!strcmp(pcvs.pcv[i].type,pcvs.pcv[i-1].type)) continue;
		wchar_t *wtype = char2wchar(pcvs.pcv[i].type,MAXANT);

		list->Add(wtype);
		if(wtype!=NULL){
			free(wtype);
			wtype=NULL;
		}
	}
	RovAnt->Items=list;
	RefAnt->Items=list;
	
	free(pcvs.pcv);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::SpeedButton1Click(TObject *Sender)
{
	KeyDialog->Flag=2;
	KeyDialog->Show();
}

/* ------------------------------------------------
* convert char to wchar
* args   : char *str
*          int length
* return : wchar_t* 
* notes  : 
*--------------------------------------------------*/
wchar_t* __fastcall TOptDialog::char2wchar(char *str, int length)
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
char* __fastcall TOptDialog::wchar2char(wchar_t *wstr, int length)
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
void __fastcall TOptDialog::PhaCycBtnClick(TObject *Sender)
{
	OpenDialog->Title="Phase Cycle Shift File";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	PhaCycFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------

void __fastcall TOptDialog::GloIfbBtnClick(TObject *Sender)
{
	OpenDialog->Title="GLONASS IFB File";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	GloIfbFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------

void __fastcall TOptDialog::ErrModBtnClick(TObject *Sender)
{
	OpenDialog->Title="Error Model File";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	ErrModFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------

void __fastcall TOptDialog::PhaCycSBtnClick(TObject *Sender)
{
	if (PhaCycFile->Text=="") return;
	TTextViewer *viewer=new TTextViewer(Application);
	viewer->Show();
	viewer->Read(PhaCycFile->Text);
}
//---------------------------------------------------------------------------

void __fastcall TOptDialog::GloIfbSBtnClick(TObject *Sender)
{
	if (GloIfbFile->Text=="") return;
	TTextViewer *viewer=new TTextViewer(Application);
	viewer->Show();
	viewer->Read(GloIfbFile->Text);
}
//---------------------------------------------------------------------------

void __fastcall TOptDialog::ErrModSBtnClick(TObject *Sender)
{
	if (ErrModFile->Text=="") return;
	TTextViewer *viewer=new TTextViewer(Application);
	viewer->Show();
	viewer->Read(ErrModFile->Text);
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnBIPMCircularTFileClick(TObject *Sender)
{
	OpenDialog->Title="BIPM Circular T File";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	BIPMCircularTFile->Text=OpenDialog->FileName;
}
void __fastcall TOptDialog::EstSatFCBChange(TObject *Sender)
{
	UpdateEnable();
}
void __fastcall TOptDialog::EstSatCloChange(TObject *Sender)
{
	UpdateEnable();
}

void __fastcall TOptDialog::SemiDCParaBtnClick(TObject *Sender)
{
	OpenDialog->Title="Semi-Dynamic Correction Parameter";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	SemiDCPara->Text=OpenDialog->FileName;
}

void __fastcall TOptDialog::SolDirBtnClick(TObject *Sender)
{
	OpenDialog->Title="Solution Dir";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	SolDirFile->Text=OpenDialog->FileName;
}

void __fastcall TOptDialog::FCBBtnClick(TObject *Sender)
{
	OpenDialog->Title="Satellite FCB";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	FCBFile->Text=OpenDialog->FileName;
}


void __fastcall TOptDialog::TemStoFileButClick(TObject *Sender)
{
	OpenDialog->Title="Temporary storage epoch parameters";
	OpenDialog->FilterIndex=1;
	if (!OpenDialog->Execute()) return;
	TemStoFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------
void __fastcall TOptDialog::BtnMaskClick(TObject *Sender)
{
	MaskOptDialog->Mask=SnrMask;
	if (MaskOptDialog->ShowModal()!=mrOk) return;
	SnrMask=MaskOptDialog->Mask;
}
//---------------------------------------------------------------------------

void __fastcall TOptDialog::GloCodePri1KeyPress(TObject *Sender, System::WideChar &Key)
{
	Key = checkCode(SYS_GLO,0,AnsiString(Key).c_str());
}
//---------------------------------------------------------------------------

void __fastcall TOptDialog::GloCodePri2KeyPress(TObject *Sender, System::WideChar &Key)
{
	Key = checkCode(SYS_GLO,1,AnsiString(Key).c_str());
}
//---------------------------------------------------------------------------


