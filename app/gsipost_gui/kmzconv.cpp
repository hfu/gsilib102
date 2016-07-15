/*------------------------------------------------------------------------------
* kmzconv.cpp
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
#include <stdio.h>
#include <string.h>
#pragma hdrstop

#include "postmain.h"
#include "kmzconv.h"
#include "rtklib.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TConvDialog *ConvDialog;
//---------------------------------------------------------------------------
static double str2dbl(UnicodeString str)
{
	double val=0.0;
	swscanf(str.c_str(),L"%lf",&val);
	return val;
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::FormShow(TObject *Sender)
{
//	GoogleEarthFile->Text=MainForm->GoogleEarthFile;
}
//---------------------------------------------------------------------------
__fastcall TConvDialog::TConvDialog(TComponent* Owner)
	: TForm(Owner)
{
	setlocale( LC_ALL, "" );
	_wsetlocale(LC_ALL, L"" );
	UpdateEnable();
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::SetInput(UnicodeString File)
{
	InputFile->Text=File;
	UpdateOutFile();
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::TimeSpanClick(TObject *Sender)
{
	UpdateEnable();	
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::AddOffsetClick(TObject *Sender)
{
	UpdateEnable();	
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::TimeIntFClick(TObject *Sender)
{
	UpdateEnable();	
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::BtnInputFileClick(TObject *Sender)
{
	OpenDialog->FileName=InputFile->Text;
	if (!OpenDialog->Execute()) return;
	InputFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::BtnGoogleClick(TObject *Sender)
{
	wchar_t cmd[1024];
	swprintf(cmd,L"\"%s\" \"%s\"",MainForm->GoogleEarthFile,OutputFile->Text.c_str());
	if (!ExecCmd(cmd)) ShowMsg(L"error : google earth execution");
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::BtnConvertClick(TObject *Sender)
{
	int stat;
	wchar_t cmd[1024],file[1024],kmlfile[1024],*p;
	double offset[3]={0},es[6]={1970,1,1},ee[6]={2038,1,1},tint=0.0;
	gtime_t ts={0},te={0};
	ShowMsg(L"");
	if (InputFile->Text==L""||OutputFile->Text==L"") return;
	ShowMsg(L"converting ...");
	if (TimeSpan->Checked) {
		swscanf(TimeY1->Text.c_str(),L"%lf/%lf/%lf",es  ,es+1,es+2);
		swscanf(TimeH1->Text.c_str(),L"%lf:%lf:%lf",es+3,es+4,es+5);
		swscanf(TimeY2->Text.c_str(),L"%lf/%lf/%lf",ee  ,ee+1,ee+2);
		swscanf(TimeH2->Text.c_str(),L"%lf:%lf:%lf",ee+3,ee+4,ee+5);
		ts=epoch2time(es);
		te=epoch2time(ee);
	}
	if (AddOffset->Checked) {
		offset[0]=str2dbl(Offset1->Text);
		offset[1]=str2dbl(Offset2->Text);
		offset[2]=str2dbl(Offset3->Text);
	}
	if (TimeIntF->Checked) tint=str2dbl(TimeInt->Text);
	wcscpy(file,InputFile->Text.c_str());
	wcscpy(kmlfile,file);
	p=wcsrchr(kmlfile,'.');
	if (p==NULL) p=kmlfile+wcslen(kmlfile);
	wcscpy(p,L".kml");

	char *cfile = wchar2char(file,1024);
	char *ckmlfile = wchar2char(kmlfile,1024);
	char *ch = wchar2char(OutputFile->Text.c_str(),-1);

	if((stat=convkml(cfile,Compress->Checked?ckmlfile:ch,
					 ts,te,tint,QFlags->ItemIndex,offset,
					 TrackColor->ItemIndex,PointColor->ItemIndex,
					 OutputAlt->ItemIndex,OutputTime->ItemIndex))<0) {
		if      (stat==-1) ShowMsg(L"error : read input file");
		else if (stat==-2) ShowMsg(L"error : input file format");
		else if (stat==-3) ShowMsg(L"error : no data in input file");
		else               ShowMsg(L"error : write kml file");
		if(cfile!=NULL){
			free(cfile);
			cfile=NULL;
		}
		if(ckmlfile!=NULL){
			free(ckmlfile);
			ckmlfile=NULL;
		}
		if(ch!=NULL){
			free(ch);
			ch=NULL;
		}
		return;
	}
	if(cfile!=NULL){
		free(cfile);
		cfile=NULL;
	}
	if(ckmlfile!=NULL){
		free(ckmlfile);
		ckmlfile=NULL;
	}
	if(ch!=NULL){
		free(ch);
		ch=NULL;
	}
	if (Compress->Checked) {
		swprintf(cmd,L"zip.exe -j -m %s %s",OutputFile->Text.c_str(),kmlfile);
		if (!ExecCmd(cmd)) {
			ShowMsg(L"error : zip execution");
			return;
		}
	}
	ShowMsg(L"done");
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::BtnCloseClick(TObject *Sender)
{
	Close();
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::CompressClick(TObject *Sender)
{
	UpdateOutFile();
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::UpdateEnable(void)
{
	Offset1->Enabled=AddOffset->Checked;
	Offset2->Enabled=AddOffset->Checked;
	Offset3->Enabled=AddOffset->Checked;
	TimeY1->Enabled=TimeSpan->Checked;
	TimeH1->Enabled=TimeSpan->Checked;
	TimeY2->Enabled=TimeSpan->Checked;
	TimeH2->Enabled=TimeSpan->Checked;
	TimeY1UD->Enabled=TimeSpan->Checked;
	TimeH1UD->Enabled=TimeSpan->Checked;
	TimeY2UD->Enabled=TimeSpan->Checked;
	TimeH2UD->Enabled=TimeSpan->Checked;
	TimeInt->Enabled=TimeIntF->Checked;
}
//---------------------------------------------------------------------------
int __fastcall TConvDialog::ExecCmd(wchar_t *cmd)
{
	STARTUPINFOW si={0};
	PROCESS_INFORMATION info;
	si.cb=sizeof(si);
	if (!CreateProcessW(NULL,cmd,NULL,NULL,false,CREATE_NO_WINDOW,NULL,NULL,&si,
					   &info)) return 0;
	CloseHandle(info.hProcess);
	CloseHandle(info.hThread);
	return 1;
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::ShowMsg(UnicodeString msg)
{
	Message->Caption=msg;
	if (wcsstr(msg.c_str(),L"error")) Message->Font->Color=clRed;
	else Message->Font->Color=clBlue;
	Application->ProcessMessages();
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::InputFileChange(TObject *Sender)
{
	UpdateOutFile();
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::UpdateOutFile(void)
{
	wchar_t file[256],*p;
	if (InputFile->Text=="") return;
	wcscpy(file,InputFile->Text.c_str());
	p=wcsrchr(file,'.');
	if (p==NULL) p=file+wcslen(file);
	wcscpy(p,Compress->Checked?L".kmz":L".kml");
	OutputFile->Text=file;
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::TimeY1UDChangingEx(TObject *Sender,
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
//---------------------------------------------------------------------------
void __fastcall TConvDialog::TimeH1UDChangingEx(TObject *Sender,
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
//---------------------------------------------------------------------------
void __fastcall TConvDialog::TimeY2UDChangingEx(TObject *Sender,
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
//---------------------------------------------------------------------------
void __fastcall TConvDialog::TimeH2UDChangingEx(TObject *Sender,
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

/* ------------------------------------------------
* convert char to wchar
* args   : char *str
*          int length
* return : wchar_t* 
* notes  : 
*--------------------------------------------------*/
wchar_t* __fastcall TConvDialog::char2wchar(char *str, int length)
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
char* __fastcall TConvDialog::wchar2char(wchar_t *wstr, int length)
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
//---------------------------------------------------------------------------
void __fastcall TConvDialog::GoogleEarthFileChange(TObject *Sender)
{
//	MainForm->GoogleEarthFile=GoogleEarthFile->Text;
}
//---------------------------------------------------------------------------
void __fastcall TConvDialog::BtnGoogleEarthFileClick(TObject *Sender)
{
	OpenDialog->Title="Google Earth Exe File";
	OpenDialog->FilterIndex=8;
	if (!OpenDialog->Execute()) return;
//	GoogleEarthFile->Text=OpenDialog->FileName;
}
//---------------------------------------------------------------------------

