/*------------------------------------------------------------------------------
* tcpoptdlg.h
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

#ifndef tcpoptdlgH
#define tcpoptdlgH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>

#define MAXHIST		10

//---------------------------------------------------------------------------
class TTcpOptDialog : public TForm
{
__published:
	TButton *BtnCancel;
	TButton *BtnOk;
	TEdit *Port;
	TLabel *LabelAddr;
	TLabel *LabelPort;
	TEdit *User;
	TEdit *Passwd;
	TLabel *LabelUser;
	TLabel *LabelPasswd;
	TLabel *LabelMntPnt;
	TEdit *Str;
	TLabel *LabelStr;
	TComboBox *Addr;
	TComboBox *MntPnt;
	TButton *BtnNtrip;
	void __fastcall FormShow(TObject *Sender);
	void __fastcall BtnOkClick(TObject *Sender);
	void __fastcall BtnNtripClick(TObject *Sender);
private:
	void __fastcall AddHist(TComboBox *list, AnsiString *hist);
	int __fastcall ExecCmd(AnsiString cmd, int show);
public:
	int Opt;
	AnsiString Path,History[MAXHIST],MntpHist[MAXHIST];
	__fastcall TTcpOptDialog(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TTcpOptDialog *TcpOptDialog;
//---------------------------------------------------------------------------
#endif
