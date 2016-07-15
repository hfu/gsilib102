/*------------------------------------------------------------------------------
* serioptdlg.h
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

#ifndef serioptdlgH
#define serioptdlgH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
//---------------------------------------------------------------------------
class TSerialOptDialog : public TForm
{
__published:
	TButton *BtnOk;
	TButton *BtnCancel;
	TComboBox *BitRate;
	TLabel *Label1;
	TLabel *Label3;
	TComboBox *Port;
	TLabel *Label2;
	TComboBox *ByteSize;
	TLabel *Label4;
	TComboBox *Parity;
	TLabel *Label5;
	TComboBox *StopBits;
	TLabel *Label8;
	TComboBox *FlowCtr;
	TButton *BtnCmd;
	void __fastcall FormShow(TObject *Sender);
	void __fastcall BtnOkClick(TObject *Sender);
	void __fastcall BtnCmdClick(TObject *Sender);
private:
	void __fastcall UpdatePortList(void);
public:
	AnsiString Path,Cmds[2];
	int Opt,CmdEna[2];
	__fastcall TSerialOptDialog(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TSerialOptDialog *SerialOptDialog;
//---------------------------------------------------------------------------
#endif
