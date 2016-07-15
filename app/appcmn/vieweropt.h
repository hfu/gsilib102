/*------------------------------------------------------------------------------
* vieweropt.h
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

#ifndef vieweroptH
#define vieweroptH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <Dialogs.hpp>
//---------------------------------------------------------------------------
class TViewerOptDialog : public TForm
{
__published:
	TButton *BtnOk;
	TButton *BtnCancel;
	TLabel *Label6;
	TPanel *Color1;
	TButton *BtnColor1;
	TLabel *Label1;
	TPanel *Color2;
	TButton *BtnColor2;
	TLabel *Label15;
	TLabel *FontLabel;
	TButton *BtnFont;
	TColorDialog *ColorDialog;
	TFontDialog *FontDialog;
	void __fastcall BtnColor1Click(TObject *Sender);
	void __fastcall BtnColor2Click(TObject *Sender);
	void __fastcall BtnFontClick(TObject *Sender);
	void __fastcall BtnOkClick(TObject *Sender);
	void __fastcall FormShow(TObject *Sender);
private:
public:
	__fastcall TViewerOptDialog(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TViewerOptDialog *ViewerOptDialog;
//---------------------------------------------------------------------------
#endif
