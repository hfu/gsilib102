/*------------------------------------------------------------------------------
* aboutdlg.cpp
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

#include "rtklib.h"
#include "aboutdlg.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TAboutDialog *AboutDialog;
//---------------------------------------------------------------------------
__fastcall TAboutDialog::TAboutDialog(TComponent* Owner)
	: TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TAboutDialog::FormShow(TObject *Sender)
{
	TImage *icon[]={Icon1,Icon2,Icon3,Icon4,Icon5,Icon6,Icon7,Icon8};
	AnsiString s;
	if (IconIndex>0) icon[IconIndex-1]->Visible=true;
	LabelAbout->Caption=About;
	LabelVer->Caption=s.sprintf("with GSILIB ver.%s",VER_RTKLIB);
	LabelCopyright->Caption=COPYRIGHT_GSILIB;
}
