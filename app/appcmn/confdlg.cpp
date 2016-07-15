/*------------------------------------------------------------------------------
* confdlg.cpp
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

#include "confdlg.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TConfDialog *ConfDialog;
//---------------------------------------------------------------------------
__fastcall TConfDialog::TConfDialog(TComponent* Owner)
	: TForm(Owner)
{
}
//---------------------------------------------------------------------------
