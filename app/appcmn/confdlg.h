/*------------------------------------------------------------------------------
* confdlg.h
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

#ifndef confdlgH
#define confdlgH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
//---------------------------------------------------------------------------
class TConfDialog : public TForm
{
__published:
	TButton *BtnOverwrite;
	TButton *BtnCancel;
	TLabel *Label1;
	TLabel *Label2;
private:
public:
	__fastcall TConfDialog(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TConfDialog *ConfDialog;
//---------------------------------------------------------------------------
#endif
