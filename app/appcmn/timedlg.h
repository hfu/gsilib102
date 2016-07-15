/*------------------------------------------------------------------------------
* timedlg.h
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

#ifndef timedlgH
#define timedlgH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include "rtklib.h"
//---------------------------------------------------------------------------
class TTimeDialog : public TForm
{
__published:
	TButton *BtnOk;
	TLabel *Message;
	void __fastcall FormShow(TObject *Sender);
private:
public:
	gtime_t Time;
	__fastcall TTimeDialog(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TTimeDialog *TimeDialog;
//---------------------------------------------------------------------------
#endif
