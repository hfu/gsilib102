/*------------------------------------------------------------------------------
* gsipost_gui.cpp
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
//---------------------------------------------------------------------------
USEFORM("..\appcmn\vieweropt.cpp", ViewerOptDialog);
USEFORM("..\appcmn\viewer.cpp", TextViewer);
USEFORM("postmain.cpp", MainForm);
USEFORM("kmzconv.cpp", ConvDialog);
USEFORM("..\appcmn\keydlg.cpp", KeyDialog);
USEFORM("..\appcmn\confdlg.cpp", ConfDialog);
USEFORM("..\appcmn\aboutdlg.cpp", AboutDialog);
USEFORM("..\appcmn\timedlg.cpp", TimeDialog);
USEFORM("..\appcmn\refdlg.cpp", RefDialog);
USEFORM("..\appcmn\maskoptdlg.cpp", MaskOptDialog);
USEFORM("postopt.cpp", OptDialog);
//---------------------------------------------------------------------------
WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
	try
	{
		Application->Initialize();
		Application->Title = "GSIPOST";
		Application->CreateForm(__classid(TMainForm), &MainForm);
		Application->CreateForm(__classid(TOptDialog), &OptDialog);
		Application->CreateForm(__classid(TConvDialog), &ConvDialog);
		Application->CreateForm(__classid(TOptDialog), &OptDialog);
		Application->CreateForm(__classid(TTextViewer), &TextViewer);
		Application->CreateForm(__classid(TViewerOptDialog), &ViewerOptDialog);
		Application->CreateForm(__classid(TRefDialog), &RefDialog);
		Application->CreateForm(__classid(TTimeDialog), &TimeDialog);
		Application->CreateForm(__classid(TConfDialog), &ConfDialog);
		Application->CreateForm(__classid(TAboutDialog), &AboutDialog);
		Application->CreateForm(__classid(TKeyDialog), &KeyDialog);
		Application->CreateForm(__classid(TMaskOptDialog), &MaskOptDialog);
		Application->Run();
	}
	catch (Exception &exception)
	{
		Application->ShowException(&exception);
	}
	catch (...)
	{
		try
		{
			throw Exception("");
		}
		catch (Exception &exception)
		{
			Application->ShowException(&exception);
		}
	}
	return 0;
}
//---------------------------------------------------------------------------
