/*------------------------------------------------------------------------------
* gsiplot.cpp
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
USEFORM("..\appcmn\timedlg.cpp", TimeDialog);
USEFORM("..\appcmn\tspandlg.cpp", SpanDialog);
USEFORM("..\appcmn\serioptdlg.cpp", SerialOptDialog);
USEFORM("..\appcmn\tcpoptdlg.cpp", TcpOptDialog);
USEFORM("conndlg.cpp", ConnectDialog);
USEFORM("fileseldlg.cpp", FileSelDialog);
USEFORM("..\appcmn\viewer.cpp", TextViewer);
USEFORM("..\appcmn\vieweropt.cpp", ViewerOptDialog);
USEFORM("..\appcmn\refdlg.cpp", RefDialog);
USEFORM("..\appcmn\confdlg.cpp", ConfDialog);
USEFORM("..\appcmn\console.cpp", Console);
USEFORM("..\appcmn\aboutdlg.cpp", AboutDialog);
USEFORM("..\appcmn\cmdoptdlg.cpp", CmdOptDialog);
USEFORM("..\appcmn\keydlg.cpp", KeyDialog);
USEFORM("..\appcmn\fileoptdlg.cpp", FileOptDialog);
USEFORM("..\appcmn\ftpoptdlg.cpp", FtpOptDialog);
USEFORM("pntdlg.cpp", PntDialog);
USEFORM("satdlg.cpp", SatDialog);
USEFORM("plotmain.cpp", Plot);
USEFORM("plotopt.cpp", PlotOptDialog);
USEFORM("geview.cpp", GoogleEarthView);
USEFORM("gmview.cpp", GoogleMapView);
USEFORM("mapdlg.cpp", MapAreaDialog);
//---------------------------------------------------------------------------
WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
	try
	{
		Application->Initialize();
		Application->Title = "GSIPLOT";
		Application->CreateForm(__classid(TPlot), &Plot);
		Application->CreateForm(__classid(TPlotOptDialog), &PlotOptDialog);
		Application->CreateForm(__classid(TSatDialog), &SatDialog);
		Application->CreateForm(__classid(TRefDialog), &RefDialog);
		Application->CreateForm(__classid(TAboutDialog), &AboutDialog);
		Application->CreateForm(__classid(TSpanDialog), &SpanDialog);
		Application->CreateForm(__classid(TTimeDialog), &TimeDialog);
		Application->CreateForm(__classid(TConnectDialog), &ConnectDialog);
		Application->CreateForm(__classid(TSerialOptDialog), &SerialOptDialog);
		Application->CreateForm(__classid(TTcpOptDialog), &TcpOptDialog);
		Application->CreateForm(__classid(TCmdOptDialog), &CmdOptDialog);
		Application->CreateForm(__classid(TFileOptDialog), &FileOptDialog);
		Application->CreateForm(__classid(TKeyDialog), &KeyDialog);
		Application->CreateForm(__classid(TTextViewer), &TextViewer);
		Application->CreateForm(__classid(TViewerOptDialog), &ViewerOptDialog);
		Application->CreateForm(__classid(TPntDialog), &PntDialog);
		Application->CreateForm(__classid(TMapAreaDialog), &MapAreaDialog);
		Application->CreateForm(__classid(TConfDialog), &ConfDialog);
		Application->CreateForm(__classid(TFileSelDialog), &FileSelDialog);
		Application->CreateForm(__classid(TGoogleEarthView), &GoogleEarthView);
		Application->CreateForm(__classid(TFtpOptDialog), &FtpOptDialog);
		Application->CreateForm(__classid(TConsole), &Console);
		Application->CreateForm(__classid(TGoogleMapView), &GoogleMapView);
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

