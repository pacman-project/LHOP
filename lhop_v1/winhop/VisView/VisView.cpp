// VisView.cpp : Defines the class behaviors for the application.
//

#include "stdafx.h"
#include "VisView.h"
#include "VisViewDlg.h"
#include "layers\initialization.h"
//#include "sinstance.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CVisViewApp

BEGIN_MESSAGE_MAP(CVisViewApp, CWinApp)
	ON_COMMAND(ID_HELP, &CWinApp::OnHelp)
END_MESSAGE_MAP()


// CVisViewApp construction

CVisViewApp::CVisViewApp()
{
	// TODO: add construction code here,
	// Place all significant initialization in InitInstance
}


// The one and only CVisViewApp object

CVisViewApp theApp;


// CVisViewApp initialization

BOOL CVisViewApp::InitInstance()
{
    //Check for the previous instance as soon as possible
    //CInstanceChecker instanceChecker(_T("VisView application")); 

    //instanceChecker.ActivateChecker();
    //if (instanceChecker.PreviousInstanceRunning()) {
    //    CCommandLineInfo cmdInfo;
    //    ParseCommandLine(cmdInfo);

    //    instanceChecker.ActivatePreviousInstance(cmdInfo.m_strFileName);
    //    return FALSE;
    //}

    // Init "recognition" specific stuff.
    init_atoms();
    init_streaming();

    // InitCommonControlsEx() is required on Windows XP if an application
	// manifest specifies use of ComCtl32.dll version 6 or later to enable
	// visual styles.  Otherwise, any window creation will fail.
	INITCOMMONCONTROLSEX InitCtrls;
	InitCtrls.dwSize = sizeof(InitCtrls);
	// Set this to include all the common control classes you want to use
	// in your application.
	InitCtrls.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&InitCtrls);

	CWinApp::InitInstance();

    // Standard initialization
	// If you are not using these features and wish to reduce the size
	// of your final executable, you should remove from the following
	// the specific initialization routines you do not need
	// Change the registry key under which our settings are stored
	// TODO: You should modify this string to be something appropriate
	// such as the name of your company or organization
	SetRegistryKey(_T("Local AppWizard-Generated Applications"));

	CVisViewDlg dlg;
	m_pMainWnd = &dlg;

    // Track this instance
    //instanceChecker.TrackFirstInstanceRunning();

    dlg.DoModal();

    return FALSE;
}
