// VideoDialog.cpp : implementation file
//

#include "stdafx.h"
#include "VisView.h"
#include "VideoDialog.h"
#include "VisViewDlg.h"
////#include "afxeditbrowsectrl.h"


// VideoDialog dialog

IMPLEMENT_DYNAMIC(VideoDialog, CDialog)

VideoDialog::VideoDialog(CD3DWnd* pD3DWindow, CWnd* pParent /*=NULL*/)
	: CDialog(VideoDialog::IDD, pParent)
	, m_pD3DWindow(NULL)
{
	m_pD3DWindow = pD3DWindow;
}

VideoDialog::~VideoDialog()
{
}

void VideoDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//DDX_Control( pDX, IDC_INPUTEDIT, m_InputEdit);
	//DDX_Control( pDX, IDC_OUTPUTEDIT, m_OutputEdit);
	
}


BEGIN_MESSAGE_MAP(VideoDialog, CDialog)
END_MESSAGE_MAP()


// VideoDialog message handlers

BOOL VideoDialog::OnInitDialog()
{
	CDialog::OnInitDialog();
	// TODO:  Add extra initialization here
	////m_InputEdit.EnableFolderBrowseButton();
	////m_OutputEdit.EnableFolderBrowseButton();


	return TRUE;  // return TRUE unless you set the focus to a control
	// EXCEPTION: OCX Property Pages should return FALSE
}
