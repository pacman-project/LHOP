#pragma once

#include "VisViewDlg.h"
////#include "afxeditbrowsectrl.h"

// VideoDialog dialog

class VideoDialog : public CDialog
{
	DECLARE_DYNAMIC(VideoDialog)

public:
	VideoDialog(CD3DWnd* pD3DWindow, CWnd* pParent = NULL);   // standard constructor
	virtual ~VideoDialog();

// Dialog Data
	enum { IDD = IDD_VIDEO_DIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	DECLARE_MESSAGE_MAP()
public:
	virtual BOOL OnInitDialog();
	CD3DWnd* m_pD3DWindow;
	////CMFCEditBrowseCtrl m_OutputEdit;
	////CMFCEditBrowseCtrl m_InputEdit;
protected:

public:
	afx_msg void OnEnChangeInputedit();
};
