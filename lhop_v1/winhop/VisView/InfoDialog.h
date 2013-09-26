#pragma once

#include "../interface/interface.h"

// CInfoDialog dialog

class CInfoDialog : public CDialog
{
	DECLARE_DYNAMIC(CInfoDialog)

protected:
    CLYInterface* m_pInterface;

    CString CharStringToWindowText(char* buf, unsigned bufSize);
public:
	CInfoDialog(CWnd* pParent = NULL);   // standard constructor
	virtual ~CInfoDialog();

// Custom methods
    void SetResult(CLYInterface* pInterface);
    CString GetGeneralInfo();
    CString GetEdgeInfoText();
    CString GetNodesInfoText();

// Dialog Data
	enum { IDD = IDD_INFO_DIALOG };

protected:
    virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
};
