// InfoDialog.cpp : implementation file
//

#include "stdafx.h"
#include "VisView.h"
#include "InfoDialog.h"


// CInfoDialog dialog

IMPLEMENT_DYNAMIC(CInfoDialog, CDialog)

CInfoDialog::CInfoDialog(CWnd* pParent /*=NULL*/)
	: CDialog(CInfoDialog::IDD, pParent)
{
    m_pInterface = NULL;
}

CInfoDialog::~CInfoDialog()
{
}

void CInfoDialog::SetResult(CLYInterface* pInterface)
{
    m_pInterface = pInterface;
}

CString CInfoDialog::GetGeneralInfo()
{
    CString result;
    char* buf;
    unsigned size;

    size = m_pInterface->GetInfo(buf);
    result = CharStringToWindowText(buf, size);
    free(buf);
    return CString(_T("Layer size info\r\n---------\r\n")) + result + _T("\r\n");
}

CString CInfoDialog::GetEdgeInfoText()
{
    CString result;
    char* buf;
    unsigned size;

    size = m_pInterface->GetEdgeInfo(buf);
    result = CharStringToWindowText(buf, size);
    free(buf);
    return CString(_T("Edge info\r\n---------\r\n")) + result;
}

CString CInfoDialog::GetNodesInfoText()
{
    CString str;

    str.Format(_T("\r\nNodes info\r\n---------------\r\nAll nodes: %d\r\nHypo nodes: %d\r\n"), 
        m_pInterface->GetAllNodesCount(), m_pInterface->GetHypoNodesCount());

    return str;
}

CString CInfoDialog::CharStringToWindowText(char* buf, unsigned bufSize)
{
    CString result(buf);

    result.Replace(_T("\n"), _T("\r\n"));
    return result;
 /*   wchar_t destBuf = malloc(bufSize, sizeof(wchar_t));

    mbstowcs(destBuf, buf, bufSize);

    CString result(destBuf)

   wcstombs(
   char *mbstr,
   const wchar_t *wcstr,
   size_t count 
   */

}


BOOL CInfoDialog::OnInitDialog()
{
    CDialog::OnInitDialog();

    CEdit* edit = (CEdit*)GetDlgItem(IDC_EDIT_INFO);
    edit->SetWindowText(GetGeneralInfo() + GetEdgeInfoText() + GetNodesInfoText());

    return TRUE;
}

void CInfoDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(CInfoDialog, CDialog)
END_MESSAGE_MAP()


// CInfoDialog message handlers
