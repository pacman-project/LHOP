// VisViewDlg.cpp : implementation file
//

#include "stdafx.h"
//#include "vld.h"
#include "VisView.h"
#include "VisViewDlg.h"
#include "InfoDialog.h"
//#include "VideoDialog.h"
#include <atlstr.h>
#include <atlpath.h>
#include "../interface/interface.h"

#define _CRTDBG_MAP_ALLOC
#define VLD_SELF_TEST
#define VLD_MAX_DATA_DUMP 0


#include <stdlib.h>
#include <crtdbg.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	enum { IDD = IDD_ABOUTBOX };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()


// CVisViewDlg dialog



CVisViewDlg::CVisViewDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CVisViewDlg::IDD, pParent)
{
	//_CrtDumpMemoryLeaks();
	m_ctrlDown = false;
	m_makingVideo = false;
    m_initialized = false;
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CVisViewDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_PARTLIST, m_partList);
	DDX_Control(pDX, IDC_COLORSEL, m_ODColors);
	DDX_Control(pDX, IDC_VIDEO, m_pMakeVideoButton);
	DDX_Control(pDX, IDC_CHECKSHOWLINES, m_showLinesCB);
	DDX_Control(pDX, IDC_CHECKSHOWCIRCLES, m_showCircles);
	DDX_Control(pDX, IDC_CHECKSHOWEll, m_showEll);
	DDX_Control(pDX, IDC_CHECKINHIB, m_inhibCB);
	DDX_Control(pDX, IDC_CHECKPICKPARTS, m_pickPartsCB);
	DDX_Control(pDX, IDC_CHECKHIDETEXTURE, m_hideTexture);
	
}

BEGIN_MESSAGE_MAP(CVisViewDlg, CDialog)
    ON_WM_DROPFILES()
    ON_WM_MOUSEWHEEL()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_WM_SYSCOMMAND()
    ON_WM_VSCROLL()
    ON_CONTROL(LBN_MOUSEITEMCHANGE, IDC_PARTLIST, &CVisViewDlg::OnLbnMouseItemChangePartlist)
	//}}AFX_MSG_MAP
	ON_BN_CLICKED(IDC_BUTTON_OPEN, &CVisViewDlg::OnBnClickedButtonOpen)
    ON_LBN_SELCHANGE(IDC_PARTLIST, &CVisViewDlg::OnLbnSelchangePartlist)
    ON_LBN_SELCHANGE(IDC_LAYERLIST, &CVisViewDlg::OnLbnSelchangeLayerlist)
    ON_BN_CLICKED(IDC_BUTTONALLPARTS, &CVisViewDlg::OnBnClickedButtonallparts)
    ON_BN_CLICKED(IDC_BUTTONNOPARTS, &CVisViewDlg::OnBnClickedButtonnoparts)
    ON_CBN_SELENDOK(IDC_COMBOFILE, &CVisViewDlg::OnCbnSelendokCombofile)
    ON_BN_CLICKED(IDC_BUTTON_OPENLIBRARY, &CVisViewDlg::OnBnClickedButtonOpenlibrary)
    ON_CBN_SELENDOK(IDC_COMBOLIBRARY, &CVisViewDlg::OnCbnSelendokCombolibrary)
    ON_BN_CLICKED(IDC_CHECKFIRSTONLY, &CVisViewDlg::OnBnClickedCheckfirstonly)
    ON_EN_UPDATE(IDC_EDITTHRESHOLD, &CVisViewDlg::OnEnUpdateEditthreshold)
    ON_WM_COPYDATA()
    ON_BN_CLICKED(IDC_CHECKSHOWALL, &CVisViewDlg::OnBnClickedCheckshowall)
    ON_BN_CLICKED(IDC_BUTTONBWD, &CVisViewDlg::OnBnClickedButtonbwd)
    ON_BN_CLICKED(IDC_BUTTONFWD, &CVisViewDlg::OnBnClickedButtonfwd)
    ON_BN_CLICKED(IDC_CHECKRECONSTRUCTION, &CVisViewDlg::OnBnClickedCheckreconstruction)
	ON_BN_CLICKED(IDC_SAVECAM, &CVisViewDlg::OnBnClickedSavecam)
	ON_BN_CLICKED(IDC_CHECKTREE, &CVisViewDlg::OnBnClickedChecktree)
	ON_BN_CLICKED(IDC_LOADCAM, &CVisViewDlg::OnBnClickedLoadcam)
	ON_BN_CLICKED(IDC_SCREENSHOT, &CVisViewDlg::OnBnClickedScreenshot)
	ON_WM_KEYDOWN()
	ON_WM_KEYUP()
	ON_CBN_SELCHANGE(IDC_COLORSEL, &CVisViewDlg::OnCbnSelchangeColorsel)
	ON_BN_CLICKED(IDC_RESETCOLOR, &CVisViewDlg::OnBnClickedResetcolor)
	ON_BN_CLICKED(IDC_CUSTOMCOL, &CVisViewDlg::OnBnClickedCustomcol)
	ON_BN_CLICKED(IDC_VIDEO, &CVisViewDlg::OnBnClickedVideo)
//	ON_CBN_SELCHANGE(IDC_COMBOFILE, &CVisViewDlg::OnCbnSelchangeCombofile)
ON_BN_CLICKED(IDC_BUTTONOPENNEXT, &CVisViewDlg::OnBnClickedButtonopennext)
	ON_BN_CLICKED(IDOK, &CVisViewDlg::OnBnClickedOk)
	ON_MESSAGE(WM_USER_SCREENSHOT_DONE, &CVisViewDlg::OnScreenshotDone)
    //ON_BN_CLICKED(IDC_CHECKSHOWHYPONODES, &CVisViewDlg::OnBnClickedCheckshowhyponodes)
	ON_BN_CLICKED(IDC_CHECKSHOWLINES, &CVisViewDlg::OnBnClickedCheckshowlines)
	ON_BN_CLICKED(IDC_CHECKSHOWCIRCLES, &CVisViewDlg::OnBnClickedCheckshowcircles)
	ON_BN_CLICKED(IDC_CHECKSHOWEll, &CVisViewDlg::OnBnClickedCheckshowell)
	ON_BN_CLICKED(IDC_CHECKINHIB, &CVisViewDlg::OnBnClickedCheckinhib)
	ON_WM_ERASEBKGND()
	ON_WM_SETFOCUS()
//  ON_BN_CLICKED(IDC_CHECKHYPO, &CVisViewDlg::OnBnClickedCheckhypo)
//  ON_EN_UPDATE(IDC_EDITTHRESHOLD2, &CVisViewDlg::OnEnUpdateEditthreshold2)
ON_EN_UPDATE(IDC_EDITTHRESHOLD2, &CVisViewDlg::OnEnUpdateEditthreshold2)
//ON_WM_SHOWWINDOW()
ON_EN_UPDATE(IDC_EDITTHRESHOLD3, &CVisViewDlg::OnEnUpdateEditthreshold3)
ON_BN_CLICKED(IDC_CHECKPICKPARTS, &CVisViewDlg::OnBnClickedCheckpickparts)
ON_MESSAGE(WM_GEOMETRY_INIT_DONE,  &CVisViewDlg::OnGeometryInitDone)
ON_BN_CLICKED(IDC_CHECKHIDETEXTURE, &CVisViewDlg::OnBnClickedCheckhidetexture)
ON_BN_CLICKED(IDC_BUTTON_EDGEINFO, &CVisViewDlg::OnBnClickedButtonEdgeinfo)
ON_EN_UPDATE(IDC_EDITTHRESHOLD4, &CVisViewDlg::OnEnUpdateEditthreshold4)
//ON_WM_SHOWWINDOW()
ON_WM_ACTIVATE()
ON_BN_CLICKED(IDC_BUTTON2, &CVisViewDlg::OnBnClickedButton2)
ON_BN_CLICKED(IDC_BUTTON_LOADSELECT, &CVisViewDlg::OnBnClickedButtonLoadselect)
END_MESSAGE_MAP()


// CVisViewDlg message handlers

BOOL CVisViewDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// Added

    // D3DWindow
	CRect rect;
	GetDlgItem(IDC_STATIC_D3D)->SetWindowPos(&wndTop, 0, 0, 800, 600, SWP_NOMOVE);
	GetDlgItem(IDC_STATIC_D3D)->GetWindowRect(rect);
    ScreenToClient(rect);
	m_D3DWindow.D3DCreate(rect, this);

    // Controls
    m_pLayerList = (CListBox*)GetDlgItem(IDC_LAYERLIST);
    m_pPartStatic = (CStatic*)GetDlgItem(IDC_STATICPART);
    m_pInfoStatic = (CStatic*)GetDlgItem(IDC_STATICINFO);
    m_pFileCombo = (CComboBox*)GetDlgItem(IDC_COMBOFILE);
    m_pLibraryCombo = (CComboBox*)GetDlgItem(IDC_COMBOLIBRARY);
    m_pFirstOnlyCB = (CButton*)GetDlgItem(IDC_CHECKFIRSTONLY);
	//m_pHypoThresholdCB = (CButton*)GetDlgItem(IDC_CHECKHYPO);
    m_pReconstructionCB = (CButton*)GetDlgItem(IDC_CHECKRECONSTRUCTION);
    m_pShowHyponodesCB = (CButton*)GetDlgItem(IDC_CHECKSHOWHYPONODES);
	m_pTreeCB = (CButton*)GetDlgItem(IDC_CHECKTREE);
    m_pThresholdEdit = (CEdit*)GetDlgItem(IDC_EDITTHRESHOLD);
    m_pThresholdEdit2 = (CEdit*)GetDlgItem(IDC_EDITTHRESHOLD2);
    m_pThresholdEdit3 = (CEdit*)GetDlgItem(IDC_EDITTHRESHOLD3);
    m_pThresholdEdit4 = (CEdit*)GetDlgItem(IDC_EDITTHRESHOLD4);
    m_pThresholdSpin = (CSpinButtonCtrl*)GetDlgItem(IDC_SPINTHRESHOLD);
    m_pThresholdSpin2 = (CSpinButtonCtrl*)GetDlgItem(IDC_SPINTHRESHOLD2);
    m_pThresholdSpin3 = (CSpinButtonCtrl*)GetDlgItem(IDC_SPINTHRESHOLD3);
    m_pThresholdSpin4 = (CSpinButtonCtrl*)GetDlgItem(IDC_SPINTHRESHOLD4);
	//m_pHypoThresholdEdit = (CEdit*)GetDlgItem(IDC_EDITHYPOTHRESHOLD);
	//m_pHypoThresholdSpin = (CSpinButtonCtrl*)GetDlgItem(IDC_SPINHYPOTHRESHOLD);
    m_pShowAllCB = (CButton*)GetDlgItem(IDC_CHECKSHOWALL);
    m_pFwd = (CButton*)GetDlgItem(IDC_BUTTONFWD);
    m_pBwd = (CButton*)GetDlgItem(IDC_BUTTONBWD);	

    // Init some controls
    m_pThresholdSpin->SetBuddy(m_pThresholdEdit);
    m_pThresholdSpin->SetRange(0, 20);
    m_pThresholdSpin2->SetBuddy(m_pThresholdEdit2);
    m_pThresholdSpin2->SetRange(0, 20);
    m_pThresholdSpin3->SetBuddy(m_pThresholdEdit3);
    m_pThresholdSpin3->SetRange(0, 20);
    m_pThresholdSpin4->SetBuddy(m_pThresholdEdit4);
    m_pThresholdSpin4->SetRange(0, 20);
	//m_pHypoThresholdSpin->SetBuddy(m_pHypoThresholdEdit);
	//m_pHypoThresholdSpin->SetRange(0, 20);

    UpdateThresholdControls();
	int show = m_showLinesCB.GetCheck();
	m_D3DWindow.SetProperty(CD3DWnd::DisplaySetting::SHOW_LINES, show);
    m_pReconstructionCB->SetCheck(
        m_D3DWindow.GetInterface()->GetShowReconstruction() ? BST_CHECKED : BST_UNCHECKED
    );
    //m_pShowHyponodesCB->SetCheck(
    //    m_D3DWindow.GetInterface()->GetShowHyponodes() ? BST_CHECKED : BST_UNCHECKED
    //);
	m_pTreeCB->SetCheck(
		m_D3DWindow.GetInterface()->GetShowTree() ? BST_CHECKED : BST_UNCHECKED
	);
    if (m_D3DWindow.GetInterface()->GetDisplayFirst()) {
        m_pFirstOnlyCB->SetCheck(BST_CHECKED);
        m_pThresholdEdit->EnableWindow(FALSE);
        m_pThresholdSpin->EnableWindow(FALSE);
        m_pThresholdEdit2->EnableWindow(FALSE);
        m_pThresholdSpin2->EnableWindow(FALSE);
        m_pThresholdEdit3->EnableWindow(FALSE);
        m_pThresholdSpin3->EnableWindow(FALSE);
        m_pThresholdEdit4->EnableWindow(FALSE);
        m_pThresholdSpin4->EnableWindow(FALSE);
    } else {
        m_pFirstOnlyCB->SetCheck(BST_UNCHECKED);
        m_pThresholdEdit->EnableWindow(TRUE);
        m_pThresholdSpin->EnableWindow(TRUE);
        m_pThresholdEdit2->EnableWindow(TRUE);
        m_pThresholdSpin2->EnableWindow(TRUE);
        m_pThresholdEdit3->EnableWindow(TRUE);
        m_pThresholdSpin3->EnableWindow(TRUE);
        m_pThresholdEdit4->EnableWindow(TRUE);
        m_pThresholdSpin4->EnableWindow(TRUE);
    }
	//if (m_D3DWindow.GetInterface()->GetHypoThreshold()) {
	//  m_pHypoThresholdCB->SetCheck(BST_CHECKED);
	//	m_pHypoThresholdEdit->EnableWindow(TRUE);
    //  m_pHypoThresholdSpin->EnableWindow(TRUE);
    //} else {
	//	m_pHypoThresholdCB->SetCheck(BST_UNCHECKED);
	//	m_pHypoThresholdEdit->EnableWindow(FALSE);
    //  m_pHypoThresholdSpin->EnableWindow(FALSE);
    //}
    if (m_D3DWindow.GetInterface()->GetShowAll()) {
        m_pShowAllCB->SetCheck(BST_CHECKED);
        m_pFwd->EnableWindow(FALSE);
        m_pBwd->EnableWindow(FALSE);
    } else {
        m_pShowAllCB->SetCheck(BST_UNCHECKED);
        m_pFwd->EnableWindow(TRUE);
        m_pBwd->EnableWindow(TRUE);
    }
    //m_pInfoStatic->ShowWindow(SW_HIDE);

    // Parse command line arguments
    CCommandLineInfo cmdInfo;

    AfxGetApp()->ParseCommandLine(cmdInfo);
    OpenLYFile(cmdInfo.m_strFileName);

	m_ODColors.AddStdColors();
	//color dialog
	ZeroMemory(&m_chooseColor, sizeof(m_chooseColor));
	m_chooseColor.lStructSize = sizeof(m_chooseColor);
	//m_chooseColor.hwndOwner = ((CButton*)GetDlgItem(IDC_CUSTOMCOL))->m_hWnd;
	m_chooseColor.lpCustColors = (LPDWORD) m_acrCustClr;
	m_chooseColor.rgbResult = 0;
	m_chooseColor.Flags = CC_FULLOPEN | CC_RGBINIT;
	CustomColorsFromReg();
    m_initialized = true;
    return TRUE;  // return TRUE  unless you set the focus to a control
}


void CVisViewDlg::UpdateThresholdControls()
{
    CString str;
    double thresh = m_D3DWindow.GetInterface()->GetThreshold();
    double thresh2 = m_D3DWindow.GetInterface()->GetThreshold2();
    double thresh3 = m_D3DWindow.GetInterface()->GetThreshold3();
    double thresh4 = m_D3DWindow.GetInterface()->GetThreshold4();
	//double hypothresh = m_D3DWindow.GetInterface()->GetHypoThresholdVal();

    str.Format(_T("%.2f"), thresh);
    m_pThresholdEdit->SetWindowText(str);
    m_pThresholdSpin->SetPos32((int)(20*thresh));
    str.Format(_T("%.2f"), thresh2);
    m_pThresholdEdit2->SetWindowText(str);
    m_pThresholdSpin2->SetPos32((int)(20*thresh2));
    str.Format(_T("%.2f"), thresh3);
    m_pThresholdEdit3->SetWindowText(str);
    m_pThresholdSpin3->SetPos32((int)(20*thresh3));
    str.Format(_T("%.2f"), thresh4);
    m_pThresholdEdit4->SetWindowText(str);
    m_pThresholdSpin4->SetPos32((int)(4*thresh4));
	//str.Format(_T("%.2f"), hypothresh);
	//m_pHypoThresholdEdit->SetWindowText(str);
    //m_pHypoThresholdSpin->SetPos32((int)(20*hypothresh));

}

void CVisViewDlg::ResetPartList()
{
    CString str;
    int count = m_D3DWindow.GetInterface()->PartCount();

    // Part list
    m_partList.ResetContent();
    for (int i = 0; i < count; ++i) {
        str.Format(_T("part %d"), m_D3DWindow.GetInterface()->GetPart(i));
        m_partList.AddString(str);
    }
    for (int i = 0; i < count; ++i) {
        m_partList.SetSel(i, m_D3DWindow.GetInterface()->PartDisplayed(i) && m_D3DWindow.GetInterface()->PartVisable(i));
    }

}


void CVisViewDlg::ResetLayerList()
{
    CString str;
    int count = m_D3DWindow.GetInterface()->LayerCount();

    m_pLayerList->ResetContent();
    for (int i = 0; i < count; ++i) {
        str.Format(_T("layer %d"), i);
        m_pLayerList->AddString(str);
    }
    m_pLayerList->SetCurSel(m_D3DWindow.GetInterface()->GetLayer());
}

void CVisViewDlg::OpenLYFile(CString fileName)
{
    CWaitCursor wait;
    int pos;
	//m_D3DWindow.OpenLYFile(fileName);
	m_D3DWindow.GetInterface()->OpenLYFile(fileName);
	//OpenLYFile(fileName);
	m_D3DWindow.SetLyFilePath(fileName);
    //m_D3DWindow.OpenLYFile(fileName);
    if ((pos = m_pFileCombo->FindStringExact(0, fileName)) == CB_ERR) {
        m_pFileCombo->InsertString(0, fileName);
        pos = 0;
    }
    m_pFileCombo->SetCurSel(pos);
    ResetLayerList();
	m_D3DWindow.Reset();
    ResetPartList();
    m_D3DWindow.SetFocus();
    
    wait.Restore();
}

void CVisViewDlg::OpenLibraryFile(CString fileName)
{
    CWaitCursor wait;
    int pos;

    m_D3DWindow.GetInterface()->OpenLibrary(fileName);
    if ((pos = m_pLibraryCombo->FindStringExact(0, fileName)) == CB_ERR) {
        m_pLibraryCombo->InsertString(0, fileName);
        pos = 0;
    }
    m_pLibraryCombo->SetCurSel(pos);
    m_D3DWindow.SetFocus();
    
    wait.Restore();
}

void CVisViewDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CVisViewDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CVisViewDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

void CVisViewDlg::OnBnClickedButtonOpen()
{
    CFileDialog dlg(TRUE, NULL, NULL, 0, _T("Layer Files (*.ly?)|*.ly?|All Files (*.*)|*.*||"));
    if (dlg.DoModal() == IDOK) {
        OpenLYFile(dlg.GetPathName());
        //OpenLYFile(_T("c:\\work\\test\\fer_mugs_puppy_2.ly6"));        
    }
}

void CVisViewDlg::OnDropFiles(HDROP hDropInfo)
{
    CDialog::OnDropFiles(hDropInfo);

    if (DragQueryFile(hDropInfo, (UINT)-1, NULL, 0) > 0) {
        TCHAR file[_MAX_PATH];
        DragQueryFile(hDropInfo, 0, file, _MAX_PATH);

        CString ext = CString(file).Right(4);

        if (ext.CompareNoCase(_T(".plb")) == 0) 
            OpenLibraryFile(file);
        else if (ext.Left(3).CompareNoCase(_T(".ly")) == 0 || ext.Left(3).CompareNoCase(_T("ly1")) == 0) 
            OpenLYFile(file);
    }
    DragFinish(hDropInfo);  
}

void CVisViewDlg::DrawPart(int pIndex)
{
    int part = m_D3DWindow.GetInterface()->GetPart(pIndex);
    int lyr = m_pLayerList->GetCurSel();
    void* bits;
    int w, h;
    
	RECT m_pPartStaticSize;
	m_pPartStatic->GetClientRect(&m_pPartStaticSize);

	int max_width = m_pPartStaticSize.right - m_pPartStaticSize.left;
	int max_height = m_pPartStaticSize.bottom - m_pPartStaticSize.top;

    m_D3DWindow.GetInterface()->GetPartImage(bits, w, h, lyr, part, max_width, max_height);
    if (bits != NULL) {
        HDC dc = ::GetWindowDC(m_pPartStatic->m_hWnd);
        HBITMAP bmp = CreateCompatibleBitmap(dc, w, h);
        UINT lpbiSize = sizeof(BITMAPINFOHEADER) + 3 * sizeof(RGBQUAD);
        LPBITMAPINFO lpbi = (LPBITMAPINFO)(malloc(lpbiSize));
        
        memset(lpbi, 0, lpbiSize);
        lpbi->bmiHeader.biSize = sizeof(BITMAPINFOHEADER); 
        lpbi->bmiHeader.biWidth = w;
        lpbi->bmiHeader.biHeight = -h;
        lpbi->bmiHeader.biPlanes = 1;
        lpbi->bmiHeader.biBitCount = 32;
        lpbi->bmiHeader.biCompression = 3;
        lpbi->bmiHeader.biSizeImage = 4*w*h;
        lpbi->bmiColors[0].rgbRed = 0xff;
        lpbi->bmiColors[1].rgbGreen =  0xff;
        lpbi->bmiColors[2].rgbBlue = 0xff;
        SetDIBits(dc, bmp, 0, h, bits, lpbi, DIB_RGB_COLORS);
        m_pPartStatic->SetBitmap(bmp);
        //CImage im;
        //im.Attach(bmp);
        //im.Save(_T("c:\\\work\\bmp.png"));
        DeleteObject(bmp);
        free(bits);
        free(lpbi);
    }
}

void CVisViewDlg::OnLbnMouseItemChangePartlist()
{
    DrawPart(m_partList.m_mouseItem);
    
    //str.Format(_T("OUT")); 
    //else str.Format(_T("I %d"), m_partList.m_mouseItem);
    //m_pPartStatic->SetWindowText(str);
}

BOOL CVisViewDlg::OnMouseWheel(UINT nFlags, short zDelta, CPoint point)
{
	RECT partRect;
	RECT layerRect;
	RECT d3dRect;	

	CPoint m_partList_point = point;
	m_partList.GetClientRect(&partRect);
	m_partList.ScreenToClient(&m_partList_point);

	CPoint m_pLayerList_point = point;
	m_pLayerList->GetClientRect(&layerRect);
	m_pLayerList->ScreenToClient(&m_pLayerList_point);

	CPoint m_D3DWindow_point = point;
	m_D3DWindow.GetClientRect(&d3dRect);
	m_D3DWindow.ScreenToClient(&m_D3DWindow_point);

	if (partRect.left <= m_partList_point.x && partRect.right >= m_partList_point.x &&
		partRect.top <= m_partList_point.y && partRect.bottom >= m_partList_point.y) {
		const MSG* pMsg = GetCurrentMessage();
		CWnd* oldFocus = GetFocus();
		m_partList.SetFocus();
		m_partList.PostMessage(pMsg->message, pMsg->wParam, pMsg->lParam);
		oldFocus->SetFocus();
	} else if (layerRect.left <= m_pLayerList_point.x && layerRect.right >= m_pLayerList_point.x &&
		layerRect.top <= m_pLayerList_point.y && layerRect.bottom >= m_pLayerList_point.y) {
		const MSG* pMsg = GetCurrentMessage();
		CWnd* oldFocus = GetFocus();
		m_pLayerList->SetFocus();
		m_pLayerList->PostMessage(pMsg->message, pMsg->wParam, pMsg->lParam);
		oldFocus->SetFocus();
	} else if (d3dRect.left <= m_D3DWindow_point.x && d3dRect.right >= m_D3DWindow_point.x &&
		d3dRect.top <= m_D3DWindow_point.y && d3dRect.bottom >= m_D3DWindow_point.y) {
		return m_D3DWindow.OnMouseWheel(nFlags, zDelta, point, m_ctrlDown);
	} else {
		return CWnd::OnMouseWheel(nFlags, zDelta, point);
	}
}

void CVisViewDlg::OnLbnSelchangePartlist()
{
    int count = m_partList.GetCount();

	int num_selected = 0;
    for (int i = 0; i < count; ++i) {
        m_D3DWindow.GetInterface()->DisplayPart(i, m_partList.GetSel(i) > 0);
	}

	

    m_D3DWindow.GetInterface()->ResetIndex();
    m_D3DWindow.ResetPoints();
    m_D3DWindow.SetFocus();
}

void CVisViewDlg::OnLbnSelchangeLayerlist()
{
    int sel = m_pLayerList->GetCurSel();
    m_D3DWindow.GetInterface()->SetLayer(sel);
	m_D3DWindow.GetInterface()->ClearLayer0Mapping();
    m_D3DWindow.ResetPoints();
    ResetPartList();
    m_D3DWindow.SetFocus();
}

void CVisViewDlg::OnBnClickedButtonallparts()
{
    m_D3DWindow.GetInterface()->ToggleParts(true);
    ResetPartList();
    m_D3DWindow.ResetPoints();
}

void CVisViewDlg::OnBnClickedButtonnoparts()
{
    m_D3DWindow.GetInterface()->ToggleParts(false);
    ResetPartList();
    m_D3DWindow.ResetPoints();
}

void CVisViewDlg::OnCbnSelendokCombofile()
{
    int sel = m_pFileCombo->GetCurSel();

    if (sel != CB_ERR) {
        CString str;
        int len = m_pFileCombo->GetLBTextLen(sel);

        m_pFileCombo->GetLBText(sel, str.GetBuffer(len));
        str.ReleaseBuffer();
        OpenLYFile(str);
    }
}

void CVisViewDlg::OnBnClickedButtonOpenlibrary()
{
    //OpenLibraryFile(_T("c:\\work\\test\\olib-archive.plb"));        
    CFileDialog dlg(TRUE, NULL, NULL, 0, _T("Part Library Files (*.plb)|*.plb|All Files (*.*)|*.*||"));
    if (dlg.DoModal() == IDOK) {
        OpenLibraryFile(dlg.GetPathName());
    }
}

void CVisViewDlg::OnCbnSelendokCombolibrary()
{
    int sel = m_pLibraryCombo->GetCurSel();

    if (sel != CB_ERR) {
        CString str;
        int len = m_pLibraryCombo->GetLBTextLen(sel);

        m_pLibraryCombo->GetLBText(sel, str.GetBuffer(len));
        str.ReleaseBuffer();
        OpenLibraryFile(str);
    }
}

void CVisViewDlg::OnBnClickedCheckfirstonly()
{
    if (m_pFirstOnlyCB->GetCheck() == BST_CHECKED && !m_D3DWindow.GetInterface()->GetDisplayFirst()) {
        m_D3DWindow.GetInterface()->SetDisplayFirst(true);
		m_D3DWindow.GetInterface()->ClearLayer0Mapping();
        m_pThresholdEdit->EnableWindow(FALSE);
        m_pThresholdSpin->EnableWindow(FALSE);
        m_pThresholdEdit2->EnableWindow(FALSE);
        m_pThresholdSpin2->EnableWindow(FALSE);
        m_pThresholdEdit3->EnableWindow(FALSE);
        m_pThresholdSpin3->EnableWindow(FALSE);
        m_pThresholdEdit4->EnableWindow(FALSE);
        m_pThresholdSpin4->EnableWindow(FALSE);
        m_D3DWindow.ResetPoints();
    } else if (m_pFirstOnlyCB->GetCheck() == BST_UNCHECKED && m_D3DWindow.GetInterface()->GetDisplayFirst()) {
        m_D3DWindow.GetInterface()->SetDisplayFirst(false);
		m_D3DWindow.GetInterface()->ClearLayer0Mapping();
        m_pThresholdEdit->EnableWindow(TRUE);
        m_pThresholdSpin->EnableWindow(TRUE);
        m_pThresholdEdit2->EnableWindow(TRUE);
        m_pThresholdSpin2->EnableWindow(TRUE);
        m_pThresholdEdit3->EnableWindow(TRUE);
        m_pThresholdSpin3->EnableWindow(TRUE);
        m_pThresholdEdit4->EnableWindow(TRUE);
        m_pThresholdSpin4->EnableWindow(TRUE);
        m_D3DWindow.ResetPoints();
    }
}

void CVisViewDlg::OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar)
{
    if (pScrollBar == GetDlgItem(IDC_SPINTHRESHOLD)) {
        double d = m_pThresholdSpin->GetPos32()/20.0;

        m_D3DWindow.GetInterface()->SetThreshold(d);
		m_D3DWindow.GetInterface()->ClearLayer0Mapping();
        UpdateThresholdControls();
        m_D3DWindow.ResetPoints();
    }
    if (pScrollBar == GetDlgItem(IDC_SPINTHRESHOLD2)) {
        double d = m_pThresholdSpin2->GetPos32()/20.0;

        m_D3DWindow.GetInterface()->SetThreshold2(d);
		m_D3DWindow.GetInterface()->ClearLayer0Mapping();
		m_D3DWindow.GetInterface()->ClearLayer0Mapping();
        UpdateThresholdControls();
        m_D3DWindow.ResetPoints();
    }
    if (pScrollBar == GetDlgItem(IDC_SPINTHRESHOLD3)) {
        double d = m_pThresholdSpin3->GetPos32()/20.0;

        m_D3DWindow.GetInterface()->SetThreshold3(d);
		m_D3DWindow.GetInterface()->ClearLayer0Mapping();
        UpdateThresholdControls();
        m_D3DWindow.ResetPoints();
    }
    if (pScrollBar == GetDlgItem(IDC_SPINTHRESHOLD4)) {
        double d = m_pThresholdSpin4->GetPos32()/4.0;

        m_D3DWindow.GetInterface()->SetThreshold4(d);
		m_D3DWindow.GetInterface()->ClearLayer0Mapping();
        UpdateThresholdControls();
        m_D3DWindow.ResetPoints();
    }
	//if (pScrollBar == GetDlgItem(IDC_SPINHYPOTHRESHOLD)) {
    //       double d = m_pHypoThresholdSpin->GetPos32()/20.0;

    //       m_D3DWindow.GetInterface()->SetHypoThresholdVal(d);
	//	m_D3DWindow.GetInterface()->ClearLayer0Mapping();
    //       UpdateThresholdControls();
    //       m_D3DWindow.ResetPoints();
    //   }
    CDialog::OnVScroll(nSBCode, nPos, pScrollBar);
}

void CVisViewDlg::OnEnUpdateEditthreshold()
{
    CString str;
    double d;
    
    m_pThresholdEdit->GetWindowText(str.GetBuffer(10), 9);
    str.ReleaseBuffer();
    d = wcstod(str, NULL);
    d = max(0.0, d);
    d = 20.0 * min(1.0, d);
    d = ((int)d)/20.0;
    
    if (abs(m_D3DWindow.GetInterface()->GetThreshold() - d) > 0.001) {
        m_D3DWindow.GetInterface()->SetThreshold(d);
		m_D3DWindow.GetInterface()->ClearLayer0Mapping();
        UpdateThresholdControls();
        m_D3DWindow.ResetPoints();
    }
}

void CVisViewDlg::OnBnClickedButtonEdgeinfo()
{
    CInfoDialog dlgInfo;

    dlgInfo.SetResult(m_D3DWindow.GetInterface());
	dlgInfo.DoModal();
}

BOOL CVisViewDlg::OnCopyData(CWnd* pWnd, COPYDATASTRUCT* pCopyDataStruct)
{
    TCHAR* pszCmdLine = (TCHAR*)(pCopyDataStruct->lpData);

    if (pszCmdLine) OpenLYFile(pszCmdLine);
    return TRUE;
}

void CVisViewDlg::DisplayThresholdInfo()
{
    char* buf;

    m_D3DWindow.GetInterface()->GetThresholdInfo(buf);

    CString str(buf);

    m_pInfoStatic->SetWindowText(str);
    free(buf);
}

void CVisViewDlg::OnBnClickedCheckshowall()
{
    if (m_pShowAllCB->GetCheck() == BST_CHECKED && !m_D3DWindow.GetInterface()->GetShowAll()) {
        m_D3DWindow.GetInterface()->SetShowAll(true);
        m_pFwd->EnableWindow(FALSE);
        m_pBwd->EnableWindow(FALSE);
        //m_pInfoStatic->ShowWindow(SW_HIDE);
        m_D3DWindow.ResetPoints();
        m_pInfoStatic->SetWindowText(_T(""));
    } else if (m_pShowAllCB->GetCheck() == BST_UNCHECKED && m_D3DWindow.GetInterface()->GetShowAll()) {
        m_D3DWindow.GetInterface()->SetShowAll(false);
        m_pFwd->EnableWindow(TRUE);
        m_pBwd->EnableWindow(TRUE);
        //m_pInfoStatic->ShowWindow(SW_SHOW);
        m_D3DWindow.ResetPoints();
        DisplayThresholdInfo();    
    }
}

void CVisViewDlg::OnBnClickedButtonbwd()
{
    m_D3DWindow.GetInterface()->DecPartIndex();
    m_D3DWindow.ResetPoints();
    DisplayThresholdInfo();
}

void CVisViewDlg::OnBnClickedButtonfwd()
{
    m_D3DWindow.GetInterface()->IncPartIndex();
    m_D3DWindow.ResetPoints();
    DisplayThresholdInfo();
}

void CVisViewDlg::OnBnClickedCheckreconstruction()
{
    m_D3DWindow.GetInterface()->SetShowReconstruction(m_pReconstructionCB->GetCheck() == BST_CHECKED);
	m_D3DWindow.GetInterface()->ClearLayer0Mapping();
    m_D3DWindow.ResetPoints();
}

void CVisViewDlg::OnBnClickedSavecam()
{
	m_D3DWindow.cameraPosToRegistry();
}


void CVisViewDlg::OnBnClickedChecktree()
{
	 m_D3DWindow.GetInterface()->SetShowTree(m_pTreeCB->GetCheck() == BST_CHECKED);
	 m_D3DWindow.ResetPoints();
}

void CVisViewDlg::OnBnClickedLoadcam()
{
	m_D3DWindow.cameraPosFromRegistry();
}

void CVisViewDlg::OnBnClickedScreenshot()
{
	Screenshot();
}

void CVisViewDlg::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
{
	// TODO: Add your message handler code here and/or call default
	if(nChar == 17 || nChar == 11) m_ctrlDown = true;
	CDialog::OnKeyDown(nChar, nRepCnt, nFlags);
}

void CVisViewDlg::OnKeyUp(UINT nChar, UINT nRepCnt, UINT nFlags)
{
	// TODO: Add your message handler code here and/or call default
	if(nChar == 17 || nChar == 11) m_ctrlDown = false;
	CDialog::OnKeyUp(nChar, nRepCnt, nFlags);
}

void CVisViewDlg::OnCbnSelchangeColorsel()
{
	m_D3DWindow.GetInterface()->SetDisplayedColor(COLORREF2D3DCOLOR(m_ODColors.GetItemData(m_ODColors.GetCurSel())));
	m_D3DWindow.ResetPoints();
}

void CVisViewDlg::OnBnClickedResetcolor()
{
	m_D3DWindow.GetInterface()->SetColorDefault();
	m_D3DWindow.ResetPoints();
}


void CVisViewDlg::OnBnClickedCustomcol()
{
	// Show the fully opened Color dialog with red as the selected color.
	CColorDialog dlg; //(RGB(255, 0, 0), CC_FULLOPEN);

	if(dlg.DoModal() == IDOK)
	{
		COLORREF c = dlg.GetColor();
		m_ODColors.AddColorREF(c);
		m_D3DWindow.GetInterface()->SetDisplayedColor(COLORREF2D3DCOLOR(c));
		m_D3DWindow.ResetPoints();
		CustomColorsToReg();
	}
}

int CVisViewDlg::CustomColorsFromReg(void)
{
	HKEY   hkey;
	DWORD dwType, dwSize;
	if(RegOpenKeyEx(HKEY_CURRENT_USER, TEXT("Software\\FRI\\VisView"), 0, KEY_READ, &hkey)
		== ERROR_SUCCESS)
	{
		dwType = REG_BINARY;
		dwSize = sizeof(m_acrCustClr);
		RegQueryValueEx(hkey, TEXT("CustomColors"), NULL, &dwType, 
			(PBYTE)&(m_acrCustClr), &dwSize);
		RegCloseKey(hkey);
	 }
	return 0;
}

int CVisViewDlg::CustomColorsToReg(void)
{
	HKEY   hkey;
	DWORD dwType, dwSize;
	DWORD dwDisposition;
	RegCreateKeyEx(HKEY_CURRENT_USER, TEXT("Software\\FRI\\VisView"), 0, NULL, 0, 0, NULL, &hkey, &dwDisposition);
	if(RegOpenKeyEx(HKEY_CURRENT_USER, TEXT("Software\\FRI\\VisView"), 0, KEY_WRITE, &hkey)
		== ERROR_SUCCESS)
	{
		dwType = REG_BINARY;
		dwSize = sizeof(m_acrCustClr);
		RegSetValueEx(hkey, TEXT("CustomColors"), NULL, dwType, 
			(PBYTE)&(m_acrCustClr), dwSize);
		RegCloseKey(hkey);
	 }
	return 0;
}

void CVisViewDlg::OnBnClickedVideo()
{
	//VideoDialog videoDlg(&m_D3DWindow);
	//videoDlg.DoModal();
	//MakeVideo(TEXT("c:\\programs"), TEXT("c:\\programs"));
	m_makingVideo = true;
	m_D3DWindow.QueueScreenshot(D3DXIFF_PNG);
}

void CVisViewDlg::Screenshot()
{
	WIN32_FIND_DATA ffd;
	HANDLE hFind = INVALID_HANDLE_VALUE;
	int idx;
	int maxidx = 0;
	TCHAR filename[20];
	hFind = FindFirstFile(TEXT("screenshot?????.png"), &ffd);
	if(hFind == INVALID_HANDLE_VALUE)
	{
		m_D3DWindow.QueueScreenshot(D3DXIFF_PNG, TEXT("screenshot00000.png"));
	}
	else
	{
		do{
			_stscanf(ffd.cFileName, TEXT("screenshot%5d.png"), &idx);
			if(idx > maxidx) maxidx = idx;
		}
		while(FindNextFile(hFind, &ffd));
		_stprintf(filename, TEXT("screenshot%05d.png"),maxidx+1);
		m_D3DWindow.QueueScreenshot(D3DXIFF_PNG, filename);
	}
}

/*void CVisViewDlg::MakeVideo(LPCTSTR inputFile, LPCTSTR outputDir)
{
	//all frames will have the same camera and displayMap settings
	float wXRot;
    float wYRot;
    float wZRot;
    float vEyeDistance;
	float wZScale;

	int layersel = m_pLayerList->GetCurSel();;
    int nlayers = m_D3DWindow.GetInterface()->LayerCount();
	map<int, dispinfo> displayMap;

	m_D3DWindow.GetInterface()->GetLayerDisplayMap(layersel, displayMap);
	m_D3DWindow.GetCameraProperties(wXRot, wYRot, wZRot, wZScale, vEyeDistance);
	
	WIN32_FIND_DATA ffd;
	HANDLE hFind = INVALID_HANDLE_VALUE;
	
	hFind = FindFirstFile(inputFile, &ffd);
	if(hFind != INVALID_HANDLE_VALUE)
	{
		if(FILE_ATTRIBUTE_DIRECTORY & ffd.dwFileAttributes)
		{
			// directory exsists, find ly* files
			CString wildcard(inputFile);
			wildcard += TEXT("/*.ly?");
			hFind = FindFirstFile(wildcard, &ffd);
			int i=0;
			if(hFind != INVALID_HANDLE_VALUE)
			{
				
				do
				{
					CString name(inputFile);
					name.Append(TEXT("\\"));
					name.Append(ffd.cFileName);
					m_D3DWindow.OpenLYFile(name);
					m_D3DWindow.SetCameraProperties(wXRot, wYRot, wZRot, wZScale, vEyeDistance);
					m_D3DWindow.GetInterface()->SetLayerDisplayMap(layersel, displayMap);
					ResetPartList();
					ResetLayerList();
					m_D3DWindow.ForceDraw();
					Screenshot();

					i++;
					
				}while(FindNextFile(hFind, &ffd));
			}
		}
	}
}*/


void CVisViewDlg::OnBnClickedButtonopennext()
{
    int sel = m_pFileCombo->GetCurSel();

    if (sel != CB_ERR) {
        CString str, str2;
        int len = m_pFileCombo->GetLBTextLen(sel);
    	WIN32_FIND_DATA fd;
    	HANDLE hFind;

        m_pFileCombo->GetLBText(sel, str.GetBuffer(len));

        CPath fileName(str);
        CPath filePath(fileName);
        CString ext = fileName.GetExtension();
        CString result = _T("");
        
        fileName.StripPath();
        filePath.RemoveFileSpec();
        str2 = filePath + _T("\\*") + ext;

        hFind = FindFirstFile(str2, &fd);
        if (hFind != INVALID_HANDLE_VALUE) {
            do {
                if (fd.cFileName == fileName.m_strPath) {
                    if (FindNextFile(hFind, &fd)) 
                        result = filePath + _T("\\") + fd.cFileName;
                    break;
                }
            } while (FindNextFile(hFind, &fd));
        }
        FindClose(hFind);
	    
        if (!result.IsEmpty()) 
            OpenLYFile(result);

        str.ReleaseBuffer();
        
    }
}

void CVisViewDlg::OnBnClickedOk()
{
	OnOK();
}

LRESULT CVisViewDlg::OnScreenshotDone(WPARAM wParam,LPARAM lParam)
{
	if(m_makingVideo)
	{
		int sel = m_pFileCombo->GetCurSel();

		if (sel != CB_ERR) {
			CString str, str2;
			int len = m_pFileCombo->GetLBTextLen(sel);
    		WIN32_FIND_DATA fd;
    		HANDLE hFind;

			m_pFileCombo->GetLBText(sel, str.GetBuffer(len));

			CPath fileName(str);
			CPath filePath(fileName);
			CString ext = fileName.GetExtension();
			CString result = _T("");
	        
			fileName.StripPath();
			filePath.RemoveFileSpec();
			str2 = filePath + _T("\\*") + ext;

			hFind = FindFirstFile(str2, &fd);
			if (hFind != INVALID_HANDLE_VALUE) {
				do {
					if (fd.cFileName == fileName.m_strPath) {
						if (FindNextFile(hFind, &fd)) 
							result = filePath + _T("\\") + fd.cFileName;
						break;
					}
				} while (FindNextFile(hFind, &fd));
			}
			FindClose(hFind);
		    
			if (!result.IsEmpty())
			{
				OpenLYFile(result);
				m_D3DWindow.QueueScreenshot(D3DXIFF_PNG);
			}
			else
			{
				m_makingVideo = false;
			}

			str.ReleaseBuffer();
		}
    }
	return NULL;
}
//void CVisViewDlg::OnBnClickedCheckshowhyponodes()
//{
//    m_D3DWindow.GetInterface()->SetShowHyponodes(m_pShowHyponodesCB->GetCheck() == BST_CHECKED);
//	m_D3DWindow.GetInterface()->ClearLayer0Mapping();
//    m_D3DWindow.ResetPoints();
//}

void CVisViewDlg::OnBnClickedCheckshowlines()
{
	int show = m_showLinesCB.GetCheck();
	m_D3DWindow.SetProperty(CD3DWnd::DisplaySetting::SHOW_LINES, show);
	m_D3DWindow.ResetPoints();
}

void CVisViewDlg::OnBnClickedCheckshowcircles()
{
	int show = m_showCircles.GetCheck();
	m_D3DWindow.SetProperty(CD3DWnd::DisplaySetting::SHOW_CIRCLEPTS, show);
	m_D3DWindow.ResetPoints();
}

void CVisViewDlg::OnBnClickedCheckshowell()
{
	int show = m_showEll.GetCheck();
	m_D3DWindow.SetProperty(CD3DWnd::DisplaySetting::SHOW_ELL, show);
	m_D3DWindow.ResetPoints();
}

void CVisViewDlg::OnBnClickedCheckinhib()
{
	int inhib = m_inhibCB.GetCheck();
	m_D3DWindow.SetProperty(CD3DWnd::DisplaySetting::INHIB, inhib == BST_CHECKED ? 1 : 0);
	m_D3DWindow.GetInterface()->ClearLayer0Mapping();
	m_D3DWindow.ResetPoints();
}

BOOL CVisViewDlg::OnEraseBkgnd(CDC* pDC)
{
	m_D3DWindow.OnDraw(NULL);
	m_D3DWindow.Invalidate(FALSE);
    return CDialog::OnEraseBkgnd(pDC);
}

void CVisViewDlg::OnSetFocus(CWnd* pOldWnd)
{
	CDialog::OnSetFocus(pOldWnd);
    m_D3DWindow.SetFocus();
    m_D3DWindow.OnDraw(NULL);
}

//void CVisViewDlg::OnBnClickedCheckhypo()
//{
//	 if (m_pHypoThresholdCB->GetCheck() == BST_CHECKED) {
//        m_D3DWindow.GetInterface()->SetHypoThreshold(true);
//		m_D3DWindow.GetInterface()->ClearLayer0Mapping();
//		m_pHypoThresholdEdit->EnableWindow(TRUE);
//        m_pHypoThresholdSpin->EnableWindow(TRUE);
//        m_D3DWindow.ResetPoints();
//    } else if (m_pHypoThresholdCB->GetCheck() == BST_UNCHECKED) {
//		CWaitCursor wait;
//        m_D3DWindow.GetInterface()->SetHypoThreshold(false);
//		m_D3DWindow.GetInterface()->ClearLayer0Mapping();
//		m_pHypoThresholdEdit->EnableWindow(FALSE);
//        m_pHypoThresholdSpin->EnableWindow(FALSE);
//        m_D3DWindow.ResetPoints();
//    }
//}

void CVisViewDlg::OnEnUpdateEditthreshold2()
{
    CString str;
    double d;
    
    m_pThresholdEdit2->GetWindowText(str.GetBuffer(10), 9);
    str.ReleaseBuffer();
    d = wcstod(str, NULL);
    d = max(0.0, d);
    d = 20.0 * min(1.0, d);
    d = ((int)d)/20.0;
    
    if (abs(m_D3DWindow.GetInterface()->GetThreshold2() - d) > 0.001) {
        m_D3DWindow.GetInterface()->SetThreshold2(d);
		m_D3DWindow.GetInterface()->ClearLayer0Mapping();
        UpdateThresholdControls();
        m_D3DWindow.ResetPoints();
    }
}

void CVisViewDlg::OnEnUpdateEditthreshold3()
{
    CString str;
    double d;
    
    m_pThresholdEdit3->GetWindowText(str.GetBuffer(10), 9);
    str.ReleaseBuffer();
    d = wcstod(str, NULL);
    d = max(0.0, d);
    d = 20.0 * min(1.0, d);
    d = ((int)d)/20.0;
    
    if (abs(m_D3DWindow.GetInterface()->GetThreshold3() - d) > 0.001) {
        m_D3DWindow.GetInterface()->SetThreshold3(d);
		m_D3DWindow.GetInterface()->ClearLayer0Mapping();
        UpdateThresholdControls();
        m_D3DWindow.ResetPoints();
    }
}

void CVisViewDlg::OnBnClickedCheckpickparts()
{
	m_D3DWindow.SetProperty(CD3DWnd::DisplaySetting::PICK_PARTS, m_pickPartsCB.GetCheck() == BST_CHECKED);
}

LRESULT CVisViewDlg::OnGeometryInitDone(WPARAM wParam, LPARAM lParam) {
	// we only need to check for status of part types and update m_partList
	
	int count = m_partList.GetCount();
	int num_selected = 0;
	bool shouldUpdate = false;

	for (int i = 0; i < count; ++i) {
		int oldValue = m_partList.GetItemData(i);
		int newValue = m_D3DWindow.GetInterface()->PartVisable(i) == true ? 1 : 0;
		m_partList.SetItemData(i, newValue);

		shouldUpdate = shouldUpdate | oldValue != newValue;
		
		if (newValue != 0 && m_partList.GetSel(i) > 0)
			num_selected++;
    }

	// if there is only one part selected then enable "Show All" checkbox otherwise disable it
	if (count > 0 && num_selected <= 1) {
		m_pShowAllCB->EnableWindow(true);		
	} else {
		m_pShowAllCB->EnableWindow(false);
		
		// also set "ShowAll" back to default (checked)
		m_pShowAllCB->SetCheck(BST_CHECKED);
		OnBnClickedCheckshowall();
	}

	if (shouldUpdate == true) {
		m_partList.Invalidate();
		m_partList.UpdateWindow();		
	}

	return S_OK;
}
void CVisViewDlg::OnBnClickedCheckhidetexture()
{
	m_D3DWindow.SetProperty(CD3DWnd::DisplaySetting::HIDE_TEXTURE, m_hideTexture.GetCheck());
}

void CVisViewDlg::OnEnUpdateEditthreshold4()
{
    CString str;
    double d;
    
    m_pThresholdEdit4->GetWindowText(str.GetBuffer(10), 9);
    str.ReleaseBuffer();
    d = wcstod(str, NULL);
    d = max(0.0, d);
    d = 4.0 * min(10.0, d);
    d = ((int)d)/4.0;
    
    if (abs(m_D3DWindow.GetInterface()->GetThreshold4() - d) > 0.001) {
        m_D3DWindow.GetInterface()->SetThreshold4(d);
		m_D3DWindow.GetInterface()->ClearLayer0Mapping();
        UpdateThresholdControls();
        m_D3DWindow.ResetPoints();
    }
}


//void CVisViewDlg::OnShowWindow(BOOL bShow, UINT nStatus)
//{
//    CDialog::OnShowWindow(bShow, nStatus);
//
//    MessageBox(_T("Event!"));
//    // TODO: Add your message handler code here
//}


void CVisViewDlg::OnActivate(UINT nState, CWnd* pWndOther, BOOL bMinimized)
{
	/*
    CDialog::OnActivate(nState, pWndOther, bMinimized);
    if (m_initialized && !bMinimized && nState == WA_ACTIVE)
        m_D3DWindow.D3DReInitialize();
    //    MessageBox(_T("Event!"));
    // TODO: Add your message handler code here
	*/
}


void CVisViewDlg::OnBnClickedButton2()
{
    //m_D3DWindow.D3DReInitialize();
    // TODO: Add your control notification handler code here
}



void CVisViewDlg::OnBnClickedButtonLoadselect()
{
	
	CFileDialog dlg(TRUE, NULL, NULL, 0, _T("All Files (*.*)|*.*||"));
    if (dlg.DoModal() == IDOK) {
		m_D3DWindow.GetInterface()->LoadPartSelection(dlg.GetPathName());
        m_D3DWindow.Reset();
		ResetPartList();
		m_D3DWindow.SetFocus();
		m_D3DWindow.cameraPosFromRegistry();
    }
}
