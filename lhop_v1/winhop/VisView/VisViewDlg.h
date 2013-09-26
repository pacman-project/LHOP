// VisViewDlg.h : header file
//

#pragma once

#include "D3D.h"
#include "controls.h"
#include "odcombobox.h"
#include "afxwin.h"

#define WM_USER_SCREENSHOT_DONE (WM_USER+0x101)

// CVisViewDlg dialog
class CVisViewDlg : public CDialog
{
// Construction
public:
	CVisViewDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	enum { IDD = IDD_VISVIEW_DIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support

// Added members
protected:
	bool m_ctrlDown;
	CD3DWnd m_D3DWindow;
    CPartListBox m_partList;
    CListBox* m_pLayerList;
    CStatic* m_pPartStatic;
    CStatic* m_pInfoStatic;
    CComboBox* m_pFileCombo;
    CComboBox* m_pLibraryCombo;
    CButton* m_pReconstructionCB;
    CButton* m_pShowHyponodesCB;
    CButton* m_pFirstOnlyCB;
	CButton* m_pHypoThresholdCB;
	CButton* m_pTreeCB;
    CEdit* m_pThresholdEdit;
    CEdit* m_pThresholdEdit2;
    CEdit* m_pThresholdEdit3;
    CEdit* m_pThresholdEdit4;
    CSpinButtonCtrl* m_pThresholdSpin;
    CSpinButtonCtrl* m_pThresholdSpin2;
    CSpinButtonCtrl* m_pThresholdSpin3;
    CSpinButtonCtrl* m_pThresholdSpin4;
	//CEdit* m_pHypoThresholdEdit;
    //CSpinButtonCtrl* m_pHypoThresholdSpin;
    CButton* m_pShowAllCB;
    CButton* m_pFwd;
    CButton* m_pBwd;
	CODComboBox m_ODColors;
	CHOOSECOLOR m_chooseColor;
	COLORREF m_acrCustClr[16]; // array of custom colors
	bool m_makingVideo;
	bool m_pickingParts;
    bool m_initialized;


    void OpenLYFile(CString fileName);
    void OpenLibraryFile(CString fileName);
    void DrawPart(int pIndex);
    void ResetPartList();
    void ResetLayerList();
    void UpdateThresholdControls();
	int CustomColorsFromReg(void);
	int CustomColorsToReg(void);
	void Screenshot();
	void MakeVideo(LPCTSTR inputFile, LPCTSTR outputDir);
    void DisplayThresholdInfo();

    afx_msg void OnDropFiles(HDROP hDropInfo);
    afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint point);
    afx_msg void OnLbnMouseItemChangePartlist();
    afx_msg void OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
 
// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedButtonOpen();
    afx_msg void OnLbnSelchangePartlist();
    afx_msg void OnLbnSelchangeLayerlist();
    afx_msg void OnBnClickedButtonallparts();
    afx_msg void OnBnClickedButtonnoparts();
    afx_msg void OnCbnSelendokCombofile();
    afx_msg void OnBnClickedButtonOpenlibrary();
    afx_msg void OnCbnSelendokCombolibrary();
    afx_msg void OnBnClickedCheckfirstonly();
    afx_msg void OnEnUpdateEditthreshold();
//    afx_msg void OnBnClickedButtonEdgeinfo();
    afx_msg BOOL OnCopyData(CWnd* pWnd, COPYDATASTRUCT* pCopyDataStruct);
    afx_msg void OnBnClickedCheckshowall();
    afx_msg void OnBnClickedButtonbwd();
    afx_msg void OnBnClickedButtonfwd();
    afx_msg void OnBnClickedCheckreconstruction();
	afx_msg void OnBnClickedSavecam();
	afx_msg void OnBnClickedChecktree();
	afx_msg void OnBnClickedLoadcam();
	afx_msg void OnBnClickedScreenshot();
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnKeyUp(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnBnClickedColor();
	afx_msg void OnCbnSelchangeColorsel();
	afx_msg void OnBnClickedResetcolor();
	afx_msg void OnBnClickedCustomcol();
	afx_msg void OnBnClickedVideo();
//	afx_msg void OnCbnSelchangeCombofile();
    afx_msg void OnBnClickedButtonopennext();
	afx_msg void OnBnClickedOk();
	afx_msg LRESULT OnScreenshotDone(WPARAM wParam,LPARAM lParam);
	CButton m_pMakeVideoButton;
    //afx_msg void OnBnClickedCheckshowhyponodes();
	CButton m_showLinesCB;
	afx_msg void OnBnClickedCheckshowlines();
	CButton m_showCircles;
	afx_msg void OnBnClickedCheckshowcircles();
	afx_msg void OnBnClickedCheckshowell();
	CButton m_showEll;
	CButton m_inhibCB;
	afx_msg void OnBnClickedCheckinhib();
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	afx_msg void OnSetFocus(CWnd* pOldWnd);
//	afx_msg void OnBnClickedCheckhypo();
//    afx_msg void OnEnUpdateEditthreshold2();
    afx_msg void OnEnUpdateEditthreshold2();
//    afx_msg void OnShowWindow(BOOL bShow, UINT nStatus);
    afx_msg void OnEnUpdateEditthreshold3();
	afx_msg void OnBnClickedCheckpickparts();
	CButton m_pickPartsCB;
	LRESULT OnGeometryInitDone(WPARAM wParam, LPARAM lParam);
	CButton m_hideTexture;
	afx_msg void OnBnClickedCheckhidetexture();
    afx_msg void OnBnClickedButtonEdgeinfo();
    afx_msg void OnEnUpdateEditthreshold4();
//    afx_msg void OnShowWindow(BOOL bShow, UINT nStatus);
    afx_msg void OnActivate(UINT nState, CWnd* pWndOther, BOOL bMinimized);
    afx_msg void OnBnClickedButton2();
	afx_msg void OnBnClickedBtnLoadselect();
	afx_msg void OnBnClickedButtonOpen2();
	afx_msg void OnBnClickedButtonLoadselect();
};
