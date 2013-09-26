// Various controls used by VisView program
//

#include "stdafx.h"
#include "controls.h"

// Method definitions 
///////////////////////////////////////////////////////////////////////////////

// CPartListBox
/////////////////

BEGIN_MESSAGE_MAP(CPartListBox, CListBox)
	ON_WM_MOUSEMOVE()
	//{{AFX_MSG_MAP(CPartListBox)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

CPartListBox::CPartListBox() 
{
    CListBox::CListBox();
    m_mouseItem = -1;
}

CPartListBox::~CPartListBox()
{
    CListBox::~CListBox();
}

void CPartListBox::OnMouseMove(UINT nFlags, CPoint point)
{
    BOOL outside;
    UINT item = ItemFromPoint(point, outside);

    if ((outside && m_mouseItem != -1) || (!outside && item != m_mouseItem)) {
        m_mouseItem = outside ? -1 : item;
        GetParent()->SendMessage(WM_COMMAND, 
            MAKEWPARAM(GetDlgCtrlID(), LBN_MOUSEITEMCHANGE), (LPARAM)m_hWnd);
    }
    CListBox::OnMouseMove(nFlags, point);
}



// default implementation
BOOL CPartListBox::IsItemEnabled(UINT nIndex) const
{
	if (nIndex>=(DWORD)GetCount())
		return TRUE;	// whatever

	DWORD uData=GetItemData(nIndex);

	return (uData&1);
}

/////////////////////////////////////////////////////////////////////////////
// CPartListBox message handlers

void CPartListBox::PreSubclassWindow() 
{
	// TODO: Add your specialized code here and/or call the base class
	ASSERT(GetStyle() & LBS_OWNERDRAWFIXED);
	ModifyStyle(0, LBS_HASSTRINGS | LBS_WANTKEYBOARDINPUT, 0);
	
	CListBox::PreSubclassWindow();
}

void CPartListBox::DrawItem(LPDRAWITEMSTRUCT lpDrawItemStruct) 
{
	// TODO: Add your code to draw the specified item
	ASSERT((GetStyle() & (LBS_OWNERDRAWFIXED | LBS_HASSTRINGS | LBS_WANTKEYBOARDINPUT)) ==
		(LBS_OWNERDRAWFIXED | LBS_HASSTRINGS | LBS_WANTKEYBOARDINPUT));

	CDC* pDC = CDC::FromHandle (lpDrawItemStruct->hDC);

	if (((LONG)(lpDrawItemStruct->itemID) >= 0) &&
		(lpDrawItemStruct->itemAction & (ODA_DRAWENTIRE | ODA_SELECT)))
	{
		BOOL fDisabled = !IsWindowEnabled () || !IsItemEnabled(lpDrawItemStruct->itemID);

		COLORREF newTextColor = fDisabled ?
			RGB(0x80, 0x80, 0x80) : GetSysColor (COLOR_WINDOWTEXT);  // light gray

		COLORREF oldTextColor = pDC->SetTextColor (newTextColor);

		COLORREF newBkColor = GetSysColor (COLOR_WINDOW);
		COLORREF oldBkColor = pDC->SetBkColor (newBkColor);

		if (newTextColor == newBkColor)
			newTextColor = RGB(0xC0, 0xC0, 0xC0);   // dark gray

		if (!fDisabled && ((lpDrawItemStruct->itemState & ODS_SELECTED) != 0))
		{
			pDC->SetTextColor (GetSysColor (COLOR_HIGHLIGHTTEXT));
			pDC->SetBkColor (GetSysColor (COLOR_HIGHLIGHT));
		}

		CString strText;
		GetText (lpDrawItemStruct->itemID, strText);

		const RECT &rc=lpDrawItemStruct->rcItem;

		pDC->ExtTextOut(rc.left + 2,
				  rc.top + 2,// + max(0, (cyItem - m_cyText) / 2),
				  ETO_OPAQUE, &rc,
				  strText, strText.GetLength (), NULL);

		pDC->SetTextColor (oldTextColor);
		pDC->SetBkColor (oldBkColor);
	}

	if ((lpDrawItemStruct->itemAction & ODA_FOCUS) != 0)
		pDC->DrawFocusRect(&lpDrawItemStruct->rcItem);
}

void CPartListBox::MeasureItem(LPMEASUREITEMSTRUCT lpMeasureItemStruct) 
{
	// TODO: Add your code to determine the size of specified item
	UNREFERENCED_PARAMETER(lpMeasureItemStruct);
}