// ODComboBox.cpp : implementation file
//

#include "stdafx.h"
#include "VisView.h"
#include "ODComboBox.h"

// CODComboBox

IMPLEMENT_DYNAMIC(CODComboBox, CComboBox)

CODComboBox::CODComboBox()
{

}

CODComboBox::~CODComboBox()
{
}


BEGIN_MESSAGE_MAP(CODComboBox, CComboBox)
END_MESSAGE_MAP()



// CODComboBox message handlers

void CODComboBox::AddStdColors()
{
	AddColorREF(0xFF0000);
	AddColorREF(0x00FF00);
	AddColorREF(0x0000FF);
	AddColorREF(0xFFFF00);
	AddColorREF(0xFF00FF);
	AddColorREF(0x00FFFF);
}

void CODComboBox::DrawItem(LPDRAWITEMSTRUCT lpDrawItemStruct)
{
		// TODO:  Add your code to draw the specified item
	CDC *pDC = CDC::FromHandle(lpDrawItemStruct->hDC);

	switch (lpDrawItemStruct->itemAction)
	{
		case ODA_DRAWENTIRE:
		case ODA_SELECT:
			{
				CBrush brush;
				CBrush *pOldBrush;
				CRect rc;
				int nMargin;

				if (lpDrawItemStruct->itemState & ODS_DISABLED)
				{
					brush.CreateSolidBrush(GetSysColor(COLOR_BTNFACE));
				}
				else if (lpDrawItemStruct->itemState & ODS_SELECTED)
				{
					brush.CreateSolidBrush(GetSysColor(COLOR_HIGHLIGHT));
				}
				else
				{
					brush.CreateSolidBrush(GetSysColor(COLOR_WINDOW));
				}

				// Paint background
				pDC->FillRect(&lpDrawItemStruct->rcItem, &brush);
				brush.DeleteObject();

				// Print color rectangle
				if (lpDrawItemStruct->itemID != -1)
				{
					/*if(lpDrawItemStruct->itemData == 0)
					{
						//add custom color dialog
						rc = lpDrawItemStruct->rcItem;
						nMargin = (rc.bottom - rc.top) / 5;
						rc.left += nMargin;
						rc.top += nMargin;
						rc.bottom -= nMargin;
						rc.right -= nMargin;

						brush.CreateSolidBrush(GetSysColor(COLOR_WINDOW));
						pOldBrush = pDC->SelectObject(&brush);
						pDC->Rectangle(&rc);
						pDC->SelectObject(pOldBrush);

						//LPCTSTR t = LPCTSTR(TEXT("Add..."));
						//pDC->DrawText(t, 6, &lpDrawItemStruct->rcItem, DT_CENTER|DT_SINGLELINE|DT_VCENTER);

					}
					else*/
					{
						rc = lpDrawItemStruct->rcItem;
						nMargin = (rc.bottom - rc.top) / 5;
						rc.left += nMargin;
						rc.top += nMargin;
						rc.bottom -= nMargin;
						rc.right -= nMargin;

						brush.CreateSolidBrush(lpDrawItemStruct->itemData);
						pOldBrush = pDC->SelectObject(&brush);
						pDC->Rectangle(&rc);
						pDC->SelectObject(pOldBrush);
					}
				}

			}

			if (!(lpDrawItemStruct->itemState & ODS_FOCUS))
				break;

			// Else fall through

		case ODA_FOCUS:
			// Draw focus rectangle
			pDC->DrawFocusRect(&lpDrawItemStruct->rcItem);
			break;
	}
	
}

int CODComboBox::AddColorD3D(D3DCOLOR color)
{
	int nItem = AddString(_T(""));
	if (nItem != CB_ERR)
	{
		SetItemData(nItem, D3DCOLOR2COLORREF(color));
	}
	return nItem;
}

int CODComboBox::AddColorREF(COLORREF color)
{
	int nItem = AddString(_T(""));
	if (nItem != CB_ERR)
	{
		SetItemData(nItem, color);
	}
	return nItem;
}
