// Various controls used by VisView program
//

#pragma once

#include "stdafx.h"

// messages

#define LBN_MOUSEITEMCHANGE 100

// CPartListBox
///////////////////////////////////////////////////////////////////////////////

class CPartListBox : public CListBox {
public:
    int m_mouseItem;

    CPartListBox();
    virtual ~CPartListBox();

	// default implementation uses LSB of item data; override to get other behaviour
	virtual BOOL IsItemEnabled(UINT) const;


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPartListBox)
	public:
	virtual void DrawItem(LPDRAWITEMSTRUCT lpDrawItemStruct);
	virtual void MeasureItem(LPMEASUREITEMSTRUCT lpMeasureItemStruct);
	protected:
	virtual void PreSubclassWindow();
	//}}AFX_VIRTUAL


protected:
    afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	// Generated message map functions
	//{{AFX_MSG(CPartListBox)
	//}}AFX_MSG

    DECLARE_MESSAGE_MAP()
};

