#pragma once


// CODComboBox
#define D3DCOLOR2COLORREF(color) ((COLORREF)( (((color)&0xff)<<16) | (((color)&0xff00)) | (((color)&0xff0000)>>16)))
#define COLORREF2D3DCOLOR(color) ((D3DCOLOR)( 0xff000000 | (((color)&0xff)<<16) | (((color)&0xff00)) | (((color)&0xff0000)>>16) ))

class CODComboBox : public CComboBox
{
	DECLARE_DYNAMIC(CODComboBox)

public:
	CODComboBox();
	virtual ~CODComboBox();
	void AddStdColors();

protected:
	DECLARE_MESSAGE_MAP()
public:
	virtual void DrawItem(LPDRAWITEMSTRUCT /*lpDrawItemStruct*/);
	int AddColorD3D(D3DCOLOR color);
	int AddColorREF(COLORREF color);
};

