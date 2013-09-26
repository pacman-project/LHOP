#pragma once
#include "stdafx.h"
#include "../interface/interface.h"

#define WM_USER_SCREENSHOT_DONE (WM_USER+0x101)
#define WM_GEOMETRY_INIT_DONE (WM_USER+0x102)

// #define's
///////////////////////////////////////////////////////////////////////////////

// Global functions
///////////////////////////////////////////////////////////////////////////////

// Direct3D
/////////////


struct CUSTOMVERTEXL
{
    D3DXVECTOR3 position; // The position
    D3DCOLOR color;    // The color
    FLOAT tu, tv;   // The texture coordinates
};

struct CUSTOMVERTEX
{
    D3DXVECTOR3 position; // The position
    D3DCOLOR color;    // The color
};

struct CUSTOMVERTEXN
{
    D3DXVECTOR3 position; // The position
	D3DXVECTOR3 normal;	// The normal
    D3DCOLOR color;    // The color
};

struct MESHVERTEX
{
    D3DXVECTOR3 position; // The position
	D3DXVECTOR3 normal;	// The normal
};


#define D3DFVF_CUSTOMVERTEXL (D3DFVF_XYZ | D3DFVF_DIFFUSE | D3DFVF_TEX1)
#define D3DFVF_CUSTOMVERTEX (D3DFVF_XYZ | D3DFVF_DIFFUSE)
#define D3DFVF_CUSTOMVERTEXN (D3DFVF_XYZ | D3DFVF_NORMAL | D3DFVF_DIFFUSE)

// DEFINITIONS
///////////////////////////////////////////////////////////////////////////////




// classes
///////////////////////////////////////////////////////////////////////////////

// Camera
///////////

class CCamera {
public:
    // World rotation
    float wXRot, wXRotMem;
    float wYRot, wYRotMem;
    float wZRot, wZRotMem;
	// World scaling
	float wZScale;

    // Viewpoint position
	float vXEyePos, vXEyePosMem;
	float vYEyePos, vYEyePosMem;
    float vEyeDistance;

	// all calculated matrixes 
	D3DXMATRIX m_mView;
	D3DXMATRIX m_mProj;
	D3DXMATRIX m_mWorld;

    // Methods
    CCamera();

    void CalculateWorldMatrix(float height);
	void CalculateWorldMatrixSphere();
    void CalculateViewMatrix();
    void CalculateProjMatrix();
    void MemorizeRotAndPos();
    void Reset();

	const D3DXMATRIX* GetViewMatrix() const { return &m_mView; }
    const D3DXMATRIX* GetProjMatrix() const { return &m_mProj; }
	const D3DXMATRIX* GetWorldMatrix() const { return &m_mWorld; }		
};


// C3D3Wnd
////////////

class CD3DWnd : public CWnd {
private:
	class GLayer
	{
		private:
			LPDIRECT3DDEVICE9       g_pd3dDevice;
			LPDIRECT3DVERTEXBUFFER9 g_pVBL;         // Buffer to hold vertices of the layer plane
			LPDIRECT3DVERTEXBUFFER9 g_pVBTL;			// Buffer to hold vertices of the picture texture layer
			LPDIRECT3DVERTEXBUFFER9 g_pVB;          // Buffer to hold point vertices 
			LPDIRECT3DVERTEXBUFFER9 g_pVBN;         // Buffer to hold line vertices
			LPDIRECT3DVERTEXBUFFER9 g_pVBC;         // Buffer to hold circle vertices

			LPDIRECT3DVERTEXBUFFER9 g_pVBELL;         // Buffer to hold ellipse vertices
			LPDIRECT3DINDEXBUFFER9	g_pIBELL;			// Buffer to hold ellipse indices

			LPDIRECT3DTEXTURE9      g_pTexture;     // layer texture
			
			CUSTOMVERTEX*			g_pPlaneVertices; // actual layer plane vertexes

			unsigned int m_npoints;
			unsigned int m_ncircles;
			unsigned int m_nell;
			unsigned int m_nlines;
			float m_dimx;
			float m_dimy;
			float m_z;			
		public:			
			GLayer(LPDIRECT3DDEVICE9 pD3D, float z);
			~GLayer();
			HRESULT InitPoints(vector<COLOREDPOINT>& pts);
			HRESULT InitLines(vector<COLOREDPOINT>& pts);
			HRESULT InitCircles(vector<COLOREDPOINT>& pts);
			HRESULT InitEll(vector<COLOREDPOINT>& pts);
			HRESULT InitPlane(float maxx, float minx, float maxy, float miny, bool display = true);
			HRESULT InitTexturePlane(float maxx, float minx, float maxy, float miny, LPDIRECT3DTEXTURE9 pTexture);
			HRESULT DrawLayer(bool fromAbove, bool m_hideTexture);
			//void GetLayerBB(D3DXVECTOR3& max, D3DXVECTOR3& min);
			CUSTOMVERTEX* GetPlaneVertices() { return g_pPlaneVertices; }
	};
	
	/* vector of CLYInterface objects for displaying more scales */
	vector<CLYInterface*> m_vpinterface;

	bool m_forbidWndRefresh;
	bool m_showLines;
	bool m_showTree;
	bool m_showCirclePts;
	bool m_showEll;
	bool m_hideTexture;
	int m_z;

	bool m_pickingPartsByMouse;

	/* making OpenLY* functions thread safe */
	HANDLE m_hKill;
	HANDLE m_hInterfaceSemaphore;	// interface vector semaphore
	HANDLE m_hInitDone;
	int HoldInterfaceOrReturn();
	int HoldInterface();
	void ReleaseInterface();
	/* InitGeometry must not be executed concurrently */
	HANDLE m_hSignaledEvent;

	bool m_reInitGeometry;
	bool m_dispalyWaitCursor;
    // D3D
	LPDIRECT3D9             g_pD3D;         // Used to create the D3DDevice
	LPDIRECT3DDEVICE9       g_pd3dDevice;   // Our rendering device
	LPDIRECT3DTEXTURE9      g_pTexture;     // Our texture
	LPDIRECT3DTEXTURE9		g_pPictureTexture;	// layer -1 image
	vector<GLayer*>			g_layers;
	
	//screenshot
	CString m_screenshotFilename;
	D3DXIMAGE_FILEFORMAT m_screenshotFormat;
	bool m_makeScreenshot;

    // Wnd
    CRect m_origSize;
    CPoint m_mousePos;
    CCamera m_camera;

    // Data
	CString m_lyFile;
	int m_treeLayerCount;


	void D3DReleaseGeometry();
	HRESULT D3DInitialize();
	void SafeInitGeometry();
	HRESULT InitGeometry();
	HRESULT InitTexture();
	void D3DDrawScene();
	void SetupMatrices();
	void SetProjMatToBB();
	void SetInterfaceProperties();	//sets display tree, etc.  for all interfaces  according to D3D's settings
	/* functions before and after adding new CLYInterface-s, making it thread safe */
	void BeforeOpenLY();
	void AfterOpenLY();
	/* private display settings */
	void SetShowTree(bool d);
	void SetShowLayer(int z);
	void SetInhib(bool b);
	/* */
//	LPDIRECT3DTEXTURE9 CreateCircleTexture(const vector<COLOREDPOINT>& pts, const rectangle2<float>& brect, const point2<int> res) const;

public:

	CD3DWnd();
	virtual ~CD3DWnd();
    // D3D
    void D3DReInitialize();
	CLYInterface* GetInterface(int n=0);	//return n-th interface object in m_vinterface
	void ForceDraw();
	void D3DCreate(CRect rect, CWnd *parent);
    void OpenLYFile(CString str);
	void OpenLYObj(boost::shared_ptr<layer1_result> res);
	void OpenLYObj(vector<boost::shared_ptr<layer1_result>>& res);
	void SetLyFilePath(CString& s) { m_lyFile = s; }
    void Reset();
    void ResetPoints();
	void cameraPosToRegistry();
	void cameraPosFromRegistry();
	void GetCameraProperties(float& wXRot, float& wYRot,float& wZRot,float& wZScale, float& vEyeDistance);
	void SetCameraProperties(float wXRot, float wYRot,float wZRot,float wZScale, float vEyeDistance);
	HANDLE GetInitDoneEvent();
	void Scroll(double d);
	void PickPartByMouseLoc(CPoint mouseLoc);
	bool IntersectTriangle(const D3DXVECTOR3& orig, const D3DXVECTOR3& dir,
							D3DXVECTOR3& v0, D3DXVECTOR3& v1, D3DXVECTOR3& v2,
							FLOAT* t, FLOAT* u, FLOAT* v );
	void SetDisplayCursor(bool display) { m_dispalyWaitCursor = display; }
	/* funtion for setting display properties, must be thread safe */
	enum DisplaySetting
	{
		SHOW_LINES,
		SHOW_TREE,
		SHOW_CIRCLEPTS,
		SHOW_ELL,
		LAYER_N,
		LAYER_N_NOINIT,
		INHIB,
		PICK_PARTS,
		HIDE_TEXTURE
	};
	void make_opaqueBGRA(char* p)
	{
		/*float val = (float)*p + (float)*(p+1) + (float)*(p+2);
		val /= 3.0f;
		*p = (char)255;
		*(p+1) = (char)255;
		*(p+2) = (char)255;
		*(p+3) = (char)val;*/
		//float val = (float)*p + (float)*(p+1) + (float)*(p+2);
		//val /= 3.0f;
		//if(val < 20.0f)
		if((unsigned int)*p < (unsigned int)25)
		{
			/**p = (char)255;
			*(p+1) = (char)255;
			*(p+2) = (char)255;
			*///*(p+3) = (char)val;// * 10;
			*(p+3) = (char)(10*(unsigned int)*p);
		}
		else *(p+3) = (char)255;
		//else if(*(p+1) == 0 && *(p+2) == 0 && *(p) == 0) *(p+3) = 0;
	}

	void make_transpBGRA(char* p, const char alpha)
	{
		// 1-(1-avg)^2
		/*float val = 1.0f-(((float)*p + (float)*(p+1) + (float)*(p+2)) / (255.0f * 3.0f));
		val = 254.0f * (1.0f-val*val);*/
		float val = (float)*p + (float)*(p+1) + (float)*(p+2);
		val /= 3.0f;
		*(p+3) = alpha;		
	}

	HRESULT ChangeTextureAlpha(LPDIRECT3DTEXTURE9 ptex, const char alpha);
	HRESULT ChangeTextureAlphaGabor(LPDIRECT3DTEXTURE9 ptex);
	void SetProperty(DisplaySetting s, int param = 0);
	HRESULT MakeScreenshot();
	void QueueScreenshot(D3DXIMAGE_FILEFORMAT format, LPCTSTR filename = NULL);
	//void GetBB(D3DXVECTOR3& max, D3DXVECTOR3& min);

	// Message response methods:
	afx_msg void OnPaint();
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg	void OnDraw(CDC* pDC);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg int  OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnDestroy( );
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
//    afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
    BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint point, bool ctrlPressed = false); // called by parent!
	DECLARE_MESSAGE_MAP()

};




