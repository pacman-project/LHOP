//TEST
#include "stdafx.h"
#include "D3D.h"
#include "D3D9Types.h"
#include <mmsystem.h>
#pragma warning( disable : 4996 ) // disable deprecated warning 
#include <strsafe.h>
#pragma warning( default : 4996 )

#define getA(c) (((c)&0xff000000)>>24)
#define getR(c) (((c)&0x00ff0000)>>16)
#define getG(c) (((c)&0x0000ff00)>>8)
#define getB(c) ((c)&0x000000ff)

#define LAYER_ALPHA 0xB0
#define TEXTURE_ALPHA 0xb0
//#define LINE_ALPHA 0x05
#define LINE_ALPHA 0xFF
#define LAYER_COLOR 0xC0
#define CIRCLE_SIZE 0.01f
#define CIRCLE_COLOR_FACTOR 1.2f
#define CIRCLE_COLOR_ADD 0.2f
#define TEXTURE_POWER 0.5f

const int N_CIRCLE_TRIANGLES = 5;
const int N_ELL_POINTS = 12;

// Global variables
///////////////////////////////////////////////////////////////////////////////


// Local functions declarations
///////////////////////////////////////////////////////////////////////////////

// CCamera
///////////////////////////////////////////////////////////////////////////////

CCamera::CCamera()
{
    wXRot = wXRotMem = 0.0f;
    wYRot = wYRotMem = 0.0f;
    wZRot = wZRotMem = 0.0f;
	vXEyePos = 0.0f;
	vYEyePos = 0.0f;
    vEyeDistance = -5.0f;
	wZScale = 1.0f;
}

void CCamera::CalculateWorldMatrix(float height)
{
    D3DXQUATERNION q;
	D3DXMATRIXA16 rot;
	D3DXMATRIXA16 scale;
	D3DXMATRIXA16 trans;
	D3DXMATRIXA16 tmp;
	
	D3DXQuaternionRotationYawPitchRoll(&q, wYRot, wXRot, wZRot);
	D3DXMatrixRotationQuaternion(&rot, &q);
	D3DXMatrixScaling(&scale, 1, 1, wZScale);
	D3DXMatrixTranslation(&trans, 0, 0,  height*wZScale);
  	
	D3DXMatrixMultiply(&tmp, &scale, &trans);
	D3DXMatrixMultiply(&m_mWorld, &tmp, &rot);
}

void CCamera::CalculateViewMatrix()
{
    D3DXVECTOR3 vEyePt(vXEyePos, -vYEyePos, vEyeDistance);
    D3DXVECTOR3 vLookatPt(vXEyePos, -vYEyePos, 0.0f); 
    D3DXVECTOR3 vUpVec(0.0f, 1.0f, 0.0f);

	D3DXMatrixLookAtLH(&m_mView, &vEyePt, &vLookatPt, &vUpVec);
}

void CCamera::CalculateProjMatrix()
{
	// WARNING: ratio between front (0.1f) and back (100.0f) MUST not be over 1000 or 10000 due
	//			to problems with z-buffer which will not be able to differentiate between close object in far view 
    //D3DXMatrixPerspectiveFovLH(&result, D3DX_PI / 4, 1.0f, 1.0f, 100.0f);
	D3DXMatrixPerspectiveFovLH(&m_mProj, D3DX_PI / 4, 1.0f, 0.1f, 100.0f);
	
}

void CCamera::MemorizeRotAndPos()
{
    wXRotMem = wXRot;
    wYRotMem = wYRot;
    wZRotMem = wZRot;

	vXEyePosMem = vXEyePos;
	vYEyePosMem = vYEyePos;
}

void CCamera::Reset()
{
    wXRot = wXRotMem = 0.0f;
    wYRot = wYRotMem = 0.0f;
    wZRot = wZRotMem = 0.0f;
    vEyeDistance = -5.0f;
}

// CD3DWnd
///////////////////////////////////////////////////////////////////////////////

CD3DWnd::CD3DWnd() :
    m_camera(),
    m_mousePos(0, 0),
	m_forbidWndRefresh(false),
	g_pD3D(NULL),
	g_pd3dDevice(NULL),
	g_pTexture(NULL),
	g_pPictureTexture(NULL),
	m_makeScreenshot(false),
	m_showLines(false),
	m_showTree(false),
	m_showCirclePts(false),
	m_showEll(false),
	m_hideTexture(false),
	m_z(0),
	m_hSignaledEvent(0),
	m_reInitGeometry(false),
	m_pickingPartsByMouse(false),
	m_dispalyWaitCursor(true)
{
	m_hInterfaceSemaphore = CreateSemaphore(NULL, 1, 1, NULL);
	m_hKill = CreateEvent( NULL, TRUE, FALSE, TEXT("D3DWndKill")); 
	m_vpinterface.push_back(new CLYInterface());
	m_hSignaledEvent = CreateEvent( NULL, TRUE, TRUE, TEXT("SignaledEvent"));
	m_hInitDone = CreateEvent( NULL, TRUE, FALSE, TEXT("InitDone"));
}

CD3DWnd::~CD3DWnd()
{
	SetEvent(m_hKill);
	HoldInterface();
	for(vector<CLYInterface*>::iterator iter = m_vpinterface.begin(); iter != m_vpinterface.end(); ++iter)
	{
		delete *iter;
	}
	if (g_pTexture != NULL){ g_pTexture->Release(); g_pTexture = NULL; }
	if (g_pPictureTexture != NULL){ g_pPictureTexture->Release(); g_pPictureTexture = NULL; }
	D3DReleaseGeometry();
    if (g_pd3dDevice != NULL) g_pd3dDevice->Release();
    if (g_pD3D != NULL) g_pD3D->Release();
	ReleaseInterface();
	CloseHandle(m_hInterfaceSemaphore);
	CloseHandle(m_hKill);
	CloseHandle(m_hSignaledEvent);
	CloseHandle(m_hInitDone);
}

BEGIN_MESSAGE_MAP(CD3DWnd, CWnd)
	ON_WM_PAINT()
	ON_WM_SIZE()
	ON_WM_CREATE()
	ON_WM_MOUSEMOVE()
    ON_WM_LBUTTONDOWN()
	ON_WM_RBUTTONDOWN()
    ON_WM_LBUTTONUP()
END_MESSAGE_MAP()

void CD3DWnd::OnPaint()
{
   CWnd::OnPaint();
   OnDraw(NULL);
}

void CD3DWnd::OnSize(UINT nType, int cx, int cy)
{
	CWnd::OnSize(nType, cx, cy);
}

int CD3DWnd::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CWnd::OnCreate(lpCreateStruct) == -1) return -1;

	cameraPosFromRegistry();
	D3DInitialize();
	InitGeometry();
	
	CWnd::OnCreate(lpCreateStruct);
	return 0;
}

void CD3DWnd::OnDestroy( )
{
	CWnd::OnDestroy();
}

void CD3DWnd::OnDraw(CDC *pDC)
{
	D3DDrawScene();
}

void CD3DWnd::OnLButtonDown(UINT nFlags, CPoint point)
{
	SetFocus();
    m_mousePos = point;
	m_camera.MemorizeRotAndPos();
	CWnd::OnLButtonDown(nFlags, point);
}
void CD3DWnd::OnRButtonDown(UINT nFlags, CPoint point)
{
	if (m_pickingPartsByMouse == true) {
		SetFocus();
		PickPartByMouseLoc(point);		
	}	
	CWnd::OnRButtonDown(nFlags, point);
}

void CD3DWnd::OnMouseMove(UINT nFlags, CPoint point)
{
	
    if (nFlags & MK_LBUTTON && nFlags & MK_CONTROL) {
		CPoint diff = m_mousePos - point;

		m_camera.vYEyePos = m_camera.vYEyePosMem + 0.001f * (float)diff.y * ::abs(m_camera.vEyeDistance);
        m_camera.vXEyePos = m_camera.vXEyePosMem + 0.001f * (float)diff.x * ::abs(m_camera.vEyeDistance);
		OnDraw(NULL);        
	} else if (nFlags & MK_LBUTTON) {
		CPoint diff = m_mousePos - point;

        m_camera.wXRot = m_camera.wXRotMem + 0.01f * (float)diff.y;
        m_camera.wYRot = m_camera.wYRotMem + 0.01f * (float)diff.x;
		OnDraw(NULL);
	}
	CWnd::OnMouseMove(nFlags, point);
}

BOOL CD3DWnd::OnMouseWheel(UINT nFlags, short zDelta, CPoint point, bool ctrlPressed)
{
	if(GetAsyncKeyState(VK_CONTROL))
	{
		m_camera.wZScale += 0.25f*(float)zDelta/WHEEL_DELTA;
		OnDraw(NULL);
		return TRUE;
	}
    CRect rect;
	
    ScreenToClient(&point);
    GetClientRect(&rect);
	

    if (point.x < rect.left || point.y > rect.right ||
        point.y < rect.top || point.y > rect.bottom)
	{
		OnDraw(NULL);
		return FALSE;
		}

    m_camera.vEyeDistance += 0.5f*(float)zDelta/WHEEL_DELTA;
	OnDraw(NULL);
    return TRUE;
}

void CD3DWnd::D3DCreate(CRect rect, CWnd *parent)
{
	CString className = AfxRegisterWndClass(CS_HREDRAW | CS_VREDRAW | CS_OWNDC, 
        NULL, (HBRUSH)GetStockObject(BLACK_BRUSH), NULL);

	CreateEx(0, className, TEXT("D3D"), WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | WS_CLIPCHILDREN, 
        rect, parent, 0);
}

HRESULT CD3DWnd::D3DInitialize(void)
{
	LPDIRECT3D9 g_pD3D;

    // Create the D3D object.
    if ((g_pD3D = Direct3DCreate9(D3D_SDK_VERSION)) == NULL)
        return E_FAIL;

    // Set up the structure used to create the D3DDevice. Since we are now
    // using more complex geometry, we will create a device with a zbuffer.
    D3DPRESENT_PARAMETERS d3dpp;
    ZeroMemory(&d3dpp, sizeof(d3dpp));
    d3dpp.Windowed = TRUE;
    d3dpp.SwapEffect = D3DSWAPEFFECT_DISCARD;
    d3dpp.BackBufferFormat = D3DFMT_UNKNOWN;
    d3dpp.EnableAutoDepthStencil = TRUE;
    d3dpp.AutoDepthStencilFormat = D3DFMT_D24S8;
	
	d3dpp.hDeviceWindow = m_hWnd;

    // Create the D3DDevice
	//D3DCREATE_HARDWARE_VERTEXPROCESSING?
    /*if (FAILED(g_pD3D->CreateDevice(D3DADAPTER_DEFAULT, D3DDEVTYPE_HAL, m_hWnd,
                                    D3DCREATE_SOFTWARE_VERTEXPROCESSING,
                                    &d3dpp, &g_pd3dDevice)))*/
	if (FAILED(g_pD3D->CreateDevice(D3DADAPTER_DEFAULT, D3DDEVTYPE_HAL, m_hWnd,
                                    D3DCREATE_HARDWARE_VERTEXPROCESSING,
                                    &d3dpp, &g_pd3dDevice)))
    {
        return E_FAIL;
    }

    // Turn off culling
    g_pd3dDevice->SetRenderState(D3DRS_CULLMODE, D3DCULL_NONE);
	//g_pd3dDevice->SetRenderState(D3DRS_CULLMODE, D3DCULL_CCW);

    // Turn off D3D lighting

	g_pd3dDevice->SetRenderState(D3DRS_LIGHTING, FALSE);

	// Turn off the zbuffer
    g_pd3dDevice->SetRenderState(D3DRS_ZENABLE, TRUE);


	if (FAILED(g_pd3dDevice->SetRenderState(D3DRS_BLENDOP, D3DBLENDOP_ADD))) {
			cout << "Unable to SetRenderState D3DRS_BLENDOP in DX3D Device" << endl;
	}
	if (FAILED(g_pd3dDevice->SetRenderState(D3DRS_SRCBLEND, D3DBLEND_SRCALPHA))) {
			cout << "Unable to SetRenderState D3DRS_SRCBLEND in DX3D Device" << endl;
	}
	if (FAILED(g_pd3dDevice->SetRenderState(D3DRS_DESTBLEND, D3DBLEND_INVSRCALPHA))) {
			cout << "Unable to SetRenderState D3DRS_DESTBLEND in DX3D Device" << endl;
	}
	if (FAILED(g_pd3dDevice->SetRenderState( D3DRS_ALPHABLENDENABLE, TRUE ))) {
			cout << "Unable to SetRenderState D3DRS_ALPHABLENDENABLE in DX3D Device" << endl;
	}

	return S_OK;
}

VOID CD3DWnd::SetupMatrices()
{
    D3DXMATRIXA16 matWorld, matView, matProj;
	float ytranslation = 0;
	if(m_vpinterface[0]->GetShowTree()) 
	{
		float zspacing;
		m_vpinterface[0]->GetLayerZSpacing(zspacing);
		ytranslation = - zspacing * m_vpinterface[0]->GetLayer() / 2;
	}
    m_camera.CalculateWorldMatrix(ytranslation);
	g_pd3dDevice->SetTransform(D3DTS_WORLD, m_camera.GetWorldMatrix());

    m_camera.CalculateViewMatrix();
	g_pd3dDevice->SetTransform(D3DTS_VIEW, m_camera.GetViewMatrix());

	m_camera.CalculateProjMatrix();
	g_pd3dDevice->SetTransform(D3DTS_PROJECTION, m_camera.GetProjMatrix());
}

/*void CD3DWnd::SetProjMatToBB()
{
		//----
		//In the perspective transform, the limits of the x- and y-directions are -1 and 1.
		//The limits of the z-direction are 0 for the front plane and 1 for the back plane.
	if(m_bb_min != D3DXVECTOR3(0,0,1) && m_bb_max != D3DXVECTOR3(0,0,1))
	{
		D3DXMATRIXA16 matWorld, matView, matProj, matAll;
		D3DXMATRIXA16 transl;
		D3DXMATRIXA16 scale;
		D3DXMATRIXA16 tmp;
		D3DXVECTOR4 bb_min_t, bb_max_t, bb_min_ts, bb_max_ts;

		g_pd3dDevice->GetTransform(D3DTS_PROJECTION, &matProj);
		g_pd3dDevice->SetTransform(D3DTS_WORLD, &matWorld);
		g_pd3dDevice->SetTransform(D3DTS_VIEW, &matView);
		D3DXMatrixMultiply(&tmp, &matWorld, &matView);
		D3DXMatrixMultiply(&matAll, &tmp, &matProj);
		D3DXVec3Transform(&bb_min_t, &m_bb_min, &matAll);
		D3DXVec3Transform(&bb_max_t, &m_bb_max, &matAll);
		D3DXVec4Scale(&bb_max_ts, &bb_max_t, 1/bb_max_t.w);
		D3DXVec4Scale(&bb_min_ts, &bb_min_t, 1/bb_min_t.w);

		D3DXVECTOR4 d = bb_max_ts - bb_min_ts;
		if(d.x > d.z && d.x > d.y)
		{
			D3DXMatrixScaling(&scale, 1/d.x, 1/d.x, 1/d.x);	
		}
		else if(d.y > d.z)
		{
			D3DXMatrixScaling(&scale, 1/d.y, 1/d.y, 1/d.y);
		}
		else
		{
			D3DXMatrixScaling(&scale, 1/d.z, 1/d.z, 1/d.z);
		}
		D3DXMatrixTranslation(&transl, -(m_bb_max.x + m_bb_min.x) / 2.0f, -(m_bb_max.y + m_bb_min.y) / 2.0f, 0);
		D3DXMatrixMultiply(&tmp ,&matProj, &transl); 
		D3DXMatrixMultiply(&matProj, &tmp, &scale);
		g_pd3dDevice->SetTransform(D3DTS_PROJECTION, &matProj);
	}
	
}*/

HRESULT CD3DWnd::InitTexture()
{
	if (g_pTexture != NULL) { g_pTexture->Release(); g_pTexture = NULL; }
	if (g_pPictureTexture != NULL) { g_pPictureTexture->Release(); g_pPictureTexture = NULL; }
	if (g_pPictureTexture != NULL) g_pPictureTexture->Release();
	if(m_vpinterface.empty()) return S_OK;

    unsigned textureFileSize;
    void* textureMemFile;
    int xDim, yDim;
	HRESULT hr;
    m_vpinterface[0]->GetLayerDimensions(xDim, yDim, 0);
    m_vpinterface[0]->GetImage(textureMemFile, textureFileSize);
    if (textureFileSize > 0) {
        //D3DXCreateTextureFromFile(g_pd3dDevice, L"c:\\text.bmp", &g_pTexture);
        //D3DXCreateTextureFromFileInMemory(g_pd3dDevice, textureMemFile, textureFileSize, &g_pTexture);
		/*D3DXCreateTextureFromFileEx(g_pd3dDevice, TEXT("c:\\text.bmp"),
			D3DX_DEFAULT, D3DX_DEFAULT, D3DX_DEFAULT, 0, D3DFMT_UNKNOWN, D3DPOOL_DEFAULT, D3DX_DEFAULT,
			D3DX_DEFAULT, 0xFF000000, NULL, NULL, &g_pTexture);*/
		D3DXCreateTextureFromFileInMemoryEx(g_pd3dDevice, textureMemFile, textureFileSize,
			D3DX_DEFAULT, D3DX_DEFAULT, D3DX_DEFAULT, 0, D3DFMT_A8R8G8B8, D3DPOOL_MANAGED, D3DX_DEFAULT,
			D3DX_DEFAULT, 0, NULL, NULL, &g_pTexture);
        free(textureMemFile);
    }
	ChangeTextureAlphaGabor(g_pTexture);
	//read picture from file
	if(!m_lyFile.IsEmpty())
	{
		WIN32_FIND_DATA ffd;
		HANDLE hFind = INVALID_HANDLE_VALUE;
		int loc = m_lyFile.ReverseFind(TEXT('.'));
		CString wildcard = m_lyFile.Left(loc+1) + TEXT("*");
		hFind = FindFirstFile(wildcard, &ffd);
		loc = m_lyFile.ReverseFind(TEXT('\\'));
		if(hFind != INVALID_HANDLE_VALUE)
		{
			do
			{
				CString filename;
				if(loc != -1)
				{
					 filename = m_lyFile.Left(loc+1) + ffd.cFileName;
				}
				else
				{
					filename = ffd.cFileName;
				}
				D3DCAPS9 caps;
				g_pd3dDevice->GetDeviceCaps(&caps);
				//hr = D3DXCreateTextureFromFile(g_pd3dDevice, filename, &g_pPictureTexture);
				hr = D3DXCreateTextureFromFileEx(g_pd3dDevice, filename, D3DX_DEFAULT, D3DX_DEFAULT, D3DX_DEFAULT, D3DUSAGE_DYNAMIC , D3DFMT_A8R8G8B8, D3DPOOL_DEFAULT, D3DX_DEFAULT, D3DX_DEFAULT, 0, NULL, NULL, &g_pPictureTexture);
			}
			//D3DPOOL_MANAGED
			while(hr != D3D_OK && FindNextFile(hFind, &ffd));
			if(hr == D3D_OK)
			{
				////add alpha
				//D3DLOCKED_RECT lockedRect;
				//D3DSURFACE_DESC surfDesc;
				//g_pPictureTexture->GetLevelDesc(0, &surfDesc);
				//if(SUCCEEDED(g_pPictureTexture->LockRect(0, &lockedRect, NULL, 0)))
				//{
				//	long* p;
				//	long* end;
				//	for(unsigned int y=0; y<surfDesc.Height; y++)
				//	{
				//		//Pitch is in bytes
				//		p = (long*)lockedRect.pBits + y*lockedRect.Pitch / 4;
				//		end = p + surfDesc.Width;
				//		long lalpha = 0x02<< 24L;
				//		for(; p < end; p++)
				//		{
				//			////long r = (*p & 0x00FFFFFF);
				//			//if(alpha != -1.0f)
				//			//{
				//				//float r = (float)getR(*p) / 255.0f;
				//				//float g = (float)getG(*p) / 255.0f;
				//				//float b = (float)getB(*p) / 255.0f;
				//				//r = pow(r,TEXTURE_POWER);
				//				//b = pow(b,TEXTURE_POWER);
				//				//g = pow(g,TEXTURE_POWER);
				//				//D3DCOLOR c = D3DCOLOR_COLORVALUE((int)(r*255.0f), (int)(g*255.0f), (int)(b*255.0f),TEXTURE_ALPHA) ;
				//				//// *p = TEXTURE_ALPHA<<24L | (*p & 0x00FFFFFFL);
				//				//*p = TEXTURE_ALPHA<<24L | c;
				//			//}
				//			//else
				//			//{
				//				*p = lalpha | *p;
				//			//}
				//		}
				//	}	
				//g_pPictureTexture->UnlockRect(0);
				//}
				ChangeTextureAlpha(g_pPictureTexture, 0xA0);
			}
		}
	}
	return S_OK;
}

void CD3DWnd::OpenLYFile(CString str)
{	
	BeforeOpenLY();
	m_vpinterface.push_back(new CLYInterface());
	wchar_t* c = (wchar_t*)str.GetBuffer(str.GetLength() + 2);
	m_vpinterface[0]->OpenLYFile(c);
	SetInterfaceProperties();
	str.ReleaseBuffer();
	AfterOpenLY();
}

void CD3DWnd::OpenLYObj(vector<boost::shared_ptr<layer1_result>>& res)
{	
	if(res.empty()) return;
	OpenLYObj(res[0]);
}

void CD3DWnd::OpenLYObj(boost::shared_ptr<layer1_result> res)
{	
	BeforeOpenLY();
	m_vpinterface.push_back(new CLYInterface());
	m_vpinterface[0]->OpenLYObj(res);
	SetInterfaceProperties();
	AfterOpenLY();
}
//
//void CD3DWnd::OpenLYObj(vector<boost::shared_ptr<layer1_result>>& res)
//{
//	if(res->empty()) return;
//	BeforeOpenLY();
//	int nnodes = 0;
//	int maxlayer = 0;
//	vector<layer1_result*>::iterator ibest = res->end();
//	for(vector<layer1_result*>::iterator iter = res->begin(); iter != res->end(); ++iter)
//	{
//		if(!*iter) continue;
//		int ml = (*iter)->layer_count()-1;
//		if(ml > maxlayer)
//		{
//			maxlayer = ml;
//			nnodes = (*iter)->get_nnodes(ml);
//			ibest = iter;
//			
//		}
//		else if(ml == maxlayer)
//		{
//			int n = (*iter)->get_nnodes(ml);
//			if(nnodes < n)
//			{
//				ibest = iter;
//				nnodes = n;
//			}
//		}
//	}
//	for(vector<layer1_result*>::iterator iter = res->begin(); iter != res->end(); ++iter)
//	{
//		if(ibest != iter)
//		{
//			delete *iter;
//		}
//		else
//		{
//			m_vpinterface.push_back(new CLYInterface());
//			m_vpinterface.back()->OpenLYObj(*ibest);
//		}
//	}
//	SetInterfaceProperties();
//	AfterOpenLY();
//}



void CD3DWnd::BeforeOpenLY()
{
	HoldInterface();
	for(vector<CLYInterface*>::iterator iter = m_vpinterface.begin(); iter != m_vpinterface.end(); ++iter)
	{
		delete *iter;
	}
	m_vpinterface.clear();
	m_forbidWndRefresh = true;
	if (g_pTexture != NULL){ g_pTexture->Release(); g_pTexture = NULL; }
	if (g_pPictureTexture != NULL){ g_pPictureTexture->Release(); g_pPictureTexture = NULL; }
	D3DReleaseGeometry();
}
void CD3DWnd::AfterOpenLY()
{
	InitTexture();
	InitGeometry();
	m_forbidWndRefresh = false;
	ReleaseInterface();
	OnDraw(0);
}

void CD3DWnd::Reset()
{
	InitTexture();
	SafeInitGeometry();	
}

void CD3DWnd::ResetPoints()
{
	SafeInitGeometry();	
}

void CD3DWnd::cameraPosFromRegistry()
{
	HKEY   hkey;
	DWORD dwType, dwSize;
//	DWORD dwDisposition;

	//if(RegCreateKeyEx(HKEY_CURRENT_USER, TEXT("Software\\FRI\\VisView"), 0, NULL, 0, 0, NULL, &hkey, &dwDisposition)
	//	== ERROR_SUCCESS)
	if(RegOpenKeyEx(HKEY_CURRENT_USER, TEXT("Software\\FRI\\VisView"), 0, KEY_READ, &hkey)
		== ERROR_SUCCESS)
	
	{
		//if(dwDisposition == REG_OPENED_EXISTING_KEY)
		//{
			dwType = REG_BINARY;
			dwSize = sizeof(float);
			RegQueryValueEx(hkey, TEXT("CameraRotX"), NULL, &dwType, 
				(PBYTE)&(m_camera.wXRot), &dwSize);
			m_camera.wXRotMem = m_camera.wXRot;
			RegQueryValueEx(hkey, TEXT("CameraRotY"), NULL, &dwType, 
				(PBYTE)&(m_camera.wYRot), &dwSize);
			m_camera.wYRotMem = m_camera.wYRot;
			RegQueryValueEx(hkey, TEXT("CameraRotZ"), NULL, &dwType, 
				(PBYTE)&(m_camera.wZRot), &dwSize);
			m_camera.wZRotMem = m_camera.wZRot;
			RegQueryValueEx(hkey, TEXT("CameraViewX"), NULL, &dwType, 
				(PBYTE)&(m_camera.vXEyePos), &dwSize);
			m_camera.vXEyePosMem = m_camera.vXEyePos;
			RegQueryValueEx(hkey, TEXT("CameraViewY"), NULL, &dwType, 
				(PBYTE)&(m_camera.vYEyePos), &dwSize);
			m_camera.vYEyePosMem = m_camera.vYEyePos;
			RegQueryValueEx(hkey, TEXT("CameraViewDist"), NULL, &dwType, 
				(PBYTE)&(m_camera.vEyeDistance), &dwSize);
			RegQueryValueEx(hkey, TEXT("CameraZScale"), NULL, &dwType, 
				(PBYTE)&(m_camera.wZScale), &dwSize);
			
		//}
		RegCloseKey(hkey);
	 }
}

void CD3DWnd::cameraPosToRegistry()
{
	HKEY   hkey;
	DWORD dwType, dwSize;
	DWORD dwDisposition;
	RegCreateKeyEx(HKEY_CURRENT_USER, TEXT("Software\\FRI\\VisView"), 0, NULL, 0, 0, NULL, &hkey, &dwDisposition);
	if(RegOpenKeyEx(HKEY_CURRENT_USER, TEXT("Software\\FRI\\VisView"), 0, KEY_WRITE, &hkey)
		== ERROR_SUCCESS)
	{
		dwType = REG_BINARY;
		dwSize = sizeof(float);
		RegSetValueEx(hkey, TEXT("CameraRotX"), NULL, dwType, 
			(PBYTE)&(m_camera.wXRot), dwSize);
		RegSetValueEx(hkey, TEXT("CameraRotY"), NULL, dwType, 
			(PBYTE)&(m_camera.wYRot), dwSize);
		RegSetValueEx(hkey, TEXT("CameraRotZ"), NULL, dwType, 
			(PBYTE)&(m_camera.wZRot), dwSize);
		RegSetValueEx(hkey, TEXT("CameraViewX"), NULL, dwType, 
			(PBYTE)&(m_camera.vXEyePos), dwSize);
		RegSetValueEx(hkey, TEXT("CameraViewY"), NULL, dwType, 
			(PBYTE)&(m_camera.vYEyePos), dwSize);
		RegSetValueEx(hkey, TEXT("CameraViewDist"), NULL, dwType, 
			(PBYTE)&(m_camera.vEyeDistance), dwSize);
		RegSetValueEx(hkey, TEXT("CameraZScale"), NULL, dwType, 
			(PBYTE)&(m_camera.wZScale), dwSize);
		RegCloseKey(hkey);
	 }
}


void CD3DWnd::QueueScreenshot(D3DXIMAGE_FILEFORMAT format, LPCTSTR filename)
{
	InitGeometry();
	if(filename == NULL)
	{
		int loc = m_lyFile.ReverseFind(TEXT('.'));
		m_screenshotFilename = m_lyFile.Left(loc) + TEXT("_vid.png");
	}
	else
	{
		m_screenshotFilename = filename;
	}
	m_screenshotFormat = format;
	m_makeScreenshot = true;
}

HRESULT CD3DWnd::MakeScreenshot()
{
	if(m_makeScreenshot)
	{
		LPDIRECT3DSURFACE9 pbackbuffer;
		g_pd3dDevice->GetBackBuffer(0, 0, D3DBACKBUFFER_TYPE_MONO, &pbackbuffer);	
		HRESULT hr = D3DXSaveSurfaceToFile(m_screenshotFilename, m_screenshotFormat, pbackbuffer, NULL, NULL);
		pbackbuffer->Release();
		//send msg to visviewdlg window
		CWnd* p_pwnd = GetParent();
		m_makeScreenshot = false;
		p_pwnd->SendMessage(WM_USER_SCREENSHOT_DONE);	
		return hr;	
	}
	return NULL;
}

void CD3DWnd::GetCameraProperties(float& wXRot, float& wYRot,float& wZRot,float& wZScale, float& vEyeDistance)
{
	HoldInterface();
	wXRot = m_camera.wXRot;
	wYRot = m_camera.wYRot;
	wZRot = m_camera.wZRot;
	wZScale = m_camera.wZScale;
	vEyeDistance = m_camera.vEyeDistance;
	ReleaseInterface();
}

void CD3DWnd::SetCameraProperties(float wXRot, float wYRot,float wZRot,float wZScale, float vEyeDistance)
{
	HoldInterface();
	m_camera.wXRot = wXRot;
	m_camera.wYRot = wYRot;
	m_camera.wZRot = wZRot;
	m_camera.wZScale = wZScale;
	m_camera.vEyeDistance = vEyeDistance;
	ReleaseInterface();
}

void CD3DWnd::ForceDraw()
{
	D3DDrawScene();
}

void CD3DWnd::D3DReInitialize()
{
    D3DInitialize();
    Reset();
}

// GLayer
CD3DWnd::GLayer::GLayer(LPDIRECT3DDEVICE9 pD3D, float z):
	g_pVBL(NULL),
	g_pVB(NULL),
	g_pVBN(NULL),
	g_pVBC(NULL),
	g_pVBELL(NULL),
	g_pIBELL(NULL),
	g_pTexture(NULL),
	g_pVBTL(NULL),
	m_npoints(0),
	m_nlines(0),
	m_ncircles(0),
	m_nell(0),
	m_dimx(0.0f),
	m_dimy(0.0f),
	g_pd3dDevice(pD3D),
	m_z(z),
	g_pPlaneVertices(NULL)
{
	}

CD3DWnd::GLayer::~GLayer()
{
	
	if (g_pVB != NULL) { g_pVB->Release(); g_pVB = NULL; }
	if (g_pVBL != NULL) { g_pVBL->Release(); g_pVBL = NULL; }
	if (g_pVBN != NULL) { g_pVBN->Release(); g_pVBN = NULL; }
	if (g_pVBC != NULL) { g_pVBC->Release(); g_pVBC = NULL; }
	if (g_pVBELL != NULL) { g_pVBELL->Release(); g_pVBELL = NULL; }
	if (g_pVBTL != NULL) { g_pVBTL->Release(); g_pVBTL = NULL; }
	if (g_pPlaneVertices != NULL) { free(g_pPlaneVertices); g_pPlaneVertices = NULL; }
}

HRESULT CD3DWnd::GLayer::InitPoints(vector<COLOREDPOINT>& pts)
{

	HRESULT hr;

	if (g_pVB != NULL) {int s = g_pVB->Release(); g_pVB = NULL; cout << "Relasing g_pVB and ref count left: " << s << endl;}
	m_npoints = (int)pts.size();
	if (m_npoints > 0) {
		// MUST use D3DPOOL_MANAGED otherwise some graphic cards will not allocate loockable memory and Lock() would always fail
        if (FAILED(g_pd3dDevice->CreateVertexBuffer(m_npoints * sizeof(CUSTOMVERTEX),
                0, D3DFVF_CUSTOMVERTEX, D3DPOOL_MANAGED, &g_pVB, NULL))) {
			cout << "Unable to create vertex buffer" << endl;
            return E_FAIL;
        }
		CUSTOMVERTEX* pVerticesV;

		if (FAILED(hr = g_pVB->Lock(0, 0, (void**)&pVerticesV, 0))) {
			cout << "Unable to lock buffer" << endl;
			return hr;
		}            
        
		for(int i=0; i<(int)m_npoints; i++)
		{
			pVerticesV[i].color = pts[i].color;
			pVerticesV[i].position.x = pts[i].x;
			pVerticesV[i].position.y = pts[i].y;
			pVerticesV[i].position.z = pts[i].z;
		}		
		
		if (FAILED(hr = g_pVB->Unlock())) {
			cout << "Unable to unlock buffer" << endl;
			return hr;
		}
		
    }
	return D3D_OK;
}

HRESULT CD3DWnd::GLayer::InitEll(vector<COLOREDPOINT>& pts)
{

	if (g_pVBELL != NULL) { g_pVBELL->Release(); g_pVBELL = NULL; }
	if (g_pIBELL != NULL) { g_pIBELL->Release(); g_pIBELL = NULL; }
	
	m_nell = (int)pts.size() / N_ELL_POINTS;
	int ni = m_nell * (N_ELL_POINTS+1);
	if (m_nell > 0) {
		// MUST use D3DPOOL_MANAGED otherwise some graphic cards will not allocate loockable memory and Lock() would always fail
        if (
			FAILED(g_pd3dDevice->CreateVertexBuffer(ni * sizeof(CUSTOMVERTEX),
                0, D3DFVF_CUSTOMVERTEX, D3DPOOL_MANAGED, &g_pVBELL, NULL))
			||
			FAILED(g_pd3dDevice->CreateIndexBuffer(ni * sizeof(D3DFMT_INDEX16),
                0, D3DFMT_INDEX16, D3DPOOL_MANAGED, &g_pIBELL, NULL))
				) {
            return E_FAIL;
        }
		CUSTOMVERTEX* pVerticesV;
	
        if (FAILED(g_pVBELL->Lock(0, 0, (void**)&pVerticesV, 0)))
            return E_FAIL;
        CUSTOMVERTEX* pv = pVerticesV;
		/*for(vector<COLOREDPOINT>::const_iterator iter = pts.begin(); iter != pts.end(); ++iter)
		{
			pv->color = iter->color;
			pv->position.x = iter->x;
			pv->position.y = iter->y;
			pv->position.z = iter->z;
			pv++;
		}*/
		int i = 0;
		for(vector<COLOREDPOINT>::const_iterator iter = pts.begin(); iter != pts.end(); ++iter) {
			pv->color = iter->color;
			pv->position.x = iter->x;
			pv->position.y = iter->y;
			pv->position.z = iter->z;
            if (++i == N_ELL_POINTS) {
                ++pv;
                *pv = *(pv - N_ELL_POINTS);
                i = 0;
            }
			++pv;
		}
        g_pIBELL->Unlock();
		WORD* pIdx;
		if (FAILED(g_pIBELL->Lock(0, 0, (void**)&pIdx, 0)))
            return E_FAIL;
        
		/*int idx = 0;
		WORD* p = pIdx;
		for(int i=0; i<m_nell; i++)
		{
			for(int j=0; j<N_ELL_POINTS; j++)
			{
				*p = idx;
				idx++;
				p++;
			}
			*p = idx - N_ELL_POINTS;
			p++;
		}*/

        int imax = (N_ELL_POINTS + 1) * m_nell;
        for (int i = 0; i < imax; ++i) {
            pIdx[i] = i;
        }
        g_pIBELL->Unlock();
    }
	return D3D_OK;
}


HRESULT CD3DWnd::GLayer::InitCircles(vector<COLOREDPOINT>& pts)
{
	if (g_pVBC != NULL) { g_pVBC->Release(); g_pVBC = NULL; }
	m_ncircles = (int)pts.size();
	if(m_ncircles > 0)
	{
		//circles in triangle fans
		D3DXVECTOR3* pos = new D3DXVECTOR3[N_CIRCLE_TRIANGLES+1];
		pos[0].x = 0;
		pos[0].y = 0;
		pos[0].z = 0;
		for(int i=1; i<=N_CIRCLE_TRIANGLES; i++)
		{
			float angle = 2 * 3.141592653f * i/ N_CIRCLE_TRIANGLES;
			pos[i].z = 0;
			pos[i].x = CIRCLE_SIZE * sin(angle);
			pos[i].y = CIRCLE_SIZE * cos(angle);
		}

		// MUST use D3DPOOL_MANAGED otherwise some graphic cards will not allocate loockable memory and Lock() would always fail
		if (FAILED(g_pd3dDevice->CreateVertexBuffer(m_ncircles * (N_CIRCLE_TRIANGLES + 2) * sizeof(CUSTOMVERTEX),
				0, D3DFVF_CUSTOMVERTEX, D3DPOOL_MANAGED, &g_pVBC, NULL))) {
			return E_FAIL;
		}
		CUSTOMVERTEX* pVBC;
		if (FAILED(g_pVBC->Lock(0, 0, (void**)&pVBC, 0)))
        return E_FAIL;
		int idx;
		int stride = N_CIRCLE_TRIANGLES+2;
		float cadd = CIRCLE_COLOR_ADD * 255.0f;
		for(int i=0; i<(int)m_ncircles; i++)
		{
			for(int p=0; p < stride-1; p++)
			{
				idx = i * stride + p;
				pVBC[idx].position.x = pos[p].x + pts[i].x;
				pVBC[idx].position.y = pos[p].y + pts[i].y;
				pVBC[idx].position.z = pos[p].z + pts[i].z;
				// make color darker or brighter
				
				//float r = min(getR(pts[i].color)+ cadd * CIRCLE_COLOR_FACTOR, 255.0f);
				//float g = min(getG(pts[i].color)+ cadd * CIRCLE_COLOR_FACTOR, 255.0f);
				//float b = min(getB(pts[i].color)+ cadd * CIRCLE_COLOR_FACTOR , 255.0f);
				////D3DCOLOR c = D3DCOLOR_ARGB(getA(pts[i].color), (int)r, (int)g, (int)b);
				//D3DCOLOR c = D3DCOLOR_ARGB(0xff, (int)r, (int)g, (int)b);
				//pVBC[idx].color = D3DXCOLOR((int)((float)pts[i].color * (float)CIRCLE_COLOR_FACTOR));
				//pVBC[idx].color = c;
				pVBC[idx].color = pts[i].color;
			}
			idx = (i+1) * stride - 1;
			pVBC[idx] = pVBC[idx - stride + 2];
		}
		g_pVBC->Unlock();
	}
	return D3D_OK;
}

HRESULT CD3DWnd::GLayer::InitLines(vector<COLOREDPOINT>& pts)
{
	if (g_pVBN != NULL) { g_pVBN->Release(); g_pVBN = NULL; }

	m_nlines = (int)pts.size() / 2;
	if(m_nlines > 0)
	{
		// MUST use D3DPOOL_MANAGED otherwise some graphic cards will not allocate loockable memory and Lock() would always fail
		if (FAILED(g_pd3dDevice->CreateVertexBuffer((m_nlines*2)  * sizeof(CUSTOMVERTEX),
				0, D3DFVF_CUSTOMVERTEX, D3DPOOL_MANAGED, &g_pVBN, NULL))) {
			return E_FAIL;
		}
		CUSTOMVERTEX* pVerticesN;
		if (FAILED(g_pVBN->Lock(0, 0, (void**)&pVerticesN, 0))) {
			cout << "Unable to lock vertex buffer in InitLines" << endl;
			return E_FAIL;
		}
		for(int i=0; i<2*(int)m_nlines; i++)
		{
			pVerticesN[i].position.x = pts[i].x;
			pVerticesN[i].position.y = pts[i].y;
			pVerticesN[i].position.z = pts[i].z;
			pVerticesN[i].color = LINE_ALPHA << 24 | pts[i].color & 0x00ffffff;

		}
		g_pVBN->Unlock();
	}
	return D3D_OK;
}

HRESULT CD3DWnd::GLayer::InitPlane(float maxx, float minx, float maxy, float miny, bool display)
{

	if (g_pPlaneVertices == NULL)
		g_pPlaneVertices = new CUSTOMVERTEX[4];
	
	int LAYER_COLOR_ = LAYER_COLOR;

	// create plane even if we do not show it
	g_pPlaneVertices[0].position = D3DXVECTOR3(minx, miny, m_z);
	g_pPlaneVertices[0].color = D3DCOLOR_ARGB(LAYER_ALPHA, LAYER_COLOR_, LAYER_COLOR_, LAYER_COLOR_);
	g_pPlaneVertices[1].position = D3DXVECTOR3(minx, maxy, m_z);
	g_pPlaneVertices[1].color = D3DCOLOR_ARGB(LAYER_ALPHA, LAYER_COLOR_, LAYER_COLOR_, LAYER_COLOR_);
	g_pPlaneVertices[2].position = D3DXVECTOR3(maxx, miny, m_z);
	g_pPlaneVertices[2].color = D3DCOLOR_ARGB(LAYER_ALPHA, LAYER_COLOR_, LAYER_COLOR_, LAYER_COLOR_);
	g_pPlaneVertices[3].position = D3DXVECTOR3(maxx, maxy, m_z);
	g_pPlaneVertices[3].color = D3DCOLOR_ARGB(LAYER_ALPHA, LAYER_COLOR_, LAYER_COLOR_, LAYER_COLOR_);

	if (display == true) {
		if (g_pVBL != NULL) { g_pVBL->Release(); g_pVBL = NULL; }

		// MUST use D3DPOOL_MANAGED otherwise some graphic cards will not allocate loockable memory and Lock() would always fail
		if (FAILED(g_pd3dDevice->CreateVertexBuffer(4 * sizeof(CUSTOMVERTEX),
				0, D3DFVF_CUSTOMVERTEX, D3DPOOL_MANAGED, &g_pVBL, NULL))) {
			return E_FAIL;
		}

		CUSTOMVERTEX* pVertices;
		if (FAILED(g_pVBL->Lock(0, 0, (void**)&pVertices, 0))) {
			cout << "Unable to lock vertex buffer in InitPlane" << endl;
			return E_FAIL;
		}

		// copy to buffer that we will need for ray-tracing when picking only one part		
		memcpy(pVertices, g_pPlaneVertices, 4 * sizeof(CUSTOMVERTEX));

		g_pVBL->Unlock();
	}
    return S_OK;
}

HRESULT CD3DWnd::GLayer::InitTexturePlane(float maxx, float minx, float maxy, float miny, LPDIRECT3DTEXTURE9 pTexture)
{
	g_pTexture = pTexture;
	if (g_pVBTL != NULL) { g_pVBTL->Release(); g_pVBTL = NULL; }
	
	// MUST use D3DPOOL_MANAGED otherwise some graphic cards will not allocate loockable memory and Lock() would always fail
	if (FAILED(g_pd3dDevice->CreateVertexBuffer(4 * sizeof(CUSTOMVERTEXL),
            0, D3DFVF_CUSTOMVERTEXL, D3DPOOL_MANAGED, &g_pVBTL, NULL))) {
        return E_FAIL;
	}

	CUSTOMVERTEXL* pVertices;
	if (FAILED(g_pVBTL->Lock(0, 0, (void**)&pVertices, 0))) {
		cout << "Unable to lock vertex buffer in InitTexturePlane" << endl;
		return E_FAIL;
	}

	pVertices[0].position = D3DXVECTOR3(minx, miny, m_z);
	pVertices[0].color = D3DCOLOR_ARGB(TEXTURE_ALPHA, LAYER_COLOR, LAYER_COLOR, LAYER_COLOR);
	pVertices[0].tu = 1.0f;
	pVertices[0].tv = 0.0f;
	pVertices[1].position = D3DXVECTOR3(minx, maxy, m_z);
	pVertices[1].color = D3DCOLOR_ARGB(TEXTURE_ALPHA, LAYER_COLOR, LAYER_COLOR, LAYER_COLOR);
	pVertices[1].tu = 1.0f;
	pVertices[1].tv = 1.0f;
	pVertices[2].position = D3DXVECTOR3(maxx, miny, m_z );
	pVertices[2].color = D3DCOLOR_ARGB(TEXTURE_ALPHA, LAYER_COLOR, LAYER_COLOR, LAYER_COLOR);
	pVertices[2].tu = 0.0f;
	pVertices[2].tv = 0.0f;
	pVertices[3].position = D3DXVECTOR3(maxx, maxy, m_z  );
	pVertices[3].color = D3DCOLOR_ARGB(TEXTURE_ALPHA, LAYER_COLOR, LAYER_COLOR, LAYER_COLOR);	
	pVertices[3].tu = 0.0f;
	pVertices[3].tv = 1.0f;

	g_pVBTL->Unlock();
    return S_OK;
}

HRESULT CD3DWnd::GLayer::DrawLayer(bool fromAbove, bool m_hideTexture)
{
	// lines if layer bottom is facing the viewer
	if (!fromAbove && g_pVBN != NULL) {
		g_pd3dDevice->SetStreamSource(0, g_pVBN, 0, sizeof(CUSTOMVERTEX));
		g_pd3dDevice->SetFVF(D3DFVF_CUSTOMVERTEX);
		g_pd3dDevice->DrawPrimitive(D3DPT_LINELIST, 0, m_nlines);
	}
	// enable stencil buffer and clear it to zeros (clearing must be done for each layer)
	g_pd3dDevice->SetRenderState(D3DRS_STENCILENABLE, TRUE);
	g_pd3dDevice->Clear(0, NULL, D3DCLEAR_STENCIL, D3DCOLOR_XRGB(0, 0, 0), 1.0f, 0);

	// Set stencil buffer to 1 for each circle, ellipse and point that will be drawn
	g_pd3dDevice->SetRenderState(D3DRS_STENCILREF, 1);
	g_pd3dDevice->SetRenderState(D3DRS_STENCILFUNC, D3DCMP_NOTEQUAL);
	g_pd3dDevice->SetRenderState(D3DRS_STENCILPASS, D3DSTENCILOP_REPLACE);

	//circles
	if(g_pVBC != NULL)
	{
		g_pd3dDevice->SetStreamSource(0, g_pVBC, 0, sizeof(CUSTOMVERTEX));
		g_pd3dDevice->SetFVF(D3DFVF_CUSTOMVERTEX);
		for(int i=0; i<(int)m_ncircles; i++)
		{
			g_pd3dDevice->DrawPrimitive(D3DPT_TRIANGLEFAN, (N_CIRCLE_TRIANGLES+2)*i , N_CIRCLE_TRIANGLES);
		}
	}
	//ellipses
	if(g_pVBELL != NULL)
	{
		g_pd3dDevice->SetStreamSource(0, g_pVBELL, 0, sizeof(CUSTOMVERTEX));
		g_pd3dDevice->SetFVF(D3DFVF_CUSTOMVERTEX);
		g_pd3dDevice->SetIndices(g_pIBELL);
		//g_pd3dDevice->DrawPrimitive(D3DPT_LINELIST, 0, m_nell * N_ELL_POINTS / 2);
		for(int i = 0; i < (int)m_nell; ++i)
		//for(int i = 0; i < (int)1; ++i)
		{
			g_pd3dDevice->DrawIndexedPrimitive(D3DPT_LINESTRIP, 0, (N_ELL_POINTS + 1) * i, N_ELL_POINTS + 1, (N_ELL_POINTS + 1) * i, N_ELL_POINTS);
		}
		//g_pd3dDevice->DrawPrimitive(D3DPT_LINESTRIP, 0, m_nell * N_ELL_POINTS);
	}
	//points
	if (g_pVB != NULL)
	{
		g_pd3dDevice->SetStreamSource(0, g_pVB, 0, sizeof(CUSTOMVERTEX));
		g_pd3dDevice->SetFVF(D3DFVF_CUSTOMVERTEX);
		g_pd3dDevice->DrawPrimitive(D3DPT_POINTLIST, 0, m_npoints);   
	}

	// Display new pixel only when stencil buffer value is 0 (i.e. at this location there is no circle, point or ellipse)
	g_pd3dDevice->SetRenderState(D3DRS_STENCILREF, 0);
	g_pd3dDevice->SetRenderState(D3DRS_STENCILFUNC, D3DCMP_EQUAL);
	g_pd3dDevice->SetRenderState(D3DRS_STENCILPASS, D3DSTENCILOP_KEEP);

	//plane
	if(g_pVBL != NULL)
	{
		g_pd3dDevice->SetStreamSource(0, g_pVBL, 0, sizeof(CUSTOMVERTEX));
		g_pd3dDevice->SetFVF(D3DFVF_CUSTOMVERTEX);
		g_pd3dDevice->DrawPrimitive(D3DPT_TRIANGLESTRIP, 0, 2);
	}
	//texture plane
	if (g_pVBTL != NULL && g_pTexture != NULL && m_hideTexture == false)
	{
		g_pd3dDevice->SetTexture(0, g_pTexture);
		g_pd3dDevice->SetStreamSource(0, g_pVBTL, 0, sizeof(CUSTOMVERTEXL));
		g_pd3dDevice->SetFVF(D3DFVF_CUSTOMVERTEXL);
		g_pd3dDevice->DrawPrimitive(D3DPT_TRIANGLESTRIP, 0, 2);
		g_pd3dDevice->SetTexture(0, NULL);
	}

	g_pd3dDevice->SetRenderState(D3DRS_STENCILENABLE, FALSE);
	//lines if layer top is facing the viewer
	if (fromAbove && g_pVBN != NULL)
	{
		g_pd3dDevice->SetStreamSource(0, g_pVBN, 0, sizeof(CUSTOMVERTEX));
		g_pd3dDevice->SetFVF(D3DFVF_CUSTOMVERTEX);
		g_pd3dDevice->DrawPrimitive(D3DPT_LINELIST, 0, m_nlines);
	}
	return D3D_OK;
}

void CD3DWnd::SafeInitGeometry()
{
	if(HoldInterfaceOrReturn() == 0)
	{
		InitGeometry();
		ReleaseInterface();
		OnDraw(NULL);
	}	
}
HRESULT CD3DWnd::InitGeometry()
{
	CWaitCursor* wait;
	if (m_dispalyWaitCursor == true)
		wait = new CWaitCursor();
	
	for(vector<GLayer*>::iterator i=g_layers.begin(); i != g_layers.end(); i++)
	{
		delete *i;
	}
	if(m_vpinterface.empty()) return S_OK;
	ResetEvent(m_hInitDone);
	g_layers.clear();
	m_treeLayerCount = m_vpinterface[0]->GetLayer();
	bool showTree  = m_vpinterface[0]->GetShowTree();
	float zspacing;
	m_vpinterface[0]->GetLayerZSpacing(zspacing);
	float maxx, minx, maxy, miny;
	if(showTree && g_pPictureTexture != NULL)
	{
		//picture layer (-1)
		m_vpinterface[0]->GetLayerBB(minx, maxx, miny, maxy, 0);
		GLayer* l = new GLayer(g_pd3dDevice, zspacing*(-1));
		l->InitTexturePlane(minx, maxx, miny, maxy, g_pPictureTexture);
		g_layers.push_back(l);

	}
	//layer0
	GLayer* l = new GLayer(g_pd3dDevice, 0);

	m_vpinterface[0]->GetLayerBB(minx, maxx, miny, maxy, 0);
	if(g_pTexture != NULL) l->InitTexturePlane(minx, maxx, miny, maxy, g_pTexture);
	//l->InitPlane(1.5f * minx, 1.5f * maxx, 1.5f * miny, 1.5f * maxy);
	g_layers.push_back(l);
	
	if(showTree)
	{
		g_layers.back()->InitPlane(1.5f * minx, 1.5f * maxx, 1.5f * miny, 1.5f * maxy);
		vector<vector<COLOREDPOINT>> lines;
		vector<vector<COLOREDPOINT>> pts;
		vector<COLOREDPOINT> layer;
		//m_vpinterface[0]->GetLayeredTree2(lines, pts);
		m_vpinterface[0]->GetMixedTree(lines, pts);
		//MergeTrees(lines, pts);
		for(int i=0; i<=m_treeLayerCount; i++)
		{
			if(i != 0)
			{
				m_vpinterface[0]->GetLayerBB(minx, maxx, miny, maxy, i);
				GLayer* l = new GLayer(g_pd3dDevice, zspacing*(i));
				g_layers.push_back(l);

			}
			if(m_showLines && i!=0) g_layers.back()->InitLines(lines[i]);
			g_layers.back()->InitPlane(minx, maxx, miny, maxy);
			if(m_showCirclePts) g_layers.back()->InitCircles(pts[i]);
			else g_layers.back()->InitPoints(pts[i]);
		}
	}
	else
	{
	
		vector<COLOREDPOINT> points;		
		m_vpinterface[0]->GetPoints(points);
		g_layers.back()->InitPlane(1.5f * minx, 1.5f * maxx, 1.5f * miny, 1.5f * maxy, false);
		if (m_showCirclePts) g_layers.back()->InitCircles(points);
		else g_layers.back()->InitPoints(points);
		if (m_showEll) {
			vector<COLOREDPOINT> ell;
			m_vpinterface[0]->GetEllipses(N_ELL_POINTS, ell);
			g_layers.back()->InitEll(ell);
		}		
	}
	SetEvent(m_hInitDone);
	// notify parent that init geometry is done
	GetParent()->PostMessage(WM_GEOMETRY_INIT_DONE, NULL, NULL);
	// restore pointer from waiting to normal
	if (m_dispalyWaitCursor == true) {
		wait->Restore();
		delete wait;
	}
    return S_OK;
}

void CD3DWnd::D3DDrawScene(void)
{	
	if(HoldInterfaceOrReturn() == 1) return;
	if(g_layers.empty())
	{
		ReleaseInterface();
		return;
	}
	//HoldInterface();
    if (FAILED(g_pd3dDevice->Clear(0, NULL, D3DCLEAR_TARGET, D3DCOLOR_XRGB(0, 0, 0), 1.0f, 0))) {
			cout << "Unable to clear D3DCLEAR_TARGET in DX3D Device" << endl;
	}
	if (FAILED(g_pd3dDevice->Clear(0, NULL, D3DCLEAR_ZBUFFER, D3DCOLOR_XRGB(0, 0, 0), 1.0f, 0))) {
			cout << "Unable to clear D3DCLEAR_ZBUFFER in DX3D Device" << endl;
	}

	//SetProjMatToBB();


    // Begin the scene
    if (SUCCEEDED(g_pd3dDevice->BeginScene())) {
        // Setup the world, view, and projection matrices
        SetupMatrices();
		
		for(int i = (int)g_layers.size(); i > 0; i--)
		{
			g_layers[i-1]->DrawLayer(true, m_hideTexture);
		}
		if (FAILED(g_pd3dDevice->EndScene())) {
			cout << "Unable to end scene in DX3D Device" << endl;
		}
    }
	if(m_makeScreenshot) MakeScreenshot();
    // Present the backbuffer contents to the display
    if (FAILED(g_pd3dDevice->Present(NULL, NULL, NULL, NULL))) {
			cout << "Unable to Present scene in DX3D Device" << endl;
	}
	ReleaseInterface();	
}

void CD3DWnd::D3DReleaseGeometry()
{
	for(vector<GLayer*>::iterator i=g_layers.begin(); i != g_layers.end(); i++)
	{
		delete *i;
	}
	g_layers.clear();
}

CLYInterface* CD3DWnd::GetInterface(int n)
{
	if((int)m_vpinterface.size() > n) return m_vpinterface[n];
	return NULL;
}

int CD3DWnd::HoldInterface()
{
	HANDLE hObjOrTerm[2];
	hObjOrTerm[0] = m_hKill;
	hObjOrTerm[1] = m_hInterfaceSemaphore;
	DWORD dwWaitResult = WaitForMultipleObjects((DWORD)2, hObjOrTerm, FALSE, INFINITE);
	switch(dwWaitResult)
	{
	case WAIT_OBJECT_0:
		return 2;
	case WAIT_OBJECT_0 + 1:
		return 0;
	}
	return 1;
}
//
//int CD3DWnd::HoldInterfaceOrReturn()
//{
//	HANDLE handles[3];
//	handles[0] = m_hKill;
//	handles[1] = m_hInterfaceSemaphore;
//	handles[2] = m_hSignaledEvent;
//	DWORD dwWaitResult = WaitForMultipleObjects((DWORD)2, handles, FALSE, INFINITE);
//	switch(dwWaitResult)
//	{
//	case WAIT_OBJECT_0:
//		return 2;
//		break;
//	case WAIT_OBJECT_0 + 1:
//		/* success */
//		return 0;
//		break;
//	case WAIT_OBJECT_0 + 2:
//		return 1;
//	}
//	return 1;
//}

int CD3DWnd::HoldInterfaceOrReturn()
{
	DWORD dwWaitResult = WaitForSingleObject(m_hInterfaceSemaphore, 0);
	switch(dwWaitResult)
	{
	case WAIT_OBJECT_0:		
		return 0;
	default:
		return 1;
	}
}

void CD3DWnd::ReleaseInterface()
{
	if (!ReleaseSemaphore(m_hInterfaceSemaphore, 1, NULL))
	{
		printf("ReleaseSemaphore error: %d\n", GetLastError());
	}
}

void CD3DWnd::SetShowTree(bool d)
{
	for(vector<CLYInterface*>::iterator iter = m_vpinterface.begin(); iter != m_vpinterface.end(); ++iter)
	{
		(*iter)->SetShowTree(d);
	}
	m_showTree = d;
}

void CD3DWnd::SetShowLayer(int z)
{
	m_z = z;
	for(vector<CLYInterface*>::iterator iter = m_vpinterface.begin(); iter != m_vpinterface.end(); ++iter)
	{
		(*iter)->SetLayer(m_z);
	}
}

void CD3DWnd::SetInhib(bool b)
{
	for(vector<CLYInterface*>::iterator iter = m_vpinterface.begin(); iter != m_vpinterface.end(); ++iter)
	{
		(*iter)->SetInhib(b);
	}
}

void CD3DWnd::SetInterfaceProperties()
{
	SetShowTree(m_showTree);
	SetShowLayer(m_z);
}


void CD3DWnd::SetProperty(DisplaySetting s, int param)
{
	HoldInterface();
    switch(s) {
	    case SHOW_LINES:
		    if (param == 0) m_showLines = false; else m_showLines = true;
		    break;
	    case SHOW_TREE:
		    if (param == 0) SetShowTree(false); else SetShowTree(true);
		    break;
	    case LAYER_N:
		    SetShowLayer(param);
		    InitGeometry();
		    break;
	    case LAYER_N_NOINIT:
		    SetShowLayer(param);
		    break;
	    case SHOW_CIRCLEPTS:
		    if (param == 0) m_showCirclePts = false; else m_showCirclePts = true;
		    break;
	    case SHOW_ELL:
		    if (param == 0) m_showEll = false; else m_showEll = true;
		    break;
	    case INHIB:
		    SetInhib(param ? true : false);
		    InitGeometry();
            break;
		case PICK_PARTS:
			m_pickingPartsByMouse = param == 0 ? false : true;
			if (m_pickingPartsByMouse == false) {
				m_vpinterface[0]->SetPickPartsLoc(-1,-1,-1);
				InitGeometry();
			}
			break;
		case HIDE_TEXTURE:
		    if (param == 0) m_hideTexture = false; else m_hideTexture = true;
			break;
	}
	ReleaseInterface();
	OnDraw(NULL);
}

void CD3DWnd::Scroll(double d)
{
	HoldInterface();
	for(vector<CLYInterface*>::iterator iter = m_vpinterface.begin(); iter != m_vpinterface.end(); ++iter)
	{
		(*iter)->SetThreshold(d);
	}
	ResetPoints();
	ReleaseInterface();
	OnDraw(NULL);
}

HANDLE CD3DWnd::GetInitDoneEvent()
{
	return m_hInitDone;
}

HRESULT CD3DWnd::ChangeTextureAlphaGabor(LPDIRECT3DTEXTURE9 ptex)
{
	D3DLOCKED_RECT d3drc;
	if(!ptex) return S_OK;
    if (FAILED(ptex->LockRect(0, &d3drc, NULL, 0)))
        return S_FALSE;   
	D3DSURFACE_DESC surfd;
	g_pTexture->GetLevelDesc(0, &surfd);

	//make image opaque
	for(int y = 0; y < (int)surfd.Height; ++y)
	{
		for(int x = 0; x < (int)surfd.Width; ++x)
		{
			//BGRA, pitch is in bytes
		/*	*((char*)d3drc.pBits + d3drc.Pitch * y + x*4 + 3) = 	
				(char)(((float)*((char*)d3drc.pBits + d3drc.Pitch * y + x*4) + (float)*((char*)d3drc.pBits + d3drc.Pitch * y + x*4 + 1)
				+ (float)*((char*)d3drc.pBits + d3drc.Pitch * y + x*4 + 2) )/ 3);*/
			//make_opaqueBGRA((char*)d3drc.pBits + d3drc.Pitch * y + x*4);
			make_opaqueBGRA((char*)d3drc.pBits + d3drc.Pitch * y + x*4);
		}
	}

	if (FAILED(ptex->UnlockRect(0)))
        return S_FALSE;   
	return S_OK;
}

HRESULT CD3DWnd::ChangeTextureAlpha(LPDIRECT3DTEXTURE9 ptex, const char alpha)
{
	D3DLOCKED_RECT d3drc;
	if(!ptex) return S_OK;
    if (FAILED(ptex->LockRect(0, &d3drc, NULL, 0)))
        return S_FALSE;   
	D3DSURFACE_DESC surfd;
	g_pTexture->GetLevelDesc(0, &surfd);

	//make image opaque
	for(int y = 0; y < (int)surfd.Height; ++y)
	{
		for(int x = 0; x < (int)surfd.Width; ++x)
		{
			//BGRA, pitch is in bytes
		/*	*((char*)d3drc.pBits + d3drc.Pitch * y + x*4 + 3) = 	
				(char)(((float)*((char*)d3drc.pBits + d3drc.Pitch * y + x*4) + (float)*((char*)d3drc.pBits + d3drc.Pitch * y + x*4 + 1)
				+ (float)*((char*)d3drc.pBits + d3drc.Pitch * y + x*4 + 2) )/ 3);*/
			//make_opaqueBGRA((char*)d3drc.pBits + d3drc.Pitch * y + x*4);
			make_transpBGRA((char*)d3drc.pBits + d3drc.Pitch * y + x*4, alpha);
		}
	}

	if (FAILED(ptex->UnlockRect(0)))
        return S_FALSE;   
	return S_OK;
}
//
//LPDIRECT3DTEXTURE9 CD3DWnd::CreateCircleTexture(const vector<COLOREDPOINT>& pts, const rectangle2<float>& brect, const point2<int> res) const
//{
//	HRESULT CreateTexture(res.x, res.y, 1, D3DX_DEFAULT    
//
//		g_pd3dDevice, textureMemFile, textureFileSize,
//			D3DX_DEFAULT, D3DX_DEFAULT, D3DX_DEFAULT, 0, D3DFMT_A8R8G8B8, D3DPOOL_MANAGED, D3DX_DEFAULT,
//			D3DX_DEFAULT, 0, NULL, NULL, &g_pTexture);
//
//    UINT Width,
//    UINT Height,
//    UINT Levels,
//    DWORD Usage,
//    D3DFORMAT Format,
//    D3DPOOL Pool,
//    IDirect3DTexture9** ppTexture,
//    HANDLE* pHandle
//);
//}


void CD3DWnd::PickPartByMouseLoc(CPoint mouseLoc) {

	CRect rect;
    GetClientRect(&rect);

	// make appropriate inverse projections from view space to world space to get ray
    D3DXVECTOR3 vPickRayDir;
    D3DXVECTOR3 vPickRayOrig;

	const D3DXMATRIX* pmatProj = m_camera.GetProjMatrix();

    // Compute the vector of the pick ray in screen space
    D3DXVECTOR3 v;
    v.x = ( ( ( 2.0f * mouseLoc.x ) / rect.Width() ) - 1 ) / pmatProj->_11;
    v.y = -( ( ( 2.0f * mouseLoc.y ) / rect.Height() ) - 1 ) / pmatProj->_22;
    v.z = 1.0f;

    // Get the inverse view matrix
	const D3DXMATRIX matView = *m_camera.GetViewMatrix();
    const D3DXMATRIX matWorld = *m_camera.GetWorldMatrix();
    D3DXMATRIX mWorldView = matWorld * matView;
    D3DXMATRIX m;
    D3DXMatrixInverse( &m, NULL, &mWorldView );

    // Transform the screen space pick ray into 3D space
    vPickRayDir.x = v.x * m._11 + v.y * m._21 + v.z * m._31;
    vPickRayDir.y = v.x * m._12 + v.y * m._22 + v.z * m._32;
    vPickRayDir.z = v.x * m._13 + v.y * m._23 + v.z * m._33;
    vPickRayOrig.x = m._41;
    vPickRayOrig.y = m._42;
    vPickRayOrig.z = m._43;

	// now find closest layer plane that intersects with ray

	FLOAT bestDist = FLT_MAX;
	int intersectPointLayerIndex = -1;

	for (int i = 0; i < g_layers.size(); i++) {
		FLOAT fBary1, fBary2;
		FLOAT fDist;
		// get intersection with layer plane
		CUSTOMVERTEX* planeVertices = g_layers[i]->GetPlaneVertices();
		if (IntersectTriangle(vPickRayOrig, vPickRayDir, 
								planeVertices[0].position, planeVertices[1].position, planeVertices[2].position,
								&fDist, &fBary1, &fBary2) || 
			IntersectTriangle(vPickRayOrig, vPickRayDir, 
								planeVertices[1].position, planeVertices[2].position, planeVertices[3].position,
								&fDist, &fBary1, &fBary2)) {
			// if ray does intersect then check if is closest and keep it
			if (intersectPointLayerIndex < 0 || fDist < bestDist) {
				bestDist = fDist;
				intersectPointLayerIndex = i;
			}
		}
		
	}

	if (intersectPointLayerIndex >= 0) {
		// calculate actualy 3D point in world
		D3DXVECTOR3 intersectPoint = vPickRayOrig + bestDist * vPickRayDir;

		// set location for filtering
		m_vpinterface[0]->SetPickPartsLoc(intersectPoint.x, intersectPoint.y, intersectPointLayerIndex);

		// re-initialize geometry
		SafeInitGeometry();
	}

	return;
}

bool CD3DWnd::IntersectTriangle(const D3DXVECTOR3& orig, const D3DXVECTOR3& dir,
								D3DXVECTOR3& v0, D3DXVECTOR3& v1, D3DXVECTOR3& v2,
								FLOAT* t, FLOAT* u, FLOAT* v ) {
	
	
    // Find vectors for two edges sharing vert0
    D3DXVECTOR3 edge1 = v1 - v0;
    D3DXVECTOR3 edge2 = v2 - v0;

    // Begin calculating determinant - also used to calculate U parameter
    D3DXVECTOR3 pvec;
    D3DXVec3Cross( &pvec, &dir, &edge2 );

    // If determinant is near zero, ray lies in plane of triangle
    FLOAT det = D3DXVec3Dot( &edge1, &pvec );

    D3DXVECTOR3 tvec;
    if( det > 0 )
    {
        tvec = orig - v0;
    }
    else
    {
        tvec = v0 - orig;
        det = -det;
    }

    if( det < 0.0001f )
        return FALSE;

    // Calculate U parameter and test bounds
    *u = D3DXVec3Dot( &tvec, &pvec );
    if( *u < 0.0f || *u > det )
        return FALSE;

    // Prepare to test V parameter
    D3DXVECTOR3 qvec;
    D3DXVec3Cross( &qvec, &tvec, &edge1 );

    // Calculate V parameter and test bounds
    *v = D3DXVec3Dot( &dir, &qvec );
    if( *v < 0.0f || *u + *v > det )
        return FALSE;

    // Calculate t, scale parameters, ray intersects triangle
    *t = D3DXVec3Dot( &edge2, &qvec );
    FLOAT fInvDet = 1.0f / det;
    *t *= fInvDet;
    *u *= fInvDet;
    *v *= fInvDet;

    return TRUE;
}