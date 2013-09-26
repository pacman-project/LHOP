// Interface ro layer1_result(_export) declarations

#pragma once;

//#include "stdafx.h"
#include "layers/layer_1.h"
#include <d3dx9.h>
#include <boost/shared_ptr.hpp>

#define BGR2RGB(X) (((X & 0x000000FF)<<16) | (X & 0x0000FF00) | ((X & 0x00FF0000)>>16))
#define RGB2BGR(X) BGR2RGB(X)

struct COLOREDPOINT {
    float x, y, z;
    D3DCOLOR color;
	COLOREDPOINT(float _x = 0.0f, float _y = 0.0f, float _z = 0.0f, D3DCOLOR _color = 0x00000000) : x(_x),y(_y),z(_z),color(_color) {}
};

struct POINT2
{
	float x,y;
};

//mixcolor class, used for color mixing wiht optimizied tree display
class mixcolor
{
public:	
	
	set<node*> backConnection;

	unsigned int nmixes;
	unsigned int a;
	unsigned int r;
	unsigned int g;
	unsigned int b;

	mixcolor( node* backConnValue = NULL):a(0),r(0),g(0),b(0),nmixes(0) {
		if (backConnValue != NULL)
			backConnection.insert(backConnValue);
	}
	mixcolor(const mixcolor& c):a(c.a),r(c.r),g(c.g),b(c.b),nmixes(c.nmixes),backConnection(c.backConnection) {
	}
	mixcolor(const D3DCOLOR c, node* backConnValue = NULL):a(getA(c)),r(getR(c)),g(getG(c)),b(getB(c)),nmixes(1) {
		if (backConnValue != NULL)
			backConnection.insert(backConnValue);
	}
	
	D3DCOLOR getD3D() const
	{
		D3DCOLOR ret = 0;
		if(nmixes == 0) return ret;
		ret += (a / nmixes)<<24;
		ret += (r / nmixes)<<16;
		ret += (g / nmixes)<<8;
		ret += (b / nmixes);
		return ret;
	}

	mixcolor& operator+=(const mixcolor& c)
	{
		a += c.a;
		r += c.r;
		g += c.g;
		b += c.b;
		nmixes += c.nmixes;
		return *this;
	}
private:
	unsigned int getA(D3DCOLOR c) { return (((c)&0xff000000)>>24); }
	unsigned int getR(D3DCOLOR c) { return (((c)&0x00ff0000)>>16); }
	unsigned int getG(D3DCOLOR c) { return (((c)&0x0000ff00)>>8); }
	unsigned int getB(D3DCOLOR c) { return ((c)&0x000000ff); }
};

// CLYInterface
////////////////////////////////////////

struct dispinfo
{
	D3DCOLOR color;
	bool selected; // if part type is manualy selected by user in part list
	bool visable; // if part type is being displayed (i.e. it may be filtered by mouse picking or thershold setting)
};

class CLYInterface {
protected:
    struct layer_data {
        vector<int> m_parts;            // vector of parts
        map<int, dispinfo> m_displayMap;    // is part displayed? could be a set<int>, but 
                                        // it can be a map to a specific color!
    };

protected:
    boost::shared_ptr<layer1_result> m_res;
    part_lib* m_library;
    int m_z;
    bool m_default_display;

    layer_data m_layers[MAX_LAYER_NUMBER];
	map<node*, pair<int,int>> m_hyporatio;	//number of hypo nodes / all nodes
    map<int, D3DCOLOR> m_colorMap;
    bool m_firstOnly;
    double m_responseThresh;  // r-response
    double m_responseThresh2; // g-response
    double m_responseThresh3; // rr-response
    double m_responseThresh4; // s-response
    double m_responseThresh5; // x-response
	bool m_hyporatioThres;
	double m_hyporatioThresVal;
    bool m_showAll;
    bool m_showReconstruction;
    bool m_showHyponodes;
	bool m_showTree;
	bool m_inhib;
    int m_partIndex;
    int m_indexCount;
    int m_infoHypoNodes;
    int m_infoAllNodes;
	float m_zspacing;
	float m_treeBBmaxx;
	float m_treeBBmaxy;
	float m_treeBBminx;
	float m_treeBBminy;

	dpoint2 m_pickParts;
	dpoint2 m_pickParts_non_reconstruct;
	bool m_isNewPickPartsLocation;
	map<node*, mixcolor> layer0_mapping;
    node* m_pSelectedNode; // pointer to the selected node or NULL if no node is selected
                           // or multiple nodes are selected
public:	
    CLYInterface();
    ~CLYInterface();

    void OpenLYFile(const wchar_t* str);
	void OpenLYObj(boost::shared_ptr<layer1_result> res);
    void OpenLibrary(const wchar_t* str);
    int PartCount() { return (int)m_layers[m_z].m_parts.size(); }
    int GetPart(int index);
    void GetPoints(vector<COLOREDPOINT>& pts);
	void GetLines(vector<COLOREDPOINT>& lines);
	void GetMixedTree(vector<vector<COLOREDPOINT>>& lines, vector<vector<COLOREDPOINT>>& pts);
	float GetTreeBBMaxX() { return m_treeBBmaxx; }
	float GetTreeBBMaxY() { return m_treeBBmaxy; }
	float GetTreeBBMinX() { return m_treeBBminx; }
	float GetTreeBBMinY() { return m_treeBBminy; }
	void GetLayerBB(float& minx, float& maxx, float& miny, float& maxy, int layer);
	//colors
	void SetDisplayedColor(D3DCOLOR color);
	void SetColorDefault();
	//
    void GetImage(void*& mem, unsigned int& size);
    void GetPartImage(void*& mem, int& width, int& height, int layer, int part, int max_width = -1, int max_height = -1);
    unsigned GetEdgeInfo(char*& buf);
    unsigned GetInfo(char*& buf);
    unsigned GetThresholdInfo(char*& buf);
	bool PartVisable(int index);
    bool PartDisplayed(int index);
    void DisplayPart(int index, bool display);
    void ToggleParts(bool select);
    int GetLayer() { return m_z; }
	int GetMaxLayer() { return  m_res->max_layer_index(); }
    int LayerCount() { return (m_res == NULL) ? 0 : m_res->layer_count(); }
    void SetLayer(int z);
    void GetLayerDimensions(int& x, int& y, int z);
	//void GetPartsRect(vector<rectangle2<float>>& partsRect);
	void GetPointsRect(vector<rectangle2<float>>& pointsRect);
	void GetLayerZSpacing(float& zspacing);
    void GetCenter(int& cx, int& cy);
    void SetDisplayFirst(bool first) { m_firstOnly = first; }
    bool GetDisplayFirst() { return m_firstOnly; }
    void SetShowHyponodes(bool hypo) { m_showHyponodes = hypo; }
    bool GetShowHyponodes() { return m_showHyponodes; }
	void SetThreshold(double thr) { m_responseThresh = thr; }
    double GetThreshold() { return m_responseThresh; }
	void SetThreshold2(double thr) { m_responseThresh2 = thr; }
    double GetThreshold2() { return m_responseThresh2; }
	void SetThreshold3(double thr) { m_responseThresh3 = thr; }
    double GetThreshold3() { return m_responseThresh3; }
	void SetThreshold4(double thr) { m_responseThresh4 = thr; }
    double GetThreshold4() { return m_responseThresh4; }
	void SetThreshold5(double thr) { m_responseThresh5 = thr; }
    double GetThreshold5() { return m_responseThresh5; }
	void SetHypoThreshold(const bool t);
	bool GetHypoThreshold() const { return m_hyporatioThres; }
	void SetHypoThresholdVal(const double t) { m_hyporatioThresVal = t; }
	double GetHypoThresholdVal() const { return  m_hyporatioThresVal; }
    void SetShowAll(bool show) { m_showAll = show; }
    bool GetShowAll() { return m_showAll; }
    void SetShowReconstruction(bool show) { m_showReconstruction = show; }
	void SetShowTree(bool show) { m_showTree = show; }
    bool GetShowReconstruction() { return m_showReconstruction; }
	bool GetShowTree() { return m_showTree; }
	void SetInhib(bool inhib) { m_inhib = inhib; }
	bool GetInhib() { return m_inhib; }
    void IncPartIndex();
    void DecPartIndex();
    int GetPartIndex() { return m_partIndex; }
    void ResetIndex() { m_partIndex = 0; }
    int GetHypoNodesCount() { return m_infoHypoNodes; }
    int GetAllNodesCount() { return m_infoAllNodes; }
	void GetLayerDisplayMap(int layer, map<int, dispinfo>& displayMap);
	void SetLayerDisplayMap(int layer, map<int, dispinfo>& displayMap);
	void DebugRelease();
	void GetEllipses(const int ptsPerEll, vector<COLOREDPOINT>& pts);
	void GetMixedReconstructionPoints(vector<COLOREDPOINT>& pts);
	void GetColorMap(map<int, COLORREF>& cmap);

	void BuildOtherLayers(map<int,set<node*>>& topnodes, vector<map<node*, mixcolor>>& lvlnodes, vector<int>& parts, vector<D3DCOLOR>& colors, int& attr_count, int& node_count);
	void SetPickPartsLoc(double x, double y, int z);
	void ClearLayer0Mapping() { layer0_mapping.clear(); }
	void LoadPartSelection(const wchar_t* str);
protected:
    void Release();
    void ReleaseLibrary();
    //void GetReconstructionPoints(vector<COLOREDPOINT>& pts);
    void GetLayerPoints(vector<COLOREDPOINT>& pts);
	//void InhibitPoints(map<int,set<node*>>& pts, const map<int,set<node*>>& allpts, const double thresh);
	bool GetEllipse(double& cx, double& cy, double& a, double& b, double& angle, node* n);
	void Ellipse2Pts(const vector<pair<node*, D3DCOLOR>>& nodes, const int npts, vector<COLOREDPOINT>& pts);
	void GetMixedColorTree(vector<map<node*, mixcolor>>& lvlnodes);
	void GetMixedReconstruction(map<node*, mixcolor>& nodes);
	void SelectIndexedNodes(map<int,set<node*>>& nodes);
	void GetHypoRatio();
	void SelectLayerPoints(const vector<int>& parts, map<int, set<node*>>& nodes, bool ignore_pick_parts = false);
};
