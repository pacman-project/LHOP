// Interface ro layer1_result(_export) declarations

#pragma once;

#include "stdafx.h"
#include "../layer_1/layer_1.h"


struct COLOREDPOINT {
    float x, y, z;
    D3DCOLOR color;
};

// CLYInterface
////////////////////////////////////////

struct dispinfo
{
	D3DCOLOR color;
	bool displayed;
};

class CLYInterface {
protected:
    struct layer_data {
        vector<int> m_parts;            // vector of parts
        map<int, dispinfo> m_displayMap;    // is part displayed? could be a set<int>, but 
                                        // it can be a map to a specific color!
    };


protected:
    layer1_result* m_res;
    part_lib* m_library;
    int m_z;

    layer_data m_layers[MAX_LAYER_NUMBER];
    set<D3DCOLOR> m_colorSet;
    bool m_firstOnly;
    double m_responseThresh;
    bool m_showAll;
    bool m_showReconstruction;
	bool m_showTree;
    int m_partIndex;
    int m_indexCount;
    int m_infoHypoNodes;
    int m_infoAllNodes;
	float m_zspacing;
	float m_treeBBmaxx;
	float m_treeBBmaxy;
	float m_treeBBminx;
	float m_treeBBminy;

public:
    CLYInterface();
    ~CLYInterface();

    void OpenLYFile(const wchar_t* str);
    void OpenLibrary(const wchar_t* str);
    int PartCount() { return (int)m_layers[m_z].m_parts.size(); }
    int GetPart(int index);
    void GetPoints(vector<COLOREDPOINT>& pts);
	void GetLines(vector<COLOREDPOINT>& lines);
	//tree
	void GetTree(vector<COLOREDPOINT>& pts, vector<COLOREDPOINT>& sphere_pts);
	void GetLayeredTree(vector<vector<COLOREDPOINT>>& pts, vector<COLOREDPOINT>& toppts);
	float GetTreeBBMaxX() { return m_treeBBmaxx; }
	float GetTreeBBMaxY() { return m_treeBBmaxy; }
	float GetTreeBBMinX() { return m_treeBBminx; }
	float GetTreeBBMinY() { return m_treeBBminy; }
	void GetLayerBB(float& minx, float& maxx, float& miny, float& maxy, int layer);
	//colors
	//D3DCOLOR GetPartColor(int index);
	void SetDisplayedColor(D3DCOLOR color);
	void SetColorDefault();
	//
    void GetImage(void*& mem, unsigned int& size);
    void GetPartImage(void*& mem, int& width, int& height, int layer, int part);
    unsigned GetEdgeInfo(char*& buf);
    bool PartDisplayed(int index);
    void DisplayPart(int index, bool display);
    void ToggleParts(bool select);
    int GetLayer() { return m_z; }
    int LayerCount() { return (m_res == NULL) ? 0 : m_res->layer_count(); }
    void SetLayer(int z);
    void GetLayerDimensions(int& x, int& y, int z);
	void GetLayerZSpacing(float& zspacing);
    void GetCenter(int& cx, int& cy);
    void SetDisplayFirst(bool first) { m_firstOnly = first; }
    bool GetDisplayFirst() { return m_firstOnly; }
    void SetThreshold(double thr) { m_responseThresh = thr; }
    double GetThreshold() { return m_responseThresh; }
    void SetShowAll(bool show) { m_showAll = show; }
    bool GetShowAll() { return m_showAll; }
    void SetShowReconstruction(bool show) { m_showReconstruction = show; }
	void SetShowTree(bool show) { m_showTree = show; }
    bool GetShowReconstruction() { return m_showReconstruction; }
	bool GetShowTree() { return m_showTree; }
    void IncPartIndex();
    void DecPartIndex();
    int GetPartIndex() { return m_partIndex; }
    void ResetIndex() { m_partIndex = 0; }
    int GetHypoNodesCount() { return m_infoHypoNodes; }
    int GetAllNodesCount() { return m_infoAllNodes; }
	void GetLayerDisplayMap(int layer, map<int, dispinfo>& displayMap);
	void SetLayerDisplayMap(int layer, map<int, dispinfo>& displayMap);
	void DebugRelease();

protected:
    void Release();
    void ReleaseLibrary();
    void GetReconstructionPoints(vector<COLOREDPOINT>& pts);
    void GetLayerPoints(vector<COLOREDPOINT>& pts);
};