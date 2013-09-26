// Interface to layer1_result(_export) definitions
//

#include "stdafx.h"
#include "interface.h"

// layer1_result_type_export
//////////////////////////////

class layer1_result_type_export {
public:
    typedef map<node*, int> container_t;
    typedef container_t::iterator iterator_t;
    typedef container_t::const_iterator const_iterator_t;

    container_t src, result;
    int attr_count;
    int node_count;
    unsigned attr;

    layer1_result_type_export(unsigned a = 0) : 
        src(), result(), attr_count(0), node_count(0), attr(a) { }

    void reset() { src.clear(); result.clear(); }

    void add_src_node(node* n) 
    {
        layer1_data* nd = (layer1_data*)n->data;
        if (n->is_attr_set2(attr)) ++attr_count;
        src.insert(container_t::value_type(n, nd->m));
        ++node_count;
    }

    int size() { return (int)src.size(); }

    void keep_index(int index) 
    {
        if (index < 0 || src.empty()) return;
        if (index >= size()) index = size() - 1;

        int t, count = 0;
        node* n;

        for (iterator_t iter = src.begin(); iter != src.end(); ++iter) {
            if (count++ == index) { n = iter->first; t = iter->second; }
        }
        src.clear();
        src.insert(container_t::value_type(n, t));
    }

    node* get_node(const container_t::value_type& iter) { return iter.first; }

    bool check_neighbor(container_t& neighbors, const container_t::value_type& n, node* nn)
    {
        neighbors.insert(container_t::value_type(nn, n.second));
        if (nn->is_attr_set2(attr)) ++attr_count;
        ++node_count;
        return true;
    }

    void insert_leaf(const container_t::value_type& n) { result.insert(n); }
};

// CLYInterface definitions
/////////////////////////////

CLYInterface::CLYInterface()
{
    m_res = NULL;
    m_library = NULL;
    m_z = 0;
    m_firstOnly = true;
    m_responseThresh = 0.3;
    m_showAll = true;
    m_showReconstruction = true;
	m_showTree = false;
    m_partIndex = 0;
    m_indexCount = 0;
	m_zspacing = -0.1f;

	m_colorSet.insert(D3DCOLOR_ARGB(0xFF, 0xFF, 0x00, 0x00));
	m_colorSet.insert(D3DCOLOR_ARGB(0xFF, 0x00, 0xFF, 0x00));
	m_colorSet.insert(D3DCOLOR_ARGB(0xFF, 0x00, 0x00, 0xFF));
	m_colorSet.insert(D3DCOLOR_ARGB(0xFF, 0xFF, 0xFF, 0x00));
	m_colorSet.insert(D3DCOLOR_ARGB(0xFF, 0xFF, 0x00, 0xFF));
	m_colorSet.insert(D3DCOLOR_ARGB(0xFF, 0x00, 0xFF, 0xFF));

    m_infoHypoNodes = 0;
    m_infoAllNodes = 0;
}

CLYInterface::~CLYInterface()
{
    Release();
    ReleaseLibrary();
}

void CLYInterface::Release()
{
    if (m_res != NULL) delete m_res;
    m_res = NULL;
}

void CLYInterface::DebugRelease()
{
	if (m_res != NULL)
	{
		delete m_res;
		m_res = NULL;
	}
	if (m_library != NULL)
	{
		delete m_library;
		m_library = NULL;
	}
    for(int i=0; i<MAX_LAYER_NUMBER; i++)
	{
		m_layers[i].m_parts.clear();
		m_layers[i].m_displayMap.clear();
	}
	m_colorSet.clear();
}

void CLYInterface::ReleaseLibrary()
{
    if (m_library != NULL) delete m_library;
    m_library = NULL;
}

void CLYInterface::OpenLibrary(const wchar_t* str)
{
    size_t slen = wcslen(str);
    char* fname = new char[slen + 1];
    wcstombs(fname, str, slen + 1);

    ReleaseLibrary();
    try {
        read_library(fname, m_library);
    } catch (...) {
        m_library = NULL;
    }
    delete fname;
}

void CLYInterface::OpenLYFile(const wchar_t* str)
{
	size_t slen = wcslen(str);
    char* fname = new char[slen + 1];
    wcstombs(fname, str, slen + 1);

    Release();
    try {
        read_layer1_result(m_res, fname);
    } catch (...) {
        m_res = NULL;
    }
    if (m_res == NULL) {
        m_z = 0;
        m_layers[0].m_parts.clear();
    } else {
		set<D3DCOLOR>::iterator citer;
        int maxz = m_res->max_layer_index();

        if (m_z > maxz) m_z = maxz;
		dispinfo di;
		di.displayed = true;
        for (int z = 0; z <= maxz; ++z) {
            set<int> partSet;
            int maxi = (int)m_res->shape_nodes[z].size() - 1;
            m_layers[z].m_parts.clear();
            m_res->get_parts(partSet, z);
			citer = m_colorSet.begin();
            for (set<int>::iterator iter = partSet.begin(); iter != partSet.end(); ++iter) {
                int p = *iter;
                m_layers[z].m_parts.push_back(p);
				di.color = *citer;
				citer++;
				if(citer == m_colorSet.end()) citer = m_colorSet.begin();
				m_layers[z].m_displayMap.insert(pair<int, dispinfo>(p, di));
            }
        }
    }
    delete fname;
}


int CLYInterface::GetPart(int index)
{
    if (index >= 0 && index < (int)m_layers[m_z].m_parts.size()) return m_layers[m_z].m_parts[index];
    else return -1;
}

bool CLYInterface::PartDisplayed(int index)
{
    if (index < 0 || index >= (int)m_layers[m_z].m_parts.size()) return false;
    else {
        map<int, dispinfo>::iterator iter = m_layers[m_z].m_displayMap.find(m_layers[m_z].m_parts[index]);
		if (iter == m_layers[m_z].m_displayMap.end()) return false;
		return iter->second.displayed;
    }
}

void CLYInterface::DisplayPart(int index, bool display)
{
    if (index < 0 || index >= (int)m_layers[m_z].m_parts.size()) return;
    else { 
        map<int, dispinfo>::iterator iter = m_layers[m_z].m_displayMap.find(m_layers[m_z].m_parts[index]);
       // if (iter != m_layers[m_z].m_displayMap.end()) iter->second = display;
		if (iter != m_layers[m_z].m_displayMap.end()) iter->second.displayed = display;
    }
}

void CLYInterface::ToggleParts(bool select)
{
    map<int, dispinfo>& dmap = m_layers[m_z].m_displayMap;

    for (map<int, dispinfo>::iterator iter = dmap.begin(); iter != dmap.end(); ++iter) {
		iter->second.displayed = select;
    }
}

void CLYInterface::GetReconstructionPoints(vector<COLOREDPOINT>& pts)
{
    pts.clear();
    if (m_res == NULL) return;

	D3DCOLOR color;
    layer1_result_type_export exporter(HYPO_NODE_ATTR);
    set<node*> nodes;
    vector<int> srcParts;
    float factor = 2.0f/m_res->x_size(0);
    float yOffset = (factor * m_res->y_size(0))/2.0f;
    COLOREDPOINT pt;
    int colorCount = (int)m_colorSet.size();
    double thresh = m_firstOnly ? -1.0 : m_responseThresh;
    int index = m_showAll ? -1 : m_partIndex;

    for (map<int, dispinfo>::iterator iter = m_layers[m_z].m_displayMap.begin(); 
            iter != m_layers[m_z].m_displayMap.end(); ++iter)
	{
        if (iter->second.displayed) srcParts.push_back(iter->first);
		else continue;
		color = iter->second.color;
		m_indexCount = m_res->get_reconstruction_nodes(exporter, m_z, srcParts, thresh, index);
		if (index >= 0 && m_partIndex >= m_indexCount) m_partIndex = m_indexCount - 1;
		for (layer1_result_type_export::iterator_t iter = exporter.result.begin(); 
				iter != exporter.result.end(); ++iter) {
			layer1_data* nd = (layer1_data*)(iter->first)->data;

			if (nd->z != 0) continue;
			pt.x = factor * nd->x - 1.0f;
			pt.y = yOffset - factor * nd->y;
			pt.z = (float)nd->z - 0.001f;
			pt.color = color;
			pts.push_back(pt);
		}
		nodes.clear();
		srcParts.clear();
	}
    m_infoHypoNodes = exporter.attr_count;
    m_infoAllNodes = exporter.node_count;
}

void CLYInterface::GetLayerPoints(vector<COLOREDPOINT>& pts)
{
    pts.clear();
    if (m_res == NULL) return;
	D3DCOLOR color;
    set<node*> nodes;
    vector<int> srcParts;
    float factor = 2.0f/m_res->x_size(m_z);
    float yOffset = (factor * m_res->y_size(m_z))/2.0f;
    COLOREDPOINT pt;
    int colorCount = (int)m_colorSet.size();
    double thresh = m_firstOnly ? -1.0 : m_responseThresh;
    
    for (map<int, dispinfo>::iterator iter = m_layers[m_z].m_displayMap.begin(); 
            iter != m_layers[m_z].m_displayMap.end(); ++iter)
	{
		if (iter->second.displayed) srcParts.push_back(iter->first);
		else continue;

		color = iter->second.color;
		m_res->get_layer_nodes(nodes, m_z, srcParts, thresh);
		for (set<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
			layer1_data* nd = (layer1_data*)(*iter)->data;

			pt.x = factor * nd->x - 1.0f;
			pt.y = yOffset - factor * nd->y;
			pt.z = /*(float)nd->z*/ - 0.001f;
			pt.color = color;
			pts.push_back(pt);
		}
		nodes.clear();
		srcParts.clear();
	}   
    m_partIndex = 0;
    m_indexCount = 1;
    m_infoHypoNodes = 0;
    m_infoAllNodes = (int)nodes.size();
}

void CLYInterface::GetTree(vector<COLOREDPOINT>& pts, vector<COLOREDPOINT>& sphere_pts)
{
	pts.clear();
	sphere_pts.clear();
    if (m_res == NULL) return;

    set<node*> nodes;
	//==
	set<node*> nodes0;
	set<node*> nodes1;
	set<node*>* neighborNodes = &nodes0;
	set<node*>* levelNodes = &nodes1;
	//==
    vector<int> srcParts;
	vector<short> color;

    float factor;
    float yOffset;
	float nfactor;
	float nyOffset;
	float xOffset;
	float nxOffset;
    COLOREDPOINT pt;
	COLOREDPOINT ptNeighbor;

	m_treeBBmaxx = FLT_MIN;
	m_treeBBmaxy = FLT_MIN;
	m_treeBBminx = FLT_MAX;
	m_treeBBminy = FLT_MAX;

    int colorCount = (int)m_colorSet.size();
    double thresh = m_firstOnly ? -1.0 : m_responseThresh;
	D3DCOLOR partColor;
    
    for (map<int, dispinfo>::iterator iterd = m_layers[m_z].m_displayMap.begin(); 
            iterd != m_layers[m_z].m_displayMap.end(); ++iterd)
	{
		if(iterd->second.displayed) srcParts.push_back(iterd->first);
		else continue;

       	
		partColor = iterd->second.color;
		// for each part
		m_res->get_layer_nodes(*levelNodes, m_z, srcParts, thresh);
		//for (set<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
		//	node* nd = *iter;
		//	levelNodes->push_back(nd);
		//}	

		//add lower levels
		for(int i=m_z; levelNodes->size() > 0; --i)
		{
			factor = 2.0f/m_res->x_size(i) * m_res->x_size(i) / m_res->x_size(0);
			yOffset = (factor * m_res->y_size(i))/2.0f;
			xOffset = -(float)m_res->x_size(i) / (float)m_res->x_size(0);

			nfactor = 2.0f /m_res->x_size(i-1) * m_res->x_size(i-1) / m_res->x_size(0);
			nyOffset = (factor * m_res->y_size(i-1))/2.0f;
			nxOffset = -(float)m_res->x_size(i-1) / (float)m_res->x_size(0);
			
			for(set<node*>::iterator iter = levelNodes->begin(); iter != levelNodes->end(); ++iter)
			{
				layer1_data* nd = (layer1_data*)(*iter)->data;
				pt.x = factor * nd->x  + xOffset;
				pt.y = yOffset - factor * nd->y;
				pt.z = (float)nd->z * m_zspacing - 0.001f;
				pt.color = partColor;
				
				if(pt.x > m_treeBBmaxx) m_treeBBmaxx = pt.x;
				if(pt.y > m_treeBBmaxy) m_treeBBmaxy = pt.y;
				if(pt.x < m_treeBBminx) m_treeBBminx = pt.x;
				if(pt.y < m_treeBBminy) m_treeBBminy = pt.y;

				forall_neighbors(*iter, iter1)
				{
					//add to pts list
					//neighborNodes->push_back(neighbor_node(iter1));
                    neighborNodes->insert(neighbor_node(iter1));
					
					layer1_data* nn = (layer1_data*)(neighbor_node(iter1))->data;
					ptNeighbor.x = nfactor * nn->x + nxOffset;
					ptNeighbor.y = nyOffset - nfactor * nn->y;
					ptNeighbor.z = (float)nn->z * m_zspacing - 0.001f;
					ptNeighbor.color = partColor;
					pts.push_back(ptNeighbor);
					pts.push_back(pt);
					if(i == 1)
					{
						sphere_pts.push_back(ptNeighbor);
					}
				}
				sphere_pts.push_back(pt);
			}
			levelNodes->clear();
			set<node*>* tmp = levelNodes;
			levelNodes = neighborNodes;
			neighborNodes = tmp;
		}
		levelNodes->clear();
		nodes.clear();
		srcParts.clear();
	}
}

void CLYInterface::GetLayeredTree(vector<vector<COLOREDPOINT>>& pts, vector<COLOREDPOINT>& toppts)
{
	pts.clear();
	toppts.clear();
	pts.reserve(m_z);
	for(int i=0; i<m_z; i++)
	{
		pts.push_back ( vector<COLOREDPOINT>());
	}

    if (m_res == NULL) return;

    set<node*> nodes;
	//==
	set<node*> nodes0;
	set<node*> nodes1;
	set<node*>* neighborNodes = &nodes0;
	set<node*>* levelNodes = &nodes1;
	//==
    vector<int> srcParts;
	vector<short> color;

    float factor;
    float yOffset;
	float nfactor;
	float nyOffset;
	float xOffset;
	float nxOffset;
    COLOREDPOINT pt;
	COLOREDPOINT ptNeighbor;

	m_treeBBmaxx = FLT_MIN;
	m_treeBBmaxy = FLT_MIN;
	m_treeBBminx = FLT_MAX;
	m_treeBBminy = FLT_MAX;

    int colorCount = (int)m_colorSet.size();
    double thresh = m_firstOnly ? -1.0 : m_responseThresh;
	D3DCOLOR partColor;
    

    for (map<int, dispinfo>::iterator iterd = m_layers[m_z].m_displayMap.begin(); 
            iterd != m_layers[m_z].m_displayMap.end(); ++iterd)
	{
		if(iterd->second.displayed) srcParts.push_back(iterd->first);
		else continue;

       	
		partColor = iterd->second.color;
		// for each part
		m_res->get_layer_nodes(*levelNodes, m_z, srcParts, thresh);
		//for (set<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
		//	node* nd = *iter;
		//	levelNodes->push_back(nd);
		//}	

		//add lower levels
		for(int i=m_z; levelNodes->size() > 0; --i)
		{
			factor = 2.0f/m_res->x_size(i) * m_res->x_size(i) / m_res->x_size(0);
			yOffset = (factor * m_res->y_size(i))/2.0f;
			xOffset = -(float)m_res->x_size(i) / (float)m_res->x_size(0);

			nfactor = 2.0f /m_res->x_size(i-1) * m_res->x_size(i-1) / m_res->x_size(0);
			nyOffset = (factor * m_res->y_size(i-1))/2.0f;
			nxOffset = -(float)m_res->x_size(i-1) / (float)m_res->x_size(0);
			
			for(set<node*>::iterator iter = levelNodes->begin(); iter != levelNodes->end(); ++iter)
			{
				layer1_data* nd = (layer1_data*)(*iter)->data;
				pt.x = factor * nd->x  + xOffset;
				pt.y = yOffset - factor * nd->y;
				pt.z = (float)nd->z * m_zspacing - 0.001f;
				pt.color = partColor;
				if(i == m_z) toppts.push_back(pt);
				
				if(pt.x > m_treeBBmaxx) m_treeBBmaxx = pt.x;
				if(pt.y > m_treeBBmaxy) m_treeBBmaxy = pt.y;
				if(pt.x < m_treeBBminx) m_treeBBminx = pt.x;
				if(pt.y < m_treeBBminy) m_treeBBminy = pt.y;

				forall_neighbors(*iter, iter1)
				{
					//add to pts list
					//neighborNodes->push_back(neighbor_node(iter1));
                    neighborNodes->insert(neighbor_node(iter1));
					
					layer1_data* nn = (layer1_data*)(neighbor_node(iter1))->data;
					ptNeighbor.x = nfactor * nn->x + nxOffset;
					ptNeighbor.y = nyOffset - nfactor * nn->y;
					ptNeighbor.z = (float)nn->z * m_zspacing - 0.001f;
					ptNeighbor.color = partColor;
					pts[i-1].push_back(pt);
					pts[i-1].push_back(ptNeighbor);
				}
			}
			levelNodes->clear();
			set<node*>* tmp = levelNodes;
			levelNodes = neighborNodes;
			neighborNodes = tmp;
		}
		levelNodes->clear();
		nodes.clear();
		srcParts.clear();
	}
}

void CLYInterface::GetLayerBB(float& minx, float& maxx, float& miny, float& maxy, int layer)
{
	//-1 .. 1
	if(m_res == NULL)
	{
		minx = 0;
		maxx = 0;
		miny = 0;
		maxy = 0;
		return;
	}
	if(layer < GetLayer() || layer == 0)
	{
		float yFactor = (float)m_res->y_size(layer)/(float)m_res->x_size(layer);
		float ratiox = (float)m_res->x_size(layer) / (float)m_res->x_size(0);
		float ratioy = (float)m_res->y_size(layer) / (float)m_res->y_size(0);
		minx = -ratiox;
		maxx = ratiox;
		miny = -ratioy * yFactor;
		maxy = ratioy * yFactor;
	}
	
}

void CLYInterface::GetPoints(vector<COLOREDPOINT>& pts)
{
    if (m_showReconstruction) GetReconstructionPoints(pts);
    else GetLayerPoints(pts);
}

void CLYInterface::GetImage(void*& mem, unsigned int& size)
{
    if (m_res == NULL) size = 0;
    else {
        img* image = m_res->get_image_reconstructed(0, 0, vector<int>());
        //image->save("c:\\text.bmp");
        size = image->save(mem);
		//+++
		delete image;
		//+++
    }
}

void CLYInterface::GetPartImage(void*& mem, int& width, int& height, int layer, int part)
{
    mem = NULL;
    if (m_library == NULL) return;
    if (layer < 0 || layer >= (int)m_library->parts.size()) return;
    if (part < 0 || part >= (int)m_library->parts[layer].size()) return;

    node* n = m_library->parts[layer][part];
    img im = ((part_data*)(n->data))->get_image(n);

    width = im.width;
    height = im.height;
    im.get_bits(mem, width, height, img::BIT_TYPE_32);
}

unsigned CLYInterface::GetEdgeInfo(char*& buf)
{
    if (m_res == NULL) {
        buf = NULL;
        return 0;
    } else {
        strstream os;

        m_res->write_edge_info(os);
        os << '\0';

        char* str = os.str();
        size_t size = os.pcount();

        buf = (char*)malloc(sizeof(char) * size);
        memcpy(buf, str, sizeof(char) * size);
        return size;
    }
}

void CLYInterface::GetLayerDimensions(int& x, int& y, int z)
{
    if (m_res == NULL) {
        x = 0; y = 0;
    } else {
        x = m_res->x_size(z);
        y = m_res->y_size(z);
    }
}

void CLYInterface::GetLayerZSpacing(float& zspacing)
{
	zspacing = m_zspacing;
}

void CLYInterface::SetLayer(int z)
{
    if (z < 0 || z > m_res->max_layer_index()) m_z = m_res->max_layer_index();
    else m_z = z;
}

void CLYInterface::IncPartIndex()
{
    if (m_indexCount == 0) m_partIndex = 0;
    else m_partIndex = (m_partIndex + 1) % m_indexCount;
}

void CLYInterface::DecPartIndex()
{
    if (m_indexCount == 0) m_partIndex = 0;
    else if (m_partIndex == 0) m_partIndex = m_indexCount - 1;
    else --m_partIndex;
}

void CLYInterface::SetDisplayedColor(D3DCOLOR color)
{
	for (map<int, dispinfo>::iterator iter = m_layers[m_z].m_displayMap.begin(); 
            iter != m_layers[m_z].m_displayMap.end(); ++iter)
	{
		if(iter->second.displayed)
		{
			iter->second.color = color;
		}
	}
}

void CLYInterface::SetColorDefault()
{
	set<D3DCOLOR>::iterator citer = m_colorSet.begin();
	for(int cur_z = 0; cur_z < LayerCount(); ++cur_z)
	{
		for (map<int, dispinfo>::iterator iter = m_layers[cur_z].m_displayMap.begin(); 
				iter != m_layers[cur_z].m_displayMap.end(); ++iter)
		{			
			iter->second.color = *citer;
			citer++;
			if(citer == m_colorSet.end()) citer = m_colorSet.begin();
		}
	}
}

void CLYInterface::GetLayerDisplayMap(int layer, map<int, dispinfo>& displayMap)
{
	if(layer >= LayerCount()) displayMap.clear();
	else
	{
		displayMap = m_layers[layer].m_displayMap;
	}
}

void CLYInterface::SetLayerDisplayMap(int layer, map<int, dispinfo>& displayMap)
{
	if(layer < LayerCount())
	{
		map<int, dispinfo>::iterator miter = m_layers[layer].m_displayMap.begin();
		for(
			map<int, dispinfo>::iterator iter = displayMap.begin();
			iter != displayMap.end() && miter != m_layers[layer].m_displayMap.end();
			++iter, ++miter)
		{
			miter->second = iter->second;
		}
		//if there are more parts, hide them
		for(;miter != m_layers[layer].m_displayMap.end(); ++miter)
		{
			miter->second.displayed = false;
		}
	}
	for(int i=0; i<LayerCount(); i++)
	{
		if(i != layer)
		{
			
			for(map<int, dispinfo>::iterator miter = m_layers[layer].m_displayMap.begin();
				miter != m_layers[layer].m_displayMap.end();
				++miter)
			{
				miter->second.displayed = false;
			}
		}
	}
}