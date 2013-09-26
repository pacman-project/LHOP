// Interface to layer1_result(_export) definitions
//

//#include "stdafx.h"
#include "interface.h"
#include "utils/structures.h"
#include <strstream>

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

    bool check_neighbor(container_t& neighbors, const container_t::value_type& n, edge_data*, node* nn)
    {
        neighbors.insert(container_t::value_type(nn, n.second));
        if (nn->is_attr_set2(attr)) ++attr_count;
        ++node_count;
        return true;
    }

    void insert_leaf(const container_t::value_type& n) { result.insert(n); }
};

//

int get_index_count(const map<int,set<node*>>& partnodes, const map<int, dispinfo> partdisplay)
{
	int retval = 0;
	for(map<int,set<node*>>::const_iterator iter = partnodes.begin();
		iter != partnodes.end(); ++iter)
	{
		//retval += (int)iter->second.size();
		map<int, dispinfo>::const_iterator partdisp_iter = partdisplay.find(iter->first);
		if (partdisp_iter != partdisplay.end() && partdisp_iter->second.selected == true)
			retval = max((int)iter->second.size(),retval);
	}
	return retval;
}

// CLYInterface definitions
/////////////////////////////

CLYInterface::CLYInterface():
	m_library(0),
    m_z(0),
    m_firstOnly(true),
    m_responseThresh(0.3),
    m_responseThresh2(0.3),
    m_responseThresh3(0.5),
    m_responseThresh4(2.5),
    m_responseThresh5(0.5),
	m_hyporatioThres(false),
	m_hyporatioThresVal(0.3),
    m_showAll(true),
    m_showReconstruction(true),
    m_showHyponodes(true),
	m_showTree(false),
    m_partIndex(0),
    m_indexCount(0),
	m_zspacing(0.1f),
    m_default_display(true),
	m_inhib(false),
	m_infoHypoNodes(0),
    m_infoAllNodes(0),
	m_isNewPickPartsLocation(false),
    m_pSelectedNode(0)
{
	m_pickParts.x = -1;
	m_pickParts.y = -1;
	m_pickParts_non_reconstruct.x = -1;
	m_pickParts_non_reconstruct.y = -1;
	m_colorMap.insert(pair<int, D3DCOLOR>(0, D3DCOLOR_ARGB(0xFF, 0xFF, 0x00, 0x00)));
	m_colorMap.insert(pair<int, D3DCOLOR>(1, D3DCOLOR_ARGB(0xFF, 0x00, 0xFF, 0x00)));
	m_colorMap.insert(pair<int, D3DCOLOR>(2, D3DCOLOR_ARGB(0xFF, 0x00, 0x00, 0xFF)));
	m_colorMap.insert(pair<int, D3DCOLOR>(3, D3DCOLOR_ARGB(0xFF, 0xFF, 0xFF, 0x00)));
	m_colorMap.insert(pair<int, D3DCOLOR>(4, D3DCOLOR_ARGB(0xFF, 0xFF, 0x00, 0xFF)));
	m_colorMap.insert(pair<int, D3DCOLOR>(5, D3DCOLOR_ARGB(0xFF, 0x00, 0xFF, 0xFF)));
	m_colorMap.insert(pair<int, D3DCOLOR>(6, D3DCOLOR_ARGB(0xFF, 0xFF, 0xFF/2, 0x00)));
	m_colorMap.insert(pair<int, D3DCOLOR>(7, D3DCOLOR_ARGB(0xFF, 0xFF, 0x00, 0xFF/2)));
}

CLYInterface::~CLYInterface()
{
    ReleaseLibrary();
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

	boost::shared_ptr<layer1_result> res;
    try {
        res.reset(read_layer1_result(fname));
    } catch (...) {
        res = boost::shared_ptr<layer1_result>();
    }
    OpenLYObj(res);
    delete fname;
}

void CLYInterface::OpenLYObj(boost::shared_ptr<layer1_result> res)
{
	m_res = res;
	if (m_res == NULL) {
        m_z = 0;
        m_layers[0].m_parts.clear();
    } else {
		set<D3DCOLOR>::iterator citer;
        int maxz = m_res->max_layer_index();

        if (m_z > maxz) m_z = maxz;
		dispinfo di;
		di.selected = m_default_display;
		di.visable = true;
        for (int z = 0; z <= maxz; ++z) {
            set<int> partSet;
            int maxi = (int)m_res->shape_nodes[z].size() - 1;
            m_layers[z].m_parts.clear();
            m_res->get_parts(partSet, z);
            for (set<int>::iterator iter = partSet.begin(); iter != partSet.end(); ++iter) {
                int p = *iter;
                m_layers[z].m_parts.push_back(p);
				di.color = m_colorMap.find( p %(int)m_colorMap.size())->second;
				m_layers[z].m_displayMap.insert(pair<int, dispinfo>(p, di));
            }
        }
		//compute hypo node ratio
		if(GetHypoThreshold()) GetHypoRatio();
	}
}

void CLYInterface::GetHypoRatio()
{
	m_hyporatio.clear();
	if(m_res == NULL || m_res->shape_nodes.empty()) return;
	////insert nodes of layer 0
	//for(vector<node*>::const_iterator iter = m_res->shape_nodes[0].begin();
	//	iter != m_res->shape_nodes[0].end(); ++iter)
	//{
	//	if((*iter)->is_attr_set(HYPO_NODE_ATTR))
	//	{
	//		m_hyporatio.insert(pair<node*,pair<int,int>>(*iter, pair<int,int>(1,1)));
	//	}
	//	else
	//	{
	//		m_hyporatio.insert(pair<node*,pair<int,int>>(*iter, pair<int,int>(0,1)));
	//	}
	//}
	////sum over neighbours for layer 1,2 etc etc
	//if(m_res->shape_nodes.size() > 1)
	//{
	//	for(vector<vector<node*>>::const_iterator iter = m_res->shape_nodes.begin() + 1;
	//		iter != m_res->shape_nodes.end(); ++iter)
	//	{
	//		for(vector<node*>::const_iterator inode = iter->begin(); inode != iter->end(); ++inode)
	//		{
	//			int hn = 0;	//hypo
	//			int an = 0;	//all nodes
	//			int name = atom("toPrevLayer");
	//			//foreach neighbour
	//			foreach_neighbor(*inode, name, ineigh)
	//			{
	//				map<node*,pair<int,int>>::iterator f = m_hyporatio.find(neighbor_node(ineigh));
	//				if(f != m_hyporatio.end())
	//				{
	//					hn += f->second.first;
	//					an += f->second.second;
	//				}
	//			}
	//			if((*inode)->is_attr_set(HYPO_NODE_ATTR)) ++hn;
	//			++an;
	//			m_hyporatio.insert(pair<node*,pair<int,int>>(*inode, pair<int,int>(hn,an)));
	//		}
	//	}
	//}

	//precise
	for(vector<vector<node*>>::const_iterator iter = m_res->shape_nodes.begin();
			iter != m_res->shape_nodes.end(); ++iter)
	{
		for(vector<node*>::const_iterator inode = iter->begin(); inode != iter->end(); ++inode)
		{
			node* n = *inode;
			while(n != NULL)
			{
				layer1_result_type_export exporter(HYPO_NODE_ATTR);
				m_res->get_hypo_node_count(exporter, *inode);
				m_hyporatio.insert(pair<node*,pair<int,int>>(n,
					pair<int,int>(exporter.attr_count, exporter.node_count)));
				n = ((layer1_data*)(n->data))->next;
			}
		}
	}
}



int CLYInterface::GetPart(int index)
{
    if (index >= 0 && index < (int)m_layers[m_z].m_parts.size()) return m_layers[m_z].m_parts[index];
    else return -1;
}

bool CLYInterface::PartVisable(int index) {
	if (index < 0 || index >= (int)m_layers[m_z].m_parts.size()) return false;
    else {
        map<int, dispinfo>::iterator iter = m_layers[m_z].m_displayMap.find(m_layers[m_z].m_parts[index]);
		if (iter == m_layers[m_z].m_displayMap.end()) return false;
		return iter->second.visable;
    }
}
bool CLYInterface::PartDisplayed(int index)
{
    if (index < 0 || index >= (int)m_layers[m_z].m_parts.size()) return false;
    else {
        map<int, dispinfo>::iterator iter = m_layers[m_z].m_displayMap.find(m_layers[m_z].m_parts[index]);
		if (iter == m_layers[m_z].m_displayMap.end()) return false;
		return iter->second.selected;
    }
}

void CLYInterface::DisplayPart(int index, bool selected)
{
    if (index < 0 || index >= (int)m_layers[m_z].m_parts.size()) return;
    else { 
        map<int, dispinfo>::iterator iter = m_layers[m_z].m_displayMap.find(m_layers[m_z].m_parts[index]);
       // if (iter != m_layers[m_z].m_displayMap.end()) iter->second = display;
		if (iter != m_layers[m_z].m_displayMap.end()) iter->second.selected = selected;
    }
}
void CLYInterface::ToggleParts(bool select)
{
    map<int, dispinfo>& dmap = m_layers[m_z].m_displayMap;

    m_default_display = select;
    for (map<int, dispinfo>::iterator iter = dmap.begin(); iter != dmap.end(); ++iter) {
		iter->second.selected = select;
    }
}

#include <iostream>
#include <fstream>
#include <string>

void CLYInterface::LoadPartSelection(const wchar_t* file_str) {
	// open file 
	ifstream in_file( file_str, ifstream::in );
	
	string line;
	while (in_file.eof() == false) {
		getline(in_file,line);

		if (line.size() <= 0)
			continue;

		stringstream ss(line);
		
		int layer = -1;
		int num_selected_parts = -1;
		int selected_part;

		// for each line get its layer number
		ss >> layer >> num_selected_parts;
		
		if (layer < 0 || LayerCount() <= layer)
			continue;

		// first delected all parts in this layer
		for (map<int, dispinfo>::iterator iter = m_layers[layer].m_displayMap.begin(); iter != m_layers[layer].m_displayMap.end(); iter++) {
			iter->second.selected = false;
		}

		// for each part found in line we select it from m_displayMap
		while (--num_selected_parts >= 0) {
			selected_part = -1;
			ss >> selected_part;
			map<int, dispinfo>::iterator iter;
			if (selected_part < 0 || (iter=m_layers[layer].m_displayMap.find(selected_part)) == m_layers[layer].m_displayMap.end())
				continue;

			iter->second.selected = true;
		}
	
	}
}

void CLYInterface::GetMixedReconstructionPoints(vector<COLOREDPOINT>& pts)
{
	pts.clear();
	if(m_res == NULL) return;
	vector<map<node*, mixcolor>> lvlnodes;

	float factor = 2.0f/m_res->x_size(0);
    float yOffset = (factor * m_res->y_size(0))/2.0f;

	if (m_isNewPickPartsLocation == true) {
		// adjust position based on factor and yOffset
		m_pickParts.x = (m_pickParts.x + 1.0f)/factor;
		m_pickParts.y = (yOffset - m_pickParts.y) / factor;

		m_isNewPickPartsLocation = false;
	}

	GetMixedColorTree(lvlnodes);

	if(lvlnodes.empty()) return;
	
	for(map<node*, mixcolor>::const_iterator iter = lvlnodes[0].begin(); iter != lvlnodes[0].end(); ++iter)
	{
		layer1_data* nd = (layer1_data*)(iter->first->data);
		if (nd->z != 0) continue;
		pts.push_back(COLOREDPOINT(factor * nd->x - 1.0f, 
			yOffset - factor * nd->y,
			(float)nd->z,
			iter->second.getD3D()));
	}
}
/*
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
    int colorCount = (int)m_colorMap.size();
    double thresh = m_firstOnly ? -1.0 : m_responseThresh;
    double thresh2 = m_firstOnly ? -1.0 : m_responseThresh2;
    double thresh3 = m_firstOnly ? -1.0 : m_responseThresh3;
    int index = m_showAll ? -1 : m_partIndex;

    for (map<int, dispinfo>::iterator iter = m_layers[m_z].m_displayMap.begin(); 
            iter != m_layers[m_z].m_displayMap.end(); ++iter)
	{
        if (iter->second.selected) srcParts.push_back(iter->first);
		else continue;
		color = iter->second.color;
		m_indexCount = m_res->get_reconstruction_nodes(exporter, m_z, srcParts, thresh, thresh2, thresh3, index);
		if (index >= 0 && m_partIndex >= m_indexCount) m_partIndex = m_indexCount - 1;
		for (layer1_result_type_export::iterator_t iter = exporter.result.begin(); 
				iter != exporter.result.end(); ++iter) {
			layer1_data* nd = (layer1_data*)(iter->first)->data;

			if (nd->z != 0) continue;
			pt.x = factor * nd->x - 1.0f;
			pt.y = yOffset - factor * nd->y;
			pt.z = (float)nd->z;// - 0.001f;
			pt.color = color;
			pts.push_back(pt);
		}
		nodes.clear();
		srcParts.clear();
	}
    m_infoHypoNodes = exporter.attr_count;
    m_infoAllNodes = exporter.node_count;
}
*/
void CLYInterface::GetEllipses(const int ptsPerEll, vector<COLOREDPOINT>& pts)
{
	pts.clear();
    if (m_res == NULL) return;
	D3DCOLOR color;
    set<node*> nodes;
    vector<int> srcParts;
    int colorCount = (int)m_colorMap.size();
    double thresh = m_firstOnly ? -1.0 : m_responseThresh;
    double thresh2 = m_firstOnly ? -1.0 : m_responseThresh2;
    double thresh3 = m_firstOnly ? -1.0 : m_responseThresh3;
    double thresh4 = m_firstOnly ? -1.0 : m_responseThresh4;
    double thresh5 = m_firstOnly ? -1.0 : m_responseThresh5;
    int index = m_showAll ? -1 : m_partIndex;
	vector<pair<node*, D3DCOLOR>> nodescolor;

    for (map<int, dispinfo>::iterator iter = m_layers[m_z].m_displayMap.begin(); 
            iter != m_layers[m_z].m_displayMap.end(); ++iter)
	{
		if (iter->second.selected) srcParts.push_back(iter->first);
		else continue;

		color = iter->second.color;
		m_res->get_layer_nodes(nodes, m_z, srcParts, thresh, thresh2, thresh3, thresh4);

		int i = 0;
		
		for (set<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
            if (((*iter)->is_attr_set(HYPO_NODE_ATTR) && !m_showHyponodes) || (index != -1 && i != index)) {
				i++;
				continue;
			}	
			nodescolor.push_back(pair<node*, D3DCOLOR>(*iter, color));
			i++;
		}
		nodes.clear();
		srcParts.clear();
	}
	Ellipse2Pts(nodescolor, ptsPerEll, pts);
   /* m_partIndex = 0;
    m_indexCount = 1;
    m_infoHypoNodes = 0;
    m_infoAllNodes = (int)nodes.size();*/
}

void CLYInterface::Ellipse2Pts(const vector<pair<node*, D3DCOLOR>>& nodes, const int npts, vector<COLOREDPOINT>& pts)
{

	//pt.x = factor * nd->x - 1.0f;
	//pt.y = yOffset - factor * nd->y;
	////pt.z = /*(float)nd->z*/ - 0.001f;
	//pt.z = 0.0f;
	
	vector<dpoint2> circle;
	double factor = 2.0f/m_res->x_size(0);
	double yOffset = (factor * m_res->y_size(0))/2.0f;
	
    circle.reserve(npts);
	for (int i = 0; i < npts; ++i) {
		double t = (double)i * 2.0 * M_PI/(double)npts;
		circle.push_back(dpoint2(cos(t), sin(t)));
	}
	
	for(vector<pair<node*, D3DCOLOR>>::const_iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
		double cx, cy, a, b, angle;

		GetEllipse(cx, cy, a, b, angle, iter->first);
		a *= 2.5;
        b *= 2.5;

		for(vector<dpoint2>::const_iterator ic = circle.begin(); ic != circle.end(); ++ic) {
			COLOREDPOINT p;

			p.x = (float)(ic->x * a * cos(angle) + ic->y * b * sin(angle) + cx);
			p.y = (float)(-ic->y * b * cos(angle) + ic->x * a * sin(angle) + cy);
			p.x = factor * p.x - 1.0f;
			p.y = yOffset - factor * p.y;
			p.z = 0.0f;
			p.color = iter->second;
			pts.push_back(p);
		}
	}	
}

bool CLYInterface::GetEllipse(double& cx, double& cy, double& a, double& b, double& angle, node* n)
{
	layer_predicate pred(0);
	int to_prev = atom("toPrevLayer").get_index();
	set<node*> end_nodes;
	vector<dpoint2> pts;

    m_res->recurse_from_node(n, to_prev, pred, end_nodes);
	pts.reserve(end_nodes.size());
	for(set<node*>::iterator iter = end_nodes.begin(); iter != end_nodes.end(); ++iter) {
		layer1_data* nd = (layer1_data*)(*iter)->data;
		if (nd->z == 0)
			pts.push_back(dpoint2((double)nd->x, (double)nd->y));
	}
    if (!pts.empty()) {
        fit_ellipse(cx, cy, a, b, angle, pts);
		return true;
	}
	return false;
}

/*void CLYInterface::InhibitPoints(map<int,set<node*>>& pts, const map<int,set<node*>>& allpts, const double thresh)
{
	//gather all nodes in one set, assumes sets in pts are disjunct
	set<node*> allnodes;
	for(map<int, set<node*>>::const_iterator iter = allpts.begin(); iter != allpts.end(); ++iter)
	{
		for(set<node*>::const_iterator inode = iter->second.begin(); inode != iter->second.end(); ++inode)
		allnodes.insert(allnodes.end(), *inode);
	}
	pts.clear();
	//included factor(A) = max(|A intersect B| / |A|) for all |B| < |A|
	vector<double> included_factor;
	m_res->get_included_factor(included_factor, allnodes);
	//for each part set, keep those that are above a included threshold
	//for(map<int,set<node*>>::const_iterator iter = allpts.begin(); iter != allpts.end(); ++iter)
	//{
	//	pts.insert(pair<int,set<node*>>(iter->first, set<node*>()));
	//}
	//set<node*>::const_iterator ialln = allnodes.begin();
	//vector<double>::const_iterator iincf = included_factor.begin(); 
	//for(map<int,set<node*>>::const_iterator iter = allpts.begin(); iter != allpts.end(); ++iter)
	//{
	//	for(set<node*>::const_iterator inode = iter->second.begin(); inode != iter->second.end(); ++inode)
	//	{
	//		//if current node 
	//		while(ialln != allnodes.end() && *ialln != *inode)
	//		{
	//			++iincf;
	//			++ialln;
	//		}
	//		if(ialln == allnodes.end())
	//		{
	//			ialln = allnodes.begin();
	//			iincf = included_factor.begin();
	//			break;
	//		}
	//		if(*iincf < thresh)
	//		{
	//			map<int,set<node*>>::iterator el = pts.find(iter->first);
	//			el->second.insert(*inode);
	//		}
	//	}
	//}

	//v2
	for(map<int,set<node*>>::const_iterator iter = allpts.begin(); iter != allpts.end(); ++iter)
	{
		pts.insert(pair<int,set<node*>>(iter->first, set<node*>()));
	}
	map<node*, double> nfac;
	vector<double>::const_iterator ifac = included_factor.begin();
	for(set<node*>::const_iterator iter = allnodes.begin(); iter != allnodes.end(); ++iter)
	{
		nfac.insert(nfac.end(),pair<node*, double>(*iter, *ifac));
		++ifac;
	}
	for(map<int,set<node*>>::const_iterator iter = allpts.begin(); iter != allpts.end(); ++iter)
	{
		for(set<node*>::const_iterator inode = iter->second.begin(); inode != iter->second.end(); ++inode)
		{
			map<node*, double>::iterator f = nfac.find(*inode);
			if(f != nfac.end() && f->second < thresh)
			{
				map<int,set<node*>>::iterator el = pts.find(iter->first);
				el->second.insert(*inode);
			}
		}
	}
}*/

void CLYInterface::GetLayerPoints(vector<COLOREDPOINT>& pts)
{ 
	pts.clear();
    if (m_res == NULL) return;
	D3DCOLOR color;
    float factor = 2.0f/m_res->x_size(m_z);
    float yOffset = (factor * m_res->y_size(m_z))/2.0f;
    COLOREDPOINT pt;
    int colorCount = (int)m_colorMap.size();
    int index = m_showAll ? -1 : m_partIndex;
	vector<int> parts;
	vector<int> allParts;
	vector<D3DCOLOR> colors;

	for (map<int, dispinfo>::iterator iter = m_layers[m_z].m_displayMap.begin(); 
            iter != m_layers[m_z].m_displayMap.end(); ++iter)
	{
		allParts.push_back(iter->first);
		if (iter->second.selected)
		{
			parts.push_back(iter->first);	//parts are oredered
			colors.push_back(iter->second.color);
		}
		// set visable to false by default
		iter->second.visable = false;
	}
	map<int, set<node*>> nodes;

	if (m_isNewPickPartsLocation == true) {
		// adjust position based on factor and yOffset
		m_pickParts_non_reconstruct.x = (m_pickParts_non_reconstruct.x + 1.0f)/factor;
		m_pickParts_non_reconstruct.y = (yOffset - m_pickParts_non_reconstruct.y) / factor;

		m_isNewPickPartsLocation = false;
	}

	SelectLayerPoints(allParts, nodes);
	// get max number of parts per type
	m_indexCount = get_index_count(nodes, m_layers[m_z].m_displayMap);
    m_pSelectedNode = NULL;

	// set visablity of all parts that SelectLayerPoints returned
	for(map<int, set<node*>>::iterator ip = nodes.begin(); ip != nodes.end(); ++ip) {
		// set this part type to visable
		if (ip->second.size() > 0)
			m_layers[m_z].m_displayMap[ip->first].visable = true;
	}

	//foreach part
	vector<D3DCOLOR>::const_iterator ic = colors.begin();
	int attrib_nodes = 0;
	int node_count = 0;
	//for(map<int, set<node*>>::iterator ip = nodes.begin(); ip != nodes.end(); ++ip)
	for(vector<int>::iterator ip = parts.begin(); ip != parts.end(); ++ip)
	{
		int ind = *ip;
		int i = 0;
		for (set<node*>::iterator iter = nodes[ind].begin(); iter != nodes[ind].end(); ++iter) {
            if (((*iter)->is_attr_set(HYPO_NODE_ATTR) && !m_showHyponodes) 
				|| (index != -1 && i != index) )
			{
				++i;
				continue;
			}
			if((*iter)->is_attr_set(HYPO_NODE_ATTR)) attrib_nodes++;
			node_count++;
			layer1_data* nd = (layer1_data*)(*iter)->data;

			pt.x = factor * nd->x - 1.0f;
			pt.y = yOffset - factor * nd->y;
			//pt.z = /*(float)nd->z*/ - 0.001f;
			pt.z = 0.0f;
			pt.color = *ic;
			pts.push_back(pt);

            if (index != -1) {
                m_pSelectedNode = *iter;
            }

			++i;
		}
		++ic;
	} 
    m_infoHypoNodes = attrib_nodes;
    m_infoAllNodes = node_count;

}

void CLYInterface::SelectLayerPoints(const vector<int>& parts, map<int, set<node*>>& nodes, bool ignore_pick_parts)
{
	if(GetInhib() == false && m_firstOnly && GetHypoThreshold() == false)
	{
		//m_res->get_layer_nodes(nodes, m_z, parts, -1, -1, -1);
		if (ignore_pick_parts || m_pickParts_non_reconstruct.x < 0 || m_pickParts_non_reconstruct.y < 0)
			m_res->get_layer_nodes(nodes, m_z, parts, -1, -1, -1, 100);
		else
			m_res->get_layer_nodes(nodes, m_z, parts, -1, -1, -1, 100, &m_pickParts_non_reconstruct);
		return;
	}
	map<int, set<node*>> allnodes;
	double thresh = m_firstOnly ? -1.0 : m_responseThresh;
    double thresh2 = m_firstOnly ? -1.0 : m_responseThresh2;
    double thresh3 = m_firstOnly ? -1.0 : m_responseThresh3;
    double thresh4 = m_firstOnly ? -1.0 : m_responseThresh4;
    double thresh5 = m_firstOnly ? -1.0 : m_responseThresh5;
	//m_res->get_layer_nodes(allnodes, m_z, parts, thresh, thresh2, thresh3);
	if (ignore_pick_parts || m_pickParts_non_reconstruct.x < 0 || m_pickParts_non_reconstruct.y < 0)
		m_res->get_layer_nodes(allnodes, m_z, parts, thresh, thresh2, thresh3, thresh4);
	else
		m_res->get_layer_nodes(allnodes, m_z, parts, thresh, thresh2, thresh3, thresh4, &m_pickParts_non_reconstruct);
	
	//remove those with too high hyponode ratio
	if(GetHypoThreshold())
	{
		for(map<int, set<node*>>::iterator i1 = allnodes.begin(); i1 != allnodes.end(); )
		{
			vector<double> vals;
			for(set<node*>::iterator i2 = i1->second.begin(); i2 != i1->second.end(); )
			{
				map<node*, pair<int,int>>::iterator f = m_hyporatio.find(*i2);
				if(f != m_hyporatio.end())
				{
					double thr = (double)f->second.first / (double)f->second.second;
					double t2 = GetHypoThresholdVal();
					vals.push_back(thr);
					if( thr > t2)
					{
						i1->second.erase(i2++);
						if(i1->second.empty())
						{
							break;
						}
					}
					else ++i2;
				}
				else ++i2;
			}
			if(i1->second.empty())
			{
				allnodes.erase(i1++);
			}
			else ++i1;
		}
	}
    nodes = allnodes;	
}

void CLYInterface::SetPickPartsLoc(double x, double y, int z) {

	// only set if selected on same layer as we are currently showing
	if (0 == z) {
		if (m_showReconstruction) {
			m_pickParts = dpoint2(x, y);
			m_pickParts_non_reconstruct = dpoint2(-1, -1);
		} else {
			m_pickParts_non_reconstruct = dpoint2(x, y);
			m_pickParts = dpoint2(-1, -1);
		}
		m_isNewPickPartsLocation = true;
	} else if (z < 0) {
		m_pickParts = dpoint2(-1, -1);
		m_pickParts_non_reconstruct = dpoint2(-1, -1);
		m_isNewPickPartsLocation = false;	
	}
}

void CLYInterface::GetMixedTree(vector<vector<COLOREDPOINT>>& lines, vector<vector<COLOREDPOINT>>& pts)
{
	//trasform tree into d3dpoints
	pts.clear();
	lines.clear();

	pts.resize(m_z+1, vector<COLOREDPOINT>());
	lines.resize(m_z+1, vector<COLOREDPOINT>());
	if (m_res == NULL) return;

    float factor;
    float yOffset;
	float nfactor;
	float nyOffset;
	float xOffset;
	float nxOffset;
    COLOREDPOINT pt;
	COLOREDPOINT ptNeighbor;
	vector<map<node*, mixcolor>> lvlnodes;

	if (m_isNewPickPartsLocation == true) {
		// adjust position based on factor, yOffset and xOffset
		factor = 2.0f/m_res->x_size(00) * m_res->x_size(00) / m_res->x_size(0);
		yOffset = (factor * m_res->y_size(00))/2.0f;
		xOffset = -(float)m_res->x_size(00) / (float)m_res->x_size(0);
		
		m_pickParts.x = (m_pickParts.x - xOffset)/factor;
		m_pickParts.y = (yOffset - m_pickParts.y) / factor;

		m_isNewPickPartsLocation = false;
	}

	GetMixedColorTree(lvlnodes);
	//foreach lvl
	//for(vector<map<node*,mixcolor>>::const_iterator ilvl = lvlnodes.begin(); ilvl != lvlnodes.end(); ++ilvl)
	for(int i = m_z; i >= 0; --i)
	{
		factor = 2.0f/m_res->x_size(i) * m_res->x_size(i) / m_res->x_size(0);
		yOffset = (factor * m_res->y_size(i))/2.0f;
		xOffset = -(float)m_res->x_size(i) / (float)m_res->x_size(0);

		nfactor = 2.0f /m_res->x_size(i-1) * m_res->x_size(i-1) / m_res->x_size(0);
		nyOffset = (factor * m_res->y_size(i-1))/2.0f;
		nxOffset = -(float)m_res->x_size(i-1) / (float)m_res->x_size(0);

		//foreach node
		//for(map<node*,mixcolor>::const_iterator inode = ilvl->begin(); inode != ilvl->end(); ++inode)
		for(map<node*,mixcolor>::const_iterator inode = lvlnodes[i].begin(); inode != lvlnodes[i].end(); ++inode)
		{
			layer1_data* nd = (layer1_data*)(inode->first)->data;
			pt.x = factor * nd->x  + xOffset;
			pt.y = yOffset - factor * nd->y;
			pt.z = (float)nd->z * m_zspacing;// - 0.001f;
			D3DCOLOR color = inode->second.getD3D();
			pt.color = color;
			//if(i == m_z) pts.back().push_back(pt);
			pts.at(i).push_back(pt);
			//int e = 0;
            int name = atom("toPrevLayer");
			//foreach neighbour
			foreach_neighbor(inode->first, name, iter1)
			{
				//add to pts list
				layer1_data* nn = (layer1_data*)(neighbor_node(iter1))->data;
				ptNeighbor.x = nfactor * nn->x + nxOffset;
				ptNeighbor.y = nyOffset - nfactor * nn->y;
				ptNeighbor.z = (float)nn->z * m_zspacing - 0.001f;
				ptNeighbor.color = color;
				lines.at(i).push_back(pt);
				lines.at(i).push_back(ptNeighbor);
			}
		}
	}
}

void CLYInterface::SelectIndexedNodes(map<int,set<node*>>& nodes)
{
	for(map<int,set<node*>>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter)
	{
		int i = 0;
		for(set<node*>::iterator inode = iter->second.begin(); inode != iter->second.end();)
		{
			set<node*>::iterator nnode = inode;
			++nnode;
			if(i != m_partIndex)
			{
				iter->second.erase(inode);
			}
			inode = nnode;
			++i;
		}
	}
}

void CLYInterface::GetMixedColorTree(vector<map<node*, mixcolor>>& lvlnodes)
{
	//top level
	vector<int> allParts;
	vector<int> parts;	
	vector<D3DCOLOR> colors;
	for (map<int, dispinfo>::iterator iterd = m_layers[m_z].m_displayMap.begin(); 
            iterd != m_layers[m_z].m_displayMap.end(); ++iterd)
	{
		// by default use all parts when filtering with SelectLayerPoints
		allParts.push_back(iterd->first);
		if(iterd->second.selected)
		{
			parts.push_back(iterd->first);
			colors.push_back(iterd->second.color);
		}
		// set visable to false by default
		iterd->second.visable = false;
	}
	
	map<int,set<node*>> topnodes;
	
	// if filtering by mouse pick then we need to find all parts on m_z layer that are constructed from selected point
	if (m_pickParts.x >= 0 && m_pickParts.y >= 0) {

		// first create mapping from layer0 parts to m_z layer parts so we can find proper part based on position of layer0 parts
		// we create this map only when is empty (i.e. does not exists or some property changed)
		if (layer0_mapping.empty()) {

			// get nodes for top layer only
			SelectLayerPoints(allParts, topnodes, false);

			// fill lvlnodes with nodes from topnodes on all layers
			BuildOtherLayers(topnodes, lvlnodes, parts, colors, m_infoHypoNodes, m_infoAllNodes);

			// only save lowest layer and reject others
			layer0_mapping = lvlnodes[0];
		}
		// clear
		topnodes.clear();
		lvlnodes.clear();

		for(vector<int>::const_iterator iter = allParts.begin(); iter != allParts.end(); ++iter) {
			topnodes.insert(pair<int, set<node*> >(*iter, set<node*>()));
		}

		// find all parts at location picked by mouse
		for (map<node*, mixcolor>::iterator iter_ly0 = layer0_mapping.begin();
				iter_ly0 != layer0_mapping.end(); iter_ly0++) {
			
			layer1_data* nd = (layer1_data*)(iter_ly0->first)->data;
			if (dpoint2::distance2(m_pickParts, dpoint2(nd->x,nd->y)) < 10) {
				// from ones that are closest get all parts that include this location as subpart
				// (we already have this values saved)
				for (set<node*>::iterator iter_highlayer = iter_ly0->second.backConnection.begin();
						iter_highlayer != iter_ly0->second.backConnection.end(); iter_highlayer++) {
					// sort all this parts into topnodes
					layer1_data* n_top = (layer1_data*)(*iter_highlayer)->data;						
					
					map<int,set<node*> >::iterator el = topnodes.find(n_top->m);
					if(el != topnodes.end()) {						
						el->second.insert(*iter_highlayer);
					}
				}
			}
		}
	} else {
		// get nodes for top layer only
		SelectLayerPoints(allParts, topnodes);
	}

	// set visablitiy of part types to all based on topnodes
	for (map<int,set<node*>>::const_iterator iter = topnodes.begin(); iter != topnodes.end(); ++iter) 
	{
		// set this part type to visable if any parts found for this type
		if (iter->second.size() > 0)
			m_layers[m_z].m_displayMap[iter->first].visable = true;
	}

	// get max number of parts for one type
	m_indexCount = get_index_count(topnodes, m_layers[m_z].m_displayMap);

	// if not showing all then select only indexed one
	if(!m_showAll) SelectIndexedNodes(topnodes);

	// fill lvlnodes with nodes from topnodes on all layers
	BuildOtherLayers(topnodes, lvlnodes, parts, colors, m_infoHypoNodes, m_infoAllNodes);
	
}

void CLYInterface::BuildOtherLayers(map<int,set<node*>>& topnodes, vector<map<node*, mixcolor>>& lvlnodes, vector<int>& parts, vector<D3DCOLOR>& colors, int& attr_count, int& node_count) 
{
    int ename = atom("toPrevLayer");

	lvlnodes.clear();
	lvlnodes.resize(m_z+1, map<node*, mixcolor>());

	vector<D3DCOLOR>::const_iterator icolor = colors.begin();
	
	attr_count = 0;
	//for(map<int,set<node*>>::const_iterator iter = topnodes.begin(); iter != topnodes.end(); ++iter)
	for(vector<int>::const_iterator iter = parts.begin(); iter != parts.end(); ++iter)
	{
		int ind = *iter;
		for(set<node*>::const_iterator inode = topnodes[ind].begin(); inode != topnodes[ind].end(); ++inode)
		{
			//counting certain nodes
			if ((*inode)->is_attr_set(HYPO_NODE_ATTR)) ++attr_count;
            if (!m_showAll) m_pSelectedNode = *inode;

			lvlnodes.back().insert(lvlnodes.back().end(), pair<node*, mixcolor>(*inode, mixcolor(*icolor, *inode)));
		}
		++icolor;
	}
	
	//build rest of the tree, top-down
	for(vector<map<node*, mixcolor>>::reverse_iterator ilvl = lvlnodes.rbegin(); ilvl != lvlnodes.rend(); ilvl++)
	{
		//for each node, add its neigh.
		for(map<node*,mixcolor>::iterator inode = ilvl->begin(); inode != ilvl->end(); ++inode)
		{
			//node* nd = inode->first;
			foreach_neighbor(inode->first, ename, ineigh)
			{
				layer1_data* n = (layer1_data*)(neighbor_node(ineigh))->data;
				if ((neighbor_node(ineigh))->is_attr_set(HYPO_NODE_ATTR)) ++attr_count;
				map<node*,mixcolor>::iterator nn = lvlnodes[n->z].find(neighbor_node(ineigh));
				if(nn != lvlnodes[n->z].end())
				{
					nn->second += inode->second;					
				}
				else
				{
					lvlnodes[n->z].insert(pair<node*, mixcolor>(neighbor_node(ineigh), mixcolor(inode->second.getD3D())));
				}
				map<node*,mixcolor>::iterator ff = lvlnodes[n->z].find(neighbor_node(ineigh));
				ff->second.backConnection.insert(inode->second.backConnection.begin(),inode->second.backConnection.end());
			}
		}
	}
	node_count = 0;
	for(vector<map<node*, mixcolor>>::const_iterator iter = lvlnodes.begin(); iter != lvlnodes.end(); ++iter)
	{
		node_count += iter->size();
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
	if(layer <= GetLayer() || layer == 0)
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
    //if (m_showReconstruction) GetReconstructionPoints(pts);
	if (m_showReconstruction) GetMixedReconstructionPoints(pts);
    else GetLayerPoints(pts);
}

void CLYInterface::GetImage(void*& mem, unsigned int& size)
{
	CStopWatch clock1,clock2, clock3;
    if (m_res == NULL) size = 0;
    else {
		clock1.startTimer();
        img* image = m_res->get_image_reconstructed(0, 0, vector<int>());
		clock1.stopTimer();
		
		clock2.startTimer();
        size = image->save(mem);
		clock2.stopTimer();

		clock3.startTimer();
		delete image;
		clock3.stopTimer();

		cout << "time reconstructed: " << clock1.getElapsedTime() << ", time copy " << clock2.getElapsedTime() << ", time delete " << clock3.getElapsedTime() << endl;
    }
}

void CLYInterface::GetPartImage(void*& mem, int& width, int& height, int layer, int part, int max_width, int max_height)
{
    mem = NULL;
    if (m_library == NULL) return;

    img im = m_library->get_image(layer, part);

	width = im.width;
    height = im.height;

	if (max_width > 0 && width > max_width) {
		height = height * (max_width / (float)width);
		width = max_width;
	}
	if (max_height > 0 && height >= max_height) {
		width = width * (max_height / (float)height);
		height = max_height;
	}
    
	if (width != im.width || height != im.height) {
		im = im.get_resized(width, height);
		width = im.width;
		height = im.height;
	}
    im.get_bits(mem, width, height, img::BIT_TYPE_32);
}

void response_to_string(strstream& os, layer1_data* nd, int response)
{
    try {
        os << nd->r(response);
    } catch (libhop_exception&) {
        os << '/';
    }
}

unsigned CLYInterface::GetThresholdInfo(char*& buf)
{
    if (m_res == NULL || m_pSelectedNode == NULL) {
        buf = NULL;
        return 0;
    }
    strstream os;
    layer1_data* nd = (layer1_data*)m_pSelectedNode->data;

    os << "  ";
    os << "x: " << nd->x << ",  ";
    os << "y: " << nd->y << ",  ";
    os << "r-response: ";
    response_to_string(os, nd, R_RESPONSE);
    os << ",  ";
    os << "g-response: ";
    response_to_string(os, nd, G_RESPONSE);
    os << ",  ";
    os << "rr-response: ";
    response_to_string(os, nd, RR_RESPONSE);
    os << ",  ";
    os << "s-response: ";
    response_to_string(os, nd, S_RESPONSE);
    os << ",  ";
    os << "x-response: ";
    response_to_string(os, nd, X_RESPONSE);
    os << '\0';

    char* str = os.str();
    size_t size = os.pcount();

    buf = (char*)malloc(sizeof(char) * size);
    memcpy(buf, str, sizeof(char) * size);
    return size;
}

unsigned CLYInterface::GetInfo(char*& buf)
{
    if (m_res == NULL) {
        buf = NULL;
        return 0;
    }

    strstream os;

    for (int l = 0; l <= m_res->max_layer_index(); ++l) {
        os << "Layer " << l << ": " << m_res->x_size(l) << " x " << m_res->y_size(l) << endl;
    }
    os << '\0';

    char* str = os.str();
    size_t size = os.pcount();

    buf = (char*)malloc(sizeof(char) * size);
    memcpy(buf, str, sizeof(char) * size);
    return size;
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
/*
void CLYInterface::GetPartsRect(vector<rectangle2<float>>& partsRect)
{
	if (m_res == NULL) return;
	partsRect.clear();
	float xmin;
	float xmax;
	float ymin;
	float ymax;

    layer1_result_type_export exporter(HYPO_NODE_ATTR);
    set<node*> nodes;
    vector<int> srcParts;
    float factor = 2.0f/m_res->x_size(0);
    float yOffset = (factor * m_res->y_size(0))/2.0f;
    COLOREDPOINT pt;
    double thresh = m_firstOnly ? -1.0 : m_responseThresh;
    double thresh2 = m_firstOnly ? -1.0 : m_responseThresh2;
    double thresh3 = m_firstOnly ? -1.0 : m_responseThresh3;
    int index = m_showAll ? -1 : m_partIndex;

    for (map<int, dispinfo>::iterator iter = m_layers[m_z].m_displayMap.begin(); 
            iter != m_layers[m_z].m_displayMap.end(); ++iter)
	{
        if (iter->second.selected) srcParts.push_back(iter->first);
		else continue;
		xmin = FLT_MAX;
		xmax = -FLT_MAX;
		ymin = FLT_MAX;
		ymax = -FLT_MAX;
	    m_indexCount = m_res->get_reconstruction_nodes(exporter, m_z, srcParts, thresh, thresh2, thresh3, index);
		if (index >= 0 && m_partIndex >= m_indexCount) m_partIndex = m_indexCount - 1;
		for (layer1_result_type_export::iterator_t iter = exporter.result.begin(); 
				iter != exporter.result.end(); ++iter) {
			layer1_data* nd = (layer1_data*)(iter->first)->data;

			if (nd->z != 0) continue;
			pt.x = factor * nd->x - 1.0f;
			// orig // pt.y = yOffset - factor * nd->y;
			pt.y = (yOffset - factor * nd->y);
			//pt.x = nd->x;
			//pt.y = nd->y;
			if(xmax < pt.x) xmax = pt.x;
			if(ymax < pt.y) ymax = pt.y;
			if(xmin > pt.x) xmin = pt.x;
			if(ymin > pt.y) ymin = pt.y;
		}
		//partsRect.push_back(rectangle2<float>(xmin, ymin, xmax, ymax));
		partsRect.push_back(rectangle2<float>(xmin, ymax, xmax, ymin));
		nodes.clear();
		srcParts.clear();
	}
    m_infoHypoNodes = exporter.attr_count;
    m_infoAllNodes = exporter.node_count;
}
*/
void CLYInterface::GetLayerZSpacing(float& zspacing)
{
	zspacing = m_zspacing;
}

void CLYInterface::SetLayer(int z)
{
	if(m_res == NULL) return;
    if (z < 0 || z > m_res->max_layer_index()) m_z = m_res->max_layer_index();
    else m_z = z;
	// reset mouse picking
	SetPickPartsLoc(-1,-1,-1);
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
		if(iter->second.selected)
		{
			iter->second.color = color;
		}
	}
}

void CLYInterface::SetColorDefault()
{
	map<int, D3DCOLOR>::iterator citer = m_colorMap.begin();
	for(int cur_z = 0; cur_z < LayerCount(); ++cur_z)
	{
		for (map<int, dispinfo>::iterator iter = m_layers[cur_z].m_displayMap.begin(); 
				iter != m_layers[cur_z].m_displayMap.end(); ++iter)
		{			
			iter->second.color = citer->second;
			citer++;
			if(citer == m_colorMap.end()) citer = m_colorMap.begin();
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
			miter->second.selected = false;
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
				miter->second.selected = false;
			}
		}
	}
}

void CLYInterface::GetColorMap(map<int, COLORREF>& cmap)
{
	cmap.clear();
	for(map<int, D3DCOLOR>::const_iterator iter = m_colorMap.begin();
		 iter != m_colorMap.end(); ++iter)
	{
		cmap.insert(cmap.end(), pair<int, COLORREF>(iter->first, RGB2BGR(iter->second)));
	}
}

void CLYInterface::SetHypoThreshold(const bool t)
{
	m_hyporatioThres = t;
	if(m_hyporatioThres == true && m_hyporatio.empty())
	{
		GetHypoRatio();
	}
}
