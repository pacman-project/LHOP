/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// hop-interface.cpp : Defines the exported functions for the DLL application.
//

#ifdef WIN32
#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN // Exclude rarely-used stuff from Windows headers
#include <windows.h>
#endif

#include <stdio.h>

#include "hop.h"
#include "hopstorage.h"

#include "layers/layer_n_creators.h"
#include "layers/layer_1_creators.h"
#include "layers/layer_1_result.h"

#define DISPOSE_HOP_STRUCT(type) \
    if (cptr) {\
        if (--cptr->count == 0) {\
            delete (type)data();\
            delete cptr;\
        }\
    }

#define DISPOSE_HOP_STRUCT_EXPR(expr) \
    if (cptr) {\
        if (--cptr->count == 0) {\
            expr;\
            delete cptr;\
        }\
    }


// hop classes
///////////////////////////////////////////////////////////////////////////////

// hop_struct
///////////////

hop_struct::hop_struct(void* data)
{
    if (!data) cptr = 0; else cptr = new counted_ptr(data);
}

void hop_struct::copy(const hop_struct& hd)
{
    if (cptr == hd.cptr) return;
    dispose();
    cptr = hd.cptr;
    if (cptr) ++cptr->count;
}

hop_blob hop_struct::to_blob()
{
    if (!data()) return hop_blob(0, 0);

    ostrstreamer os(Z_DEFAULT_COMPRESSION);

    write_to_stream(&os);

    string s = os.rdbuf()->str();
    int size = (int)s.size();

    void* blobdata = malloc(sizeof(char)*size);

    memcpy(blobdata, s.data(), sizeof(char)*size);
    return hop_blob(blobdata, size);
}

void hop_struct::to_file(const char* fname, int flags)
{
    if (data()) 
		((streamable*)data())->save(fname, Z_DEFAULT_COMPRESSION);
}


void hop_struct::write_to_stream(void* os)
{
    throw hop_exception("Method 'write_to_stream' is not implemented for this class.");
}

void hop_struct::read_from_stream(void* is)
{
    throw hop_exception("Method 'read_from_stream' is not implemented for this class.");
}

// hop_blob - arbitrary data (pointer to malloc-ated memory)
//////////////////////////////////////////////////////////////

hop_blob::hop_blob(void* data, int s) : hop_struct(data), size(s) { }

hop_blob::~hop_blob()
{
    DISPOSE_HOP_STRUCT_EXPR(free(data()));
}

void hop_blob::dispose()
{
    DISPOSE_HOP_STRUCT_EXPR(free(data()));
}

// hop_streamable - arbitrary streamable data
/////////////////////////////////////////////////////////////

hop_streamable::hop_streamable(hop_blob& blob) {
	stringbuf buf;
	buf.pubsetbuf((char*)blob.data(), blob.get_size());
	istrstreamer is(&buf);
	streamable* obj = is.read_structure();

	if (obj) cptr = new counted_ptr(obj);
}

void hop_streamable::dispose() {
		DISPOSE_HOP_STRUCT(streamable*);
}

void hop_streamable::write_to_stream(void* os) {
	((ostreamer*)os)->write_structure((streamable*)data());
}


// hop_image - wrapper for img class
//////////////////////////////////////

hop_image::hop_image(void* pimg) : hop_streamable(pimg) { }

hop_image::hop_image(hop_blob& blob) : hop_streamable(0)
{
    string str((char*)blob.data(), blob.get_size());
    stringbuf buf(str);
    istrstreamer is(&buf);
    img* im = new img();

    im->read_from_stream(is);
    cptr = new counted_ptr(im);
}

hop_image::~hop_image()
{
    DISPOSE_HOP_STRUCT(img*);
}

void hop_image::dispose()
{
    DISPOSE_HOP_STRUCT(img*);
}

int hop_image::get_width()
{
    if (data()) return ((img*)data())->width; else return -1;
}

int hop_image::get_height()
{
    if (data()) return ((img*)data())->height; else return -1;
}

double* hop_image::get_data_ptr()
{
    if (!data()) return 0; else return ((img*)data())->ptr(0, 0);
}
/*
void hop_image::write_to_stream(void* os)
{
    img* im = (img*)data();

    im->write_to_stream(*((ostreamer*)os));
}*/

hop_image hop_image::from_bytes_rgb32(int width, int height, void* bytes) {
	// create img class from byte
	unsigned long colorMask[] = {0xFF0000,0xFF00, 0xFF};
	img* im = new img();
	img::from_bits_grayscale(*im, bytes, width, height, colorMask, img::BIT_TYPE_32);

	return hop_image(im);
}

//void hop_image::read_from_stream(void* is)
//{
//}


// hop_library - wrapper for part_lib class
/////////////////////////////////////////////

hop_library_part hop_library::invalid_part(-1, -1);

hop_library::hop_library(void* plibrary) : hop_streamable(plibrary) { }

hop_library::hop_library(hop_blob& blob) : hop_streamable(blob)
{
	// !!!!! WARNING !!!!!!
	// anything added here MUST also be added to read_library(const string& name, part_lib*& library) constructor in layer_1.cpp
	// !!!!! WARNING !!!!!!    
	void* lib = data();
    if (lib)
		((part_lib*)lib)->update_var_fields();
}

hop_library::~hop_library()
{
    DISPOSE_HOP_STRUCT(part_lib*)
}

void hop_library::dispose()
{
    DISPOSE_HOP_STRUCT(part_lib*);
}

int hop_library::get_layer_count()
{
    return ((part_lib*)data())->layer_count;
}

int hop_library::get_part_count(int l)
{
    part_lib* library = (part_lib*)data();

    if (!library || l < 0 || l >= library->layer_count) return 0;
    else return (int)library->parts[l].size();
}

hop_library_part hop_library::get_part(int l, int i)
{
    part_lib* library = (part_lib*)data();

    if (!library || l < 0 || l >= library->layer_count || i < 0 || i >= (int)library->parts[l].size())
        return invalid_part;
    else {
        part_data* pd = (part_data*)(library->parts[l][i])->data;

        return hop_library_part(pd->layer, pd->type);
    }
}
void hop_library::load_svm_models(char* filename) {
	part_lib* library = (part_lib*)data();
	
	if (library == nullptr)
		cout << "Unable to load svm models storage. Invalid library data" << endl;
	else {
		//library->svm_models_storage = cvOpenFileStorage(filename, nullptr, CV_STORAGE_READ);
		cv::FileStorage fs(filename, cv::FileStorage::READ);
        if (fs.isOpened()) {
			cout << "reading cv storable data" << endl;
            library->read_cv_storable_data(*fs);
            fs.release();
        }
	}
}

void hop_library::print_layer_info(int layer) {
	part_lib* library = (part_lib*)data();
	
	if (library == nullptr)
		cout << "Unable to display layer info. Invalid library data" << endl;
	else
		library->display_layer_info(layer);
}
/*
void hop_library::write_to_stream(void* os)
{
    ((ostreamer*)os)->write_structure((part_lib*)data());
}
*/
// hop_result - wrapper for layer1_data class
///////////////////////////////////////////////

hop_node hop_result::invalid_node(-1, -1, -1, 0, 0.0, 0.0, 0.0, 0.0);

hop_result::hop_result(void* presult) : hop_streamable(presult) { }

hop_result::hop_result(hop_blob& blob) : hop_streamable(blob) { }

hop_result::~hop_result()
{
    DISPOSE_HOP_STRUCT(layer1_result*);
}

void hop_result::dispose()
{
    DISPOSE_HOP_STRUCT(layer1_result*);
}

int hop_result::get_original_image_width() const 
{
	return ((layer1_result*)data())->original_width;
}


int hop_result::get_border() const 
{
	return ((layer1_result*)data())->border;
}

int hop_result::get_height(int l) const
{
    return ((layer1_result*)data())->y_size(l);
}

int hop_result::get_width(int l) const
{
    return ((layer1_result*)data())->x_size(l);
}


int hop_result::get_layer_count() const
{
    return ((layer1_result*)data())->layer_count();
}

int hop_result::get_node_count(int l) const {
    layer1_result* result = (layer1_result*)data();

    if (!result || l < 0 || l > result->max_layer_index())
        return -1;

    return  (int)result->shape_nodes[l].size();
}

hop_node hop_result::get_node_at(int l, int x, int y) const
{
    layer1_result* result = (layer1_result*)data();
    
    if (result && l >= 0 && l <= result->max_layer_index()) {
        if (!result->grid(l)) result->init_grid(l);

        node* n = result->node_at(x, y, l);

        if (n) {
            layer1_data* nd = (layer1_data*)n->data;

            return hop_node(nd->z, nd->x, nd->y, nd->m, nd->r(R_RESPONSE), nd->r(G_RESPONSE), nd->r(RR_RESPONSE), nd->r(S_RESPONSE));
        }
    }
    return invalid_node;
}

hop_node hop_result::get_node(int l, int i) const
{
    layer1_result* result = (layer1_result*)data();

    if (!result || l < 0 || l > result->max_layer_index() || i < 0 || i >= (int)result->shape_nodes[l].size()) 
        return invalid_node;
    else {
        node* n = result->shape_nodes[l][i];
        layer1_data* nd = (layer1_data*)n->data;

        return hop_node(l, nd->x, nd->y, nd->m, nd->r(R_RESPONSE), nd->r(G_RESPONSE), nd->r(RR_RESPONSE), nd->r(S_RESPONSE));
    }
}

hop_nodes hop_result::get_child_nodes(const hop_node& pn) const
{
    layer1_result* result = (layer1_result*)data();
    hop_nodes children;
    
    if (!result || pn == invalid_node) return children;

    if (!result->grid(pn.layer)) result->init_grid(pn.layer);

    node* n = result->node_at(pn.x, pn.y, pn.layer);
    
    if (!n) return children;

    layer1_data* nd;
    int edgename = atom("toPrevLayer").get_index();

    while (n && (nd = ((layer1_data*)n->data))->m != pn.type) 
        n = nd->next;

    if (!n) return children;

    vector<hop_node> cvector;
    
    foreach_neighbor(n, edgename, iter) {
        layer1_data* nnd = (layer1_data*)neighbor_node_data(iter);
        cvector.push_back(hop_node(nnd->z, nnd->x, nnd->y, nnd->m, nnd->r(R_RESPONSE), nnd->r(G_RESPONSE), nnd->r(RR_RESPONSE), nnd->r(S_RESPONSE)));
    }
    if (cvector.empty()) 
        return children;
    children.items = new hop_node[cvector.size()];
    children.size = (int)cvector.size();
    memcpy(children.items, &cvector.at(0), sizeof(hop_node)*cvector.size());
    return children;
}

hop_indices hop_result::get_child_nodes(int l, int i) const
{
    typedef map<node*, int> map_t;

    layer1_result* result = (layer1_result*)data();
    hop_indices children;
    
    if (!result || l < 1 || l > result->max_layer_index() || i < 0 || i >= (int)result->shape_nodes[l].size())
        return children;
    
    map_t indexmap;
    node* n = result->shape_nodes[l][i];
    int edgename = atom("toPrevLayer").get_index();
    vector<int> cvector;
    vector<node*>& snodes = result->shape_nodes[l - 1];
    
    for (int i = 0; i < (int)snodes.size(); ++i) {
        indexmap.insert(map_t::value_type(snodes[i], i));
    }
    foreach_neighbor(n, edgename, iter) {
        map_t::iterator miter = indexmap.find(neighbor_node(iter));

        if (miter != indexmap.end()) cvector.push_back(miter->second);
    }
    if (cvector.empty()) 
        return children;
    children.items = new int[cvector.size()];
    children.size = (int)cvector.size();
    memcpy(children.items, &cvector.at(0), sizeof(int)*cvector.size());
    return children;
}
/*
void hop_result::write_to_stream(void* os)
{
    ((ostreamer*)os)->write_structure((layer1_result*)data());
}
*/
void hop_result::to_file(const char* fname, int flags)
{
    if (!data()) 
	return;

    if (flags & HOP_FORMAT_PROTOBUF) {
        _hop_store_result_protobuf(this, fname, (flags & HOP_FORMAT_COMPRESS) != 0);
    } else {
        _hop_store_result_blob(this, fname, (flags & HOP_FORMAT_COMPRESS) != 0);
    }
}

// hop functions
///////////////////////////////////////////////////////////////////////////////

HOPINTERFACE_API void hop_force_unbuffer_stdout() {
	setvbuf (stdout, nullptr, _IONBF, BUFSIZ);
}

HOPINTERFACE_API void hop_redirect_stdout(const char* filename) {
	freopen( filename, "a", stdout );
	cout << "stdout redirected to " << filename << endl;
	fflush(stdout);
}
HOPINTERFACE_API void hop_redirect_stderr(const char* filename) {
	freopen( filename, "a", stderr );
	cout << "stderr redirected to " << filename << endl;
	fflush(stdout);
}
HOPINTERFACE_API hop_image hop_read_image(const char* fname)
{ 
    img* i = new img(string(fname));
    if (i->empty()) {
       throw hop_exception("Unable to read image");
    }
    return hop_image(i);
}

HOPINTERFACE_API hop_image hop_image_from_array(int width, int height, double* arr)
{
    img* rimg = new img(width, height, true);
    memcpy(rimg->ptr(0, 0), arr, width * height * sizeof(double));
    return hop_image(rimg);
}

HOPINTERFACE_API hop_library hop_read_library(const char* fname)
{
    return hop_library(streamable::read(fname));
}

HOPINTERFACE_API hop_result hop_read_result(const char* fname)
{
    return hop_result(streamable::read(fname));
}

HOPINTERFACE_API int hop_inference(hop_result*& result, hop_image& image, const std::list<irectangle2> mask_regions, const char* params)
{
    config_dictionary cfg;
	cfg.from_string(params);
	cfg.update_namespace_references();

	config_dictionary cfgi;
	cfgi.from_namespace_priority(cfg, 1, "inference.ly1");

    vector<layer1_result*> results;    
	vector<img> srcimg(1);
	srcimg[0] = *((img*)image.data());
	
	vector<irectangle2> mask_regions_vector(mask_regions.begin(), mask_regions.end());
	
	lay1create::create_layer1(results, cfgi, srcimg, mask_regions_vector);

    result = new hop_result[results.size()];

    for (int i = 0; i < (int)results.size(); ++i) 
        result[i] = hop_result(results[i]);
    return (int)results.size();

}

HOPINTERFACE_API int hop_inference(hop_result*& result, hop_image& image, const char* params)
{
    config_dictionary cfg;
	cfg.from_string(params);
	cfg.update_namespace_references();

	config_dictionary cfgi;
	cfgi.from_namespace_priority(cfg, 1, "inference.ly1");

    vector<layer1_result*> results;    
	vector<img> srcimg(1);
	srcimg[0] = *((img*)image.data());
	
	lay1create::create_layer1(results, cfgi, srcimg, img());

    result = new hop_result[results.size()];

    for (int i = 0; i < (int)results.size(); ++i) 
        result[i] = hop_result(results[i]);
    return (int)results.size();

}

HOPINTERFACE_API void hop_inference(hop_result& result, hop_library& lib, const char* params)
{
	vector<layer1_result*> vresult;
    config_dictionary cfg;
    
	cfg.from_string(params);
	cfg.update_namespace_references();

	config_dictionary cfgi;
	cfgi.from_namespace_priority(cfg, 1, "inference");

    vresult.push_back((layer1_result*)result.data());
    layncreate::create_layern(vresult, cfgi, (part_lib*)lib.data());
}

HOPINTERFACE_API void hop_save_image(hop_image& im, const char* fname)
{
    img* i = (img*)im.data();
    if (i) i->save_normalized(fname);
}

HOPINTERFACE_API hop_image hop_display(hop_result& result, const char* params)
{
    config_dictionary cfg;
    vector<int> parts, icolors;
    vector<color> colors;
    color defcolor, uncovered;
    bool paintuncov;

    cfg.from_string(params);
	cfg.update_namespace_references();

    cfg.get_value(parts, "parts");
    cfg.get_value(icolors, "part_colors");
    color::to_color_vector(colors, icolors, 0);

    icolors.clear();
    cfg.get_value(icolors, "default_color");
    color::to_color(defcolor, icolors, 255);

    icolors.clear();
    cfg.get_value(icolors, "uncovered_color");
    paintuncov = (icolors.size() != 0);
    color::to_color(uncovered, icolors, 20);

    layer1_result* res = (layer1_result*)result.data();
    img* im = res->get_image(cfg.get_value_int("layer", -1),
        cfg.get_value_int("end_layer", 0),
        paintuncov, uncovered, defcolor, parts, colors); 
    return hop_image(im);
}

HOPINTERFACE_API const char* hop_time_stamp() {
	return "( with DLL timestamp: " __DATE__ " / " __TIME__ ")";
}

HOPINTERFACE_API void hop_enable_opencl(const char* c_kernel_paths, const char* c_used_devices) {
#ifdef OPENCL
	const static char* seperator = ",;";

	vector<string> kernel_paths_list;
	vector<string> used_devices_list;

	string kernel_paths(c_kernel_paths);

	size_t last_pos, pos;

	pos = last_pos= 0;
	while ((pos = kernel_paths.find_first_of(seperator, last_pos)) != string::npos) {
		pos++;
		kernel_paths_list.push_back(kernel_paths.substr(last_pos+1, pos - last_pos));
		last_pos = pos;
	}
	kernel_paths_list.push_back(kernel_paths.substr(last_pos));

	string used_devices(c_used_devices);

	pos = last_pos= 0;
	while ((pos = used_devices.find_first_of(seperator, last_pos)) != string::npos) {
		pos++;
		used_devices_list.push_back(used_devices.substr(last_pos, pos - last_pos));
		last_pos = pos;
	}	
	used_devices_list.push_back(used_devices.substr(last_pos));

	OpenCL::initializeContexts(kernel_paths_list, used_devices_list);
#else
	throw custom_libhop_exception(opencl_exception, "Unable to enable OpenCL. Binary not build with OpenCL support");
#endif
}

HOPINTERFACE_API int hop_create_histograms(std::vector<hop_histogram_descriptor>*& result, hop_result& res_wraper, std::list<irectangle2> export_windows, hop_histogram_generator& hoc_generator_wraper) {
	
	using namespace hoc;

	layer1_result* res = (layer1_result*)res_wraper.data();
	hoc_histogram_generator* hoc_generator = (hoc_histogram_generator*)hoc_generator_wraper.data();

	result = (std::vector<hop_histogram_descriptor>*)hoc_generator->generate_descriptors(res, export_windows);
    
    return result->size();
}

HOPINTERFACE_API hop_result hop_merge_scales(std::vector<hop_result>& result, const char* params) {
	layer1_result* merged_res = (layer1_result*)((layer1_result*)result.front().data())->get_copy_s();
	
	for (std::vector<hop_result>::iterator iter = result.begin(); iter != result.end(); iter++) {
		layer1_result* res = (layer1_result*)iter->data();
		merged_res->merge(res, res->border);
	}

	//printf("number of nodes on merged result: last ly %d, last-1 ly %d\n", merged_res->get_nnodes(merged_res->max_layer_index()), merged_res->get_nnodes(merged_res->max_layer_index()-1));
	//fflush(stdout);

	return hop_result(merged_res);
}

HOPINTERFACE_API int hop_get_detections(std::list<hop_detection>*& result, hop_result& res_wraper, hop_library& lib_wraper, const int only_category /*= -1*/) {
	layer1_result* res = (layer1_result*)res_wraper.data();
	part_lib* lib = (part_lib*)lib_wraper.data();

	int last_layer = res->max_layer_index() - 1;
	
	int border = res->border;
	double scale_factor = res->original_width / ((double)(res->x_size(0) - 2*border));	
	double contraction_factor = (double)res->x_size(0) / (double)res->x_size(last_layer);

	//printf("border %d, scale factor %f, contraction factor %f, last layer %d\n", border, scale_factor, contraction_factor, last_layer);

	result = new std::list<hop_detection>();
	
	//vector<irectangle2> lib_part_boxes_predictions = lib->get_predicted_boxes(last_layer);
	
	vector<node*> lnodes;
	res->get_layer_nodes(lnodes, last_layer);

	map<node*, vector<node*> > cluster_nodes;

	cluster_detections_ms(cluster_nodes, res,  list<node*>(lnodes.begin(), lnodes.end()), 2.0, 100);

	// iterate through the list of all the nodes on the last layer
	//for (vector<node*>::iterator iter = lnodes.begin(); iter != lnodes.end(); iter++) {	
		//node* n = *iter;		
        //layer1_data* nd = (layer1_data*)n->data;

		// get bounding box for this detection
		//irectangle2 bbox = res->get_box_with_cached_link(n);
		//irectangle2 bbox = lib_part_boxes_predictions[nd->m];
	for (map<node*, vector<node*> >::iterator iter = cluster_nodes.begin(); iter != cluster_nodes.end(); iter++) {
		node* n = iter->first;

		layer1_data* nd = (layer1_data*)n->data;

		// skip invalid categories if requested so
		if (only_category >= 0 && nd->m != only_category)
			continue;

		irectangle2 bbox = res->get_box_with_cached_link(n);

		for (vector<node*>::iterator iter_cluster = iter->second.begin(); iter_cluster != iter->second.end(); iter_cluster++) {
			 irectangle2 bbox_cluster = res->get_box_with_cached_link(*iter_cluster);
			 bbox.eat(bbox_cluster.ll);
			 bbox.eat(bbox_cluster.ur);
		}

		// translate bbox relative to detected part (adjust location of part to ly1 coordinates i.e. apply contraction factors)
		//bbox += ipoint2(nd->x , nd->y) * contraction_factor;
				
		// each bounding box sould grow by 20% since some parts will not cover whole object
		bbox.grow(1.2);
				
		// adjust by removing border and scale from location
		bbox.ll.x = (bbox.ll.x - border) * scale_factor;
		bbox.ll.y = (bbox.ll.y - border) * scale_factor;
		bbox.ur.x = (bbox.ur.x - border) * scale_factor;
		bbox.ur.y = (bbox.ur.y - border) * scale_factor;

		// search for any similar bbox in current result list
		bool found_similar = false;
		for (std::list<hop_detection>::iterator iter_bbox = result->begin(); iter_bbox != result->end(); iter_bbox++) {
			if (abs(iter_bbox->box.ll.x - bbox.ll.x) <= 3 &&
				abs(iter_bbox->box.ll.y - bbox.ll.y) <= 3 && 
				abs(iter_bbox->box.ur.x - bbox.ur.x) <= 3 &&
				abs(iter_bbox->box.ur.y - bbox.ur.y) <= 3) {
					// found similar example - do not add it to list
					found_similar = true;
					break;
			}
		}
		// add to result list
		if (found_similar == false && bbox.x_dim() > 0 && bbox.y_dim() > 0) {
			hop_detection det;
			det.box = bbox;
			det.category_id = nd->m;
			det.responses = vector<float>(5);
			det.responses[0] = nd->r(R_RESPONSE);
			det.responses[1] = nd->r(G_RESPONSE);
			det.responses[2] = nd->r(RR_RESPONSE);
			det.responses[3] = reverse_s_response(nd->r(S_RESPONSE));
			det.responses[4] = 0;

			// get X_RESPONSE i.e. get SVMt response

			node* p = lib->parts[nd->z][nd->m];
			vs_part_data* pd = dynamic_cast<vs_part_data*>(p->data);
			double svmdfval = 0.0;

            bool ok = p->is_attr_set(NODE_DELETED_ATTR) == false;
            if (ok && pd != nullptr && pd->svmt.svm != nullptr) {
                int nchildren = p->count_neighbors(atom("lyrSrc"));

                det.responses[4] = check_svmt(pd->svmt, nchildren, n, true);
            }

			result->push_back(det);
		}

	}
	fflush(stdout);

	return result->size();
}
//// This is an example of an exported variable
//HOPINTERFACE_API int nhopinterface=0;
//
//// This is an example of an exported function.
//HOPINTERFACE_API int fnhopinterface(void)
//{
//	return 42;
//}
//
//// This is the constructor of a class that has been exported.
//// see hop-interface.h for the class definition
//Chopinterface::Chopinterface()
//{
//	return;
//}
