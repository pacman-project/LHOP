/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// utils_toolset.cpp 

#include "main_toolset.h"

#include "utils/utils.h"

using namespace std;

// groundtruth functions 
// (TODO: maybe move them into static_utils library - create new class to handle it)
///////////////////////////////////////////////////////////////////////////////

void read_groundtruth(list<irectangle2>& rectangles, const string& fname, 
					  const string& cat_name_only, const string& gt_extension)
{
	list<std::pair<irectangle2,string> > rectangles_with_cat;
	// read with category names
	read_groundtruth(rectangles_with_cat, fname, cat_name_only, gt_extension);
	// copy only rectangles and delete strings
	for (list<std::pair<irectangle2,string> >::iterator iter = rectangles_with_cat.begin(); 
		iter != rectangles_with_cat.end(); iter++) {
		std::pair<irectangle2,string> p = *iter;
		rectangles.push_back(p.first);
	}
	
}

void read_groundtruth(list<std::pair<irectangle2,string> >& rectangles, const string& fname, 
					  const string& cat_name_only, const string& gt_extension, bool display_warnings )
{
	const string file = change_extension(fname, gt_extension);
    ifstream is(file.c_str());

	rectangles.clear();
	bool warning_displayed = !display_warnings;
    if (is.is_open()) {
        while (true) {
			string line;
			// read and parse whole line as stream of string
			getline(is,line);
			stringstream sline(line);

			rectangle2<double> rect;
			string cat_name;

			// first read rectangle values
			sline >> rect;
            if (is.fail()) break;
			
			// then read category name and trim it
			sline >> cat_name;
			cat_name = trimString(cat_name);
			
			if (cat_name.length() <= 0) {
				// category name is missing in groundtruth file
				// this may be old version so only ignore it and use all values
				if (warning_displayed == false) {
					// diplay warning only once for one file
					cout << "\nWARNING: category name is missing in groundtruth data (file: '" << file << "') - will use all bounding boxes in file" << endl;
					warning_displayed = true;
				}
			} else if (cat_name_only.length() > 0 && cat_name != cat_name_only) {
				// skip all ground truth values that do not belog to category cat_name_only
				continue;
			}
			rectangles.push_back(std::pair<irectangle2,string>(irectangle2((int)rect.ll.x, (int)rect.ll.y, (int)rect.ur.x, (int)rect.ur.y), cat_name));
        }
	} else {
		cout << "\nWARNING: Unable to open groundtruth file '" << file << "'" << endl;
	}
}

void save_groundtruth(const list<std::pair<irectangle2,string> >& rectangles, const string& outdir, const string& srcname, 
    const string& prefix, const string& suffix, bool changeext )
{
    string fstr; 

    if (changeext) fstr = outdir + prefix + change_extension(srcname, suffix + ".groundtruth");
    else fstr = outdir + prefix + srcname + suffix;

    ofstream os(fstr.c_str());

    for (list<std::pair<irectangle2,string> >::const_iterator iter = rectangles.begin(); iter != rectangles.end(); ++iter) {
		const std::pair<irectangle2,string> p = *iter;
        const irectangle2& r = p.first;
		const string& cat_name = p.second;
        
        os << r.ll.x << ' ' << r.ll.y << ' ' << r.ur.x << ' ' << r.ur.y << ' ' << cat_name << endl;
    }
    os.close();
}