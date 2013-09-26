// Atom class
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _ATOM_H_
#define _ATOM_H_

#include <stdio.h>
#include <string>
#include <vector>
#include <iostream>
#include <map>

using namespace std;


class atom {
private: 
    int index;

    static map<string, int> _name_dict;
    static vector<string> _name_list;

    static map<string, int>& get_name_dict() { static map<string, int> _name_dict; return _name_dict; }
    static vector<string>& get_name_list() { static vector<string> _name_list; return _name_list; }

public:
    atom(const atom& a) : index(a.index) {  }
    atom(string s);

    int get_index() const { return index; }
    string get_name() const { return get_name_list()[index]; }
    operator int() const { return index; }
    static string get_name(int i);
    static void print_all_atoms();
};

#endif /* _ATOM_H_ */
