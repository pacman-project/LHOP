// atiom.h

#include "atom.h"
/*
map<string, int> atom::name_dict;
vector<string> atom::name_list;
*/
atom::atom(string s)
{
	pair<map<string, int>::iterator, bool> pr;
	//printf("%s %d\n", s.c_str(), get_name_dict().size());
	pr = get_name_dict().insert(pair<string, int>(s, (int)get_name_list().size()));
	if (pr.second) {
	    index = (int)get_name_list().size();
	    get_name_list().push_back(s);
	} else index = (pr.first)->second;
}

string atom::get_name(int i) { return (i >= 0 && i < (int)get_name_list().size()) ? get_name_list()[i] : "nullptr"; }

void atom::print_all_atoms() 
{
for (map<string, int>::const_iterator iter = get_name_dict().begin(); iter != get_name_dict().end(); ++iter) {
    cout << iter->first << '\t' << iter->second << endl;
}
}
