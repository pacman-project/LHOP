// layserver.cpp : Defines the entry point for the console application.
//

#include <ctime>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <list>

#include "../../interface/hop.h"

/**


lay1create:
	layers.cpp::create_layer1()

possible commands:
	- create layer "cfg string in BASE64" input_filename
	- learn "cfg string in BASE64" input_filename
*/

using namespace std;

list<string> parse_line(const string& line) {

	list<string> result;

	bool find_end_of_quote = false;
	int sub_start = 0;
	
	while (true) {
		// find next delimiter
		size_t space_pos = line.find(' ', sub_start);
		size_t quote_pos = line.find('"', sub_start);
		
		size_t pos = 0;
		if (find_end_of_quote)
			pos = quote_pos;
		else
			pos = min<size_t>(space_pos,quote_pos);

		// get data from sub_start to pos - 1 and put it to result list
		if (pos - sub_start > 0) {
			result.push_back(line.substr(sub_start, pos - sub_start));
		}

		if (pos == string::npos)
			break;

		// set sub_start for next string based on type of delimiter
		if (line[pos] == '"') {			
			find_end_of_quote = !find_end_of_quote;
		} else {		
			find_end_of_quote = false;
		}
		sub_start = pos + 1;
	}
	return result;
}

string process_command(const list<string> args, bool &stop) {

	// first we need to redirect all output to a file (cout, cerr and clog)
//	ofstream cout_file, cerr_file, clog_file;
//	cout_file.open("std_output.txt");
//	cerr_file.open("std_error.txt");
//	clog_file.open("std_log.txt");
	stringstream cout_redir, cerr_redir, clog_redir;

	// save original outputs
	streambuf* cout_sbuf = cout.rdbuf();
	streambuf* cerr_sbuf = cerr.rdbuf();
	streambuf* clog_sbuf = clog.rdbuf();

	// redirect outputs to files
	cout.rdbuf(cout_redir.rdbuf());
	cerr.rdbuf(cerr_redir.rdbuf());
	clog.rdbuf(clog_redir.rdbuf());

	list<string>::const_iterator args_it = args.begin();

	stringstream result;

	// first value MUST be command
	const string& cmd = *(args_it++);
	bool ok;
	if (cmd.compare("create_layer_1") == 0) {
		
		const string& cfg_values = *(args_it++);
		const string& filename = *(args_it++);
		
		ok = hop_1_inference("", filename.c_str(), cfg_values.c_str());		

	} else if (cmd.compare("create_layer_n") == 0) {

		const string& cfg_values = *(args_it++);
		const string& filename = *(args_it++);
		
		ok = hop_n_inference("", filename.c_str(), cfg_values.c_str());

	} else if (cmd.compare("learn_layer") == 0) {
		const string& cfg_values = *(args_it++);
		const string& filename = *(args_it++);
		
		ok = hop_learning("", filename.c_str(), cfg_values.c_str());

	} else if (cmd.compare("exit") == 0) {
		stop = true;
		ok = true;
	} else if (cmd.compare("ping") == 0) {
		ok = true;
	} else {
		ok = false;
		cout << "ERROR: INVALID COMMAND\ntest line 1\ntest line 2"; // this will go directly to std_output file
	}
	// restore original outputs
	cout.rdbuf(cout_sbuf);
	cerr.rdbuf(cerr_sbuf);
	clog.rdbuf(clog_sbuf);

	// construct result string based on error and std_output
	if (ok)
		result << "STATUS: ok" << endl;
	else
		result << "STATUS: error" << endl;

	streambuf* result_rdbuf = result.rdbuf();
	while (cout_redir) {
		string line;
		getline(cout_redir, line);
		if (line.length() > 0)
			result << "std_out: " << line << endl;
	}

	while (cerr_redir) {
		string line;
		getline(cerr_redir, line);
		if (line.length() > 0)
			result << "std_err: " << line << endl;
	}

	return result.str();
}

int main(int argc, char* argv[])
{
	std::cout << "hopserver " << hop_time_stamp() << std::endl;

	bool stop = false;
	
	while (!stop) {

		// read command from console input and parse it into list of values
		string raw_line;
		getline(cin, raw_line);

		list<string> values = parse_line(raw_line);

		// skip if empty list
		if (values.size() <= 0)
			continue;
		
		// process command and print result
		cout << process_command(values, stop) << endl;
	}
	return 0;
}