// Atom class
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _MAPREDUCE_H_
#define _MAPREDUCE_H_

#include <stdio.h>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <list>
#include "utils.h"
#include "streaming.h"

using namespace std;

/**
 * A simple wrapper arount list of streamable* representing main class for holding input items for map/reduce deployer.
 * This class can be sent to different computers therefore it MIST be streamable and MUST implement make_instance(), read_from_stream() and write_to_stream() functions.
 */

class mapreduce_items : public streamable {
public:
	list<streamable*> items;

	// streamable implementations
	virtual streamable* make_instance() const { return new mapreduce_items(); }
	virtual void read_from_stream(istreamer& is) {
		int size;
		is.read(size);
		while (size-- > 0) {
			streamable* obj;
			is.read(obj);
			items.push_back(obj);
		}
	}
	virtual void write_to_stream(ostreamer& os) {
		os.write((int)items.size());
		for (auto iter = items.begin(); iter != items.end(); ++iter)
			os.write(*iter);
	}
};

/**
 * A simple wrapper arount list of streamable* representing main class for holding map/reduce results.
 * This class can be sent to different computers therefore it MIST be streamable and MUST implement make_instance(), read_from_stream() and write_to_stream() functions.
 */
class mapreduce_result : public streamable {
	list<streamable*> results;
	
public:
	// define if each added result should be wrapped using streamed_pointer*
	bool wrap_results;
	
	mapreduce_result() : results(), wrap_results(0) {}
	mapreduce_result(bool w) : results(), wrap_results(w) {}

	// streamable implementations
	virtual streamable* make_instance() const { return new mapreduce_result(); }
	virtual void read_from_stream(istreamer& is) {
		is.read(wrap_results);
		int size;
		is.read(size);
		while (size-- > 0) {
			streamable* obj;
			is.read(obj);
			results.push_back(obj);
		}
	}
	virtual void write_to_stream(ostreamer& os) {
		os.write(wrap_results);
		os.write((int)results.size());
		for (auto iter = results.begin(); iter != results.end(); ++iter)
			os.write(*iter);
	}

	list<streamable*> get_result() { return results; }

	void add(streamable* obj) { 
		streamable* wrapped_obj = obj;
		// check if whe should wrap streamable* with streamed_pointer*
		if (wrap_results) {
			// check if this is already streamed_pointer*
			wrapped_obj = dynamic_cast<streamed_pointer*>(obj);
			// create streamed_pointer* otherwise
			if (wrapped_obj == nullptr)
				wrapped_obj = new streamed_pointer(obj);
		} 
		results.push_back(wrapped_obj); 

	}
	
};
/**
 * Common base class for every function implemented as map/reduce.
 *
 * This method must implement map() and reduce() function that transforms streamable* input items into streamable* results.
 * This class can be sent to different computers therefore it MIST be streamable and MUST implement make_instance(), read_from_stream() and write_to_stream() functions.
 *
 * A typical implementation of this class are map_learning_mapreduce, part_learning_mapreduce, layern_creator_mapreduce
 */
class base_mapreduce : public streamable {
public:
	virtual streamable* map(streamable* item) = 0;
	virtual streamable* reduce(list<streamable*> &item_list) = 0;
	virtual string map_get_key(streamable* item) = 0;

	virtual bool has_reduce() { return true; }
	virtual bool should_wrap_results() { return false; }

	virtual mapreduce_result* initilaize_results() { return new mapreduce_result(should_wrap_results()); }
};

/**
 * Common base class for all mapreduce deployer implementations. 
 *
 * Each deployer must implement submit method that recieves list of input items (for parallel execution) and a function that implements map and reduce methods.
 * The implementing submit method must handle map/reduce calls by first calling map on each input item and then calling reduce on all returned grouped items.
 *
 * This version already has simple non-parallel implementation of submit method. 
 * Other implementations are hadoop_deployer_mapreduce (from rosette/leoparts project) that extends base_deployer_mapreudce (jni wrapper for this class).
 */
class base_deployer_mapreudce {

public:
	virtual mapreduce_result* submit(mapreduce_items* parallel_data, base_mapreduce* mapreduce_func) {

		mapreduce_result* result = mapreduce_func->initilaize_results();

		if (mapreduce_func->has_reduce()) {
			// run map for each item
			list<streamable*> map_results;
			for (auto iter = parallel_data->items.begin(); iter != parallel_data->items.end(); ++iter)
				map_results.push_back(mapreduce_func->map(*iter));
		
			// run reduce on returned results
			result->add(mapreduce_func->reduce(map_results));

			// delete intermediate map results
			for (auto iter = map_results.begin(); iter != map_results.end(); ++iter)
				delete *iter;
		} else {
			// run map for each item and output it to result (delete intermediate map results only if wrapping is being used)
			for (auto iter = parallel_data->items.begin(); iter != parallel_data->items.end(); ++iter) { 
				streamable* map_item = mapreduce_func->map(*iter);
				result->add(map_item);

				if (mapreduce_func->should_wrap_results())
					delete map_item;
			}
		}

		return result;
	}
};


class openmp_deployer_mapreudce : public base_deployer_mapreudce {

public:
	virtual mapreduce_result* submit(mapreduce_items* parallel_data, base_mapreduce* mapreduce_func) {
		int input_count = parallel_data->items.size();

		vector<streamable*> input_vector;
		
		input_vector.reserve(input_count);

		for (auto iter = parallel_data->items.begin(); iter != parallel_data->items.end(); ++iter)
			input_vector.push_back(*iter);
		
		mapreduce_result* result = mapreduce_func->initilaize_results();
		
		if (mapreduce_func->has_reduce()) {			
			// run map for each item
			list<streamable*> map_results;
			
			#pragma omp parallel  
			{
				 #pragma omp master
				{
				  printf("There are %d threads\n", omp_get_num_threads());
				}
			#pragma omp for
			for (int i = 0; i < input_count; i++) {
				streamable* map_res = mapreduce_func->map(input_vector[i]);
				
				#pragma omp critical
				map_results.push_back(map_res);

				
			}
			}
			// run reduce on returned results
			result->add(mapreduce_func->reduce(map_results));

			// delete intermediate map results
			for (auto iter = map_results.begin(); iter != map_results.end(); ++iter)
				delete *iter;
		} else {
			bool should_wrap = mapreduce_func->should_wrap_results();
			// run map for each item and output it to result (delete intermediate map results only if wrapping is being used)
			list<streamable*> map_results;

			#pragma omp parallel for
			for (int i = 0; i < input_count; i++) { 
				streamable* map_item = mapreduce_func->map(input_vector[i]);

				if (should_wrap) {
					#pragma omp critical
					map_results.push_back(new streamed_pointer(map_item));

					delete map_item;
				} else {
					#pragma omp critical
					map_results.push_back(map_item);
				}
			}
			for (auto iter = map_results.begin(); iter != map_results.end(); ++iter)
				result->add(*iter);
		}

		return result;
	}
};

#endif /* _ATOM_H_ */
