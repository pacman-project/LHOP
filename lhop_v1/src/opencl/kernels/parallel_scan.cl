// include files
#include <utils.cl>

/**
 * Normal implementation of sequential prefix sum 
 */

__kernel
void sequential_prefix_sum(__global ocl_layer1_data_coordinates* coords, int coord_size, __global int* error) {

	// only one work dimension supported 
	if (get_work_dim() != 1) {
		*error = 1;
		return;
	}
	
	int sum = 0;
	for (int i = 0; i < coord_size; i++) {
		// set coordinate from sum of previous elements
		coords[i].offset += sum;  
		// set sum to include value of this elemnet
		sum += coords[i].size;
	}
	
}



/**
 * NVIDIA implementation of parallel prefix sum 
 */

inline int scan1Inclusive(int idata, __local int *l_Data, int size){	
	//int pos = 2 * get_local_id(0) - (get_local_id(0) & (size - 1));
    int pos = 2 * get_local_id(0) - (get_local_id(0) % (size ));
    l_Data[pos] = 0;
    pos += size;
    l_Data[pos] = idata;

    for(int offset = 1; offset < size; offset <<= 1){
        barrier(CLK_LOCAL_MEM_FENCE);
        int t = l_Data[pos] + l_Data[pos - offset];
        barrier(CLK_LOCAL_MEM_FENCE);
        l_Data[pos] = t;
    }

    return l_Data[pos];
}

inline int scan1Exclusive(int idata, __local int *l_Data, int size){
    return scan1Inclusive(idata, l_Data, size) - idata;
}


//Vector scan: the array to be scanned is stored
//in work-item private memory as uint4
inline int4 scan4Inclusive(int4 data4, __local int *l_Data, int size){
    //Level-0 inclusive scan
    data4.y += data4.x;
    data4.z += data4.y;
    data4.w += data4.z;

    //Level-1 exclusive scan
    int val = scan1Inclusive(data4.w, l_Data, size / 4) - data4.w;

    return (data4 + (int4)val);
}

inline int4 scan4Exclusive(int4 data4, __local int *l_Data, int size){
    return scan4Inclusive(data4, l_Data, size) - data4;
}

////////////////////////////////////////////////////////////////////////////////
// Scan kernels
////////////////////////////////////////////////////////////////////////////////
// First step that computes prefix sum only within one work-group
__kernel
void scanExclusiveLocal1(
    __global ocl_layer1_data_coordinates* coords, int coords_size,
    __local int *l_Data,
    int local_work_gorup_size
){
	int4 idata = (int4)0;
	if (4 * get_global_id(0) < coords_size) {
		//Load data
		idata = (int4)(coords[4*get_global_id(0) + 0].size,
							coords[4*get_global_id(0) + 1].size,
							coords[4*get_global_id(0) + 2].size,
							coords[4*get_global_id(0) + 3].size);
	}
	
	//Calculate exclusive scan
	int4 odata  = scan4Exclusive(idata, l_Data, local_work_gorup_size * 4);

	if (4 * get_global_id(0) < coords_size) {
		//Write back
		coords[4*get_global_id(0) + 0].offset = odata.x;
		coords[4*get_global_id(0) + 1].offset = odata.y;
		coords[4*get_global_id(0) + 2].offset = odata.z;
		coords[4*get_global_id(0) + 3].offset = odata.w;
    }
}


// Second step in prefix sum that makes temporary buffer which 
// gets sum of each original work-gorup (in first step) and makes prefix sum on this temp buffer
__kernel
void scanExclusiveLocal2(
	__global ocl_layer1_data_coordinates* coords,
    __global ocl_layer1_data_coordinates* buffer,
    __local int *l_Data,
    int buffer_length, // number of additional prefix sums
    int local_work_gorup_processed_size // size of each work-group
){
    //Load top elements
    //Convert results of bottom-level scan back to inclusive
    //Skip loads and stores for inactive work-items of the work-group with highest index(pos >= buffer_length)
    int4 data = 0;
    if(4*get_global_id(0) < buffer_length) {
		data.x = coords[(local_work_gorup_processed_size - 1) + (local_work_gorup_processed_size) * (4*get_global_id(0) + 0)].offset + 
				coords[(local_work_gorup_processed_size - 1) + (local_work_gorup_processed_size) * (4*get_global_id(0) + 0)].size;
				
		data.y = coords[(local_work_gorup_processed_size - 1) + (local_work_gorup_processed_size) * (4*get_global_id(0) + 1)].offset + 
				coords[(local_work_gorup_processed_size - 1) + (local_work_gorup_processed_size) * (4*get_global_id(0) + 1)].size;
				
		data.z = coords[(local_work_gorup_processed_size - 1) + (local_work_gorup_processed_size) * (4*get_global_id(0) + 2)].offset + 
				coords[(local_work_gorup_processed_size - 1) + (local_work_gorup_processed_size) * (4*get_global_id(0) + 2)].size;
				
		data.w = coords[(local_work_gorup_processed_size - 1) + (local_work_gorup_processed_size) * (4*get_global_id(0) + 3)].offset + 
				coords[(local_work_gorup_processed_size - 1) + (local_work_gorup_processed_size) * (4*get_global_id(0) + 3)].size;
    }

    //Compute
    int4 odata = scan4Exclusive(data, l_Data, min(buffer_length, (int)get_local_size(0)*4));

    //Avoid out-of-bound access
    if(4*get_global_id(0) < buffer_length) {
		buffer[4*get_global_id(0) + 0].size = data.x;
        buffer[4*get_global_id(0) + 0].offset = odata.x;
        
        buffer[4*get_global_id(0) + 1].size = data.y;
        buffer[4*get_global_id(0) + 1].offset = odata.y;
        
        buffer[4*get_global_id(0) + 2].size = data.z;
        buffer[4*get_global_id(0) + 2].offset = odata.z;
        
        buffer[4*get_global_id(0) + 3].size = data.w;
        buffer[4*get_global_id(0) + 3].offset = odata.w;
	}
}

// Third step in prefix sum that updates result from step 1 with computed sums of step two.
// (each value in coords must be added by sum computed for its work-group in second step)
__kernel
void uniformUpdate(
    __global ocl_layer1_data_coordinates* coords, int coords_size,
    __global ocl_layer1_data_coordinates* buffer
){
    __local int buf[1];
	
	if (4*get_global_id(0) < coords_size) {
		
		int4 data4 = (int4)(coords[4*get_global_id(0) + 0].offset,
							coords[4*get_global_id(0) + 1].offset,
							coords[4*get_global_id(0) + 2].offset,
							coords[4*get_global_id(0) + 3].offset);
	
		if(get_local_id(0) == 0)
			buf[0] = buffer[get_group_id(0)].offset;

		barrier(CLK_LOCAL_MEM_FENCE);
		data4 += (int4)buf[0]; 
	    
		coords[4*get_global_id(0) + 0].offset = data4.x;
		coords[4*get_global_id(0) + 1].offset = data4.y;
		coords[4*get_global_id(0) + 2].offset = data4.z;
		coords[4*get_global_id(0) + 3].offset = data4.w;  
		
	} else {
		barrier(CLK_LOCAL_MEM_FENCE);
	}
}
