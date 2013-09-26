// include files
#include <utils.cl>

/**
 * Normal implementation of parallel reduction
 */

__attribute__((vec_type_hint(float4)))
__kernel
void get_sum(__global float* input, int input_size, int offset,
			__local float* local_data, int local_data_size,
			__global float* result, __global int* error) {

	#define VECTOR_SIZE 4

	#define vector_float _vector_n(float,VECTOR_SIZE)	// this should produce float16, float8, float4 or float2 - based on value of VECTOR_SIZE
	#define vector_load _vector_n(vload,VECTOR_SIZE)

	if (get_work_dim() != 1) {
		*error = 1;
		return;
	}

	// local_data_size must be at least of size VECTOR_SIZE
	if (local_data_size < VECTOR_SIZE) {
		*error = 2;
		return;
	}

	if (*error == 0) {

		vector_float sum_value = (vector_float)0;	
		
		__global float* input_offset = input + offset;

		for (int i = 0; i < input_size/VECTOR_SIZE; i+=local_data_size) {
			if ((int)get_global_id(0) + i >= 0 && (int)get_global_id(0) + i  < input_size) {
				sum_value += vector_load(get_global_id(0) + i, input_offset);
			}
		}
		
		// additionaly check for all values outside of input_size/VECTOR_SIZE (in case input_size is not dividable by VECTOR_SIZE)
		for (int i = VECTOR_SIZE * (input_size/VECTOR_SIZE); i < input_size; i++) {
			sum_value.s0 += input_offset[i];
		}
		
		sum_value.s0 += sum_value.s1;
		sum_value.s0 += sum_value.s2;
		sum_value.s0 += sum_value.s3;
#if VECTOR_SIZE > 4
		sum_value.s0 += sum_value.s4;
		sum_value.s0 += sum_value.s5;
		sum_value.s0 += sum_value.s6;
		sum_value.s0 += sum_value.s7;
	#if VECTOR_SIZE > 8	
		sum_value.s0 += sum_value.s8;
		sum_value.s0 += sum_value.s9;
		sum_value.s0 += sum_value.sA;
		sum_value.s0 += sum_value.sB;
		sum_value.s0 += sum_value.sC;
		sum_value.s0 += sum_value.sD;
		sum_value.s0 += sum_value.sE;
	#endif
#endif
		
		local_data[get_local_id(0)] = sum_value.s0;
		
		barrier(CLK_LOCAL_MEM_FENCE);
		
		if (get_local_id(0) == get_local_size(0) -1) {
			
		
			sum_value = (vector_float)0;
			
			for (int i = 0; i < local_data_size/VECTOR_SIZE; i++) {
				sum_value += vector_load(i, local_data);
			}
			
			// additionaly check for all values outside of local_data_size/VECTOR_SIZE (in case local_data_size is not dividable by VECTOR_SIZE)
			for (int i = VECTOR_SIZE * (local_data_size/VECTOR_SIZE); i < local_data_size; i++) {
				sum_value.s0 += local_data[i];
			}
					
			sum_value.s0 += sum_value.s1;
			sum_value.s0 += sum_value.s2;
			sum_value.s0 += sum_value.s3;
#if VECTOR_SIZE > 4
			sum_value.s0 += sum_value.s4;
			sum_value.s0 += sum_value.s5;
			sum_value.s0 += sum_value.s6;
			sum_value.s0 += sum_value.s7;
	#if VECTOR_SIZE > 8	
			sum_value.s0 += sum_value.s8;
			sum_value.s0 += sum_value.s9;
			sum_value.s0 += sum_value.sA;
			sum_value.s0 += sum_value.sB;
			sum_value.s0 += sum_value.sC;
			sum_value.s0 += sum_value.sD;
			sum_value.s0 += sum_value.sE;
	#endif
#endif

			*result = sum_value.s0;
		}
	}
	#undef vector_float
	#undef vector_load
	#undef VECTOR_SIZE
}

