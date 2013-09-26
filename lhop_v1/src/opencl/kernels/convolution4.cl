// include files
#include <utils.cl>

#if __DEVICE_LOCAL_MEM_TYPE__ > 0
	#define HAS_FAST_LOCAL_MEM 1
#endif

__attribute__((vec_type_hint(float4)))
__kernel
void image_convolve_vector4(__global float* image, int2 image_size,
							__constant float* convol_mask, int4 convol_mask_size,
							__local float* local_data, int2 local_data_size,
							uint save_max, __global float* max_image, __global float* max_image_index,
							__global float* output, uint layer1_creator_type, __global int* error) {
	// define vector size 4 and include body of actualy function for convolution
	#define VECTOR_SIZE 4

	#include <core_convolution.cl>

	#undef 	VECTOR_SIZE
}
