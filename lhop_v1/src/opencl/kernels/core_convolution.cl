

#ifndef VECTOR_SIZE
	#error Unable to include 'core_convolution.cl' without VECTOR_SIZE definition. (VECTOR_SIZE is missing)
#endif

#define vector_float _vector_n(float,VECTOR_SIZE)	// this should produce float16, float8, float4 or float2 - based on value of VECTOR_SIZE
#define vector_float_convert _vector_n(convert_float,VECTOR_SIZE)
#define vector_int _vector_n(int,VECTOR_SIZE)
#define vector_load _vector_n(vload,VECTOR_SIZE)
#define vector_store _vector_n(vstore,VECTOR_SIZE)


/**
 * This file should be included only within function since it is without its own headers.
 */

{
	// only two dimensions are allowed
	if (get_work_dim() != 2) {
		*error = 1;
		return;
	}


	int2 convol_mask_size_half = convol_mask_size.xy/(int2)2;

	int get_global_id_1 = convol_mask_size_half.y + get_global_id(1);
	int get_global_id_0 = convol_mask_size_half.x + get_global_id(0) * VECTOR_SIZE;

	if ((int)get_global_id_1 >= convol_mask_size_half.y && (int)get_global_id_1 < image_size.y - convol_mask_size_half.y &&
		(int)get_global_id_0 >= convol_mask_size_half.x && (int)get_global_id_0 < image_size.x - convol_mask_size_half.x && *error == 0) {



// if device has fast local memory we can use it (CPU uses global memory for local mem so this would only slow it down
// (HAS_FAST_LOCAL_MEM should be set during compile time in HOST run-time)
#ifdef HAS_FAST_LOCAL_MEM

		int get_local_id_1 = convol_mask_size_half.y + get_local_id(1);
		int get_local_id_0 = convol_mask_size_half.x + get_local_id(0) * VECTOR_SIZE;

		int2 get_global_id_ = (int2)(get_global_id_0, get_global_id_1);
		int2 get_local_id_ = (int2)(get_local_id_0, get_local_id_1);

		// load local data with additional borders based on size of convolution kernel
		load_local_data(image, image_size, local_data, local_data_size, get_global_id_, get_local_id_, VECTOR_SIZE, (int4)(convol_mask_size_half.xy, convol_mask_size_half.xy));

		int img_local_pos_offset_start = getImgPos( ((int)get_local_id_1 - convol_mask_size_half.y), local_data_size.x , (int)get_local_id_0 - convol_mask_size_half.x);
		int skipped_local_pos_x = max(0, local_data_size.x - convol_mask_size.x);
#else
		int img_global_pos_offset_start = getImgPos( ((int)get_global_id_1 - convol_mask_size_half.y), image_size.x, (int)get_global_id_0 - convol_mask_size_half.x);
		int skipped_global_pos_x = max(0, image_size.x - convol_mask_size.x);
#endif

		vector_float max_value = (vector_float)-MAXFLOAT;
		vector_float max_index = (vector_float)-1;

		switch (layer1_creator_type) {

			case LAYER1_TYPE_STRUCT:
				for (int convol_j = 0; convol_j < convol_mask_size.z/2; convol_j++) {
					// make convolution with i-th filter and save it to sum_a
					vector_float sum_a;
					{
						// this 'template function' that is able to handle vector_float
						int convol_i = convol_j; // this is input (must be value 'convol_i')
						#include <make_convolution.cl>
						sum_a = sum; // output of make_convolution.cl is saved in 'sum'
					}
					// then make convolution with (convol_mask_size.z/2 + i)-th filter
					vector_float sum_b;
					{
						int convol_i = convol_mask_size.z/2 + convol_j;
						#include <make_convolution.cl>
						sum_b = sum;
					}

					// calculate result from two filtered values
					vector_float result = sqrt(sum_a * sum_a + sum_b * sum_b);

					// save value back to global memory
					int pos = getImgPos( (int)get_global_id_1 + image_size.y * convol_j, image_size.x, (int)get_global_id_0);

					if ((int)get_global_id_0 + VECTOR_SIZE < image_size.x - convol_mask_size_half.x) {
						vector_store(result, 0, output + pos);
					} else {
						float tmp_storage[VECTOR_SIZE];
						private float* p_tmp_storage = tmp_storage;
						vector_store(result, 0, p_tmp_storage);

						int count = (image_size.x - convol_mask_size_half.x) - get_global_id_0;
						for (int i = 0; i < count; i++)
							output[pos + i] = tmp_storage[i];
					}
					// save max value for this pixel positions
					if (save_max != 0) {

						// when using operator > or < on vector it returnes true as -1
						vector_float is_max = vector_float_convert((result > max_value) * -1);
						vector_float is_not_max = vector_float_convert((result <= max_value ) * -1);

						max_value = max_value * is_not_max + result * is_max;
						max_index = max_index * is_not_max + ((vector_float)convol_j) * is_max;
					}
				}
				break;
			case LAYER1_TYPE_APP:
				for (int convol_j = 0; convol_j < convol_mask_size.z/2; convol_j++) {
					// make convolution with i-th filter

					// this 'template function' that is able to handle vector_float
					vector_float result;
					{
						int convol_i = convol_j; // this is input (must be value 'convol_i')
						#include <make_convolution.cl>
						result = sum; // output of make_convolution.cl is saved in 'sum'
					}

					// clamp result a to min 0
					vector_float result_b = max(result, 0.0f); // ==> if (*ptr < 0.0) *ptr = 0.0;

					// for result b we can reverse sign and clamp it to 0
					vector_float result_a = max(-result, 0.0f); // ==> if (*ptr > 0.0) *ptr = 0.0; else *ptr = -(*ptr);

					// save result a back to global memory
					{
						int pos = getImgPos( (int)get_global_id_1 + image_size.y * convol_j, image_size.x, (int)get_global_id_0);

						if ((int)get_global_id_0 + VECTOR_SIZE < image_size.x - convol_mask_size_half.x) {
							vector_store(result_a, 0, output + pos);
						} else {
							float tmp_storage[VECTOR_SIZE];
							private float* p_tmp_storage = tmp_storage;
							vector_store(result_a, 0, p_tmp_storage);

							int count = (image_size.x - convol_mask_size_half.x) - get_global_id_0;
							for (int i = 0; i < count; i++)
								output[pos + i] = tmp_storage[i];
						}
					}

					// save result b back to global memory (at position convol_mask_size.z + convol_j)
					{
						int pos = getImgPos( (int)get_global_id_1 + image_size.y * (convol_mask_size.z/2 + convol_j), image_size.x, (int)get_global_id_0);

						if ((int)get_global_id_0 + VECTOR_SIZE < image_size.x - convol_mask_size_half.x) {
							vector_store(result_b, 0, output + pos);
						} else {
							float tmp_storage[VECTOR_SIZE];
							private float* p_tmp_storage = tmp_storage;
							vector_store(result_b, 0, p_tmp_storage);

							int count = (image_size.x - convol_mask_size_half.x) - get_global_id_0;
							for (int i = 0; i < count; i++)
								output[pos + i] = tmp_storage[i];
						}
					}
					// save max value for this pixel positions
					if (save_max != 0) {

						// for result a
						{
							// when using operator > or < on vector it returnes true as -1
							vector_float is_max = vector_float_convert((result_a > max_value) * -1);
							vector_float is_not_max = vector_float_convert((result_a <= max_value ) * -1);

							max_value = max_value * is_not_max + result_a * is_max;
							max_index = max_index * is_not_max + ((vector_float)convol_j) * is_max;
						}
						// for result b
						{
							// when using operator > or < on vector it returnes true as -1
							vector_float is_max = vector_float_convert((result_b > max_value) * -1);
							vector_float is_not_max = vector_float_convert((result_b <= max_value ) * -1);

							max_value = max_value * is_not_max + result_b * is_max;
							max_index = max_index * is_not_max + ((vector_float)(convol_mask_size.z/2 + convol_j)) * is_max;
						}
					}
				}
				break;
			case LAYER1_TYPE_DOG:
			case LAYER1_TYPE_LOGGABOR:
			case MULTIPLY_CONVOLUTION:
			default:
				for (int convol_i = 0; convol_i < convol_mask_size.z; convol_i++) {

					// make normal convolution
					#include <make_convolution.cl>

				#if VECTOR_SIZE == 2
					printf("%d result is %f\n",VECTOR_SIZE,sum.x);
					printf("%d result is %f\n",VECTOR_SIZE,sum.y);
				#endif

					// save value back to global memory
					int pos = getImgPos( (int)get_global_id_1 + image_size.y * convol_i, image_size.x, (int)get_global_id_0);


					if ((int)get_global_id_0 + VECTOR_SIZE < image_size.x - convol_mask_size_half.x) {
						vector_store(sum, 0, output + pos);
					} else {
						float tmp_storage[VECTOR_SIZE];
						private float* p_tmp_storage = tmp_storage;
						vector_store(sum, 0, p_tmp_storage);

						int count = (image_size.x - convol_mask_size_half.x) - get_global_id_0;
						for (int i = 0; i < count; i++)
							output[pos + i] = tmp_storage[i];
					}

					// save max value for this pixel positions
					if (save_max != 0) {

						// when using operator > or < on vector it returnes true as -1
						vector_float is_max = vector_float_convert((sum > max_value) * -1);
						vector_float is_not_max = vector_float_convert((sum <= max_value ) * -1);

						max_value = max_value * is_not_max + sum * is_max;
						max_index = max_index * is_not_max + ((vector_float)convol_i) * is_max;
					}

				}
				break;

		}

		// save max value and index if max calculated
		if (save_max != 0) {
			int pos = getImgPos( (int)get_global_id_1 , image_size.x, (int)get_global_id_0);
			if ((int)get_global_id_0 + VECTOR_SIZE < image_size.x - convol_mask_size_half.x) {
				vector_store(max_value, 0, max_image + pos);
				vector_store(max_index, 0, max_image_index + pos);
			} else {
				int count = (image_size.x - convol_mask_size_half.x) - get_global_id_0;

				float tmp_storage[VECTOR_SIZE];
				private float* p_tmp_storage = tmp_storage;

				vector_store(max_value, 0, p_tmp_storage);
				for (int i = 0; i < count; i++)
					max_image[pos + i] = tmp_storage[i];

				vector_store(max_index, 0, p_tmp_storage);
				for (int i = 0; i < count; i++)
					max_image_index[pos + i] = tmp_storage[i];
			}
		}

	} else {
#ifdef HAS_FAST_LOCAL_MEM
		// this barrier should be here since AMD implementation will not let threads finish unless all threads in same block get to barrier
		barrier(CLK_LOCAL_MEM_FENCE);
#endif
	}
}


#undef vector_float
#undef vector_float_convert
#undef vector_int
#undef vector_load
#undef vector_store


/*
 saved old version of load data just in case

		int2 copy_offset = (int2)(0,0);
		int2 copy_size = (int2)(0,0);

		copy_offset = (int2)(0,0);
		copy_size = (int2)(VECTOR_SIZE, 1);
		copy_to_local_data(image, image_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);

		// if upper border
		if (get_local_id(1) == 0) {
			copy_offset = (int2)(0, -convol_mask_size_half.y);
			copy_size = (int2)(VECTOR_SIZE, convol_mask_size_half.y);
			copy_to_local_data(image, image_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);

			// if upper left corrner
			if (get_local_id(0) == 0) {
				copy_offset = (int2)(- convol_mask_size_half.x, -convol_mask_size_half.y);
				copy_size = (int2)(convol_mask_size_half.x, convol_mask_size_half.y);
				copy_to_local_data(image, image_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);
			}

			// if upper right corrner
			if (get_local_id(0) == get_local_size(0) - 1 || get_global_id_0 + VECTOR_SIZE >= image_size.x - convol_mask_size_half.x) {
				copy_offset = (int2)(VECTOR_SIZE, -convol_mask_size_half.y);
				copy_size = (int2)(convol_mask_size_half.x, convol_mask_size_half.y);
				copy_to_local_data(image, image_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);
			}
		}

		// if lower border
		if (get_local_id(1) == get_local_size(1) - 1 || get_global_id_1 + VECTOR_SIZE >= image_size.y - convol_mask_size_half.y) {
			copy_offset = (int2)(0,1);
			copy_size = (int2)(VECTOR_SIZE, convol_mask_size_half.y);
			copy_to_local_data(image, image_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);

			// if lower left corrner
			if (get_local_id(0) == 0) {
				copy_offset = (int2)(- convol_mask_size_half.x,1);
				copy_size = (int2)(convol_mask_size_half.x, convol_mask_size_half.y);
				copy_to_local_data(image, image_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);
			}

			// if lower right corrner
			if (get_local_id(0) == get_local_size(0) - 1 || get_global_id_0 + VECTOR_SIZE >= image_size.x - convol_mask_size_half.x) {
				copy_offset = (int2)(VECTOR_SIZE,1);
				copy_size = (int2)(convol_mask_size_half.x, convol_mask_size_half.y);
				copy_to_local_data(image, image_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);
			}
		}

		// if left border
		if (get_local_id(0) == 0) {
			copy_offset = (int2)(- convol_mask_size_half.x, 0);
			copy_size = (int2)(convol_mask_size_half.x, 1);
			copy_to_local_data(image, image_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);
		}

		// if right border
		if (get_local_id(0) == get_local_size(0) - 1 || get_global_id_0 + VECTOR_SIZE >= image_size.x - convol_mask_size_half.x) {
			copy_offset = (int2)(VECTOR_SIZE, 0);
			copy_size = (int2)(convol_mask_size_half.x, 1);
			copy_to_local_data(image, image_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);
		}

		barrier(CLK_LOCAL_MEM_FENCE);
	*/