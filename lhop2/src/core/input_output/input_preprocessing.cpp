
#include "input_preprocessing.h"

std::vector<AbstractInputObject*> ResizeImagePreprocessing::doPreprocessing(const std::vector<AbstractInputObject*> input_list) const  {
	std::vector<AbstractInputObject*> result;
		
	if (resize_value == -100) {
		result = input_list;
	} else {
		for (auto iter = input_list.begin(); iter != input_list.end(); ++iter) {
			// make sure we have images as input
			AbstractImageInputObject* input = castToImageObject(*iter);

			// retrive image
			std::shared_ptr<img> image = input->getImage();
			
			// resize image
			std::shared_ptr<img> resized_image(new img());
			image->get_resized(*resized_image, resize_value);

			// also normalize the image
			resized_image->normalize(0,1);

			// also apply resizing to groundtruths
			std::vector<GroundtruthObject> groundtruths = input->getGroundtruths();

			double scale_factor = resized_image->width/(double)image->width;

			std::vector<GroundtruthObject> resized_groundtruths;
			for (auto gt_iter = groundtruths.begin(); gt_iter != groundtruths.end(); ++gt_iter) {
				resized_groundtruths.push_back(GroundtruthObject());
				resized_groundtruths.back().label = gt_iter->label;
				resized_groundtruths.back().region = gt_iter->region.get_scaled(scale_factor); 
			}

			// save to result as MemoryImageInputObject
			result.push_back(new MemoryImageInputObject(resized_image, resized_groundtruths, input->getGroupMap()));

			// DO NOT DELETE input as this must be taken care by the caller
		}
	}

	return result;
}

std::vector<AbstractInputObject*> BlurImagePreprocessing::doPreprocessing(const std::vector<AbstractInputObject*> input_list) const {
	std::vector<AbstractInputObject*> result;
		
	if (scale_mask_size < 3 || scale_sigma <= 0.0) {
		result = input_list;
	} else {
		for (auto iter = input_list.begin(); iter != input_list.end(); ++iter) {
			// make sure we have images as input
			AbstractImageInputObject* input = castToImageObject(*iter);

			// retrive image
			std::shared_ptr<img> image = input->getImage();
			
			// blur image
			image->blur(scale_mask_size, scale_sigma);

			// save to result as MemoryImageInputObject (we can optimize and skip new alocation in case input is already type of MemoryImageInputObject)
			if (dynamic_cast<MemoryImageInputObject*>(input) == nullptr)
				result.push_back(new MemoryImageInputObject(image, input->getGroundtruths(), input->getGroupMap()));

			// DO NOT DELETE input as this must be taken care by the caller
		}
	}

	return result;
}


std::vector<AbstractInputObject*> ScaleImagePreprocessing::doPreprocessing(const std::vector<AbstractInputObject*> input_list) const {
	// use "scale" identifier as grouping for image scaling
	static GroupableObject::Type group_type;
	group_type.name = "scale";

	std::vector<AbstractInputObject*> result;

	for (auto iter = input_list.begin(); iter != input_list.end(); ++iter) {
		// make sure we have images as input
		AbstractImageInputObject* input = castToImageObject(*iter);

		// retrive image
		std::shared_ptr<img> image = input->getImage();
		std::vector<GroundtruthObject> groundtruths = input->getGroundtruths();
		std::multimap<GroupableObject::Type, GroupableObject::Member> groups = input->getGroupMap();

		int image_width = image->width;
		int image_height = image->height;

		// prepare group member struct;
		GroupableObject::Member group; 
		group.group_id = groups.count(group_type); // just count how many times we have already used this pre-processing

		// perform scaling until limit is reached
		int count = 0;
		bool process = true;

		while (process) {
			// add image (should also add the original unscaled image in first loop)
			MemoryImageInputObject* resized_input = new MemoryImageInputObject(image, input->getGroundtruths(), groups);
				
			// add group identifier for this scaling
			group.group_member_id = count;
			resized_input->addGroup(group_type, group);

			result.push_back(resized_input);
				
			image_width = (int)(image_width*scale_factor);
			image_height = (int)(image_height*scale_factor);
			++count;
			
			if (image_width < scale_limit && image_height < scale_limit || count >= max_scales) {
				process = false;
			} else {
				img* resized_image = new img();

				// resize image
				*resized_image = image->get_resized(image_width, image_height);				

				// continue with resized image
				image = std::shared_ptr<img>(resized_image);

				// also resize groundtruths
				for (auto gt_iter = groundtruths.begin(); gt_iter != groundtruths.end(); ++gt_iter) {
					gt_iter->region = gt_iter->region.get_scaled(scale_factor);
				}
			}
		}

	}

	return result;
}


std::vector<AbstractInputObject*> SeperateImageColorsPreprocessing::doPreprocessing(const std::vector<AbstractInputObject*> input_list) const {
	// use "scale" identifier as grouping for image scaling
	static GroupableObject::Type group_type;
	group_type.name = "rgb";

	std::vector<AbstractInputObject*> result;

	for (auto iter = input_list.begin(); iter != input_list.end(); ++iter) {
		// make sure we have images as input
		AbstractImageInputObject* input = castToImageObject(*iter);

		// retrive image
		std::shared_ptr<img> image = input->getImage();
		std::vector<GroundtruthObject> groundtruths = input->getGroundtruths();
		std::multimap<GroupableObject::Type, GroupableObject::Member> groups = input->getGroupMap();

		int image_width = image->width;
		int image_height = image->height;

		// prepare group member struct;
		GroupableObject::Member group; 
		group.group_id = groups.count(group_type); // just count how many times we have already used this pre-processing

		// split into R, G and B channels
		std::shared_ptr<img> R(new img());
		std::shared_ptr<img> G(new img());
		std::shared_ptr<img> B(new img());

		image->get_colors(*R, *G, *B);

		// create input from each color channel
		{
			// R channel (group_member_id == GROUP_CHANNEL_R)
			group.group_member_id = GROUP_CHANNEL_R; 
			MemoryImageInputObject* resized_input = new MemoryImageInputObject(R, input->getGroundtruths(), groups); resized_input->addGroup(group_type, group);
			result.push_back(resized_input);
		}
		{
			// G channel (group_member_id == GROUP_CHANNEL_G)
			group.group_member_id = GROUP_CHANNEL_G; 
			MemoryImageInputObject* resized_input = new MemoryImageInputObject(G, input->getGroundtruths(), groups); resized_input->addGroup(group_type, group);
			result.push_back(resized_input);
		}
		{
			// B channel (group_member_id == GROUP_CHANNEL_B)
			group.group_member_id = GROUP_CHANNEL_B; 
			MemoryImageInputObject* resized_input = new MemoryImageInputObject(B, input->getGroundtruths(), groups); resized_input->addGroup(group_type, group);
			result.push_back(resized_input);
		}
	}

	return result;
}


std::vector<AbstractInputObject*> FlipImagePreprocessing::doPreprocessing(const std::vector<AbstractInputObject*> input_list) const {
	std::vector<AbstractInputObject*> result;

	for (auto iter = input_list.begin(); iter != input_list.end(); ++iter) {
		// make sure we have images as input
		AbstractImageInputObject* input = castToImageObject(*iter);

		// retrive image
		std::shared_ptr<img> image = input->getImage();
		std::vector<GroundtruthObject> groundtruths = input->getGroundtruths();
		std::multimap<GroupableObject::Type, GroupableObject::Member> groups = input->getGroupMap();
		
		// flip image 
		image->flip_horizontal();
		
		// also flip groundtruths
		for (auto iter = groundtruths.begin(); iter != groundtruths.end(); ++iter) {
			irectangle2& r = iter->region;
			int tmp = r.ll.x;
	        
			r.ll.x = image->width - 1 - r.ur.x;
			r.ur.x = image->width - 1 - tmp;
		}
		
		// push new input object into the results
		result.push_back(new MemoryImageInputObject(image, groundtruths, groups));
	}

	return result;
}

std::vector<AbstractInputObject*> ExtractImageGroundtruthPreprocessing::doPreprocessing(const std::vector<AbstractInputObject*> input_list) const {
	// use "scale" identifier as grouping for image scaling
	static GroupableObject::Type group_type;
	group_type.name = "region";

	std::vector<AbstractInputObject*> result;

	for (auto iter = input_list.begin(); iter != input_list.end(); ++iter) {
		// make sure we have images as input
		AbstractImageInputObject* input = castToImageObject(*iter);

		// retrive image
		std::shared_ptr<img> image = input->getImage();
		std::vector<GroundtruthObject> groundtruths = input->getGroundtruths();
		std::multimap<GroupableObject::Type, GroupableObject::Member> groups = input->getGroupMap();

		int image_width = image->width;
		int image_height = image->height;

		// prepare group member struct;
		GroupableObject::Member group; 
		group.group_id = groups.count(group_type); // just count how many times we have already used this pre-processing

		int count = 0;
		for (auto gt_iter = groundtruths.begin(); gt_iter != groundtruths.end(); ++gt_iter) {
			// extract only regions with specific label (if label not empty)
			if (only_label.empty() == false && only_label.compare(gt_iter->label) == 0)
				continue;

			std::shared_ptr<img> cropped_image(new img());

			// crop image based on groundtruth
			*cropped_image = image->cut(gt_iter->region);

			// create only one new groundtruth based on cropped image
			std::vector<GroundtruthObject> cropped_groundtruth_list(1);
			cropped_groundtruth_list[0].label = gt_iter->label;
			cropped_groundtruth_list[0].region = irectangle2(0,0, gt_iter->region.x_dim(), gt_iter->region.y_dim());

			MemoryImageInputObject* resized_input = new MemoryImageInputObject(cropped_image, cropped_groundtruth_list, groups); 
			
			group.group_member_id = count; 
			
			resized_input->addGroup(group_type, group);
			
			result.push_back(resized_input);
			
			count++;
		}
	}

	return result;
}


std::vector<AbstractInputObject*> SplitImagePreprocessing::doPreprocessing(const std::vector<AbstractInputObject*> input_list) const {
	// use "scale" identifier as grouping for image scaling
	static GroupableObject::Type group_type;
	group_type.name = "tile";

	std::vector<AbstractInputObject*> result;

	for (auto iter = input_list.begin(); iter != input_list.end(); ++iter) {
		// make sure we have images as input
		AbstractImageInputObject* input = castToImageObject(*iter);

		// retrive image
		std::shared_ptr<img> image = input->getImage();
		std::vector<GroundtruthObject> groundtruths = input->getGroundtruths();
		std::multimap<GroupableObject::Type, GroupableObject::Member> groups = input->getGroupMap();

		int image_width = image->width;
		int image_height = image->height;

		// prepare group member struct;
		GroupableObject::Member group; 
		group.group_id = groups.count(group_type); // just count how many times we have already used this pre-processing

		if (n <= 1) {			
			// skip spliting if n <= 1
			result.push_back(input);
		} else {
			int count = 0;

			int img_width = image->width;
			int img_height = image->height;
			int overlap = overlap / 2;
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					int xmin = max<int>(0, i*img_width/n - overlap);
					int xmax = min<int>(img_width, (i + 1)*img_width/n + overlap);
					int ymin = max<int>(0, j*img_height/n - overlap);
					int ymax = min<int>(img_height, (j + 1)*img_height/n + overlap);
					irectangle2 rect(xmin, ymin, xmax, ymax);
					
					// cut the image into tile
					std::shared_ptr<img> tile_image(new img());
					*tile_image = image->cut(rect);

					// also copy only groundtruths that belong to this tile (i.e. they must intersect with the tile at least by a factor of rthresh)
					std::vector<GroundtruthObject> tile_groundtruths;
					for (auto iter = groundtruths.begin(); iter != groundtruths.end(); ++iter) {
						irectangle2 isection = rect.intersection(iter->region);

						// add only if PASCAL intersection factor >= rthresh
						if (!isection.invalid() && (double)isection.area()/iter->region.area() >= rthresh) {
							tile_groundtruths.push_back(GroundtruthObject());
							tile_groundtruths.back().label = iter->label;
							tile_groundtruths.back().region = isection - ipoint2(xmin, ymin);
						}
					}

					MemoryImageInputObject* resized_input = new MemoryImageInputObject(tile_image, tile_groundtruths, groups); 
			
					group.group_member_id = count;
			
					resized_input->addGroup(group_type, group);
			
					result.push_back(resized_input);

					count++;
				}
			}
		}
	}

	return result;
}