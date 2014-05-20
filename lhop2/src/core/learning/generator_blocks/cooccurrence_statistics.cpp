
// cooccurrence_statistics
///////////////////////////////////////////////////////////////////////////////

#include "cooccurrence_statistics.h"

#include "core/structures/subparts/support_extension.h"

void CooccurrenceStatistics::updateStatistics(const InferenceLayer& inference_layer) {
	
	InferenceSupportAccess support_parts_access = inference_layer.getInferenceTree().getAccess<InferenceSupportAccess>();

	// whole size of the receptive field is twice of its radius
	cv::Size receptive_field_size(2 * this->receptive_field_radius, 2 * this->receptive_field_radius);

	// iterate over all location of the inference_tree 
	for (auto loc_iter = inference_layer.beginIterator(); loc_iter != inference_layer.endIterator(); ++loc_iter) {
		// get list of parts at that location
		auto location_parts = inference_layer.getPartsAt(loc_iter);

		// and iterate over all parts at that location
		// !!! old libhop version collected stat for the best part at each location only (one with the highest R_RESPONSE score) !!!
		for (auto part_iter = location_parts.begin(); part_iter != location_parts.end(); ++part_iter) {

			const InferredPart& center_part = **part_iter;
			
			// use part only if its vval >= center_part_threshold: TODO: change vval into something meaningfull
			if (center_part.getVVal() < center_part_threshold)
				continue;

			cv::Rect neighborhood_region(center_part.getLocation() - cv::Point2i(receptive_field_radius,receptive_field_radius), receptive_field_size);

			// find all neighborhood parts around center_part within specific radius
			std::vector<InferredPart*> neighborhood_parts = inference_layer.getPartsInRegion(neighborhood_region);

			CentralPart& central_part_statistics = central_parts[center_part.getCorrespondingVocabularyPartUUID()];

			std::set<InferredPart*> central_part_support = support_parts_access.getOnlyDifferentInitialLayerSupportParts(center_part);

			// and collect statistic for each one
			for (auto neighbor_iter = neighborhood_parts.begin(); neighbor_iter != neighborhood_parts.begin(); ++neighbor_iter) {
				const InferredPart& surrounding_part = **neighbor_iter;

				// TODO: verify if we should include this neighbor part into statistics:
				// 1. check if surrounding_part->vval() above threshold
				// 2. check if distance to center_part not to close 
				

				// check if intersection between surrounding_part and center_part is between thresholds
				// current old version uses vocabulary part region map (i.e. now activation points) to get area of parts: this must also be implemented in the calculatePartsIntersection
				std::set<InferredPart*> surrounding_part_support = support_parts_access.getOnlyDifferentInitialLayerSupportParts(surrounding_part);

				int intersection_value = InferenceSupportExtension::calculatePartsIntersection(central_part_support, surrounding_part_support);

				if (intersection_value < min_intersection_threhsold || intersection_value > max_intersection_threhsold) 
					continue;

				// get reference to the storage of pairwise statistics between center_part and surrounding_part
				PairwisePart &pairwise_statistics = central_part_statistics.sourunding_parts[surrounding_part.getCorrespondingVocabularyPartUUID()];
				
				// do initialization if empty
				if (pairwise_statistics.distribution_array.rows <= 0) {
					pairwise_statistics.distribution_array = cv::Mat::zeros(receptive_field_size, CV_32F);
				}
				
				// collect statistics at the appropriate offset from the center_part
				pairwise_statistics.distribution_array.at<float>(surrounding_part.getLocation() - neighborhood_region.tl()) += surrounding_part.getVVal();
			}
		}
	}
}

#include <queue>

std::vector<CooccurrenceStatistics::LocalMaximaDistribution>* CooccurrenceStatistics::findBestLocalMaximas(const cv::Mat& distribution_array, int max_neighbor_radius, int blur_size, float blur_sigma,	float intensity_threshold, int max_count_of_maximas, int local_distribution_map_size) {
	auto local_maximas = new std::vector<CooccurrenceStatistics::LocalMaximaDistribution>();

	// finding maxima
	//const int max_count_of_maximas = 4; //max number of maxima
	//const double max_val_threshold = 0.01; //set all map values <= ~*max to 0 (min_update_count_percent)
	//const bool individual_max = false; //max (see map_val_threshold) is set for each map individually
	//const int max_nbhood_mask = 5; //size of "distribution" (max_nbhood_mask)
	//const double max_sigma = 0.0; //if <= 0, then distribution is cut from the map
	//const int max_radius = 2; //min distance between two maxima and between each maximum and the center

	//typedef pair<double, ipoint2> queue_item_t;
	//typedef priority_queue<queue_item_t> queue_t;
	
	// convolve with gaussian mask and save to distribution_map
	cv::Mat smooth_distribution_map;
	cv::GaussianBlur(distribution_array, smooth_distribution_map, cv::Size2i(blur_size,blur_size), blur_sigma);

	// set all values below intensity_threshold to 0 (as noise reduction)
	smooth_distribution_map *= (smooth_distribution_map > intensity_threshold)/255; // divide by 255 since cv::compare will return 255 as true

	cv::Mat local_maxima_morph_filter = cv::Mat::ones(cv::Size(2 * max_neighbor_radius +1, 2 * max_neighbor_radius +1), smooth_distribution_map.type());

	cv::Mat local_maxima_matrix, non_plateau_mask;

	// find local maxima using cv::dilate (i.e. dilate is defined as max operation over all elements within the kernel)	
	cv::dilate(smooth_distribution_map, local_maxima_matrix, local_maxima_morph_filter);
		
	// filter out pixels that are equal to the local minimum ('plateaus')	
	cv::erode(smooth_distribution_map, non_plateau_mask, local_maxima_morph_filter);

	// local maxima is where each element is higher or the same as dilated and higher or the same as eroded value
	// TODO: maybe we can add eps to compares in order to migrate numerical instability ?
	cv::Mat local_maxima_mask = smooth_distribution_map >= local_maxima_matrix & smooth_distribution_map >= non_plateau_mask;
	
	struct LocalMaximaElement {
		cv::Point2i offset; 
		float value;
		LocalMaximaElement(cv::Point2i offset, float value) : offset(offset), value() {}

		bool operator<(const LocalMaximaElement& obj) const { return value < obj.value; }
	};
	// get local maxima locations from mask and add them to the 
	// priority queue for automatic sorting (by the local maxima value)
	std::priority_queue<LocalMaximaElement> sorted_local_maxima_points;
	for (int i = 0; i < local_maxima_mask.rows; ++i) {
		for (int j = 0; j < local_maxima_mask.cols; ++j) {
			if (local_maxima_mask.at<float>(i,j) > 0) {
				// j == x (column number) and i == y (row number)
				sorted_local_maxima_points.push(LocalMaximaElement(cv::Point2i(j,i), smooth_distribution_map.at<float>(i,j))); 
			}
		}
	}

	cv::Point2i local_distribution_map_radius(local_distribution_map_size / 2, local_distribution_map_size / 2);
	cv::Point2i distribution_center_point(smooth_distribution_map.rows / 2, smooth_distribution_map.cols / 2);

	// extract only max_number_of_maximas best locations and their local distribution maps
	while (sorted_local_maxima_points.empty() && (int)sorted_local_maxima_points.size() < max_count_of_maximas) {
		const LocalMaximaElement& local_maxima = sorted_local_maxima_points.top();

		cv::Point2i local_maxima_offset = local_maxima.offset - distribution_center_point;

		// ignore local maximas within center
		if (std::abs(cv::norm(local_maxima_offset)) >= max_neighbor_radius) {
		
			// extract distribution map around local maxima
			cv::Mat local_maxima_distributon = smooth_distribution_map(cv::Rect(local_maxima.offset - local_distribution_map_radius, 
																				local_maxima.offset + local_distribution_map_radius));

			// add to list of local maximas
			local_maximas->push_back(LocalMaximaDistribution(local_maxima_offset, local_maxima_distributon));
		}

		sorted_local_maxima_points.pop();
	}

	return local_maximas;
}
