
// cooccurrence_statistics
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _CORE_LEARNING_COOCCURRENCE_STAT_
#define _CORE_LEARNING_COOCCURRENCE_STAT_

#include <vector>
#include <unordered_map>

#include <opencv2/opencv.hpp>

#include "core/structures/parse_tree.h"

/// @addtogroup core
/// @{
/// @addtogroup learning
/// @{
/// @addtogroup learning_building_blocks
/// @{
/// @addtogroup learning_building_blocks_cooccurrence
/// @{



class CooccurrenceStatistics {
public:

	struct LocalMaximaDistribution {		
		cv::Point2i offset; // TODO: changing offset will affect uniqueness of vocabulary part candidate definition 
		cv::Mat distribution_matrix;

		LocalMaximaDistribution(const cv::Point2i& offset, const cv::Mat& distribution_matrix) : offset(offset), distribution_matrix(distribution_matrix) {}
	};

	/**
	* Pairwise statistics collected around a central part type for specific
	* surrounding vocabulary part type. 
	*/
	class PairwisePart {
		// list of local maxima points in the distribution_array
		std::vector<LocalMaximaDistribution>* local_maximas;
	public:
		// spatial distribution map (2D) of occurrence of the surrounding_part around the central part
		cv::Mat distribution_array;
	
		/**
		 * Returns list of local maximas of distribution array or 
		 * computes it using CooccurrenceStatistics::computeLocalMaximas(..) if empty.
		 */
		std::vector<LocalMaximaDistribution>& getLocalMaximas() {
			if (local_maximas == nullptr) {
				// TODO: read from configuration (HOW ?? if we do not know about parameters of CooccurrenceStatistics !!!)
				int max_neighbor_radius = 2;
				int blur_size = 5;
				float blur_sigma = 0.75;
				float intensity_threshold = 0; // TODO: obtained by looking at all possible parts, but current code does not communicate between pairwise parts; but not suited for incremental computation; replace with simple noise pruning
				int max_count_of_maximas = 4;
				int local_distribution_map_size = 5;

				local_maximas = CooccurrenceStatistics::findBestLocalMaximas(distribution_array, max_neighbor_radius, blur_size, blur_sigma, intensity_threshold, max_count_of_maximas, local_distribution_map_size);
			}
			return *local_maximas;
		}
	};
	/**
	* All the statistics collected around one type of vocabulary part. We collect 
	* statistics for each vocabulary part separately in the sourunding_parts map.
	*/
	struct CentralPart {
		// we collect statistics for each corresponding vocabulary part separately 
		// key = UUIDType of surrounding part
		// value = PairwisePartStatistics statistics with collected distribution map
		std::unordered_map<UUIDType, PairwisePart> sourunding_parts;
	};

private:
	// parameters:
	const double center_part_threshold;
	const int receptive_field_radius;

	const double min_intersection_threhsold;
	const double max_intersection_threhsold;

	// for each vocabulary part we collect statistics of parts found around it
	// key = UUIDType of center part
	// value = SinglePartStatistics with list of PairwisePartStatistics stat for each surrounding part
	std::unordered_map<UUIDType, CentralPart> central_parts;
public:

	/**
	 * Update co-occurrence statistics from the specific layer of input parse tree.
	 */
	void updateStatistics(const InferenceLayer& inference_tree);


	CentralPart& getCentralPartStatistics(const UUIDType& vocabulary_part_uuid) {
		return central_parts[vocabulary_part_uuid];
	}

private:
	/**
	 * Finds locations of best local maximas within input distribution_array (intensity map of CV_32F) using the following steps:
	 *  1. Smooths the input array using gaussian filter (parameter blur_size = 5 and blur_sigma = 0.75)
	 *  2. Thresholds values below parameter intensity_threshold by setting them to zero (as noise reduction)
	 *  3. Computes local maxima matrix using cv::dilate and cv::erode with kernle size of 2 * max_neighbor_radius+1 (max_neighbor_radius = 2).
	 *	4. Extracts local maximas from matrix and sorts them by their intensity value (i.e. smoothed value from input matrix)
	 *  5. Takes only N number of best maximas (N is defined by parameter max_count_of_maximas)
	 *     and extracts local map distribution for each maxima (matrix centered at local maxima with size max_nbhood_mask)
	 *  
	 *  Returns a list of best local maximas (max N number of them) with their offset from center and local map distribution
	 *  of size local_distribution_map_size.
	 * @param distribution_array: 2D input matrix of CV_32F
	 * @param max_neighbor_radius: max neighbor radius that are checked for best value
	 * @param blur_size and blur_sigma: gaussian blur parameters
	 * @param intensity_threshold: removes values with blur(input) < intensity_threshold by setting them to zero
	 * @param max_count_of_maximas: output parameter that defines how many local maximas can we return
	 * @param local_distribution_map_size: output parameter that defines the size of extracted local distribution map around each local maxima
	 */
	static std::vector<LocalMaximaDistribution>* findBestLocalMaximas(const cv::Mat& distribution_array,
																			int max_neighbor_radius,			// local maxima specific parameter
																			int blur_size, float blur_sigma,	// initial blurring parameters
																			float intensity_threshold,			// additional thresholding parameter
																			// output specific parameters
																			int max_count_of_maximas, int local_distribution_map_size);
};

/// @}
/// @}
/// @}
/// @}


#endif /* _CORE_LEARNING_COOCCURRENCE_STAT_ */
