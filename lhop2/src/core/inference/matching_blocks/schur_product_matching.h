// schur_product_matching
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _CORE_SCHUR_PRODUCT_MATCHING_
#define _CORE_SCHUR_PRODUCT_MATCHING_

#include "core/inference/abstract_inference.h"

#include "core/structures/subparts/subparts_extension.h"
#include "core/structures/shape_context/appearance_extension.h"

/// @addtogroup core
/// @{
/// @addtogroup inference
/// @{
/// @addtogroup inferrence_building_blocks
/// @{
/// @addtogroup schur_product
/// @{


// forward definitions
template <class T1, class T2, class T3> 
class SchurProductMatchingBlock;

template <class T1> 
class SchurProductBlock;
class SchurProductResult;

class ISchurProductFunction;
class IResponseCalculationBlock;
class ICandidateValidationBlock;

/////////////////////////////////////////////////////////////////
// Candidate matching implementations

template <class TSchurProductFunction,
		class TResponseCalculationBlock, // implements IResponseCalculationBlock
		class TCandidateValidationBlock> // implements ICandidateValidationBlock
class SchurProductMatchingBlock : public IInferredCandidatesFilteringBlock {
	SchurProductBlock<TSchurProductFunction> schur_product;

	// calculates the responses from the results of the schur product
	TResponseCalculationBlock response_calculation_block;
	// validates the candidates based on the matching results
	TCandidateValidationBlock candidate_validation_block;

	double convolution_link_threshold;
public:

	SchurProductMatchingBlock(SchurProductBlock<TSchurProductFunction>& schur_product,
								TResponseCalculationBlock& response_calculation_block,
								TCandidateValidationBlock& candidate_validation_block)
		: schur_product(schur_product), response_calculation_block(response_calculation_block), candidate_validation_block(candidate_validation_block) {
	}


	virtual InferredPartCandidateList* filterInferenceCandidates(InferredPartCandidateList* inference_candidates, const VocabularyTree& vocabulary);
};

///////////////////////////////////////////////////////////////
/// Schur product implementations used by candidate matching

//  TODO: can SchurProductMatchingBlock and SchurProductBlock be merged ? how much abstraction would be needed ? how would merging affect the performance (virtual calles etc) ?

/**
 * Class handles matching between a subpart from the vocabulary (with its distribution map) and inferred layer parts around
 * one specific inferred central part.
 *
 * A response value (obtained by the schur_product_function.getValue()) of the inferred part at specific location is multiplied with the
 * value from the distribution matrix (at the corresponding location) and if the computed value is above certain threshold (convolution_threshold)
 * we retain the part. We return the retained parts (and theirs schur product values) in the valid_subparts list of the SchurProductResult class.
 */ 
template <class TSchurProductFunction> // implements ISchurProductFunction
class SchurProductBlock {
	double convolution_threshold;
	// function that extracts appropriate value for schur product from the inferred part
	TSchurProductFunction schur_product_function;	

	// access to the appearance model of the vocabulary
	VocabularyAppearanceAccess* vocabulary_appearance_access;
public:
	SchurProductBlock(TSchurProductFunction& schur_product_function, double convolution_threshold) 
		: schur_product_function(schur_product_function), convolution_threshold(convolution_threshold) {
	}

	SchurProductResult* match(InferenceLayer& inferred_layer, InferredPart& inferred_center_part, VocabularySubpartData& candidate_vocabulary_part);

	void setVocabularyAppearacnceAccess(VocabularyAppearanceAccess& vocabulary_appearance_access) {
		this->vocabulary_appearance_access = &vocabulary_appearance_access;
	}
};

struct SchurProductPart {
	InferredPart* part;
	double schur_product_value;

	// default constructor
	SchurProductPart(InferredPart* part, double value) : part(part), schur_product_value(value) {}
	SchurProductPart(const SchurProductPart& obj) : part(obj.part), schur_product_value(obj.schur_product_value) {}

	// move operation needed for the std::vector and std::remove_if
	SchurProductPart& operator=(SchurProductPart&& obj) {
		part = obj.part;
		schur_product_value = obj.schur_product_value;
		return *this;
	}

};

/**
 * Result return by the SchurProductBlock that contains:
 * - vocabulary_subpart: pointer to the vocabulary subpart (with all the link informations)
 * - valid_subparts: an array of SchurProductPart (containing InferredPart and schur_product_value) parts that survived the convolution thresholding
 */ 
class SchurProductResult {	
	SchurProductPart* max_subpart;
public:
	const VocabularySubpartData& vocabulary_subpart;
	std::vector<SchurProductPart> valid_subparts;

	SchurProductResult(const VocabularySubpartData& vocabulary_subpart) : vocabulary_subpart(vocabulary_subpart), max_subpart(nullptr) {}
	SchurProductResult(const SchurProductResult& obj)
		: max_subpart(nullptr), vocabulary_subpart(obj.vocabulary_subpart), valid_subparts(obj.valid_subparts) {}
	SchurProductResult(SchurProductResult&& obj)
		: max_subpart(nullptr), vocabulary_subpart(obj.vocabulary_subpart), valid_subparts(std::move(obj.valid_subparts)) {}

#if 0
	SchurProductResult& operator=(SchurProductResult&& obj) {
		max_subpart = nullptr;
		vocabulary_subpart = std::move(obj.vocabulary_subpart);
		valid_subparts = std::move(obj.valid_subparts);
	}
#endif

	SchurProductPart* getMaxSubpart() {
		if (max_subpart == nullptr && valid_subparts.size() > 0) {
			max_subpart = &valid_subparts[0];

			for (auto iter = valid_subparts.begin(); iter != valid_subparts.end(); ++iter) {
				if (max_subpart->schur_product_value < iter->schur_product_value) {
					max_subpart = &*iter;
				}
			}
		}

		return max_subpart;
	}
};


///////////////////////////////////////////////////////////////
/// Types of schur product function used to extract one specific value from one part
/// One of this classes must be supplied to the SchurProductMatchingBlock as template parameter 

/// \name SchurProductFunction 
/// @{
struct ISchurProductFunction {
	virtual void initialize(const InferredPart& current_center_part, const VocabularySubpartData& vocabulary_subpart) {}
	virtual double getValue(const InferredPart& central_part, const InferredPart& part) const = 0;
};

/**
 * Function always returns 1
 */
struct IdentitySchurProductFunction : ISchurProductFunction {	
	virtual double getValue(const InferredPart& central_part, const InferredPart& part) const;
};

/**
 * Function returns R_RESPONSE value
 */
struct R_ResponseSchurProductFunction : ISchurProductFunction {
	virtual double getValue(const InferredPart& central_part, const InferredPart& part) const;
};

/**
 * Function returns G_RESPONSE value
 */
struct V_ResponseSchurProductFunction : ISchurProductFunction {
	virtual double getValue(const InferredPart& central_part, const InferredPart& part) const;
};

/**
 * Function returns G_RESPONSE multiplied by the Gaussian ????
 */
class G_ResponseSchurProductFunction : public ISchurProductFunction {
	double g_response_var_factor;
	// set by the initialize
	double quotient;
    normal_distribution1 dist;
public:	
	G_ResponseSchurProductFunction(double g_response_var_factor) : g_response_var_factor(g_response_var_factor) {
	}
	virtual void initialize(const InferredPart& current_center_part, const VocabularySubpartData& vocabulary_subpart);
	virtual double getValue(const InferredPart& central_part, const InferredPart& part) const;
};

/**
 * Function returns G_RESPONSE/RR_RESPONSE value
 */
struct Simple_G_ResponseSchurProductFunction : ISchurProductFunction  {
	virtual double getValue(const InferredPart& part) const;
};

/// @}

///////////////////////////////////////////////////////////////
/// Response calculations used by candidate matching
/// @{
class IResponseCalculationBlock {
public:
	virtual ResponsesArray& updateResponse(ResponsesArray& current_responses, const std::vector<SchurProductResult*>& schur_product_results, const InferredPartCandidate& part_candidate) = 0;
};

class R_ResponseCalculationBlock : public IResponseCalculationBlock  {
	double r_response_pow;
public:
	R_ResponseCalculationBlock(double r_response_pow) : r_response_pow(r_response_pow) {
	}
	virtual ResponsesArray& updateResponse(ResponsesArray& current_responses, const std::vector<SchurProductResult*>& schur_product_results, const InferredPartCandidate& part_candidate);
};

class RR_ResponseCalculationBlock : public IResponseCalculationBlock  {
public:
	virtual ResponsesArray& updateResponse(ResponsesArray& current_responses, const std::vector<SchurProductResult*>& schur_product_results, const InferredPartCandidate& part_candidate);
};


class G_ResponseCalculationBlock : public IResponseCalculationBlock  {	
	double g_response_pow;
	bool g_response_operation;

	// use RR_ResponseCalculationBlock to calculate first the RR_RESPONSE 
	// as it is required for the calculation of the G_RESPONSE
	RR_ResponseCalculationBlock rr_response_block;
public:
	G_ResponseCalculationBlock(RR_ResponseCalculationBlock& rr_response_block, double g_response_pow, bool g_response_operation) 
		: rr_response_block(rr_response_block), g_response_pow(g_response_pow), g_response_operation(g_response_operation) {
	}

	virtual ResponsesArray& updateResponse(ResponsesArray& current_responses, const std::vector<SchurProductResult*>& schur_product_results, const InferredPartCandidate& part_candidate);
};

/**
 * Calculates the R_RESPONSE, G_RESPONSE and RR_RESPONSE (implicitly when G_RESPONSE is calculated) from the schur product results. 
 */ 
class ResponseCalculationBlock : IResponseCalculationBlock  {
	R_ResponseCalculationBlock r_response_block;
	G_ResponseCalculationBlock g_response_block;
	// RR_ResponseCalculationBlock is already called by the g_response_block
public:
	ResponseCalculationBlock(R_ResponseCalculationBlock r_response_block, G_ResponseCalculationBlock g_response_block) 
		: r_response_block(r_response_block), g_response_block(g_response_block) {
	}
	virtual ResponsesArray& updateResponse(ResponsesArray& current_responses, const std::vector<SchurProductResult*>& schur_product_results, const InferredPartCandidate& part_candidate) {
		
		current_responses = r_response_block.updateResponse(current_responses, schur_product_results, part_candidate);
		current_responses = g_response_block.updateResponse(current_responses, schur_product_results, part_candidate);

		return current_responses;
	}
};
/// @}

///////////////////////////////////////////////////////////////
/// Candidate validation check based on calculated response values
/// @{
class ICandidateValidationBlock {
public:
	virtual bool isCandidateValid(const InferredPartCandidate& part_candidate) = 0;
};

/**
 * Passes validation check if calculated R_RESPONSE is higher then R_RESPONSE threshold of the candidate vocabulary part
 */ 
class R_ResponseCandidateValidationBlock : public ICandidateValidationBlock {
public:
	virtual bool isCandidateValid(const InferredPartCandidate& part_candidate);
};

/**
 * Passes validation check if calculated RR_RESPONSE is higher then RR_RESPONSE threshold of the candidate vocabulary part
 */ 
class RR_ResponseCandidateValidationBlock : public ICandidateValidationBlock {
public:
	virtual bool isCandidateValid(const InferredPartCandidate& part_candidate);
};

/**
 * Passes validation check if calculated G_RESPONSE is higher then G_RESPONSE threshold of the candidate vocabulary part
 */ 
class G_ResponseCandidateValidationBlock : public ICandidateValidationBlock {
public:
	virtual bool isCandidateValid(const InferredPartCandidate& part_candidate);
};

/**
 * Passes validation check if calculated R_RESPONSE, RR_RESPONSE and G_RESPONSE are all higher then 
 * R_RESPONSE, RR_RESPONSE and G_RESPONSE thresholds of the candidate vocabulary part
 */ 
class ResponseCandidateValidationBlock : public ICandidateValidationBlock {
	R_ResponseCandidateValidationBlock r_response_validation;
	RR_ResponseCandidateValidationBlock rr_response_validation;
	G_ResponseCandidateValidationBlock g_response_validation;
public:
	virtual bool isCandidateValid(const InferredPartCandidate& part_candidate) {
		return r_response_validation.isCandidateValid(part_candidate) &&
				rr_response_validation.isCandidateValid(part_candidate) &&
				g_response_validation.isCandidateValid(part_candidate);
	}
};
/// @}

// include template definitions
#include "core/inference/matching_blocks/schur_product_matching.hpp"

/// @}
/// @}
/// @}
/// @}

#endif /* _CORE_SCHUR_PRODUCT_MATCHING_ */
