/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// shape_context_extension

#pragma once
#ifndef _CORE_SHAPE_CONTXT_EXTENSION_
#define _CORE_SHAPE_CONTXT_EXTENSION_

#include "core/structures/parse_tree.h"
#include "core/structures/vocabulary.h"
#include "core/structures/attached_extension.h"

#include "core/structures/subparts/support_extension.h"

/// @addtogroup core
/// @{
/// @addtogroup main_structures
/// @{

/// @addtogroup ext
/// @{

/// @addtogroup ext_inference_tree
/// @{

////////////////////////////////////////////////////////////////////////////////
//// Implementation for Belongie histogram caching as extension

// forward definitions
class IBelongieHistogramAccess;
class BelongieHistogramData;

class BelongieHistogramAccess;
class BelongieHistogramModifier;
class BelongieHistogramExtension;


class IBelongieHistogramAccess : public IExtensionAccess {
public:
	virtual BelongieHistogramData* getHistogram(const IAttachableClass& part, InferenceTree& inference_tree) = 0;
};


/**
 * 
 */ 
struct BelongieHistogramData {
	cv::Mat histogram;
	cv::Point2f center_point;

	BelongieHistogramData(cv::Mat histogram, cv::Point2f center_point) 
		: histogram(histogram), center_point(center_point) {

	}
};


/**
 * 
 */
class BelongieHistogramExtension : public SimpleAbstractExtension<BelongieHistogramExtension,
																	BelongieHistogramAccess,
																	BelongieHistogramModifier,
																	InferredPart> {
public:
	ExtensionHolder& holder;

	std::unordered_map<cv::Point, BelongieHistogramData*> histogram_map;

	// default constructor
	BelongieHistogramExtension(ExtensionHolder& holder) 
		: holder(holder) {
	}
	
	// moves all data into dst_extension and clear this one
	virtual void moveTo(IExtension& dst_extension); 
	
	class Factory: public IFactory {
	public:
		virtual IExtension* newInstance(ExtensionHolder& holder) const { return new BelongieHistogramExtension(holder); }
	};
};

////////////////////////////////////////////////
/// Access and modifier classes for Belongie shape histogram

/**
 * TODO
 */ 
class BelongieHistogramAccess : public IBelongieHistogramAccess {	
	BelongieHistogramExtension& ext;
public:
	BelongieHistogramAccess(BelongieHistogramExtension& ext) : ext(ext) {
	}

	virtual BelongieHistogramData* getHistogram(const IAttachableClass& part, InferenceTree& inference_tree);
};

/**
 * Used only for deletion of data when reference is delete 
 */
class BelongieHistogramModifier : public IExtensionModifier {	
	BelongieHistogramExtension& ext;
public:
	BelongieHistogramModifier(BelongieHistogramExtension& ext) : ext(ext) {
	}

	virtual void deleteAllReferences(const std::vector<IAttachableClass*>& reference_to);
};

/// @}

/// @addtogroup ext_vocabulary_tree
/// @{

////////////////////////////////////////////////////////////////////////////////
//// 


// forward definitions
class IShapeContextAccess;
class IShapeContextModifier;
class ShapeContextAccess;
class ShapeContextModifier;
class ShapeContextExtension;
class ShapeContextGeometryData;

////////////////////////////////////////////////////////////////////////////////
//// Implementation shape context attached to the vocabulary subpart data


class IShapeContextAccess : public IExtensionAccess {
public:
	virtual ShapeContextGeometryData* getGeometry(const IAttachableClass& subpart_data) = 0;
};


class IShapeContextModifier : public IExtensionModifier {
public:
	virtual void deleteAllReferences(const std::vector<IAttachableClass*>& reference_to) = 0;

	virtual void insertGeometry(const IAttachableClass& subpart_data, ShapeContextGeometryData* geometry) = 0;
};


/**
 * 
 */ 
struct ShapeContextGeometryData {
	// save one histogram for each leaf i.e. path between root part and one if its supports
	std::map<PartSupportPath, BelongieHistogramData*> geometry;

	ShapeContextGeometryData(const std::map<PartSupportPath, BelongieHistogramData*>& geometry) 
		: geometry(geometry) {
	}
	~ShapeContextGeometryData() {
		for (auto iter = geometry.begin(); iter != geometry.end(); ++iter) delete iter->second;
	}
};


/**
 * 
 */
class ShapeContextExtension : public SimpleAbstractExtension<ShapeContextExtension,
																ShapeContextAccess,
																ShapeContextModifier,
																VocabularySubpartData>{
public:
	ExtensionHolder& holder;

	std::unordered_map<UUIDType, ShapeContextGeometryData*> geometry_data;

	// default constructor
	ShapeContextExtension(ExtensionHolder& holder) 
		: holder(holder) {
	}

	// moves all data into dst_extension and clear this one
	virtual void moveTo(IExtension& dst_extension); 
	
	virtual AbstractSerializer::IFactory* getSerializer() const;

	class Factory: public IFactory {
	public:
		virtual IExtension* newInstance(ExtensionHolder& holder) const { return new ShapeContextExtension(holder); }
	};

	// TODO: not implemented yet !!
	static void* calculateGeometryDistance(const ShapeContextGeometryData& geometry_a, const ShapeContextGeometryData& geometry_b, bool calc_benergy, int max_matching_size);
};

////////////////////////////////////////////////
/// Access and modifier classes for shape context information

/**
 * TODO
 */ 
class ShapeContextAccess : public IShapeContextAccess {	
	ShapeContextExtension& ext;
public:
	ShapeContextAccess(ShapeContextExtension& ext) : ext(ext) {
	}

	virtual ShapeContextGeometryData* getGeometry(const IAttachableClass& subpart_data) {
		return ext.geometry_data[subpart_data.getUUID()];
	}
};

class ShapeContextModifier : public IShapeContextModifier {	
	ShapeContextExtension& ext;
public:
	ShapeContextModifier(ShapeContextExtension& ext) : ext(ext) {
	}

	virtual void deleteAllReferences(const std::vector<IAttachableClass*>& reference_to);

	virtual void insertGeometry(const IAttachableClass& subpart_data, ShapeContextGeometryData* geometry)  {
		ext.geometry_data[subpart_data.getUUID()] = geometry;
	}
};


/// @}



/// @}

/// @}
/// @}

#endif /* _CORE_SHAPE_CONTXT_EXTENSION_ */

