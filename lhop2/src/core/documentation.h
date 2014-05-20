
// documentation
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _CORE_DOCUMENTATION_
#define _CORE_DOCUMENTATION_

/**
 * @defgroup core Core
 * @brief Core definitions of structures, feature extraction, inference and learning
 */ 
/**
 * @defgroup inference Inference
 * @brief Inference classes with its main classes and building blocks
 */ 
/**
 * @defgroup main_structures Main structures
 * @brief Main structures such as VocabularyTree, InferenceTree etc.
 */
/**
 * @defgroup func Attachable functionalities (extensions)  
 * @brief Attached functionalities are extensions to the VocabularyTree and InferenceTree 
 * to provide additional information such as subparts, indexing, appearance (shape context) etc.
 */
/**
 * @defgroup func_inference_tree Inference tree related
 * @brief 
 */
/**
 * @defgroup func_vocabulary_tree Vocabulary tree related
 * @brief 
 */
/**
 * @defgroup inference_tree Inference tree
 * @brief Structures specific for the parse tree
 */
/**
 * @defgroup vocabulary_tree Vocabulary tree
 * @brief Structures specific for the vocabulary 
 */
/**
 * @defgroup feature_extraction Feature extraction
 * @brief 
 */ 
/**
 * @defgroup inference_building_blocks Building blocks
 * @brief All the building blocks used in the inference process
 */
/**
 * @defgroup shape_check Shape checking
 * @brief Shape context specific blocks
 */
/**
 * @defgroup schur_product Schur product
 * @brief Main blocks for implementing matching using the schur product
 */
/**
 * @defgroup indexing Indexing
 * @brief Indexing blocks
 */
/**
 * @defgroup serialization Serialization
 * @brief Serializers for all the structures
 */
/**
 * @defgroup input_output Input/output structures
 * @brief Wrapper around structures that can be written and read from any kind of source.
 */
/**
 * @defgroup preprocessing Preprocessing 
 * @brief Preprocessing blocks for input wrappers
 */
/**
 * @defgroup postprocessing Postprocessing
 * @brief Postprocessing blocks for output wrappers
 */

/**
 * @defgroup postprocessing Postprocessing
 * @brief Postprocessing blocks for output wrappers
 */ 
/**
 * @defgroup toolset Toolset
 * @brief Toolsets provide calling mechanism for the executable to run the actual code. 
 */
/**
 * @defgroup file_input_output File input/output structures
 * @brief Wrapper around structures for reading/writing using files
 */
 /**
 * @defgroup learning Learning
 * @brief Learning related classes
 */
 
 /**
 * @defgroup learning_building_blocks Building blocks
 * @brief All building blocks for learning
 */
 
 /**
 * @defgroup learning_building_blocks_cooccurrence Cooccurrence learning
 * @brief Vocabulary learning with cooccurrence statistics
 */
 
 /**
 * @defgroup learning_building_blocks_coverage Coverage optimization
 * @brief Vocabulary optimization/pruning using coverage metrics
 */
 
 /**
 * @defgroup func_support Initial layer (ly1) support
 * @brief Layer 1 support calculation of part
 */
 
#endif /* _CORE_DOCUMENTATION_ */
