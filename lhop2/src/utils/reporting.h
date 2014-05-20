
#pragma once
#ifndef _UTILS_REPORTING_
#define _UTILS_REPORTING_

/**
 * Reporting interface
 */
class IReporting {
public:
	virtual void reportImage(const std::string& name, const img& image) = 0;
	virtual void reportCvMat(const std::string& name, const cv::Mat& mat) = 0;
	virtual void reportMatrix(const std::string& name, const matrix& mat) = 0;
	virtual void reportVector(const std::string& name, const std::vector<??>& mat) = 0;
};
 

#endif


