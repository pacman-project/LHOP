// Interface to "munkres.cpp"
//

#pragma once

#include <vector>

// Solves optimal assignment problem using ...
// 'data' points to an array of 'rows' x 'cols' elements representing cost matrix 
// (elements of matrix are written "row"-wise, i.e. first 'cols' elements represent
// the first row, etc.)
// Result is written to 'result' which should be a pointer to a chunck of memory of 
// sufficient size; it represents the assignment matrix filled with ones and zeroes,
// ones represent matching.
void best_matching(int* result, int rows, int cols, double* data);

