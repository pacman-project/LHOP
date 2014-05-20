
#include "matching.h"
#include "munkres.h"

void best_matching(int* result, int rows, int cols, double* data)
{
    MatrixMunkres<double> matrix(rows, cols);

    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c) 
            matrix(r, c) = *(data++);

    Munkres m;

    m.solve(matrix);

    for (int r = 0; r < rows; ++r) 
        for (int c = 0; c < cols; ++c)
            *(result++) = ((matrix(r, c) == 0) ? 1 : 0);
}

