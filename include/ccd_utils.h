//
// Created by Henry Hollis on 1/28/24.
//

#ifndef LEIDEN_CCD_CCD_UTILS_H
#define LEIDEN_CCD_CCD_UTILS_H
#include <vector>
#include <unordered_map>

class ccd_utils {
public:
    static std::vector<double> sliceColumns(const std::vector<double> &matrix, const std::vector<size_t> &columnsToAccess,
                                            size_t nrow, size_t ncol);

    static double calcCCDsimple(const std::vector<double> &ref, int num_ref_rows,
                                const std::vector<double> &emat, size_t emat_row, size_t emat_col,
                                bool scale);

    static std::vector<double> calcCorMat(const std::vector<double> &rect, size_t numRows, size_t numCols);
    static long choose(size_t n, int k);

    static std::vector<double> rankVector(const std::vector<double> &X);

    static double cor(const std::vector<double> &X, const std::vector<double> &Y);

    static void sumColumnsByGroup(const std::vector<double>& matrix, size_t rows, size_t cols, const std::vector<size_t>& membership, std::vector<double>& result);
    };
#endif //LEIDEN_CCD_CCD_UTILS_H
