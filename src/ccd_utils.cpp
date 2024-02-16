//
// Created by Henry Hollis on 1/28/24.
//

#include "ccd_utils.h"
#include <iostream>
#include <cmath>
#include <stdexcept>
// Definition of the constant matrix
const std::vector<Vector> ccd_utils::refCor = {
        {1.0000000,   0.77547090,  0.72492855,  0.27817942, -0.63637681, -0.60375141, -0.8614806, -0.7471112, -0.59455286, -0.8234182, -0.9146447, -0.8473980}, //ARNTL
        {0.7754709,   1.00000000,  0.63439613,  0.07402797, -0.62632300, -0.34987550, -0.7461844, -0.6450780, -0.70865725, -0.7845410, -0.7654845, -0.7983427}, //NPAS2
        {0.7249286,   0.63439613,  1.00000000,  0.06541974, -0.59727560, -0.30024636, -0.6031795, -0.6364953, -0.56958405, -0.7144612, -0.6455111, -0.7595101}, //CLOCK
        {0.2781794,   0.07402797,  0.06541974,  1.00000000, -0.01245765, -0.72253596, -0.4099044, -0.1411756,  0.25538496, -0.0252816, -0.3401805, -0.0781101}, //CRY1
        {-0.6363768, -0.62632300, -0.59727560, -0.01245765,  1.00000000,  0.28367324,  0.6234166,  0.6454257,  0.59510653,  0.6712806,  0.6618797,  0.7597038}, //CRY2
        {-0.6037514, -0.34987550, -0.30024636, -0.72253596,  0.28367324,  1.00000000,  0.6772739,  0.4242223, -0.06776682,  0.3366267,  0.6955807,  0.3810191}, //NR1D1
        {-0.8614806, -0.74618443, -0.60317949, -0.40990436,  0.62341661,  0.67727389,  1.0000000,  0.7132144,  0.52923596,  0.7673822,  0.9111478,  0.7487607}, //NR1D2
        {-0.7471112, -0.64507795, -0.63649530, -0.14117556,  0.64542570,  0.42422234,  0.7132144,  1.0000000,  0.60794410,  0.7467579,  0.7732704,  0.7756198}, //PER1
        {-0.5945529, -0.70865725, -0.56958405,  0.25538496,  0.59510653, -0.06776682,  0.5292360,  0.6079441,  1.00000000,  0.7868302,  0.5543211,  0.7530874}, //PER2
        {-0.8234182, -0.78454102, -0.71446119, -0.02528160,  0.67128060,  0.33662668,  0.7673822,  0.7467579,  0.78683019,  1.0000000,  0.8117621,  0.8738338}, //PER3
        {-0.9146447, -0.76548454, -0.64551113, -0.34018047,  0.66187971,  0.69558073,  0.9111478,  0.7732704,  0.55432112,  0.8117621,  1.0000000,  0.8443479}, //DBP
        {-0.8473980, -0.79834269, -0.75951011, -0.07811010,  0.75970381,  0.38101906,  0.7487607,  0.7756198,  0.75308740,  0.8738338,  0.8443479,  1.0000000} //TEF
};

double ccd_utils::calcCCDsimple(const std::vector<Vector> &ref,
                     const std::vector<Vector> &emat,
                      bool scale) {

    std::vector<Vector> cormat = calcCorMat(emat); //calculate cormat of expression matrix
    if (cormat.size() != ref.size() || cormat.empty() || ref[0].size() != cormat[0].size()) {
        std::cerr << "cormat size:" << cormat.size() << std::endl;
        throw std::invalid_argument("Matrices must be of the same size for calcCCDsimple");
    }
    size_t numRows = ref.size();
    double upperTriDiff = 0.0;
    //loop through triangular matrix and accumulate the difference between entries of ref and cormat
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = i; j < numRows; ++j) {
            upperTriDiff += pow(ref[i][j] - cormat[i][j], 2);
        }
    }
    double ccd = sqrt(upperTriDiff);

    if (!ref.empty() && !ref[0].empty()) {
        // Number of columns is the size of any row (assuming all rows have the same size)
        size_t numColumns = ref[0].size();
        if (scale) {
            size_t nPairs = choose(numColumns, 2);
            ccd /= static_cast<double>(nPairs);
        }
    } else
        std::cerr << "Matrix ref is empty or has empty rows." << std::endl;

    return ccd;
}

std::vector<Vector> ccd_utils::calcCorMat(const std::vector<Vector> &ref) {
   // size_t numCols = ref[0].size();
    size_t numRows = ref.size();

    // Initialize the correlation matrix with zeros
    std::vector<Vector> correlationMatrix(numRows, Vector(numRows, 0.0));

    // Calculate pairwise correlations
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = i; j < numRows; ++j) {
            if (i == j) {
                // Diagonal elements (correlation with itself) are always 1
                correlationMatrix[i][j] = 1.0;
            } else {
                // Off-diagonal elements are calculated using the cor function
                Vector rank_x = rankVector(const_cast<Vector &>(ref[i])); //rank vectors first to calculate spearman coeff
                Vector rank_y = rankVector(const_cast<Vector &>(ref[j]));

                correlationMatrix[i][j] = correlationMatrix[j][i] = cor(rank_x, rank_y);
            }
        }
    }

    return correlationMatrix;
}


long ccd_utils::choose(size_t n, int k) {
    if (0 == k)
        return 1;
    return (n * choose(n - 1, k - 1)) / k;
}

/**
 * Spearman Correlation Code
 * Inspired from https://www.geeksforgeeks.org/program-spearmans-rank-correlation/
 */
// Utility Function to print vector
void ccd_utils::printVector(const std::vector<size_t> &X) {
    for (auto i: X)
        std::cout <<i<<" ";
    std::cout << std::endl;
}

// Function returns the rank vector
// of the set of observations
Vector ccd_utils::rankVector(Vector &X) {

    int N = X.size();

    // Rank Vector
    Vector Rank_X(N);

    for(int i = 0; i < N; i++)
    {
        int r = 1, s = 1;

        // Count no of smaller elements
        // in 0 to i-1
        for(int j = 0; j < i; j++) {
            if (X[j] < X[i] ) r++;
            if (X[j] == X[i] ) s++;
        }

        // Count no of smaller elements
        // in i+1 to N-1
        for (int j = i+1; j < N; j++) {
            if (X[j] < X[i] ) r++;
            if (X[j] == X[i] ) s++;
        }

        // Use Fractional Rank formula
        // fractional_rank = r + (n-1)/2
        Rank_X[i] = r + (s-1) * 0.5;
    }

    // Return Rank Vector
    return Rank_X;
}

// function that returns
// Pearson correlation coefficient.
float ccd_utils::cor(const Vector &X, const Vector &Y) {
    int n = X.size();
    float sum_X = 0, sum_Y = 0,
            sum_XY = 0;
    float squareSum_X = 0,
            squareSum_Y = 0;

    for (int i = 0; i < n; i++)
    {
        // sum of elements of array X.
        sum_X = sum_X + X[i];

        // sum of elements of array Y.
        sum_Y = sum_Y + Y[i];

        // sum of X[i] * Y[i].
        sum_XY = sum_XY + X[i] * Y[i];

        // sum of square of array elements.
        squareSum_X = squareSum_X +
                      X[i] * X[i];
        squareSum_Y = squareSum_Y +
                      Y[i] * Y[i];
    }

    // use formula for calculating
    // correlation coefficient.
    float corr = (float)(n * sum_XY -
                         sum_X * sum_Y) /
                 sqrt((n * squareSum_X -
                       sum_X * sum_X) *
                      (n * squareSum_Y -
                       sum_Y * sum_Y));

    return corr;
}




std::vector<Vector> ccd_utils::sliceColumns(const std::vector<Vector> &matrix, const std::vector<size_t> &columns)  {
    if (matrix.empty() || columns.empty()) {
        return {}; // Return an empty matrix if either the matrix or columns are empty
    }

    std::vector<Vector> result(matrix.size(), std::vector<double>(columns.size()));

    for (size_t i = 0; i < matrix.size(); ++i) {
        std::transform(columns.begin(), columns.end(), result[i].begin(),
                       [&matrix, i](size_t col) { return matrix[i][col]; });
    }

    return result;
}