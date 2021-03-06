//
// Created by Admin on 10/15/20.
//

#include "ThomasAlgorithm.h"
#include <iostream>
#include <vector>
//#include <Kokkos_Core.hpp>

using namespace std;
//using namespace matplot;

//void ThomasAlgorithm::setTMatrix(vector <vector<double>> &matrix) {
//    ThomasAlgorithm::Tmatrix = matrix;
//}

double ThomasAlgorithm::getRandom(int seed) {
    srand(10);
    return (rand() % 10) + 1;
}

void ThomasAlgorithm::createRandomTMatrix(vector<vector<double> > &matrix,
                                          int Mdim, int Ndim, bool stdmatrix) {
    // Mdim =Ndim because we have a suqre matrix always
    for (int i = 0; i < Mdim; i++) {
        double a, b, c;
        for (int j = 0; j < Ndim; j++) {
            if (i == j) {
                matrix[i][j] = !stdmatrix ? getRandom(j) : 2.0;
            }
            if (i - j == 1) {
                matrix[i][j] = !stdmatrix ? getRandom(j) : 1.0;;
            }
            if (j - i == 1) {
                matrix[i][j] = !stdmatrix ? getRandom(j) : 1.0;;
            }

        }
    }
}
void ThomasAlgorithm::createMatrix(vector<vector<double> > &matrix,
                                          int Mdim, int Ndim, double arow,double brow,double crow) {
    // Mdim =Ndim because we have a sqare matrix always
    for (int i = 0; i < Mdim; i++) {
        double a, b, c;
        for (int j = 0; j < Ndim; j++) {
            if (i == j) {
                matrix[i][j] = brow;
            }
            if (i - j == 1) {
                matrix[i][j] = arow;;
            }
            if (j - i == 1) {
                matrix[i][j] = crow;
            }

        }
    }
}

void
ThomasAlgorithm::backSubstitution(vector<double> &upperDiagonal, vector<double> &rho, vector<double> &results,
                                  int size) {
    for (int i = size - 1; i >= 0; i--) {
        results[i] = rho[i] - (upperDiagonal[i] * (i + 1 == size ? 0.0 : results[i + 1]));
    }
}

void
ThomasAlgorithm::forwardSweep(vector<vector<double>> &matrix, vector<double> &rmatrix, vector<double> &upperDiagonal,
                              vector<double> &rho, int Mdim, int Ndim) {
    for (int i = 0; i < Mdim; i++) {
        double a, b, c;
        for (int j = 0; j < Ndim; j++) {
            if (i == j) {
                b = matrix[i][j];
            }
            if (i - j == 1) {
                a = matrix[i][j];
            }
            if (j - i == 1) {
                c = matrix[i][j];
            }
        }
        upperDiagonal[i] = computeUpperDiagonal(c, b, a, i == 0 ? 0.0 : upperDiagonal[i - 1]);
        rho[i] = computeRho(rmatrix[i], b, a, i == 0 ? 0.0 : upperDiagonal[i - 1], i == 0 ? 0.0 : rho[i - 1]);

    }
}

double ThomasAlgorithm::computeRho(double r, double b, double a, double prevUpperDiagonal, double prevRho) {
    if (b - (a * prevUpperDiagonal) == 0) throw "Instability detected. Perform pivoting or check you matrix";
    return (r - (a * prevRho)) / (b - (a * prevUpperDiagonal));
}

double ThomasAlgorithm::computeUpperDiagonal(double c, double b, double a, double prevUpperDiagonal) {
    if (b - (a * prevUpperDiagonal) == 0) throw "Instability detected. Perform pivoting or check you matrix";
    return c / (b - (a * prevUpperDiagonal));
}

void ThomasAlgorithm::computeThomas(vector<vector<double> > &matrix, vector<double> &rmatrix,
                                    vector<double> &upperDiagonal, vector<double> &rho, vector<double> &results,
                                    int Mdim, int Ndim) {
    forwardSweep(matrix, rmatrix, upperDiagonal, rho, Mdim, Ndim);
    backSubstitution(upperDiagonal, rho, results, Mdim);
}

void ThomasAlgorithm::print(vector<double> &matrix) {
    for (auto i = matrix.begin(); i != matrix.end(); ++i)
        std::cout << *i << ' ';
    std::cout << "\n";
}

void ThomasAlgorithm::print(vector<vector<double>> &matrix) {
    for (auto i = matrix.begin(); i != matrix.end(); ++i) {
        for (auto j = i->begin(); j != i->end(); ++j) {
            std::cout << *j << ' ';
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}


int ThomasAlgorithm::DETGTRI(vector<vector<double> > &matrix) {
    //difficult to implement.
}