//
// Created by Godfred on 10/15/20.
//


#ifndef HW2_CFD_THOMASALGORITHM_H
#define HW2_CFD_THOMASALGORITHM_H

#include <vector>

using namespace std;

class ThomasAlgorithm {
public:
    vector<vector<double>> *Tmatrix;

    void createRandomTMatrix(vector<vector<double>> &matrix, int Mdim, int Ndim, bool stdmatrix = true);

    void createMatrix(vector<vector<double> > &matrix,
                            int Mdim, int Ndim, double brow, double arow, double crow);

    void backSubstitution(vector<double> &upperDiagonal, vector<double> &rho, vector<double> &results, int size);

    void
    backSubstitutionParallel(vector<double> &upperDiagonal, vector<double> &rho, vector<double> &results, int size);

    void forwardSweep(vector<vector<double>> &matrix, vector<double> &rmatrix, vector<double> &upperDiagonal,
                      vector<double> &rho, int Mdim,
                      int Ndim);

    void computeThomas(vector<vector<double>> &matrix, vector<double> &rmatrix, vector<double> &upperDiagonal,
                       vector<double> &rho, vector<double> &results, int Mdim,
                       int Ndim);

    void
    forwardSweepParallel(vector<vector<double>> &matrix, vector<double> &rmatrix, vector<double> &upperDiagonal,
                         vector<double> &rho, int Mdim,
                         int Ndim);

    void setTMatrix(vector<vector<double>> &Tmatrin);

    double computeUpperDiagonal(double c, double b, double a, double prevUpperDiagonal);

    double computeRho(double r, double b, double a, double prevUpperDiagonal, double prevRho);

    double getRandom(int seed);

    void print(vector<vector<double>> &matrix);

    void print(vector<double> &matrix);

    int DETGTRI(vector<vector<double>> &matrix);

};


#endif //HW2_CFD_THOMASALGORITHM_H
