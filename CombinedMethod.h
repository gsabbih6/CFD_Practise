//
// Created by Admin on 11/5/20.
//

#ifndef HW2_CFD_COMBINEDMETHOD_H
#define HW2_CFD_COMBINEDMETHOD_H

#include <iostream>
#include <vector>

using namespace std;

class CombinedMethod {
public:
    double simpleExplicit(double current, double prev,
                          double next, double CFL);

    void crankNickolson(vector<vector<double> > &matrix, int Mdim, int Ndim,
//                        double arow, double brow, double crow,
                        vector<double> rmatrix);

    double rmsd(vector<double> residualTemp);

    double exact(double initTemp, double currT, double alpha, double space, double time);

    vector<vector<double>> computeCrankNicholson(vector<vector<double>> grid, double CFL);

    vector<vector<double>> computeSimpleExplicit(vector<vector<double>> grid, double CFL);

    vector<vector<double>> computeExact(vector<vector<double>> grid,
                                        double alpha, double deltaX, double deltaT);

    void myplot(int totalLength = 1, double deltaX = 0.05, int maxT = 3, double alpha = 0.1);

    void crankNickPlot(double timestep = 0.05, double totalLength = 1.0,
                       double deltaX = 0.05, double maxT = 3, double alpha = 0.1);

    void
    simpleExplicPlot(double CFL = 0.45, double totalLength = 1, double deltaX = 0.05, double maxT = 3, double alpha = 0.1);

    vector<vector<double>> initialize(vector<vector<double>> grid, double solStepSize, double deltaX);

    double getTimeStep(double CLF, double deltaX, double alpha);

};


#endif //HW2_CFD_COMBINEDMETHOD_H
