//
// Created by Admin on 10/15/20.
//

#ifndef HW2_CFD_FINITEDIFFERENCE_H
#define HW2_CFD_FINITEDIFFERENCE_H

#include <iostream>
#include <vector>

using namespace std;

class FiniteDifference {
public:
    double exact(double initialSolution, double x, double c, double t);

    double exact(double x);
    double LaxWendoff(double cur, double CFL, double prev, double next);
    double MacCormack(double cur, double CFL, double prev, double next);

    double FTCS(double cur, double CFL, double prev, double next);

    double getTimeStep(double CLF, double deltaX, double c);


    double HTCS(double CFL, double prev, double next, double deltaX);

    double FTBS(double cur, double CFL, double prev);

//template<typename T>
    double lax(double next, double CFL, double prev);

    vector<vector<double>> initialize(vector<vector<double>> grid, int solStepSize, double deltaX);

    vector<vector<double>> computeFTCS(vector<vector<double>> grid, double CFL);

    vector<vector<double>> computeExact(vector<vector<double>> grid, double c, double deltaX, double deltaT);

    vector<vector<double>> computeHTCS(vector<vector<double>> grid, double CFL, double deltaX);

    vector<vector<double>> computeLax(vector<vector<double>> grid, double CFL);

    vector<vector<double>> computeFTBS(vector<vector<double>> grid, double CFL);

//    void myplot();

    vector<vector<double>> computeMacCormack(vector<vector<double>> grid, double CFL);

    vector<vector<double>> computeLaxWendof(vector<vector<double>> grid, double CFL);

    void myplot(int maxX, double deltaX, int maxT, double c);
};

#endif //HW2_CFD_FINITEDIFFERENCE_H
