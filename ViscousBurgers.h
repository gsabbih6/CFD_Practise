//
// Created by Admin on 11/17/20.
//

#ifndef HW2_CFD_VISCOUSBURGERS_H
#define HW2_CFD_VISCOUSBURGERS_H

#include <vector>

using namespace std;

class ViscousBurgers {
public:
    double exact(double viscosity, double space, double time);

    vector<vector<double>>
    computeImplicitViscourBurger(vector<vector<double>> grid, double deltaT, double deltaX, double mu);

    vector<vector<double>>
    computeNewtonImplicitViscourBurger(vector<vector<double>> grid, double deltaT, double deltaX, double mu,
                                       int iterations = 1);

    vector<double> computeNewton(vector<double> Ut, vector<double> Un1, double deltaT, double tmax, double deltaX, double mu,
                       int iterations = 0);

    vector<double> init(vector<double> Ut, double deltaX = 0.04, double startX = -1.0, double endX = 3.0, double t = 0,
              double mu = 0.2);

    vector<vector<double>> plot(double timestep = 0.2, double startX = -1.0, double endX = 3.0,
                                double deltaX = 0.04, double maxT = 0.4, double mu = 0.2);

    vector<vector<double>> plotNewton(double timestep = 0.2, double startX = -1.0, double endX = 3.0,
                                      double deltaX = 0.04, double maxT = 0.4, double mu = 0.2, int iterations = 1);

    vector<vector<double>> plotExact(double timestep = 0.2, double startX = -1.0, double endX = 3.0,
                                     double deltaX = 0.04, double maxT = 0.4, double mu = 0.2);

    vector<double> x(double startX = -1.0, double endX = 3.0,
                     double deltaX = 0.04);

    int it(double timestep = 0.2, double maxT = 0.4);

    int solutionSize(double deltaX = 0.04, double startX = -1.0, double endX = 3.0);

    vector<vector<double>> initialize(vector<vector<double>> grid, vector<double> x,
                                      double maxT, double deltaT, double mu);

    vector<vector<double>> computeExact(vector<vector<double>> grid, vector<double> x,
                                        double mu, double deltaT);
    vector<double> solve(double timestep = 0.2, double startX = -1.0, double endX = 3.0,
                                      double deltaX = 0.04, double maxT = 0.4, double mu = 0.2, int iterations = 1);


};


#endif //HW2_CFD_VISCOUSBURGERS_H
