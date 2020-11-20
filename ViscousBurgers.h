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

    vector<vector<double>> computeImplicitViscourBurger(vector<vector<double>> grid, double deltaT,double deltaX,double mu);
    void plot(double timestep=0.2, double maxX=4.0,
                       double deltaX = 0.04, double maxT = 0.4, double mu = 0.2);

    vector<vector<double>> initialize(vector<vector<double>> grid,vector<double> x,
                                      double maxT,double deltaT, double mu);
    vector<vector<double>> computeExact(vector<vector<double>> grid,vector<double> x,
                                         double mu,double deltaT);

};


#endif //HW2_CFD_VISCOUSBURGERS_H
