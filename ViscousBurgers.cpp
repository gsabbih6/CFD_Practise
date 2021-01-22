//
// Created by Admin on 11/17/20.
//

#include "ViscousBurgers.h"

#include "math.h"
#include "ThomasAlgorithm.h"
#include "CombinedMethod.h"

vector<double> upperDiagv, rhov, resultsv, rmatrixv;
ThomasAlgorithm tm;

double ViscousBurgers::exact(double mu, double x, double t) {
    return 2.0 + (
                         (x - 2.0 * t) / (t + 0.1)
                 ) /
                 (
                         1.0 + mu * mu * sqrt(t + 0.1 *
                                                  exp(((x - 2 * t) * (x - 2 * t)) / (4 * mu * (t + 0.1))))
                 );
}

vector<vector<double>>
ViscousBurgers::computeImplicitViscourBurger(vector<vector<double>> grid, double deltaT, double deltaX, double mu) {
    int dim = grid.at(0).size() - 2;
    vector<vector<double>> matrix;
    matrix.resize(dim);
    for (int i = 0; i < dim; i++) {
        matrix.at(i).resize(dim);
    }
    double x = deltaT / deltaX;
    double r = (mu * x) / deltaX;
    vector<vector<double>> result1 = grid;
    cout << result1.size() << "\n";
    for (int n = 0; n < result1.size(); n++) {
        matrix.clear();
        matrix.resize(dim);
        for (int i = 0; i < dim; i++) {
            matrix.at(i).resize(dim);
        }
        resultsv.clear();
        resultsv.resize(dim);
        rmatrixv.clear();
        rmatrixv.resize(dim);
        upperDiagv.clear();
        upperDiagv.resize(dim);
        rhov.clear();
        rhov.resize(dim);

        for (int i = 1; i <= dim; i++) {

            double arow = (-0.5 * x * result1[n][i - 1]) - r;
            double brow = 1.0 + 2.0 * r;
            double crow = (0.5 * x * result1[n][i + 1]) - r;
            double d =
                    (
                            0.5 * x *
                            (
                                    (0.5 * result1[n][i + 1] * result1[n][i + 1]) -
                                    (0.5 * result1[n][i - 1] * result1[n][i - 1])
                            ) -
                            (
                                    r *
                                    (result1[n][i + 1] - (2 * result1[n][i]) + result1[n][i - 1])
                            )
                    );
            if (i - 1 == 0 && n < result1.size() - 1) {
                rmatrixv[i - 1] = -d -
                                  (

                                          arow * (result1[n + 1][i] - result1[n][i])
                                  );
            } else if (i == dim && n < result1.size() - 1) {
                rmatrixv[i - 1] = -d -
                                  (
                                          arow * (result1[n + 1][i] - result1[n][i])
                                  );
            } else {
                rmatrixv[i - 1] = isnan(-d) ? 0 : -d;;
            }
            for (int k = 1; k <= dim; ++k) {
                if (i - k == 0) {
                    matrix[i - 1][k - 1] = isnan(brow) ? 0 : brow;
                }
                if (i - k == -1) {
                    matrix[i - 1][k - 1] = isnan(crow) ? 0 : crow;;
                }
                if (i - k == 1) {
                    matrix[i - 1][k - 1] = isnan(arow) ? 0 : arow;;
                }

            }
        }

        tm.computeThomas(matrix, rmatrixv, upperDiagv, rhov, resultsv, dim, dim);

        resultsv.insert(resultsv.begin(), 0.0);
        resultsv.push_back(0.0);
        // update velocity
        for (int i = 0; i < resultsv.size(); ++i) {
            resultsv[i] = resultsv[i] + result1[n][i];
        }
        if (n < result1.size() - 1) {
            result1[n + 1] = (resultsv);
        }
        tm.print(rmatrixv);

    }

    return result1;
}

vector<vector<double>>
ViscousBurgers::initialize(vector<vector<double>> grid,
                           vector<double> x, double maxT, double deltaT, double mu) {

    vector<double> solutionSpace;
    solutionSpace.resize(grid.at(0).size());
    for (int i = 0; i < grid.at(0).size(); i++) {
        solutionSpace[i] = exact(mu, x[i], 0);
    }

    grid[0] = solutionSpace;

    for (int n = 1; n < grid.size(); n++) {
        solutionSpace[0] = exact(mu, -1.0, n * deltaT);

        solutionSpace[solutionSpace.size() - 1] = exact(mu, 3.0, n * deltaT);

        grid[n] = solutionSpace;


    }
    return grid;
}

vector<vector<double>>
ViscousBurgers::computeExact(vector<vector<double>> grid, vector<double> xs, double mu, double deltaT) {

    for (double n = 0; n < grid.size(); n++) {
        vector<double> x = grid.at(n);
        for (int i = 0; i < x.size(); ++i) {

            x[i] = exact(mu, xs[i], n * deltaT);
        }

        grid[n] = x;
    }
    return grid;
}

int ViscousBurgers::it(double timestep, double maxT) {
    return (int) ((maxT / timestep) + .5) + 1;
}

int ViscousBurgers::solutionSize(double deltaX, double startX, double endX) {

    double maxX = endX - startX;
    return (int) ((maxX / deltaX) + .5) + 1;
}

vector<double> ViscousBurgers::x(double startX, double endX, double deltaX) {

    vector<double> x;
    x.resize(solutionSize(deltaX,startX, endX));
    x[0] = startX;
    for (int i = 1; i < solutionSize(deltaX,startX, endX); i++) {
        x[i] = x[i - 1] + deltaX;
    }
    return x;
}

vector<vector<double>>
ViscousBurgers::plotExact(double deltaT, double startX, double endX, double deltaX, double maxT, double mu) {

    vector<double> mx = x(startX, endX, deltaX);
    vector<vector<double>> grid;//, exactgrid;
    grid.resize(it(deltaT, maxT));
    for (int i = 0; i < it(deltaT, maxT); i++) {
        grid.at(i).resize(solutionSize(deltaX,startX, endX));
    }
    grid = initialize(grid, mx, maxT, deltaX, mu);
    grid = computeExact(grid, mx, mu, deltaT);
    return grid;
}

vector<vector<double>>
ViscousBurgers::plot(double deltaT, double startX, double endX, double deltaX, double maxT, double mu) {

    vector<double> mx=x(startX,endX,deltaX);
    vector<vector<double>> grid;//, exactgrid;
    grid.resize(it(deltaT, maxT));
    cout<<solutionSize(deltaX,startX, endX);
    for (int i = 0; i < it(deltaT, maxT); i++) {
        grid.at(i).resize(solutionSize(deltaX,startX, endX));
    }
    grid = initialize(grid, mx, maxT, deltaX, mu);
    grid = computeImplicitViscourBurger(grid, deltaT, deltaX, mu);


    return grid;
}