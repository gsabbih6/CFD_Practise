//
// Created by Admin on 11/17/20.
//

#include "ViscousBurgers.h"

#include "math.h"
#include "ThomasAlgorithm.h"
#include "CombinedMethod.h"

#include <matplot/matplot.h>

using namespace matplot;
vector<double> upperDiagv, rhov, resultsv, rmatrixv;
ThomasAlgorithm tm;

double ViscousBurgers::exact(double mu, double x, double t) {
    return 2.0 + ((x - 2.0 * t) * (t + 0.1)) /
                 (1.0 + mu * mu * sqrt(t + 0.1 *
                                           exp(((x - 2 * t) * (x - 2 * t)) / (4 * mu * (t + 0.1)))));
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
    vector<vector<double>> result = grid;
    vector<vector<double>> result1;
    result1.resize(grid.size()+1);
    result1[0] = grid.at(0);
    cout << result.size() << "\n";
    for (int j = 0; j < grid.size(); j++) {
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
        vector<double> nextime = grid.at(j);

        for (int i = 0; i < dim; i++) { //space 1-7

            double arow = (-0.5 * x * result[j ][i + 1 - 1]) - r;
            double brow = 1.0 + 2.0 * r;
            double crow = (0.5 * x * result[j][i + 1 + 1]) - r;
            double d = (
                    0.5 * x * (
                            (0.5 * result[j][i + 1 + 1] * result[j][i + 1 + 1]) -
                            (0.5 * result[j][i + 1 - 1]
                             * result[j][i + 1 - 1]))
                    - (r * (result[j][i + 1 + 1] - (2 * result[j][i + 1] + result[j ][i + 1 - 1]))));
//            if (i == 0 || i == dim - 1) {
//                rmatrixv[i] = -d - (arow * (result[j][i] - result[j ][i]));
//            } else {
                rmatrixv[i] = -d;
//            }
            for (int k = 0; k < dim; ++k) {
                if (i == k) {
                    matrix[i][k] = brow;

                }
                if (i - k == 1) {
                    matrix[i][k] = arow;;
                }
                if (k - i == 1) {
                    matrix[i][k] = crow;
                }

            }

        }

        tm.computeThomas(matrix, rmatrixv, upperDiagv, rhov, resultsv, dim, dim);

        resultsv.insert(resultsv.begin(), 0.0);
        resultsv.push_back(0.0);
//        tm.print(matrix);

        for (int i = 0; i < resultsv.size(); ++i) {
            resultsv[i] = resultsv[i] + result[j - 1][i];
//            cout<<resultsv.size()<<" ";
        }
//        grid.at(j) = resultsv;
//        auto itPos = result.begin();
//        result.insert(itPos, resultsv);
        result1[j+1] = (resultsv);
//        cout << result1.size() << "\n";
        tm.print(resultsv);
//
//        result[j + 1] = resultsv;
    }
//    tm.print(resultsv);

    return result1;
}

vector<vector<double>>
ViscousBurgers::initialize(vector<vector<double>> grid,
                           vector<double> x, double maxT, double deltaT, double mu) {
    vector<vector<double>> result;

    for (int j = 0; j < grid.size(); j++) {
        vector<double> solutionSpace = grid.at(j);
        solutionSpace[0] = 2.0 - (((1.0 + 2 * j * deltaT) / (j * deltaT + 0.1)) /
                                  (1.0 + (mu * mu) * sqrt(j * deltaT + 0.1) *

                                         exp((-1.0 - 2.0 * j * deltaT) *
                                             (-1.0 - 2.0 * j * deltaT) /
                                             (4.0 * mu * (j * deltaT + 0.1)))));
        solutionSpace[solutionSpace.size() - 1] = 2.0 + (((3.0 + 2 * j * deltaT) / (j * deltaT + 0.1)) /
                                                         (1.0 + (mu * mu) * sqrt(j * deltaT + 0.1) *
                                                                exp((3.0 -
                                                                     2.0 * j * deltaT) *
                                                                    (3.0 -
                                                                     2.0 * j * deltaT) /
                                                                    (4.0 * mu *
                                                                     (j * deltaT +
                                                                      0.1)))));

        for (int i = 1; i < solutionSpace.size() - 1; ++i) {
            if (j == 0) {
                solutionSpace[i] = 2.0 + ((10 * x[i]) /
                                          (1.0 + (mu * mu) * sqrt(0.1) *
                                                 exp((x[i] * x[i]) / (0.4 * mu))));
            }
        }


        grid[j] = solutionSpace;


    }
    return grid;
}

vector<vector<double>>
ViscousBurgers::computeExact(vector<vector<double>> grid, vector<double> xs, double mu, double deltaT) {

    for (double j = 0; j < grid.size(); j++) {
        vector<double> x = grid.at(j);
        for (int i = 0; i < x.size() - 1; ++i) {

            x[i] = 2.0 + (((xs[i] + 2 * j * deltaT) / (j * deltaT + 0.1)) /
                          (1.0 + (mu * mu) * sqrt(j * deltaT + 0.1) *
                                 exp((xs[i] - 2.0 * j * deltaT) *
                                     (xs[i] - 2.0 * j * deltaT) /
                                     (4.0 * mu * (j * deltaT + 0.1)))));
        }

        grid[j] = x;
    }
    return grid;
}

void ViscousBurgers::plot(double deltaT, double maxX, double deltaX, double maxT, double mu) {
    int timeStepSize = (int) ((maxT / deltaT) + .5) + 1;
    vector<double> x;
    int solStepSize = (int) ((maxX / deltaX) + .5) + 1;
//    std::cout << timeStepSize << "\n";

    int a[] = {(int) ((.1 / deltaT) + .5), (int) ((.5 / deltaT) + .5), int((1.0 / deltaT) + .5),
               int((3.0 / deltaT) + .5)};
    string as[] = {"T=0.1hr ", "T=0.5hr", "T=1hr", "T=3hrs"};
    x.resize(solStepSize);
//    double xi = -1.0;
    x[0] = -1.0;
    for (int i = 1; i < solStepSize; i++) {
        x[i] = x[i - 1] + deltaX;
    }


    vector<vector<double>> grid, exactgrid;

    grid.resize(timeStepSize);
    exactgrid.resize(timeStepSize);

    for (int i = 0; i < timeStepSize; i++) {
        grid.at(i).resize(solStepSize);
        exactgrid.at(i).resize(solStepSize);

    }
    grid = initialize(grid, x, maxT, deltaX, mu);

//    std::cout << 1<<"\n";
    exactgrid = initialize(exactgrid, x, maxT, deltaX, mu);
    exactgrid = computeExact(exactgrid, x, mu, deltaT);

    tm.print(exactgrid);
    grid = computeImplicitViscourBurger(grid, deltaT, deltaX, mu);
//
    matplot::plot( x, exactgrid.at(timeStepSize - 1),x, grid.at(timeStepSize - 1));
    show();
}