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
    return
            2.0 +
            ((x - 2.0 * t) / (t + 0.1)) /
            (1.0 + mu * mu * sqrt(t + 0.1) *
                   exp(((x - 2.0 * t) * (x - 2.0 * t)) / (4.0 * mu * (t + 0.1))));
}

vector<double> ViscousBurgers::solve(double deltaT, double startX, double endX, double deltaX, double maxT, double mu,
                                     int iterations) {
    vector<double> Ut;
    Ut = init(Ut, deltaX, startX, endX, 0, mu);
//    std::cout << Ut.size();
    Ut = computeNewton(Ut, Ut, deltaT, maxT, deltaX, mu, iterations);
    return Ut;
}

vector<double> ViscousBurgers::init(vector<double> Ut, double deltaX, double startX, double endX, double t, double mu) {
    int idim = solutionSize(deltaX, startX, endX);


    vector<double> xi = x(startX, endX, deltaX);
    Ut.resize(idim);
    for (int i = 0; i < idim; i++) {
        Ut[i] = exact(mu, xi[i], 0);
    }
    return Ut;
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
        tm.print(matrix);

    }

    return result1;
}

vector<double>
ViscousBurgers::computeNewton(vector<double> Un, vector<double> Un1, double deltaT, double tmax, double deltaX,
                              double mu, int iterations) {
    int idim = Un.size(); // size(i);
    int ndim = it(deltaT, tmax);// size of timesteps

    vector<vector<double>> thomasMatrix; // thomas matrix
    thomasMatrix.resize(idim - 2);
    for (int i = 0; i < idim - 2; i++) {
        thomasMatrix.at(i).resize(idim - 2);
    }
    resultsv.clear();
    resultsv.resize(idim - 2);
    rmatrixv.clear();
    rmatrixv.resize(idim - 2);
    upperDiagv.clear();
    upperDiagv.resize(idim - 2);
    rhov.clear();
    rhov.resize(idim - 2);

    // constants
    double dxdt = deltaT / deltaX;
    double r = (mu * dxdt) / deltaX;

    for (int n = 0; n < ndim; n++) { //time domain
        for (int m = 0; m < iterations; m++) { // Newton's iterations m
            if (m == 0) { Un1 = Un; }
            for (int i = 1; i < idim - 1; i++) { //space domain
                double a = (-0.5 * dxdt * Un1[i - 1]) - r;
                double b = 1.0 + 2.0 * r;
                double c = (0.5 * dxdt * Un1[i + 1]) - r;
                double d =
                        Un1[i] - Un[i] +
                        (0.25 * dxdt * ((Un1[i + 1] * Un1[i + 1]) - (Un1[i - 1] * Un1[i - 1]))) -
                        (r * (Un1[i + 1] - (2 * Un1[i]) + Un1[i - 1]));
                if (i == 1) {
                    rmatrixv[i - 1] = -d - (a * (Un1[i - 1] - Un[i - 1]));
                } else if (i == idim - 1) {
                    rmatrixv[i - 1] = -d - (a * (Un1[i + 1] - Un[i + 1]));
                } else {
                    rmatrixv[i - 1] = isnan(-d) ? 0 : -d;
                }
                for (int k = 1; k <= idim - 2; ++k) {
                    if (i - k == 0) {
                        thomasMatrix[i - 1][k - 1] = isnan(b) ? 0 : b;
                    }
                    if (i - k == -1) {
                        thomasMatrix[i - 1][k - 1] = isnan(c) ? 0 : c;;
                    }
                    if (i - k == 1) {
                        thomasMatrix[i - 1][k - 1] = isnan(a) ? 0 : a;;
                    }

                }
            }
            tm.computeThomas(thomasMatrix, rmatrixv,
                             upperDiagv, rhov, resultsv, idim - 2, idim - 2);

            for (int i = 0; i < resultsv.size(); i++) {
                Un1[i + 1] = resultsv[i] + Un1[i + 1];
            }
            tm.print(Un1);

        }
        //update boundary conditions
        Un1[0] = exact(mu, -1.0, (n + 1) * deltaT);
        Un1[Un1.size()] = exact(mu, 3.0, (n + 1) * deltaT);
        Un = Un1;
    }

    return Un;
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
    x.resize(solutionSize(deltaX, startX, endX));
    x[0] = startX;
    for (int i = 1; i < solutionSize(deltaX, startX, endX); i++) {
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
        grid.at(i).resize(solutionSize(deltaX, startX, endX));
    }
    grid = initialize(grid, mx, maxT, deltaX, mu);
    grid = computeExact(grid, mx, mu, deltaT);
    return grid;
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
//vector<vector<double>>
//ViscousBurgers::plot(double deltaT, double startX, double endX, double deltaX, double maxT, double mu) {
//
//    vector<double> mx = x(startX, endX, deltaX);
//    vector<vector<double>> grid;//, exactgrid;
//    grid.resize(it(deltaT, maxT));
//    cout << solutionSize(deltaX, startX, endX);
//    for (int i = 0; i < it(deltaT, maxT); i++) {
//        grid.at(i).resize(solutionSize(deltaX, startX, endX));
//    }
//    grid = initialize(grid, mx, maxT, deltaX, mu);
////    grid = computeImplicitViscourBurger(grid, deltaT, deltaX, mu);
//    grid = computeImplicitViscourBurger(grid, deltaT, deltaX, mu);
//
//
//    return grid;
//}

//vector<vector<double>> ViscousBurgers::plotNewton(double deltaT, double startX, double endX, double deltaX,
//                                                  double maxT, double mu, int iterations) {
//
//    vector<double> mx = x(startX, endX, deltaX);
//    vector<vector<double>> grid;//, exactgrid;
//    grid.resize(it(deltaT, maxT));
//    cout << solutionSize(deltaX, startX, endX);
//    for (int i = 0; i < it(deltaT, maxT); i++) {
//        grid.at(i).resize(solutionSize(deltaX, startX, endX));
//    }
//    grid = initialize(grid, mx, maxT, deltaX, mu);
////    grid = computeImplicitViscourBurger(grid, deltaT, deltaX, mu);
//    grid = computeNewtonImplicitViscourBurger(grid, deltaT, deltaX, mu, iterations);
//
//
//    return grid;
//}
//vector<vector<double>>
//ViscousBurgers::computeNewtonImplicitViscourBurger(vector<vector<double>> grid, double deltaT, double deltaX,
//                                                   double mu,
//                                                   int iterations) {
//    int dim = grid.at(0).size() - 2;
//    vector<vector<double>> matrix;
//    matrix.resize(dim);
//    for (int i = 0; i < dim; i++) {
//        matrix.at(i).resize(dim);
//    }
//
//    // constants
//    double x = deltaT / deltaX;
//    double r = (mu * x) / deltaX;
//    vector<vector<double>> result1 = grid; // time grid
//    cout << result1.size() << "\n";
//    for (int n = 0; n < result1.size() - 1; n++) {
//        vector<double> mans = result1.at(n);
//        for (int m = 0; m < iterations; m++) {
//            matrix.clear();
//            matrix.resize(dim);
//            for (int i = 0; i < dim; i++) {
//                matrix.at(i).resize(dim);
//            }
//            resultsv.clear();
//            resultsv.resize(dim);
//            rmatrixv.clear();
//            rmatrixv.resize(dim);
//            upperDiagv.clear();
//            upperDiagv.resize(dim);
//            rhov.clear();
//            rhov.resize(dim);
//
//            /* compute the tri-diagonal matrix system to use in the Thomas Alg*/
//            for (int i = 1; i <= dim; i++) {
//
//                double arow = (-0.5 * x * mans[i - 1]) - r;
//                double brow = 1.0 + 2.0 * r;
//                double crow = (0.5 * x * mans[i + 1]) - r;
//                double d =
//                        (0.5 * x * ((0.5 * mans[i + 1] * mans[i + 1]) -
//                                    (0.5 * mans[i - 1] * mans[i - 1])) - (r *
//                                                                          (mans[i + 1] - (2 * mans[i]) +
//                                                                           mans[i - 1])));
//                if (i - 1 == 0 && n < result1.size() - 1) {
//                    rmatrixv[i - 1] = -d - (arow * (result1[n + 1][i] - mans[i]));
//                } else if (i == dim && n < result1.size() - 1) {
//                    rmatrixv[i - 1] = -d - (arow * (result1[n + 1][i] - mans[i]));
//                } else {
//                    rmatrixv[i - 1] = isnan(-d) ? 0 : -d;
//                }
//                for (int k = 1; k <= dim; ++k) {
//                    if (i - k == 0) {
//                        matrix[i - 1][k - 1] = isnan(brow) ? 0 : brow;
//                    }
//                    if (i - k == -1) {
//                        matrix[i - 1][k - 1] = isnan(crow) ? 0 : crow;;
//                    }
//                    if (i - k == 1) {
//                        matrix[i - 1][k - 1] = isnan(arow) ? 0 : arow;;
//                    }
//
//                }
//            }
//
//            tm.computeThomas(matrix, rmatrixv, upperDiagv, rhov, resultsv, dim, dim);
//
////            resultsv.insert(resultsv.begin(), 0.0);
////            resultsv.push_back(0.0);
//            // update velocity
//            for (int i = 0; i < dim; i++) {
//                mans[i + 1] = resultsv[i] + mans[i + 1];
//            }
//
////            tm.print(rmatrixv);
//        }
//        for (int i = 1; i < result1.size() - 1; i++) {
//            result1[n + 1][i] = mans[i];
//        }
//        //update boundary conditions
//        result1[n + 1][0] = exact(mu, -1, (n + 1) * deltaT);
//        result1[n + 1][result1[n].size() - 1] = exact(mu, -3, (n + 1) * deltaT);
//    }
//
//    return result1;
//}
