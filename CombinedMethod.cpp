//
// Created by Admin on 11/5/20.
//

#include "CombinedMethod.h"
#include <iostream>
#include <vector>
#include <math.h>
#include "ThomasAlgorithm.h"
#include <matplot/matplot.h>

using namespace std;using namespace matplot;
vector<double> upperDiag, rho, results, rmatrix, rmsError, residualTemp;

double CombinedMethod::simpleExplicit(double current, double prev,
                                      double next, double CFL) {
    return CFL * (next - (2.0 * current) + prev) + current;
}

double CombinedMethod::rmsd(vector<double> residualTemp) {
    double sqauredSum;
    int size = residualTemp.size();

    for (int i = 0; i < size; ++i) {
        sqauredSum += (residualTemp[i] * residualTemp[i]);
    }
    return sqrt(sqauredSum) / (size - 2);
}

vector<vector<double>> CombinedMethod::computeSimpleExplicit(vector<vector<double>> grid, double CFL) {
    vector<double> prevtime = grid.at(0);
    vector<vector<double>> result;
    residualTemp.clear();
    residualTemp.resize(grid.at(0).size());

    rmsError.resize(grid.size());

    result.push_back(prevtime);
    for (int j = 1; j < grid.size(); j++) {
        vector<double> nextime = grid.at(j);

        //fill next time vector space

        for (int i = 1; i < nextime.size() - 1; i++) {

            nextime[i] = simpleExplicit(prevtime[i], prevtime[i - 1], prevtime[i + 1], CFL);
            residualTemp[i] = nextime[i] - prevtime[i];
        }

        residualTemp[0] = nextime[0] - prevtime[0];// equal to 0 always
        residualTemp[nextime.size() - 1] =
                nextime[nextime.size() - 1] - prevtime[nextime.size() - 1];// equal to 0 always

        rmsError[j] = rmsd(residualTemp);
        prevtime = nextime;
        result.push_back(nextime);

    }
    return result;
}

void CombinedMethod::crankNickolson(vector<vector<double>> &matrix, int Mdim,
                                    int Ndim,
//                                    double arow, double brow,double crow,
                                    vector<double> rmatrix
) {


    ThomasAlgorithm t;
    try {
//        t.createMatrix(matrix, Mdim, Ndim, arow, brow, crow);
//        t.print(matrix);
        t.computeThomas(matrix, rmatrix, upperDiag, rho, results, Mdim, Ndim);
//        t.print(rmatrix);
//        t.print(rho);
//        t.print(upperDiag);
        t.print(results);
    } catch (const char *e) {
        std:
        cout << e << "\n";
    }

}

vector<vector<double>> CombinedMethod::computeCrankNicholson(vector<vector<double>> grid,
                                                             double CFL) {
    int dim = grid.at(0).size() - 2;//7

    vector<vector<double>> matrix;
    ThomasAlgorithm t;
    matrix.resize(dim);
    residualTemp.clear();
    residualTemp.resize(grid.at(0).size());
    rmsError.clear();
    rmsError.resize(grid.size());
    for (int i = 0; i < dim; i++) {
        matrix.at(i).resize(dim);
    }
    double arow = -CFL / 2.0;
    double brow = CFL + 1.0;
    double crow = arow;

    t.createMatrix(matrix, dim, dim, arow, brow, crow);

    vector<double> prevtime = grid.at(0);
    vector<vector<double>> result;
    result.push_back(prevtime);
    for (int j = 1; j < grid.size(); j++) { //timesetp
        results.clear();
        results.resize(dim);
        rmatrix.clear();
        rmatrix.resize(dim);
        upperDiag.clear();
        upperDiag.resize(dim);
        rho.clear();
        rho.resize(dim);
        vector<double> nextime = grid.at(j);

        //fill next time vector space
        rmatrix[0] = ((CFL / 2.0) * prevtime[0])
                     + ((1.0 - CFL) * prevtime[0])
                     + ((CFL / 2.0) * prevtime[1]) - (arow * prevtime[0]);
        for (int i = 2; i < dim - 1; i++) { //space 1-7

//            if (i == 1) {
//
////                cout<<1;
//            } else if (i == (dim - 1)) {
//
//            } else {
            rmatrix[i - 1] = ((CFL / 2) * prevtime[i - 1])
                             + ((1.0 - CFL) * prevtime[i])
                             + ((CFL / 2.0) * prevtime[i + 1]);
//            }
//            double d = ((CFL / 2) * prevtime[i - 1])
//                       + ((1.0 - CFL) * prevtime[i])
//                       + ((CFL / 2.0) * prevtime[i + 1]);
//            if (i == 1) {
//                rmatrix[i - 1] = d - (arow * prevtime[0]);
//            } else if (i == dim) { //7-2=5
//                rmatrix[i - 1] = d - (arow * prevtime[i + 1]);
//            } else {
//                rmatrix[i - 1] = d;
//            }

        }
        rmatrix[dim - 1] = ((CFL / 2) * prevtime[dim - 2])
                           + ((1.0 - CFL) * prevtime[dim - 1])
                           + ((CFL / 2.0) * prevtime[dim]) - (arow * prevtime[dim]);

        crankNickolson(matrix, dim, dim, rmatrix);

        results.insert(results.begin(), prevtime[0]);
        results.push_back(0);

        prevtime = results;
        result.push_back(results);

        for (int k = 1; k < residualTemp.size() - 1; ++k) {
            residualTemp[k] = result.at(j)[k] - result.at(j - 1)[k];
        }
        rmsError[j] = rmsd(residualTemp);

    }
    cout << results.size();
    return result;
}

double CombinedMethod::exact(double initTemp, double currT, double alpha, double space, double time) {

    return initTemp + ((currT - initTemp) * erf(space / (2.0 * sqrt(alpha * time))));
}

vector<vector<double>> CombinedMethod::computeExact
        (vector<vector<double>> grid, double alpha, double deltaX,
         double deltaT) {
    vector<double> prevtime = grid.at(0);
    vector<vector<double>> result;
    result.push_back(prevtime);
    for (int i = 1; i < grid.size(); i++) {
        vector<double> nextime = grid.at(i);

        for (int j = 1; j < nextime.size() - 1; j++) {
            double x = j * deltaX;
            double t = i * deltaT;
            nextime[j] = exact(prevtime[0], prevtime[j], alpha, x, t);
        }

        result.push_back(nextime);
        prevtime = nextime;

    }

    return result;
}

double CombinedMethod::getTimeStep(double CFL, double deltaX, double alpha) {
    return ((deltaX * deltaX) * CFL) / alpha;
}

vector<vector<double>> CombinedMethod::initialize(vector<vector<double>> grid, double solStepSize, double deltaX) {

    vector<double> solutionSpace = grid.at(0);

    for (double j = 0; j < solStepSize; j++) {
//        double x = j * deltaX;

        if (j == 0) {
            solutionSpace[j] = 500.0;
        } else { solutionSpace[j] = 0; }


    }


    for (int j = 0; j < grid.size(); j++) {
        grid[j] = solutionSpace;
    }
//    grid.insert(grid.begin(), solutionSpace);
//    std::cout << grid.at(100)[0] << std::endl;
    return grid;
}

void CombinedMethod::simpleExplicPlot(double CFL, double totalLength, double deltaX, double maxT, double alpha) {
    vector<double> x;

    double timestep = getTimeStep(CFL, deltaX, alpha);

    int timeStepSize = (int) ((maxT / timestep) + .5) + 1;

    int solStepSize = (int) ((totalLength / deltaX) + .5) + 1;
    std::cout << solStepSize << "\n";

    int a[] = {(int) ((.1 / timestep) + .5), (int) ((.5 / timestep) + .5), int((1.0 / timestep) + .5),
               int((3.0 / timestep) + .5)};
    string as[] = {"T=0.1hr ", "T=0.5hr", "T=1hr", "T=3hrs"};

    for (int i = 0; i <= solStepSize; i++) {
        x.push_back(i * deltaX);
    }

    vector<vector<double>> grid, exactgrid;

    grid.resize(timeStepSize);
    exactgrid.resize(timeStepSize);
    for (int i = 0; i < timeStepSize; i++) {
        grid.at(i).resize(solStepSize);
        exactgrid.at(i).resize(solStepSize);

    }
    grid = initialize(grid, solStepSize, deltaX);
//    std::cout << 1<<"\n";
    exactgrid = initialize(exactgrid, solStepSize, deltaX);
    exactgrid = computeExact(exactgrid, alpha, deltaX, timestep);
    std::cout << exactgrid.size() << std::endl;
    grid = computeSimpleExplicit(grid, CFL);


    for (int j = 0; j < 4; j++) {
        plot(x, exactgrid.at(a[j]), x, grid.at(a[j]));
//        plot(x, grid.at(a[j]));
//                    legend({"Solution at " + as[j]});
                    legend({"Exact", "Solution at " + as[j]});
        title("Implicit Euler at CFL= " + to_string(CFL) + " and at Time T= " + as[j]);
        save("Implicit_Euler_at_CFL=" + to_string(CFL) + "_and_at_Time_T= " + as[j] + ".jpeg");
//                    show();
    }

    plot(rmsError);
    title("Simple Explicit RMS Error vrs time(iterations)");
    save("SimpleExplicit.png");
}

void CombinedMethod::crankNickPlot(double timestep, double totalLength,
                                   double deltaX, double maxT, double alpha) {
    vector<double> x;

//    double timestep = getTimeStep(CFL, deltaX, alpha);
    double CFL = (alpha * timestep) / (deltaX * deltaX);

//    std::cout << CFL << std::endl;

    int timeStepSize = (int) ((maxT / timestep) + .5) + 1;

    int solStepSize = (int) ((totalLength / deltaX) + .5) + 1;
    std::cout << solStepSize << "\n";

    int a[] = {(int) ((.1 / timestep) + .5), (int) ((.5 / timestep) + .5), int((1.0 / timestep) + .5),
               int((3.0 / timestep) + .5)};
    string as[] = {"T=0.1hr ", "T=0.5hr", "T=1hr", "T=3hrs"};

    for (int i = 0; i <= solStepSize; i++) {
        x.push_back(i * deltaX);
    }


    vector<vector<double>> grid, exactgrid;

    grid.resize(timeStepSize);
    exactgrid.resize(timeStepSize);

    for (int i = 0; i < timeStepSize; i++) {
        grid.at(i).resize(solStepSize);
        exactgrid.at(i).resize(solStepSize);

    }
    grid = initialize(grid, solStepSize, deltaX);

//    std::cout << 1<<"\n";
    exactgrid = initialize(exactgrid, solStepSize, deltaX);
    exactgrid = computeExact(exactgrid, alpha, deltaX, timestep);
    grid = computeCrankNicholson(grid, CFL);


    for (int j = 0; j < 4; j++) {
        plot(x, exactgrid.at(a[j]), x, grid.at(a[j]));
//        plot(x, grid.at(a[j]));
//                    legend({"Solution at " + as[j]});
                    legend({"Exact ", "Solution at " + as[j]});
        title("Crank Nicolson at timestep= " + to_string(timestep) + " and at Time T= " + as[j]);
        save("Crank_Nicolson_at_timestep=" + to_string(timestep) + "_and_at_Time_T= " + as[j] + ".jpeg");
//                    show();
    }

    plot(rmsError);
    legend("RMS Error ");
    title("Crank Nicolson RMS Error vrs time(iterations)");

    save("CrankNicolsonRMS_more" + to_string(timestep) + ".png");


}