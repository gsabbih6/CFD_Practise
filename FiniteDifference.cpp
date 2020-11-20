//
// Created by Admin on 10/15/20.
//

#include "FiniteDifference.h"
#include <iostream>
#include <cmath>
#include <matplot/matplot.h>
#include <vector>
#include <string>
#include <unistd.h>

using namespace std;
using namespace matplot;

double FiniteDifference::exact(double x) {
    if (x >= 0 && x <= 3) {
        return x;
    } else if (x >= 3 && x <= 6) {
        return 6.0 - x;
    }
    return 0.0;// c * initialSolution;//*(x-c*t);
}

double FiniteDifference::FTCS(double cur, double CFL, double prev, double next) {
    // method returns next time in current space
    return cur - ((CFL / 2.0) * (next - prev));
}

double FiniteDifference::LaxWendoff(double cur, double CFL, double prev, double next) {
    return cur - (0.5 * CFL * (next - prev)) + (0.5 * CFL * CFL * (next - (2 * cur) + prev));
}

double FiniteDifference::MacCormack(double cur, double CFL, double prev, double next) {
    double curPredictor = cur - CFL * (cur - prev);
    double nextPredictor = next - CFL * (next - cur);

    return (0.5 * (cur + curPredictor)) - (0.5 * CFL * (nextPredictor - curPredictor));
}


double FiniteDifference::getTimeStep(double CLF, double deltaX, double c) {
    return (CLF * deltaX) / c;
}

double FiniteDifference::HTCS(double CFL, double prev, double next, double deltaX) {
    // method returns next time in current space
    return ((1.0 / (2.0 * deltaX)) - (CFL / 2.0)) * (next - prev);
}

double FiniteDifference::FTBS(double cur, double CFL, double prev) {
    // method returns next time in current space
    return cur - (CFL * (cur - prev));
}

//template<typename T>
double FiniteDifference::lax(double next, double CFL, double prev) {
    // method returns next time in current space
    return (0.5 * (next + prev)) - (0.5 * CFL * (next - prev));
}

//template<typename T>
vector<vector<double>> FiniteDifference::initialize(vector<vector<double>> grid, int solStepSize, double deltaX) {
// resize the solution vector
    // init at t=0;

    vector<double> solutionSpace = grid.at(0);
    for (double j = 0; j < solStepSize; j++) {
        double x = j * deltaX;
        if (x >= 0 && x <= 3) {
            solutionSpace[j] = x;
        } else if (x >= 3 && x <= 6) {
            solutionSpace[j] = 6.0 - x;
        } else {
            solutionSpace[j] = 0.0;
        }


    }
    cout << solutionSpace.size() << endl;
//        cout << solStepSize << endl;

    grid.insert(grid.begin(), solutionSpace);
    return grid;
}

vector<vector<double>> FiniteDifference::computeFTCS(vector<vector<double>> grid, double CFL) {
    vector<double> prevtime = grid.at(0);
    vector<vector<double>> result;
    result.push_back(prevtime);
    for (int i = 1; i < grid.size(); i++) {
        vector<double> nextime = grid.at(i);

        //fill next time vector space

        for (int j = 1; j < nextime.size() - 1; j++) {

            nextime[j] = FTCS(prevtime[j], CFL, prevtime[j - 1], prevtime[j + 1]);
        }
        prevtime = nextime;
        result.push_back(nextime);

    }
    return result;
}

vector<vector<double>>
FiniteDifference::computeExact(vector<vector<double>> grid, double c, double deltaX, double deltaT) {
    vector<double> prevtime = grid.at(0);
    vector<vector<double>> result;
    result.push_back(prevtime);
    for (int i = 1; i < grid.size(); i++) {
        vector<double> nextime = grid.at(i);

        //fill next time vector space

        for (int j = 1; j < nextime.size() - 1; j++) {
            double x = j * deltaX;
            double ct = c * i * deltaT;
            nextime[j] = exact(x - ct);
        }

        result.push_back(nextime);
        prevtime = nextime;

    }

    return
            result;
}

vector<vector<double>> FiniteDifference::computeHTCS(vector<vector<double>> grid, double CFL, double deltaX) {
    vector<double> prevtime = grid.at(0);
    vector<vector<double>> result;
    result.push_back(prevtime);
    for (int i = 1; i < grid.size(); i++) {
        vector<double> nextime = grid.at(i);

        //fill next time vector space

        for (int j = 1; j < nextime.size() - 1; j++) {

            nextime[j] = HTCS(CFL, prevtime[j - 1], prevtime[j + 1], deltaX);
        }
        prevtime = nextime;
        result.push_back(nextime);

    }
    return result;
}

vector<vector<double>> FiniteDifference::computeLaxWendof(vector<vector<double>> grid, double CFL) {
    vector<double> prevtime = grid.at(0);
    vector<vector<double>> result;
    result.push_back(prevtime);
    for (int i = 1; i < grid.size(); i++) {
        vector<double> nextime = grid.at(i);

        //fill next time vector space

        for (int j = 1; j < nextime.size() - 1; j++) {

            nextime[j] = LaxWendoff(prevtime[j], CFL, prevtime[j - 1], prevtime[j + 1]);
        }
        prevtime = nextime;
        result.push_back(nextime);

    }
    return result;
}

vector<vector<double>> FiniteDifference::computeMacCormack(vector<vector<double>> grid, double CFL) {
    vector<double> prevtime = grid.at(0);
    vector<vector<double>> result;
    result.push_back(prevtime);
    for (int i = 1; i < grid.size(); i++) {
        vector<double> nextime = grid.at(i);

        //fill next time vector space

        for (int j = 1; j < nextime.size() - 1; j++) {

            nextime[j] = MacCormack(prevtime[j], CFL, prevtime[j - 1], prevtime[j + 1]);
        }
        prevtime = nextime;
        result.push_back(nextime);

    }
    return result;
}

vector<vector<double>> FiniteDifference::computeLax(vector<vector<double>> grid, double CFL) {
    vector<double> prevtime = grid.at(0);
    vector<vector<double>> result;
    result.push_back(prevtime);
    for (int i = 1; i < grid.size(); i++) {
        vector<double> nextime = grid.at(i);

        //fill next time vector space

        for (int j = 1; j < nextime.size() - 1; j++) {

            nextime[j] = lax(prevtime[j + 1], CFL, prevtime[j - 1]);
        }
        prevtime = nextime;
        result.push_back(nextime);

    }
    return result;
}

vector<vector<double>> FiniteDifference::computeFTBS(vector<vector<double>> grid, double CFL) {
    vector<double> prevtime = grid.at(0);
    vector<vector<double>> result;
    result.push_back(prevtime);
    for (int i = 1; i < grid.size(); i++) {
        vector<double> nextime = grid.at(i);

        //fill next time vector space

        for (int j = 1; j < nextime.size() - 1; j++) {

            nextime[j] = FTBS(prevtime[j], CFL, prevtime[j - 1]);
        }
        prevtime = nextime;
        result.push_back(nextime);

    }
    return result;
}

void FiniteDifference::myplot(int maxX = 25, double deltaX = 0.25, int maxT = 10, double c = 1.0) {
//    vector<vector<double>> grid;
    vector<double> x;
//    int maxX = 25;
//    double deltaX = 0.25;
//    int maxT = 10;
////    double CFL = 1.0;
//    double c = 1.0;

//    int a = 2 / deltaT;
//    int b = 5 / deltaT;
//    int ca = 10 / deltaT;
//    plot(x, grid.at(a), x, grid.at(b), x, grid.at(ca));
//    legend({"T = 2","T = 5","T = 10"});
//    plot(x, grid.at(a));
//    legend("T = 2");
//    plot(x, grid.at(b));
//    legend({"T = 5"});
//    plot(x, grid.at(ca));
//    legend({"T = 10"});
//    title(" Lax Method at CFL=1.0 and at T = 2");

//    show();

    double CFLa[] = {0.25, 0.5, 1.0, 1.25};
    for (double cfl:CFLa) {
        double deltaT = getTimeStep(cfl, deltaX, c);

        int timeStepSize = (int) ((maxT / deltaT) + .5);
        int solStepSize = (int) ((maxX / deltaX) + .5);


//        grid = initialize(grid, solStepSize, deltaX);
//    grid = computeFTCS(grid, CFL);
//    grid = computeFTBS(grid, CFL);
//        grid = computeLax(grid, cfl);

        int a[] = {(int) ((25.0 / deltaT) + .5), (int) ((50.0 / deltaT) + .5), int((75.0 / deltaT) + .5),
                   int((100.0 / deltaT) + .5)};
        string as[] = {"T=25 ", "T=50", "T=75", "T=100"};
        string met[] = {"FTBS", "Lax-Wendoff", "MacCormack"};
        for (int i = 0; i <= solStepSize; i++) {
            x.push_back(i * deltaX);
        }


        for (int i = 0; i < 4; i++) {
            vector<vector<double>> grid, exactgrid;
            grid.resize(timeStepSize);
//            exactgrid.resize(timeStepSize);
            for (int i = 0; i < timeStepSize; i++) {
                grid.at(i).resize(solStepSize);
//                exactgrid.at(i).resize(solStepSize);

            }
            grid = initialize(grid, solStepSize, deltaX);
//            exactgrid = initialize(exactgrid, solStepSize, deltaX);
//            exactgrid = computeExact(exactgrid, c, deltaX, deltaT);
            if (i == 0) { //LAX
                grid = computeFTBS(grid, cfl);
                for (int j = 0; j < 4; j++) {
//                    plot(x, exactgrid.at(a[j]), x, grid.at(a[j]));
                    plot(x, grid.at(a[j]));
//                    legend({"Solution at " + as[j]});
//                    legend({"Initial value", "Solution at " + as[j]});
                    title(met[i] + " at CFL= " + to_string(cfl) + " and at " + as[j]);
                    save(met[i] + "_at_CFL=_" + to_string(cfl) + "_and_at " + as[j] + ".jpeg");
//                    show();
                }
            }
            if (i == 1) { //FTCS
                grid = computeLaxWendof(grid, cfl);
                for (int j = 0; j < 4; j++) {
                    //                    plot(x, exactgrid.at(a[j]), x, grid.at(a[j]));
                    plot(x, grid.at(a[j]));
//                    legend({"Solution at " + as[j]});
//                    legend({"Initial value", "Solution at " + as[j]});
                    title(met[i] + " t CFL= " + to_string(cfl) + " and at " + as[j]);
                    save(met[i] + "_at_CFL=_" + to_string(cfl) + "_and_at " + as[j] + ".jpeg");
//                    show();
                }
            }
            if (i == 2) { //FTBS
                grid = computeMacCormack(grid, cfl);
                for (int j = 0; j < 4; j++) {
                    //                    plot(x, exactgrid.at(a[j]), x, grid.at(a[j]));
                    plot(x, grid.at(a[j]));
//                    legend({"Solution at " + as[j]});
//                    legend({"Initial value", "Solution at " + as[j]});
                    title(met[i] + " at CFL= " + to_string(cfl) + " and at " + as[j]);
                    save(met[i] + "_at_CFL=_" + to_string(cfl) + "_and_at " + as[j] + ".jpeg");
//                    show();
                }
            }


        }
    }

//    double CFLa[] = {0.05,0.5, 1.0, 1.25};

    for (double cfl:CFLa) {
        double deltaT = getTimeStep(cfl, deltaX, c);

        int timeStepSize = (int) ((maxT / deltaT) + .5);
        int solStepSize = (int) ((maxX / deltaX) + .5);

        int a[] = {(int) ((25.0 / deltaT) + .5), (int) ((50.0 / deltaT) + .5), int((75.0 / deltaT) + .5),
                   int((100.0 / deltaT) + .5)};
        for (int i = 0; i <= solStepSize; i++) {
            x.push_back(i * deltaX);
        }
//        vector<vector<double>> grid;
        vector<vector<double>> grid2;
        vector<vector<double>> grid3, grid4, grid5;
//        grid.resize(timeStepSize);
        grid2.resize(timeStepSize);
        grid3.resize(timeStepSize);
        grid4.resize(timeStepSize);
        grid5.resize(timeStepSize);
        for (int i = 0; i < timeStepSize; i++) {
//            grid.at(i).resize(solStepSize);
            grid2.at(i).resize(solStepSize);
            grid3.at(i).resize(solStepSize);
            grid4.at(i).resize(solStepSize);
            grid5.at(i).resize(solStepSize);


        }


//        grid = initialize(grid, solStepSize, deltaX);
        grid2 = initialize(grid2, solStepSize, deltaX);
        grid3 = initialize(grid3, solStepSize, deltaX);
        grid4 = initialize(grid4, solStepSize, deltaX);
        grid5 = initialize(grid5, solStepSize, deltaX);

//        grid = computeLax(grid, cfl);
        grid2 = computeExact(grid2, c,deltaX,deltaT);
        grid3 = computeFTBS(grid3, cfl);
        grid4 = computeLaxWendof(grid4, cfl);
        grid5 = computeMacCormack(grid5, cfl);
//        string as[] = {"T=2 ", "T=5", "T=10"};
        string as[] = {"T=25 ", "T=50", "T=75" ,"T=100"};
        string met[] = {"FTBS", "Lax-Wendoff", "MacCormack"};
        int i = 0;
        for (int aa:a) {
//            plot(x, grid5.at(aa));
//            plot(x, grid5.at(aa), x, grid.at(aa), x, grid2.at(aa), x, grid3.at(aa));
            plot(x, grid2.at(aa),x, grid3.at(aa), x, grid4.at(aa), x, grid5.at(aa));
            legend({"Exact","FTBS", "Lax-Wendoff", "MacCormack"});
//            legend({"Initial value", "Lax", "FTBS"});
            title("At CFL= " + to_string(cfl) + " and at " + as[i]);
            save("At CFL= " + to_string(cfl) + " and at " + as[i] + ".jpeg");
            i++;
        }


    }
}

