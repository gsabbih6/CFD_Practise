# Implementation of Computational Fluid Dynamics Algorithms

This repository is a compilation of implementations of several algorithms used in **CFD**. Codes are a result of
assigmenst given and course taught by **Dr. Sreenivas at UTC**. If you find them interesting, helpful, any mistakes or
just want to say Hi!, send me an email
[gsabbih5@gmail.com](mailto:gsabbih5@gmail.com)
---

## Finite Difference methods for solving PDEs

The implementation is found the FiniteDifference.cpp file which contains an application of the FTCS, FTBS and Lax
methods. This include the following methods Forward Time Central Space Forward Time Backward Space", Lax-Wendoff &
MacCormack Simply calling the method  ```myplot(int maxX, double deltaX, int maxT, double c);``` returns the sample plot
below

## Thomas Algorthm

The thomas algorithm is used to solve tridiagonal systems of equations. The implementation here creates a sparse matrix
which is not may not be memory efficient. I will update it soon, but that should be easy

### Usage

1. First download the header file for the Thomas Algorith in the directory
2. Create another a class file and use the code below
3. In your program,
   call ```computeThomas(vector<vector<double> > &matrix, vector<double> &rmatrix, vector<double> &upperDiagonal, vector<double> &rho, vector<double> &results, int Mdim, int Ndim)```
   and pass in your information. Example

### Implementation

``` c++
//
// Created by Godfred Sabbih on 10/15/20.
//

#include "ThomasAlgorithm.h"
#include <iostream>
#include <vector>
//#include <Kokkos_Core.hpp>

using namespace std;
double ThomasAlgorithm::getRandom(int seed) {
    srand(10);
    return (rand() % 10) + 1;
}

void ThomasAlgorithm::createRandomTMatrix(vector<vector<double> > &matrix,
                                          int Mdim, int Ndim, bool stdmatrix) {
    // Mdim =Ndim because we have a suqre matrix always
    for (int i = 0; i < Mdim; i++) {
        double a, b, c;
        for (int j = 0; j < Ndim; j++) {
            if (i == j) {
                matrix[i][j] = !stdmatrix ? getRandom(j) : 2.0;
            }
            if (i - j == 1) {
                matrix[i][j] = !stdmatrix ? getRandom(j) : 1.0;;
            }
            if (j - i == 1) {
                matrix[i][j] = !stdmatrix ? getRandom(j) : 1.0;;
            }

        }
    }
}
void ThomasAlgorithm::createMatrix(vector<vector<double> > &matrix,
                                          int Mdim, int Ndim, double arow,double brow,double crow) {
    // Mdim =Ndim because we have a sqare matrix always
    for (int i = 0; i < Mdim; i++) {
        double a, b, c;
        for (int j = 0; j < Ndim; j++) {
            if (i == j) {
                matrix[i][j] = brow;
            }
            if (i - j == 1) {
                matrix[i][j] = arow;;
            }
            if (j - i == 1) {
                matrix[i][j] = crow;
            }

        }
    }
}

void
ThomasAlgorithm::backSubstitution(vector<double> &upperDiagonal, vector<double> &rho, vector<double> &results,
                                  int size) {
    for (int i = size - 1; i >= 0; i--) {
        results[i] = rho[i] - (upperDiagonal[i] * (i + 1 == size ? 0.0 : results[i + 1]));
    }
}

void
ThomasAlgorithm::forwardSweep(vector<vector<double>> &matrix, vector<double> &rmatrix, vector<double> &upperDiagonal,
                              vector<double> &rho, int Mdim, int Ndim) {
    for (int i = 0; i < Mdim; i++) {
        double a, b, c;
        for (int j = 0; j < Ndim; j++) {
            if (i == j) {
                b = matrix[i][j];
            }
            if (i - j == 1) {
                a = matrix[i][j];
            }
            if (j - i == 1) {
                c = matrix[i][j];
            }
        }
        upperDiagonal[i] = computeUpperDiagonal(c, b, a, i == 0 ? 0.0 : upperDiagonal[i - 1]);
        rho[i] = computeRho(rmatrix[i], b, a, i == 0 ? 0.0 : upperDiagonal[i - 1], i == 0 ? 0.0 : rho[i - 1]);

    }
}

double ThomasAlgorithm::computeRho(double r, double b, double a, double prevUpperDiagonal, double prevRho) {
    if (b - (a * prevUpperDiagonal) == 0) throw "Instability detected. Perform pivoting or check you matrix";
    return (r - (a * prevRho)) / (b - (a * prevUpperDiagonal));
}

double ThomasAlgorithm::computeUpperDiagonal(double c, double b, double a, double prevUpperDiagonal) {
    if (b - (a * prevUpperDiagonal) == 0) throw "Instability detected. Perform pivoting or check you matrix";
    return c / (b - (a * prevUpperDiagonal));
}

void ThomasAlgorithm::computeThomas(vector<vector<double> > &matrix, vector<double> &rmatrix,
                                    vector<double> &upperDiagonal, vector<double> &rho, vector<double> &results,
                                    int Mdim, int Ndim) {
    forwardSweep(matrix, rmatrix, upperDiagonal, rho, Mdim, Ndim);
    backSubstitution(upperDiagonal, rho, results, Mdim);
}

void ThomasAlgorithm::print(vector<double> &matrix) {
    for (auto i = matrix.begin(); i != matrix.end(); ++i)
        std::cout << *i << ' ';
    std::cout << "\n";
}

void ThomasAlgorithm::print(vector<vector<double>> &matrix) {
    for (auto i = matrix.begin(); i != matrix.end(); ++i) {
        for (auto j = i->begin(); j != i->end(); ++j) {
            std::cout << *j << ' ';
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}


int ThomasAlgorithm::DETGTRI(vector<vector<double> > &matrix) {
    // to implement.
}
```

## Shock Tube Exact and Lax-Wendroff Two-Step Solution C++

The shock tube problem can be modelled using the Eulers Equations of state. The Lax-Wendroff solution to the problem is
found below

```c++
map<string, vector<double>> ShockTubeLaxWendroff::compute_lax_wendrof(tuple<double, double, double> leftstate,
tuple<double, double, double> rightstate,
tuple<double, double, double> geometry, double t,
int numtimesteps,
double gamma, int meshPoints) {
double pl, rhol, ul, pr, rhor, ur, xl, xr, xi;
double art_viscoticity(2.5);
vector<double> x, W1, W2, W3, F1, F2, F3, W1_p1,
W2_p1, W3_p1, F1_mh, F2_mh, F3_mh, F1_ph,
F2_ph, F3_ph, W1_ph, W2_ph, W3_ph, W1_mh, W2_mh, W3_mh, P;
tie(pl, rhol, ul) = leftstate;
tie(pr, rhor, ur) = rightstate;
tie(xl, xr, xi) = geometry;
double delX = (xr - xl) / meshPoints;
//    double delT = t / numtimesteps;

double delT = 0.25 * delX / sqrt(1.4 * max(700 + pr / rhor, 700 + pl / rhol));
cout << "Delta X: " << delX << endl;
cout << "Delta T: " << delT << endl;
x = xx(xl, xr, delX);
print("X", x);
W1_p1.resize(meshPoints);
W2_p1.resize(meshPoints);
W3_p1.resize(meshPoints);
W1_ph.resize(meshPoints);
W2_ph.resize(meshPoints);
W3_ph.resize(meshPoints);
W1_mh.resize(meshPoints);
W2_mh.resize(meshPoints);
W3_mh.resize(meshPoints);
F1_ph.resize(meshPoints);
F2_ph.resize(meshPoints);
F3_ph.resize(meshPoints);
F1_mh.resize(meshPoints);
F2_mh.resize(meshPoints);
F3_mh.resize(meshPoints);
P.resize(meshPoints);
tie(W1, W2, W3) = W(x, xi, leftstate, rightstate, gamma);
double nt = t / delT;
for (int n = 0; n < nt; n++) {
W1_ph[0] = W1[0];
W2_ph[0] = W2[0];
W3_ph[0] = W3[0];
W1_ph[meshPoints - 1] = W1[meshPoints - 1];
W2_ph[meshPoints - 1] = W2[meshPoints - 1];
W3_ph[meshPoints - 1] = W3[meshPoints - 1];
W1_mh[0] = W1[0];
W2_mh[0] = W2[0];
W3_mh[0] = W3[0];
W1_mh[meshPoints - 1] = W1[meshPoints - 1];
W2_mh[meshPoints - 1] = W2[meshPoints - 1];
W3_mh[meshPoints - 1] = W3[meshPoints - 1];
tie(F1, F2, F3, P) = flux(make_tuple(W1, W2, W3), gamma);

// predictor
for (int i = 1; i < meshPoints - 1; i++) {
W1_ph[i] = 0.5 * (W1[i] + W1[i + 1]) - ((0.5 * delT / delX) * (F1[i + 1] - F1[i]));
W2_ph[i] = 0.5 * (W2[i] + W2[i + 1]) - (0.5 * delT / delX) * (F2[i + 1] - F2[i]);
W3_ph[i] = 0.5 * (W3[i] + W3[i + 1]) - (0.5 * delT / delX) * (F3[i + 1] - F3[i]);

W1_mh[i] = 0.5 * (W1[i] + W1[i - 1]) - (0.5 * delT / delX) * (F1[i] - F1[i - 1]);
W2_mh[i] = 0.5 * (W2[i] + W2[i - 1]) - (0.5 * delT / delX) * (F2[i] - F2[i - 1]);
W3_mh[i] = 0.5 * (W3[i] + W3[i - 1]) - (0.5 * delT / delX) * (F3[i] - F3[i - 1]);
}

// Flux predictor
tie(F1_mh, F2_mh, F3_mh, P) = flux(make_tuple(W1_mh, W2_mh, W3_mh), gamma);
tie(F1_ph, F2_ph, F3_ph, P) = flux(make_tuple(W1_ph, W2_ph, W3_ph), gamma);


//corrector
// boundary conditions
W1_p1[0] = W1[0];
W2_p1[0] = W2[0];
W3_p1[0] = W3[0];
W1_p1[meshPoints - 1] = W1[meshPoints - 1];
W2_p1[meshPoints - 1] = W2[meshPoints - 1];
W3_p1[meshPoints - 1] = W3[meshPoints - 1];
for (int i = 1; i < meshPoints - 1; i++) {
W1_p1[i] = W1[i] - (delT / delX) * (F1_ph[i] - F1_mh[i]);


W2_p1[i] = W2[i] - (delT / delX) * (F2_ph[i] - F2_mh[i]
// adding the artificial viscosity terms
+
art_viscoticity * delX * delX * W1[i]
*
abs(((W2[i] / W1[i]) - (W2[i - 1] / W1[i - 1])) /
delX) *
((W2[i] / W1[i]) - (W2[i - 1] / W1[i - 1])) /
delX);
W3_p1[i] = W3[i] - (delT / delX) * (F3_ph[i] - F3_mh[i]
// adding the artificial viscosity terms
+
art_viscoticity * delX * delX * W1[i] * (W2[i] / W1[i])
*
abs(((W2[i] / W1[i]) - (W2[i - 1] / W1[i - 1])) /
delX) *
((W2[i] / W1[i]) - (W2[i - 1] / W1[i - 1])) /
delX);
}
W1 = W1_p1;
W2 = W2_p1;
W3 = W3_p1;

}
print("Pressure", P);
transform(W2_p1.begin(), W2_p1.end(), W1_p1.begin(), W2_p1.begin(), std::divides<>{});
transform(W3_p1.begin(), W3_p1.end(), W1_p1.begin(), W3_p1.end(), std::divides<>{});
return {
{
"x", x
},
{
"p", P
},
{
"rho", W1_p1
},
{
"u", W2_p1
},
{
"energy", W3_p1
}
//            {"APR",    temp}
};

}
```

### Usage

I am using matplot library for ploting the results as shown below

```c++
void plotting(map<string, vector<double>> exact, map<string, vector<double>> results, string key,
map<string, string> lookup) {
   plot(results.at("x"), results.at(key), exact.at("x"), exact.at(key));
   title(lookup.at(key));
   legend({"Lax-Wendroff", "Exact"});;
   yrange({0, *max_element(results.at(key).begin(), results.at(key).end()) + 0.1});
   xlabel("Distance");
   ylabel(lookup.at(key));
   save("L-W-" + lookup.at(key) + ".png");
}
int main() {
   tuple<double, double, double> left_state = make_tuple(1, 1., 0.);
   tuple<double, double, double> right_state = make_tuple(.1, 0.125, 0.);
   tuple<double, double, double> geometry = make_tuple(0., 1., 0.5);
   double t = .2;
   double gamma = 1.4;
   int npts = 500;
   int timesteps = 128;
   ShockTubeLaxWendroff laxWendroff;
   map<string, vector<double>> exact = laxWendroff.compute_exact(left_state, right_state, geometry, t, gamma, npts);
   map<string, vector<double>> results = laxWendroff.compute_lax_wendrof(left_state, right_state, geometry, t,
   timesteps,
   gamma, npts);
   map<string, string> lookup;
   lookup = {{"p",      "Pressure"},
   {"rho",    "Density"},
   {"u",      "Velocity"},
   {"energy", "Energy"}};
   
   plotting(exact, results, "p", lookup);
   plotting(exact, results, "rho", lookup);
   plotting(exact, results, "u", lookup);
   plotting(exact, results, "energy", lookup);


return 0;
}
```

![test](/cmake-build-debug/L-W-Velocity.png)
![test](/cmake-build-debug/L-W-Density.png)
![test](/cmake-build-debug/L-W-Pressure.png)
![test](/cmake-build-debug/L-W-Energy.png)

## Crank Nicholson and Simple Explicit Method

