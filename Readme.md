# Computational Fluid Dynamics Algorithms

This repositotory is a compilation of the application of several algorithms that 
are used the solution of Partial Differential equations as taught during the semester. 
## Finite Difference methods for solving PDEs
The implementation is found  the FiniteDifference.cpp file which contains an application of the FTCS, FTBS and Lax 
methods. This include the following methods
    Forward Time Central Space
    Forward Time Backward Space", 
    Lax-Wendoff &
    MacCormack

Simply calling the method  ```myplot(int maxX, double deltaX, int maxT, double c);``` returns the sample plot below

## Thomas Algorthm
The thomas algorithm is used to solve tridiagonal systems of equations
Use the method```computeThomas(vector<vector<double> > &matrix, vector<double> &rmatrix,
vector<double> &upperDiagonal, vector<double> &rho, vector<double> &results,
int Mdim, int Ndim)```

## Viscous Burgers Method

see code

## Crank Nicholson and Simple Explicit Method

