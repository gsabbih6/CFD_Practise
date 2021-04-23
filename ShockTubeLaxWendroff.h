//
// Created by Admin on 4/13/21.
//

#ifndef HW_SHOCKTUBELAXWENDROFF_H
#define HW_SHOCKTUBELAXWENDROFF_H

#include <vector>
#include <map>

using namespace std;

class ShockTubeLaxWendroff {
public:
    // region 1,3,4,5 not the rare-fraction region 2
    tuple<tuple<double, double, double>, tuple<double, double, double>, tuple<double, double, double>,
            tuple<double, double, double>, double> calculate_regions(double pl, double ul, double rhol,
                                                                     double pr, double ur, double rhor,
                                                                     double gamma, double dustFrac = 0.);

    static double
    shock_tube_function(double p4);

    static double sound_speed(double gamma, double pressure, double density, double dustFrac = 0.);

    map<string, vector<double>> compute_exact(tuple<double, double, double> leftstate,
                                              tuple<double, double, double> rightstate,
                                              tuple<double, double, double> geometry,
                                              double t, double gamma = 1.4, double meshPoints = 32,
                                              double dustFrac = 0.);

    tuple<double, double, double, double>
    calc_positions(double pl, double pr, tuple<double, double, double> region1, tuple<double, double, double> region3,
                   double w, double xi, double t, double gamma, double dustFrac = 0.);

    tuple<vector<double>, vector<double>, vector<double>, vector<double>>
    create_arrays(double pl, double pr, double xl, double xr, tuple<double, double, double, double> positions,
                  tuple<double, double, double> state1, tuple<double, double, double> state3,
                  tuple<double, double, double> state4, tuple<double, double, double> state5,
                  double npts, double gamma, double t, double xi, double dustFrac = 0.);

    double u_imhalf(double u_i, double u_im1, double h, double deltaT);

    double u_iphalf(double u_i, double u_jp1, double h, double deltaT);

    double u_leap_frog(double u_i, double u_iphalf, double u_imhalf, double h, double deltaT);

//    map<string, vector<double>> compute_lax_wendrof(tuple<double, double, double> leftstate,
//                                                    tuple<double, double, double> rightstate,
//                                                    tuple<double, double, double> geometry,
//                                                    double t, double gamma = 1.4, double meshPoints = 32,
//                                                    double dustFrac = 0.);

    vector<double> init(vector<double> x, double xshock, tuple<double, double> state);

//    vector<double> impl_lax_wendroff(vector<double> initialstate,double deltaT, double tmax, double deltaX);
    void print(string msg, vector<double> v);

    vector<double>
    impl_lax_wendroff(vector<double> W_j, vector<double> F_j, int n, int i, double deltaT, double deltaX);

    map<string, vector<double>>
    compute_lax_wendrof(tuple<double, double, double> leftstate,
                        tuple<double, double, double> rightstate,
                        tuple<double, double, double> geometry, double t,
                        int numtimesteps,
                        double gamma, int meshPoints);

    int solutionSize(double deltaX, double startX, double endX);

    tuple<vector<double>, vector<double>, vector<double>>
    W(vector<double> x, double xshock, tuple<double, double, double> leftstate,
      tuple<double, double, double> rightstate, double gamma);

    tuple<vector<double>, vector<double>, vector<double>>
    F(vector<double> x, double xshock, tuple<double, double, double>
    leftstate, tuple<double, double, double> rightstate, double gamma
    );

    tuple<vector<double>, vector<double>, vector<double>, vector<double>>
    flux(tuple<vector<double>, vector<double>, vector<double>> Ws,
         double gamma);

    vector<double> xx(double startX, double endX, double deltaX);
};


#endif //HW_SHOCKTUBELAXWENDROFF_H
