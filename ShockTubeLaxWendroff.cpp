//
// Created by Admin on 4/13/21.
//

#include "ShockTubeLaxWendroff.h"
#include <boost/math/tools/roots.hpp>
#include <map>

namespace tools = boost::math::tools;

double ShockTubeLaxWendroff::u_imhalf(double u_i, double u_im1, double h, double deltaT) {
    double f_i = u_i * u_i / 2;
    double f_im1 = u_im1 * u_im1 / 2;
    return 0.5 * (u_i + u_im1) - 0.5 * ((deltaT / h) * (f_i - f_im1));
}

double ShockTubeLaxWendroff::u_iphalf(double u_i, double u_ip1, double h, double deltaT) {
    double f_i = u_i * u_i / 2;
    double f_ip1 = u_ip1 * u_ip1 / 2;
    return 0.5 * (u_i + u_ip1) - 0.5 * ((deltaT / h) * (f_ip1 - f_i));
}

double ShockTubeLaxWendroff::u_leap_frog(double u_i, double u_iphalf, double u_imhalf, double h, double deltaT) {
    double f_iphalf = u_iphalf * u_iphalf / 2;
    double f_imhalf = u_imhalf * u_imhalf / 2;
    return u_i - ((deltaT / h) * (f_iphalf - f_imhalf));
}

void ShockTubeLaxWendroff::print(string msg, vector<double> v) {
    std::cout << msg << endl;
    for (auto i = v.begin(); i != v.end(); ++i)
        std::cout << *i << ' ';
    std::cout << "\n";
}

tuple<vector<double>, vector<double>, vector<double>>
ShockTubeLaxWendroff::W(vector<double> x, double xshock, tuple<double, double, double> leftstate,
                        tuple<double, double, double> rightstate, double gamma) {
    vector<double> W1, W2, W3;
    W1.resize(x.size());
    W2.resize(x.size());
    W3.resize(x.size());
    double pl, rhol, ul, pr, rhor, ur;
    tie(pl, rhol, ul) = leftstate;
    tie(pr, rhor, ur) = rightstate;
    vector<double> d;
    for (int i = x.size() / 2; i < x.size(); i++) { // use xshcock

        W1[i] = (rhor);
        W2[i] = (rhor * ur);
        W3[i] = (((pr / (gamma - 1.)) + (rhor * ur * ur * 0.5)));
    }
    for (int i = 0; i < x.size() / 2; i++) {
        W1[i] = (rhol);
        W2[i] = (rhol * ul);
        W3[i] = (((pl / (gamma - 1.)) + (rhol * ul * ul * 0.5)));
    }
//        elsese {
//            W1[i] = 0;
//            W2[i] = 0;
//            W3[i] = 0;
//        }
//}
//    print("W1", W1);
//    print("W2", W2);
//    print("W3", W3);
    return make_tuple(W1, W2, W3);
}


tuple<vector<double>, vector<double>, vector<double>, vector<double>>
ShockTubeLaxWendroff::flux(tuple<vector<double>, vector<double>, vector<double>> Ws,
                           double gamma) {
    vector<double> F1, F2, F3, P;
    vector<double> W1, W2, W3;
    tie(W1, W2, W3) = Ws;
    F1 = W2;
    F2.resize(W1.size());
    F3.resize(W1.size());
    P.resize(W1.size());

    for (int i = 0; i < W1.size(); i++) {
        P[i] = (gamma - 1.) * ((2. * W3[i] * W1[i] - W2[i] * W2[i]) / (2. * (W1[i])));
    }
    for (int i = 0; i < W1.size(); i++) {
        F2[i] = (W2[i] * W2[i] / W1[i]) + P[i];
    }
    for (int i = 0; i < W1.size(); i++) {
        F3[i] = (W2[i] / W1[i]) * (W3[i] + P[i]);
    }

//    print("W1", W1);
////    print("F1", F1);
//    print("W2", W2);
////    print("F2", F2);
//    print("W3", W3);
////    print("F3", F3);
    return make_tuple(F1, F2, F3, P);
}

int ShockTubeLaxWendroff::solutionSize(double deltaX, double startX, double endX) {

    double maxX = endX - startX;
    return (int) ((maxX / deltaX) + .5) + 1;
}

vector<double> ShockTubeLaxWendroff::xx(double startX, double endX, double deltaX) {

    vector<double> x;
    x.resize(solutionSize(deltaX, startX, endX));
    x[0] = startX;
    for (int i = 1; i < solutionSize(deltaX, startX, endX); i++) {
        x[i] = x[i - 1] + deltaX;
    }
    return x;
}

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
    return {{"x",      x},
            {"p",      P},
            {"rho",    W1_p1},
            {"u",      W2_p1},
            {"energy", W3_p1}
//            {"APR",    temp}
    };

}

double ShockTubeLaxWendroff::sound_speed(
        double gamma, double pressure, double density, double dustFrac) {
    double scale = sqrt(1 - dustFrac);
    return sqrt(gamma * pressure / density) * scale;
}

double p1;
double p5;
double rho1;
double rho5;
double gamm;
double dustFrac1;

double ShockTubeLaxWendroff::shock_tube_function(double p4) {

    double z = (p4 / p5 - 1.);
    double c1 = sound_speed(gamm, p1, rho1, dustFrac1);
    double c5 = sound_speed(gamm, p5, rho5, dustFrac1);

    double gm1 = gamm - 1.;
    double gp1 = gamm + 1.;
    double g2 = 2. * gamm;
    double fact = gm1 / g2 * (c5 / c1) * z / sqrt(1. + gp1 / g2 * z);
    fact = pow((1. - fact), (g2 / gm1));
    return p1 * fact - p4;
}

bool root_termination(double min, double max) {
    return abs(max - min) <= 0.000001;
}

tuple<tuple<double, double, double>, tuple<double, double, double>, tuple<double, double, double>,
        tuple<double, double, double>, double> ShockTubeLaxWendroff::calculate_regions(
        double pl, double ul, double rhol, double pr, double ur, double rhor, double gamma, double dustFrac) {
    rho1 = rhol;
    p1 = pl;
    double u1 = ul;
    rho5 = rhor;
    p5 = pr;
    gamm = gamma;
    dustFrac1 = dustFrac;
    double u5 = ur;

//# unless...
    if (pl < pr) {
        rho1 = rhor;
        p1 = pr;
        u1 = ur;
        rho5 = rhol;
        p5 = pl;
        u5 = ul;
    }


//# solve for post-shock pressure
//std::pair<double, double> result = tools::bisect(root_function, 0.01, 1.0, root_termination);
//pcrit = (result.first + result.second)/2;
//cout << pcrit << endl;

    std::pair<double, double> result = tools::bisect(shock_tube_function, pr, pl, root_termination);
    double p4 = result.first;
    cout << p4 << " " << result.second << "\n";;
    //(shock_tube_function(p4, p1, p5, rho1, rho5, gamma,
    // dustFrac);// tools::bisect(shock_tube_function, 0.01, 1.0,

//                                    cout<<p4<<"\n";
    // [](double min, double max) { abs(max - min) <= 0.000001; }).first;

//# compute post-shock density and velocity
    double z = (p4 / p5 - 1.);
    double c5 = sound_speed(gamma, p5, rho5, dustFrac);

    double gm1 = gamma - 1.;
    double gp1 = gamma + 1.;
    double gmfac1 = 0.5 * gm1 / gamma;
    double gmfac2 = 0.5 * gp1 / gamma;

    double fact = sqrt(1. + gmfac2 * z);

    double u4 = c5 * z / (gamma * fact);
    double rho4 = rho5 * (1. + gmfac2 * z) / (1. + gmfac1 * z);

//# shock speed
    double w = c5 * fact;

//# compute values at foot of rarefaction
    double p3 = p4;
    double u3 = u4;
    double rho3 = rho1 * pow((p3 / p1), (1. / gamma));
    cout << p1 << " " << p3 << " " << p4 << " " << p5 << " " << w << "\n";
    return make_tuple(make_tuple(p1, rho1, u1),
                      make_tuple(p3, rho3, u3),
                      make_tuple(p4, rho4, u4),
                      make_tuple(p5, rho5, u5), w);
}

tuple<double, double, double, double> ShockTubeLaxWendroff::calc_positions(double pl, double pr,
                                                                           tuple<double, double, double> region1,
                                                                           tuple<double, double, double> region3,
                                                                           double w, double xi, double t, double gamma,
                                                                           double dustFrac) {
    double p1;
    double rho1;
    double u1;
    double p3;
    double rho3;
    double u3;
    double xsh;
    double xcd;
    double xft;
    double xhd;
    tie(p1, rho1, u1) = region1; //# don't need velocity
    tie(p3, rho3, u3) = region3;
    double c1 = sound_speed(gamma, p1, rho1, dustFrac);
    double c3 = sound_speed(gamma, p3, rho3, dustFrac);

    if (pl > pr) {
        xsh = xi + w * t;
        xcd = xi + u3 * t;
        xft = xi + (u3 - c3) * t;
        xhd = xi - c1 * t;
    } else {
//# pr > pl
        xsh = xi - w * t;
        xcd = xi - u3 * t;
        xft = xi - (u3 - c3) * t;
        xhd = xi + c1 * t;
    }
    return make_tuple(xhd, xft, xcd, xsh);


}

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N) {
    double h = (b - a) / static_cast<double>(N - 1);
    std::vector<double> xs(N);
    std::vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
        *x = val;
    }
    return xs;
}

tuple<vector<double>, vector<double>, vector<double>, vector<double>>
ShockTubeLaxWendroff::create_arrays(double pl, double pr, double xl,
                                    double xr, tuple<double, double, double, double> positions,
                                    tuple<double, double, double> state1, tuple<double, double, double> state3,
                                    tuple<double, double, double> state4, tuple<double, double, double> state5,
                                    double npts, double gamma, double t, double xi, double dustFrac) {
    double xhd, xft, xcd, xsh, p1, rho1, u1, p3, rho3, u3, p4, rho4, u4, p5, rho5, u5, fact;

    tie(xhd, xft, xcd, xsh) = positions;
    tie(p1, rho1, u1) = state1;
    tie(p3, rho3, u3) = state3;
    tie(p4, rho4, u4) = state4;
    tie(p5, rho5, u5) = state5;
    double gm1 = gamma - 1.;
    double gp1 = gamma + 1.;

    vector<double> x_arr = LinearSpacedArray(xl, xr, npts);
    vector<double> rho = LinearSpacedArray(0, 0, npts);
    vector<double> p = LinearSpacedArray(0, 0, npts);
    vector<double> u = LinearSpacedArray(0, 0, npts);
    double c1 = sound_speed(gamma, p1, rho1, dustFrac);
    if (pl > pr) {
        for (int i = 0; i < x_arr.size(); i++) {
            if (x_arr[i] < xhd) {
                rho[i] = rho1;
                p[i] = p1;
                u[i] = u1;

            } else if (x_arr[i] < xft) {
                u[i] = 2. / gp1 * (c1 + (x_arr[i] - xi) / t);
                fact = 1. - (0.5 * gm1 * u[i] / c1);
                rho[i] = rho1 * pow(fact, (2. / gm1));
                p[i] = p1 * pow(fact, (2. * gamma / gm1));
            } else if (x_arr[i] < xcd) {
                rho[i] = rho3;
                p[i] = p3;
                u[i] = u3;
            } else if (x_arr[i] < xsh) {
                rho[i] = rho4;
                p[i] = p4;
                u[i] = u4;
            } else {
                rho[i] = rho5;
                p[i] = p5;
                u[i] = u5;
            }
        }
    } else {
        for (int i = 0; i < x_arr.size(); i++) {
            if (x_arr[i] < xsh) {
                rho[i] = rho5;
                p[i] = p5;
                u[i] = -u1;
            } else if (x_arr[i] < xcd) {
                rho[i] = rho4;
                p[i] = p4;
                u[i] = -u4;
            } else if (x_arr[i] < xft) {
                rho[i] = rho3;
                p[i] = p3;
                u[i] = -u3;
            } else if (x_arr[i] < xhd) {
                u[i] = -2. / gp1 * (c1 + (xi - x_arr[i]) / t);
                fact = 1. + 0.5 * gm1 * u[i] / c1;
                rho[i] = rho1 * pow(fact, (2. / gm1));
                p[i] = p1 * pow(fact, (2. * gamma / gm1));
            } else {
                rho[i] = rho1;
                p[i] = p1;
                u[i] = -u1;
            }
        }
    }
//    cout << rho5 << endl;
    return make_tuple(x_arr, p, rho, u);

}

map<string, vector<double>> ShockTubeLaxWendroff::compute_exact(tuple<double, double, double> leftstate,
                                                                tuple<double, double, double> rightstate,
                                                                tuple<double, double, double> geometry,
                                                                double t,
                                                                double gamma, double npts, double dustFrac) {
    double pl;
    double rhol;
    double ul;
    double pr;
    double rhor;
    double ur;
    double xl;
    double xr;
    double xi;
    vector<double> x;
    vector<double> p;
    vector<double> rho;
    vector<double> u;
    tie(pl, rhol, ul) = leftstate;
    tie(pr, rhor, ur) = rightstate;
    tie(xl, xr, xi) = geometry;

//# basic checking
    if (xl >= xr) {
        std::cout << "xl has to be less than xr!";
        exit(0);
    }
    if (xi >= xr || xi <= xl) {
        std::cout << "xi has in between xl and xr!";
        exit(0);
    }

//# calculate regions
    tuple<double, double, double> region1;
    tuple<double, double, double> region3;
    tuple<double, double, double> region4;
    tuple<double, double, double> region5;
    double w;
    tie(region1, region3, region4, region5, w) =
            calculate_regions(pl, ul, rhol, pr, ur, rhor, gamma, dustFrac);

//    regions = region_states(pl, pr, region1, region3, region4, region5);

//# calculate positions
    tuple<double, double, double, double> x_positions;
    x_positions = calc_positions(pl, pr, region1, region3, w, xi, t, gamma,
                                 dustFrac);
    tuple<string, string, string, string> pos_description;
    pos_description = make_tuple("Head of Rarefaction", "Foot of Rarefaction",
                                 "Contact Discontinuity", "Shock");
//    positions = dict(zip(pos_description, x_positions));

//# create arrays
    tie(x, p, rho, u) = create_arrays(pl, pr, xl, xr, x_positions,
                                      region1, region3, region4, region5,
                                      npts, gamma, t, xi, dustFrac);
    std::vector<double> energy, rho_total;
    energy.resize(npts);
    rho_total.resize(npts);
    for (int i = 0; i < npts; i++) {
        energy[i] = p[i] / (gamma - 1.0);
    }
    for (int i = 0; i < npts; i++) {
        rho_total[i] = rho[i] / (1.0 - dustFrac);
    }
    map<string, vector<double>> val_dict;
    val_dict = {{"x",      x},
                {"p",      p},
                {"rho",    rho},
                {"u",      u},
                {"energy", energy},
                {"APR",    rho_total}
    };

    return val_dict;

}


