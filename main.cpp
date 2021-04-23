
#include "ViscousBurgers.h"

using namespace std;

#include <matplot/matplot.h>
#include "ShockTubeLaxWendroff.h"
#include <vector>
#include <map>

using namespace matplot;
using namespace std;

void print(vector<double> &matrix) {
    for (auto i = matrix.begin(); i != matrix.end(); ++i)
        std::cout << *i << ' ';
    std::cout << "\n";
}

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
//    ViscousBurgers vb;
//    double delT = 0.2;
//    vector<double> x = vb.x(-1.0, 3, 0.04);
//    int numIterations = vb.it(delT, 0.4);
//    vector<double> exact = vb.plotExact(delT, -1, 3.0, 0.04, 0.4, 0.2)
//            .at(numIterations - 1);
//    vector<double> euler = vb.solve(delT, -1, 3.0, 0.04, 0.4, 0.2, 0
//    );
////    vector<double> euler = vb.plot(delT, -1, 3.0, 0.04, 0.4, 0.2
////    ).at(numIterations-1);
//    vector<double> euler1 = vb.solve(delT, -1, 3.0, 0.04, 0.4, 0.2, 1);
//    vector<double> euler3 = vb.solve(delT, -1, 3.0, 0.04, 0.4, 0.2, 3);
//    vector<double> euler10 = vb.solve(delT, -1, 3.0, 0.04, 0.4, 0.2,
//    10);

//    for (int i = 0; i < euler.size(); i++) {
//        euler1[i] = abs((euler1[i] - exact[i]) / exact[i]) * 100;
//        euler3[i] = abs((euler3[i] - exact[i]) / exact[i]) * 100;
//        euler10[i] = abs((euler10[i] - exact[i]) / exact[i]) * 100;
//    }
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



