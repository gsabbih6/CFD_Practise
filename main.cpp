
#include "ViscousBurgers.h"

using namespace std;

#include <matplot/matplot.h>

using namespace matplot;

int main() {
    ViscousBurgers vb;
    vector<double> x = vb.x(-1.0, 3, 0.04);
    int numIterations = vb.it(0.2, 0.4);
    vector<double> exact = vb.plotExact(0.2, -1, 3.0, 0.04, 0.4, 0.2)
            .at(numIterations - 1);
    vector<double> euler = vb.plot(0.2, -1, 3.0, 0.04, 0.4, 0.2)
            .at(numIterations - 1);

    plot(
            x, exact,
            x, euler);
    title("Solution to the Viscous Burgers Equation at different timesteps");
    legend({"Exact(0.2)", "Implicit Euler(0.2)"});
    xlabel("Distance(x)");
    ylabel("Velocity(u)");
    show();

    return 0;
}



