
#include "ViscousBurgers.h"

using namespace std;

#include <matplot/matplot.h>

using namespace matplot;

int main() {
    ViscousBurgers vb;
    double delT = 0.2;
    vector<double> x = vb.x(-1.0, 3, 0.04);
    int numIterations = vb.it(delT, 0.4);
    vector<double> exact = vb.plotExact(delT, -1, 3.0, 0.04, 0.4, 0.2)
            .at(numIterations - 1);
    vector<double> euler = vb.solve(delT, -1, 3.0, 0.04, 0.4, 0.2, 0
    );
//    vector<double> euler = vb.plot(delT, -1, 3.0, 0.04, 0.4, 0.2
//    ).at(numIterations-1);
    vector<double> euler1 = vb.solve(delT, -1, 3.0, 0.04, 0.4, 0.2, 1);
    vector<double> euler3 = vb.solve(delT, -1, 3.0, 0.04, 0.4, 0.2, 3);
    vector<double> euler10 = vb.solve(delT, -1, 3.0, 0.04, 0.4, 0.2,
                                      10);

//    for (int i = 0; i < euler.size(); i++) {
//        euler1[i] = abs((euler1[i] - exact[i]) / exact[i]) * 100;
//        euler3[i] = abs((euler3[i] - exact[i]) / exact[i]) * 100;
//        euler10[i] = abs((euler10[i] - exact[i]) / exact[i]) * 100;
//    }

    plot(
            x, exact,
//            x, euler,
            x, euler1,
            x, euler3,
            x, euler10
    );
    title("Percent Absolute Error  at  timestep=" + to_string(delT));
    legend({
                   "Exact",
//                   "NewtonMethod(m<=0)",
                   "NewtonMethod(m<=1)",
                   "NewtonMethod(m<=3)",
                   "NewtonMethod(m<=10)"
           });
//    xrange({0,3});
    xlabel("Distance(x)");
    ylabel("Velocity(u)");
    save("Sol" + to_string(delT) + ".png");







    /*
//     * Example D while loop*/
//
//    int nIterations = 0; // number of iterations
//    const double epsilon = 0.01; // epsilon, ie precision must be <= epsilo; set epsilon to a small value
//    double preSum=10; // previous upper or lower sum, you can set it to any value or 0
//    double curSum=0.0; // you set your calculated value to this variable
//    double presition=0; // the variable that hold the precision calculated values
//
//    //presision=prevAns-currentAns
//
//    do {
//
//        curSum=upperSum(); // call upper or lower sum function
//        presition=preSum-curSum; c// calculate precisiion
//
//        preSum=curSum; // set previous sum to current sum
//
//        nIterations++; // increase iteration counts by 1 during each loop
//    }
//    while (presition<=epsilon);
//
//    cout<<nIterations; // return number of iteration for the epsilon you specified

    return 0;
}



