#include <stdio.h>
#include <cmath>
#include <time.h>
#include <cassert>
#include <iostream>
#include <vector>
#include <string>

#define Boltzman 0.001987191683

using namespace std;

// Gaussian random number generator (Box-Muller)
double randn(double mu, double sigma) {
    static bool deviateAvailable = false;
    static float storedDeviate;
    double polar, rsquared, var1, var2;

    if (!deviateAvailable) {
        do {
            var1 = 2.0 * (double(rand()) / double(RAND_MAX)) - 1.0;
            var2 = 2.0 * (double(rand()) / double(RAND_MAX)) - 1.0;
            rsquared = var1 * var1 + var2 * var2;
        } while (rsquared >= 1.0 || rsquared == 0.0);

        polar = sqrt(-2.0 * log(rsquared) / rsquared);

        storedDeviate = var1 * polar;
        deviateAvailable = true;

        return var2 * polar * sigma + mu;
    } else {
        deviateAvailable = false;
        return storedDeviate * sigma + mu;
    }
}

int main(int argc, char* argv[]) {
    // Parameters
    string key;
    double temp = 300;         // temperature
    double damping = 1;        // damping coefficient
    double dt = 1e-6;          // time step
    double t0=0; 	       //release time
    double T = 1;              // total time
    double K = 0;              // harmonic bias
    double R = 1e-3;           // reporting frequency (default every 1e-3 units)
    int N = 2;                 // dimension

    // Muller potential parameters
    double aa[4] = { -1, -1, -6.5, 0.7 };
    double bb[4] = { 0, 0, 11, 0.6 };
    double cc[4] = { -10, -10, -6.5, 0.7 };
    double AA[4] = { -200, -100, -170, 15 };
    double XX[4] = { 1, 0, -0.5, -1 };
    double YY[4] = { 0, 0.5, 1.5, 1 };
    double k = 1; // potential strength

    double* x0 = new double[N];
    double* xc = new double[N];
    for (int i = 0; i < N; i++) { x0[i] = 0.0; xc[i] = 0.0; }

    // Parse parameters from stdin (or file redirected to stdin)
    while (!std::cin.eof()) {
        std::cin >> key;
        if (key == "timestep") {
            std::cin >> dt;
        } else if (key == "runtime") {
            std::cin >> T;
        } else if (key == "temperature") {
            std::cin >> temp;
        } else if (key == "damping") {
            std::cin >> damping;
        } else if (key == "strength") {
            std::cin >> k;
        } else if (key == "initial") {
            for (int i = 0; i < N; i++) std::cin >> x0[i];
        } else if (key == "center") {
            for (int i = 0; i < N; i++) std::cin >> xc[i];
        } else if (key == "report") {
            std::cin >> R;
        } else if (key == "US") {
            std::cin >> K;
        } else {
            std::cin.ignore(0, '\n');
        }
    }

    std::string outfilename = "traj.txt.gz";
    if (argc > 1) {
        outfilename = argv[1];
    }
    std::string cmd = "gzip > " + outfilename;
    FILE* out = popen(cmd.c_str(), "w");
    if (!out) {
        std::cerr << "Failed to open gzip pipe\n";
        exit(1);
    }

    // Write header
    fprintf(out, "# t");
    for (int i = 0; i < N; i++) fprintf(out, " x%d", i);
    fprintf(out, "\n");

    // Seed RNG
    srand(time(NULL) + ((argc > 2) ? atoi(argv[2]) : 0));

    double t = 0.0;
    double* x = new double[N];
    double* v = new double[N];
    double beta = 1.0 / (Boltzman * temp);
    double D = 1.0 / (beta * damping);

    for (int i = 0; i < N; i++) x[i] = x0[i];

    // Main simulation loop
    while (t < T + dt) {
        double dU[2] = { 0, 0 };
        for (int j = 0; j < 4; j++) {
            double dx = x[0] - XX[j];
            double dy = x[1] - YY[j];
            double U = k * AA[j] * exp(aa[j] * dx * dx + bb[j] * dx * dy + cc[j] * dy * dy);
            dU[0] += U * (2 * aa[j] * dx + bb[j] * dy);
            dU[1] += U * (2 * cc[j] * dy + bb[j] * dx);
        }

        for (int i = 0; i < N; i++) {
            v[i] = sqrt(2*D/dt) * randn(0,1) - beta * D * ( dU[i] + K * (x[i]-xc[i]) );
            x[i] += dt * v[i];
        }

        // Output every R units
        if (R > 0 && int(t / dt) % int(R / dt) == 0) {
            fprintf(out, "%.8f", t);
            for (int i = 0; i < N; i++)
                fprintf(out, " %.8f", x[i]);
            fprintf(out, "\n");
        }

        t += dt;
    }

    // Cleanup
    pclose(out);
    delete[] x0;
    delete[] xc;
    delete[] x;
    delete[] v;
    return 0;
}
