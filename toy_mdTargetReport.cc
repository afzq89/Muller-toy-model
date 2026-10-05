#include <stdio.h>
#include <cmath>
#include <ctime>
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

// Euclidean distance in N dimensions
double distance(double *a, double *b, int N) {
    double s = 0;
    for(int i=0;i<N;i++) s += (a[i]-b[i])*(a[i]-b[i]);
    return sqrt(s);
}

int main(int argc, char* argv[]) {
    // Parameters
    string key;
    double temp = 300;
    double damping = 1;
    double dt = 1e-6;
    double T = 1;
    double K = 0;
    int N = 2;
    double k = 1;
    int repeat = 10000;
    double R = 0.001;
    string outfilename = "traj.txt.gz";
    int seed = time(NULL);

    double* x0 = new double[N];
    double* xc = new double[N];
    double* target = new double[N];
    double min_distance = 0.1;

    for (int i = 0; i < N; i++) { x0[i] = 0.0; xc[i] = 0.0; target[i] = 0.0; }

    // Parse output filename and seed from command line
    if (argc > 1) {
        outfilename = argv[1];
    }
    if (argc > 2) {
        seed = atoi(argv[2]);
    }
    srand(seed);

    // Muller potential parameters
    double aa[4] = { -1, -1, -6.5, 0.7 };
    double bb[4] = { 0, 0, 11, 0.6 };
    double cc[4] = { -10, -10, -6.5, 0.7 };
    double AA[4] = { -200, -100, -170, 15 };
    double XX[4] = { 1, 0, -0.5, -1 };
    double YY[4] = { 0, 0.5, 1.5, 1 };

    // Parse input parameters from stdin
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
        } else if (key == "US") {
            std::cin >> K;
        } else if (key == "target") {
            for (int i = 0; i < N; i++) std::cin >> target[i];
        } else if (key == "min_distance") {
            std::cin >> min_distance;
        } else if (key == "repeat") {
            std::cin >> repeat;
        } else if (key == "report") {
            std::cin >> R;
        } else {
            std::cin.ignore(0, '\n');
        }
    }

    std::string cmd = "gzip > " + outfilename;
    FILE* out = popen(cmd.c_str(), "w");
    if (!out) {
        std::cerr << "Failed to open gzip pipe\n";
        exit(1);
    }

    // Try up to 'repeat' times
    bool success = false;
    double* x = new double[N];
    double* v = new double[N];
    double beta = 1.0 / (Boltzman * temp);
    double D = 1.0 / (beta * damping);
    double final_time = -1;
    for(int trial = 0; trial < repeat; trial++) {
        double t = 0;
        for(int i=0;i<N;i++) x[i]=x0[i];

        double next_output = 0.0; // next output time

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
                v[i] = sqrt(2 * D / dt) * randn(0, 1) - beta * D * (dU[i] + K * (x[i] - xc[i]));
                x[i] += dt * v[i];
            }
            t += dt;

            if (R > 0 && t >= next_output) {
                fprintf(out, "%.8f", t);
                for (int i = 0; i < N; i++)
                    fprintf(out, " %.8f", x[i]);
                fprintf(out, "\n");
                next_output += R;
            }

            if (distance(x, target, N) < min_distance) {
                final_time = t;
                success = true;
                break;
            }
        }
        if (success) break; // done!
    }

    pclose(out);

    // Output: final time and coordinates (to stdout)
    if (success) {
        printf("%.8f", final_time);
        for (int i=0; i<N; i++) printf(" %.8f", x[i]);
        printf("\n");
    } else {
        printf("-1");
        for (int i=0; i<N; i++) printf(" %.8f", x[i]);
        printf("\n");
    }

    delete[] x0;
    delete[] xc;
    delete[] x;
    delete[] v;
    delete[] target;
    return 0;
}

