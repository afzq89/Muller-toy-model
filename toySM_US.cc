// A 2-dimensional overdamped langevin dynamics on Muller potential surface

#include <stdio.h>
#include <cmath>
#include <time.h>
#include <cassert>
#include <iostream>

#define Boltzman 0.001987191683

using namespace std;

double randn(double mu, double sigma) {
        static bool deviateAvailable=false;
        static float storedDeviate;
        double polar, rsquared, var1, var2;
       
        //        If no deviate has been stored, the polar Box-Muller transformation is
        //        performed, producing two independent normally-distributed random
        //        deviates.  One is stored for the next round, and one is returned.
        if (!deviateAvailable) {
               
                //        choose pairs of uniformly distributed deviates, discarding those
                //        that don't fall within the unit circle
                do {
                        var1=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
                        var2=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
                        rsquared=var1*var1+var2*var2;
                } while ( rsquared>=1.0 || rsquared == 0.0);
               
                //        calculate polar tranformation for each deviate
                polar=sqrt(-2.0*log(rsquared)/rsquared);
               
                //        store first deviate and set flag
                storedDeviate=var1*polar;
                deviateAvailable=true;
               
                //        return second deviate
                return var2*polar*sigma + mu;
        }
       
        //        If a deviate is available from a previous call to this function, it is
        //        returned, and the flag is set to false.
        else {
                deviateAvailable=false;
                return storedDeviate*sigma + mu;
        }
}


int main(int argc,char *argv[]) {
        
    // reading the parameters    
    string key;

    double temp=300; // temperature
    double damping=1; // damping coefficient
    double dt=1e-6; // time step
    double T=1; // total time
    double t0=0; //release time
    double K=0; //harmonic bias
    double R=1; // reporting frequency
    int N=2; // dimension

    //     Muller Potential Parameters
    double aa [4] = { -1, -1, -6.5, 0.7 };
    double bb [4] = { 0, 0, 11, 0.6 };
    double cc [4] = { -10, -10, -6.5, 0.7 };
    double AA [4] = { -200, -100, -170, 15 };
    double XX [4] = { 1, 0, -0.5, -1 };
    double YY [4] = { 0, 0.5, 1.5, 1 };
    double k = 1; // potential strength

    double* x0= new double[N];
    double* xc= new double[N];

    while(!std::cin.eof()) {
	std::cin >> key;
	if (key=="timestep") {
	    std::cin >> dt;
	    std::cout << "# "<< key << " is " << dt << "\n";
	} else if (key=="runtime") {
	    std::cin >> T;
	    std::cout << "# "<< key << " is " << T << "\n";
	} else if (key=="temperature") {
	    std::cin >> temp;
	    std::cout << "# " << key << " is " << temp << "\n";
	} else if (key=="damping") {
	    std::cin >> damping;
	    std::cout << "# " << key << " is " << damping << "\n";
	} else if (key=="strength") {
	    std::cin >> k;
	    std::cout << "# "<< key << " is " << k << "\n";
	} else if (key=="initial") {
	    for (int i=0; i<N; i++) std::cin >> x0[i];
	    std::cout << "# "<< key << " is ( ";
	    for (int i=0; i<N; i++) std::cout << x0[i] << " ";
	    std::cout << ")\n";
	} else if (key=="center") {
	    for (int i=0; i<N; i++) std::cin >> xc[i];
	    std::cout << "# "<< key << " is ( ";
	    for (int i=0; i<N; i++) std::cout << xc[i] << " ";
	    std::cout << ")\n";
	} else if (key=="report") {
	    std::cin >> R;
	    std::cout << "# "<< key << " is " << R << "\n";
	} else if (key=="US") {
	    std::cin >> K;
	    std::cout << "# "<< key << " is " << K << "\n";
	} else if (key=="US0") {
	    std::cin >> t0;
	    std::cout << "# "<< key << " is " << t0 << "\n";
	} else std::cin.ignore(0,'\n');
    }

    // running the verlet algorithm

    srand ( time(NULL)+((argc>1)?atoi(argv[1]):0) );

    double t;
    double* x = new double[N];
    double* v = new double[N];

    double beta = 1/(Boltzman*temp);
    double D = 1/(beta*damping);

    t=0;
    for (int i=0; i<N; i++) x[i]=x0[i];

    while (t<T+dt) {

//	for (int i=0; i<N; i++) {
//	    r2+=x[i]*x[i];
//	}

	// defining the potential
    	// U=k*A*exp(a(x-X)^2+b(x-X)(y-Y)+c(y-Y)^2)
	// defining the biasing potential
	// V(x[i])=K/2 (x[i]-xc[i])^2 -> -dV/dx[i]=-K(x[i]-xc[i])
	double dU [2] = { 0, 0 };
	for (int j=0; j<4; j++) {
	    double U = k*AA[j]*exp(aa[j]*(x[0]-XX[j])*(x[0]-XX[j])+bb[j]*(x[0]-XX[j])*(x[1]-YY[j])+cc[j]*(x[1]-YY[j])*(x[1]-YY[j]));
	    dU[0] += U*(2*aa[j]*(x[0]-XX[j])+bb[j]*(x[1]-YY[j]));
	    dU[1] += U*(2*cc[j]*(x[1]-YY[j])+bb[j]*(x[0]-XX[j]));
	}
	for (int i=0; i<N; i++) {
	    v[i] = sqrt(2*D/dt) * randn(0,1) - beta * D * ( dU[i] + (t<t0?K:0) * (x[i]-xc[i]) );
	    x[i] += dt * v[i];
	    if (R>0 && int(t/dt)%int(R/dt)==0) {
		if (i==0) std::cout << t;
		std::cout << " " << x[i];
		if (i==N-1) std::cout << "\n";
	    }
	}
	t += dt;
    }
}
