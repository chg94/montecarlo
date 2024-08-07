#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

const double damping = 0.1;    
const double alpha = 0.5;      
const double beta = 1.0;       
const double omega = 0.2;      /

double d2x_dt2(double x, double dx_dt, double t) {
    return -damping * dx_dt - alpha * x + beta * sin(omega * t);
}

void monte_carlo_integration(double (*func)(double, double, double), double x0, double dx_dt0, double a, double b, int n_samples, gsl_rng* r, const char* filename) {
    std::ofstream file(filename);
    file << "Tiempo Posicion Velocidad\n";

    double delta_t = (b - a) / n_samples;

    for (int i = 0; i < n_samples; ++i) {
        double x = x0;
        double dx_dt = dx_dt0;
        double t = a;

        file << t << " " << x << " " << dx_dt << "\n";

        for (; t < b; t += delta_t) {
            double d2x = func(x, dx_dt, t);
            dx_dt += d2x * delta_t;
            x += dx_dt * delta_t;
            file << t << " " << x << " " << dx_dt << "\n";
        }
    }

    file.close();
}

void solve_differential_equation(double x0, double dx_dt0, double a, double b, int n_samples, const char* filename) {
    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, 0);

    monte_carlo_integration(d2x_dt2, x0, dx_dt0, a, b, n_samples, r, filename);

    gsl_rng_free(r);
}

int main() {
    double x_initial = 0.5;     
    double dx_dt_initial = 0.0;  
    double a = 0.0;             
    double b = 50.0;             
    int n_samples = 1000;        

    solve_differential_equation(x_initial, dx_dt_initial, a, b, n_samples, "simulation_results.txt");

    return 0;
}
