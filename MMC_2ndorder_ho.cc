#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

const double omega = 1.0;

double du_dx(double y) {
    return -omega * omega * y;
}

double dy_dx(double u) {
    return u;
}

double monte_carlo_integration(double (*func)(double), double a, double b, int n_samples, gsl_rng* r) {
    double sum = 0.0;
    for (int i = 0; i < n_samples; ++i) {
        double x = gsl_ran_flat(r, a, b);
        sum += func(x);
    }
    return (b - a) * sum / n_samples;
}

void solve_differential_equation_system(double y0, double u0, double a, double b, int n_samples, int n_generations, const char* filename) {
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    std::ofstream file(filename);

    double y = y0;
    double u = u0;
    double delta_x = (b - a) / n_generations;

    for (int i = 0; i < n_generations; ++i) {
        double x = a + i * delta_x;
        
        auto dy_func = [u](double) { return dy_dx(u); };
        auto du_func = [y](double) { return du_dx(y); };

        double dy_sum = 0.0;
        for (int j = 0; j < n_samples; ++j) {
            double x = gsl_ran_flat(r, x, x + delta_x);
            dy_sum += dy_func(x);
        }
        double dy = (delta_x) * dy_sum / n_samples;
        
        double du_sum = 0.0;
        for (int j = 0; j < n_samples; ++j) {
            double x = gsl_ran_flat(r, x, x + delta_x);
            du_sum += du_func(x);
        }
        double du = (delta_x) * du_sum / n_samples;

        y += dy;
        u += du;
        
        file << x << " " << y << " " << u << "\n";
    }

    file.close();
    gsl_rng_free(r);
}

int main() {
    double a = 0.0;
    double b = 10.0;
    int n_samples = 100000;
    int n_generations = 100;
    double y0 = 1.0;  
    double u0 = 0.0;  

    solve_differential_equation_system(y0, u0, a, b, n_samples, n_generations, "monte_carlo_differential_results_oscillator.txt");

    return 0;
}
