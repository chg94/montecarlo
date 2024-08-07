#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


double heaviside(double x) {
    return (x >= 0) ? 1.0 : 0.0;
}

double func(double y, double x) {
    return x * sqrt(fabs(y)) + pow(sin(x * M_PI / 2), 3) - 5 * heaviside(x - 2);
}

double monte_carlo_integration(double (*func)(double, double), double y0, double a, double b, int n_samples, gsl_rng* r) {
    double sum = 0.0;
    for (int i = 0; i < n_samples; ++i) {
        double x = gsl_ran_flat(r, a, b);
        sum += func(y0, x);
    }
    return (b - a) * sum / n_samples;
}

int main() {
    double a = -4.0;
    double b = 5.0; 
    int n_samples = 100000;
    int n_generations = 100;
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    std::ofstream file("monte_carlo_differential_results2.txt");

    double y0 = 4.0;
    double y = y0;
    double delta_x = (b - a) / n_generations;

    for (int i = 0; i < n_generations; ++i) {
        double x = a + i * delta_x;
        double dy = monte_carlo_integration(func, y, x, x + delta_x, n_samples, r);
        y += dy;
        file << x << " " << y << "\n";
    }

    file.close();

    gsl_rng_free(r);

    return 0;
}
