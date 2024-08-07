#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double func(double x) {
    return pow(x, 3) + 2 * pow(x, 2) - pow(3, x);
}

double monte_carlo_integration(double (*func)(double), double a, double b, int n_samples, gsl_rng* r) {
    double sum = 0.0;
    for (int i = 0; i < n_samples; ++i) {
        double x = gsl_ran_flat(r, a, b);
        sum += func(x);
    }
    return (b - a) * sum / n_samples;
}

void solve_differential_equation(double (*func)(double), double y0, double a, double b, int n_samples, int n_generations, const char* filename) {
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    std::ofstream file(filename);

    double y = y0;
    double delta_x = (b - a) / n_generations;

    for (int i = 0; i < n_generations; ++i) {
        double x = a + i * delta_x;
        double dy = monte_carlo_integration(func, x, x + delta_x, n_samples, r);
        y += dy;
        file << x << " " << y << "\n";
    }

    file.close();

    gsl_rng_free(r);
}

int main() {
    double a = -2.5;
    double b = 2.5;
    int n_samples = 100000;
    int n_generations = 100;
    double y0 = -11.2;

    solve_differential_equation(func, y0, a, b, n_samples, n_generations, "monte_carlo_differential_results1.txt");

    return 0;
}
