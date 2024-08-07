#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector> 

double integrand(double x) {
    return exp(-x * x);
}

double monte_carlo_integration(double (*func)(double), double a, double b, int n_samples, gsl_rng* r) {
    double sum = 0.0;

    for (int i = 0; i < n_samples; ++i) {
        double x = a + gsl_rng_uniform(r) * (b - a);
        sum += func(x);
    }

    return sum * (b - a) / n_samples;
}

int main() {
    double a = -5.0; 
    double b = 5.0;  

    std::vector<int> sample_sizes = {1000, 10000, 100000, 1000000};

    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, 0);

    std::ofstream file("monte_carlo_results.csv");
    file << "SampleSize,IntegralResult\n";

    for (int n_samples : sample_sizes) {
        double result = monte_carlo_integration(integrand, a, b, n_samples, r);
        file << n_samples << "," << result << "\n";
    }

    file.close();
    gsl_rng_free(r);

    return 0;
}
