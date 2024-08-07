#include <iostream>
#include <fstream>
#include <cmath>

const double C = -10.2895; 

double analytical_solution(double x) {
    return (std::pow(x, 4) / 4.0) + (2.0 * std::pow(x, 3) / 3.0) - (std::pow(3, x) / std::log(3)) + C;
}

int main() {
    double x_start = -2.5;
    double x_end = 2.5;
    double h = 0.01; 

    std::ofstream file("analytical_solution.csv");

    for (double x = x_start; x <= x_end; x += h) {
        file << x << "," << analytical_solution(x) << std::endl;
    }

    file.close();
    return 0;
}
