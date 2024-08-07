#include <iostream>
#include <fstream>
#include <cmath>

double heaviside(double x) {
    return (x >= 2) ? 1.0 : 0.0;
}

double dydx(double x, double y) {
    return x * std::sqrt(std::abs(y)) + std::pow(std::sin(M_PI * x / 2), 3) - 5 * heaviside(x);
}

void runge_kutta(double x0, double y0, double x_end, double h, const char* filename) {
    std::ofstream file(filename);
    file << "x,y\n";

    double x = x0;
    double y = y0;

    while (x < x_end) {
        file << x << "," << y << "\n";

        double k1 = h * dydx(x, y);
        double k2 = h * dydx(x + 0.5 * h, y + 0.5 * k1);
        double k3 = h * dydx(x + 0.5 * h, y + 0.5 * k2);
        double k4 = h * dydx(x + h, y + k3);

        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        x += h;
    }

    file.close();
}

int main() {
    double x_initial = -4.0; 
    double y_initial = 4.0; 
    double x_final = 5.0;   
    double step_size = 0.0001;  

    runge_kutta(x_initial, y_initial, x_final, step_size, "solution_data.csv");

    return 0;
}
