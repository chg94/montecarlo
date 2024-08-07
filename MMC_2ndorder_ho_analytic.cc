#include <iostream>
#include <fstream>
#include <cmath>

double x(double t, double x0, double v0) {
    return x0 * cos(t) + v0 * sin(t);
}

double v(double t, double x0, double v0) {
    return -x0 * sin(t) + v0 * cos(t);
}

void generate_data(double x0, double v0, double t_start, double t_end, double dt, const char* filename) {
    std::ofstream file(filename);
    file << "t,x,v\n";  

    for (double t = t_start; t <= t_end; t += dt) {
        file << t << "," << x(t, x0, v0) << "," << v(t, x0, v0) << "\n";
    }

    file.close();
}

int main() {
    double x0 = 1.0;      
    double v0 = 0.0;      
    double t_start = 0.0; 
    double t_end = 10.0;  
    double dt = 0.01;     

    generate_data(x0, v0, t_start, t_end, dt, "harmonic_oscillator_data.csv");

    std::cout << "Datos generados en 'harmonic_oscillator_data.csv'\n";

    return 0;
}
