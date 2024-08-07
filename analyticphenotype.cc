#include <iostream>
#include <fstream>
#include <cmath>


const double damping = 0.1;    
const double alpha = 0.5;      
const double beta = 1.0;       
const double omega = 0.2;      

double compute_omega0() {
    return sqrt(alpha - pow(damping / 2, 2));
}

double x_h(double t, double C1, double C2, double omega0) {
    return exp(-damping / 2 * t) * (C1 * cos(omega0 * t) + C2 * sin(omega0 * t));
}

double x_p(double t, double A, double B) {
    return A * sin(omega * t) + B * cos(omega * t);
}

void compute_AB(double &A, double &B) {
    double omega_sq = omega * omega;
    A = (beta + damping * omega * B) / (alpha + omega_sq);
    B = (damping * A * omega) / (alpha + omega_sq);
}

double v_h(double t, double C1, double C2, double omega0) {
    return exp(-damping / 2 * t) * (-C1 * omega0 * sin(omega0 * t) + C2 * omega0 * cos(omega0 * t) - damping / 2 * (C1 * cos(omega0 * t) + C2 * sin(omega0 * t)));
}

double v_p(double t, double A, double B) {
    return A * omega * cos(omega * t) - B * omega * sin(omega * t);
}


int main() {
    double C1 = 1.0;     
    double C2 = 0.0;     
    double x_initial = 0.5; 
    double dx_dt_initial = 0.0; 

    double omega0 = compute_omega0();
    
    double A, B;
    compute_AB(A, B);

    std::ofstream file("solution_dataphenotype.csv");

    file << "Time,Homogeneous.Solution,Particular.Solution,Total.Solution,Homogeneous.Velocity,Particular.Velocity,Total.Velocity\n";

    double t_start = 0.0;
    double t_end = 50.0;
    double dt = 0.1;

    for (double t = t_start; t <= t_end; t += dt) {
        double x_h_value = x_h(t, C1, C2, omega0);
        double x_p_value = x_p(t, A, B);
        double x_total = x_h_value + x_p_value;
        double v_h_value = v_h(t, C1, C2, omega0);
        double v_p_value = v_p(t, A, B);
        double v_total = v_h_value + v_p_value;

        file << t << ", " << x_h_value << ", " << x_p_value << ", " << x_total << ", " << v_h_value << ", " << v_p_value << ", " << v_total << "\n";
    }

    file.close();
    return 0;
}
