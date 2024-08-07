#include <iostream>
#include <fstream>
#include <cmath>

double analytical_integral() {
    return sqrt(M_PI);
}

int main() {
    double result = analytical_integral();
    std::ofstream file("analytical_result.txt");
    file << "Analytical Integration Result: " << result << std::endl;
    file.close();

    return 0;
}
