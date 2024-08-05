#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

// Function to simulate the option price using Monte Carlo
double simulateOptionPrice(double S0, double K, double T, double r, double sigma, 
                           std::string optionType, int numSimulations = 10000) {
    std::mt19937 gen(42); // Standard mersenne_twister_engine seeded with 42
    std::normal_distribution<> d(0, 1); // Standard normal distribution

    double dt = T / numSimulations;
    std::vector<double> pricePaths(numSimulations, S0);

    for (int t = 1; t <= numSimulations; ++t) {
        for (int i = 0; i < numSimulations; ++i) {
            double Z = d(gen);
            pricePaths[i] *= std::exp((r - 0.5 * sigma * sigma) * dt + sigma * std::sqrt(dt) * Z);
        }
    }

    double payoffSum = 0.0;
    for (double price : pricePaths) {
        if (optionType == "call") {
            payoffSum += std::max(price - K, 0.0);
        } else if (optionType == "put") {
            payoffSum += std::max(K - price, 0.0);
        }
    }

    return std::exp(-r * T) * (payoffSum / numSimulations);
}

// Function to find the implied volatility using a bisection method
double findImpliedVolatility(double S0, double K, double T, double r, double marketPrice, 
                             std::string optionType, double tol = 1e-5, int maxIter = 100) {
    double lowVol = 0.01;
    double highVol = 3.0;
    double impliedVol = (lowVol + highVol) / 2.0;

    for (int iter = 0; iter < maxIter; ++iter) {
        double simulatedPrice = simulateOptionPrice(S0, K, T, r, impliedVol, optionType);
        
        if (std::abs(simulatedPrice - marketPrice) < tol) {
            return impliedVol;
        }

        if (simulatedPrice < marketPrice) {
            lowVol = impliedVol;
        } else {
            highVol = impliedVol;
        }

        impliedVol = (lowVol + highVol) / 2.0;
    }

    return impliedVol;
}

int main() {
    double S0 = 100;  // Initial stock price
    double K = 100;   // Strike price
    double T = 1;     // Time to maturity (in years)
    double r = 0.05;  // Risk-free rate
    double marketPrice = 10; // Market price of the option

    double implied_vol = findImpliedVolatility(S0, K, T, r, marketPrice, "call");
    std::cout << "Implied Vol: " << implied_vol << std::endl;

    return 0;
}
