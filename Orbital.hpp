#ifndef ORBITAL_H
#define ORBITAL_H

#define _USE_MATH_DEFINES

#include <cmath>
#include <complex>
#include <vector>

class Orbital
{
    public:
        Orbital () = default;
        Orbital (int n, int l, int m): principalNum(n), angularNum(l), magneticNum(m) {};
        Orbital (int c): chargeNum(c) {};
        Orbital (int n, int l, int m, int c): principalNum(n), angularNum(l), magneticNum(m), chargeNum(c) {};

        int n () const { return principalNum; }
        int l () const { return angularNum; }
        int m () const { return magneticNum; }
        int charge () const { return chargeNum; }

        int& setPrincipalNum (const int &val);
        int& setAngularNum (const int &val);
        int& setMagneticNum (const int &val);

        std::vector<double> radialFunc (std::vector<double> radii);
        std::vector<std::complex<double>> sphericalFunc (std::vector<double> thetas, std::vector<double> phis);

    private:
        int principalNum = 1;
        int angularNum = 0;
        int magneticNum = 0;

        int chargeNum = 1;
        const double bohrRadius = 5.29177210544e-11; //5.29177210544×10−11

        unsigned long long int factorial(const int n);
        double binomial (const double n,  const double k);

        double generalizedLaguerre(const int n, const int m, const double x);
        double generalizedLegendre(const int l, const int m, const double x);

};

int& Orbital::setPrincipalNum (const int &val) {
    if ( val <= 0 )
        throw std::runtime_error("Principal quantum number musst be positive!");
    principalNum = val;
    return principalNum;
}

int& Orbital::setAngularNum (const int &val) {
    if ( val < 0  )
        throw std::runtime_error("Angular quantum number must not be negative!");
    if ( val >= principalNum )
        throw std::runtime_error("Angular quantum number must be smaller than the principal quantum number!");
    angularNum = val;
    return angularNum;
}

int& Orbital::setMagneticNum (const int &val) {
    if ( abs(val) < angularNum  )
        throw std::runtime_error("Absolute value of the magnetic quantum number must be smaller than the angular quantum number!");
    magneticNum = val;
    return magneticNum;
}

// Factorial function
unsigned long long int Orbital::factorial(const int n) {
    if ( n < 0  )
        throw std::runtime_error("Cannot calculate the factorial of a negative number!");
    long int f = 1;
    for (long int i=1; i<=n; ++i)
        f *= i;
        return f;
}

// Binomial coefficient
double Orbital::binomial (const double n,  const double k) {
    if ( n != static_cast<int>(n) )
        return 0;
    if ( k != static_cast<int>(k) )
        return 0;
    double result = factorial(n)/( factorial(k) * factorial(n - k) );
    return result;
}

////////////////////////////////////////////////////////////////////////////

// Radial Function
std::vector<double> Orbital::radialFunc (std::vector<double> radii) 
{
    std::vector<double> results(radii.size(), 0);
    const double factor = ( 2*charge() )/( n()*bohrRadius );
    const double preFactor = sqrt( pow(factor, 3) * factorial( n()-l()-1 )/( 2*n()*factorial( n()+l() ) ));
    for (decltype(radii.size()) i = 0; i < radii.size(); ++i) {
        const double expFactor = exp( -1*( factor/2 ) * radii[i] );
        const double powFactor = pow(factor*radii[i], l());
        const double laguerre = generalizedLaguerre(n()-l()-1, 2*l()+1, radii[i]);
        results[i] = preFactor * expFactor * powFactor * laguerre;
    }
    return results;
}

// Generalized Laguerre polynomial function L_n^m (x)
double Orbital::generalizedLaguerre(const int n, const int m, const double x) 
{
    double result = 0.0;
    for (int i = 0; i <= n; ++i) {
        auto parityVal = [](int m) { if ( m%2 == 0 ) return 1; else return -1; };
        double binomialCoefficient = binomial(n + m, n - i);
        double value = pow(x, i) / factorial(i);
        result += parityVal(i) * binomialCoefficient * value;
    }
    return result;
}

// Spherical Harmonics
std::vector<std::complex<double>> Orbital::sphericalFunc (std::vector<double> thetas, std::vector<double> phis) 
{
    std::vector<std::complex<double>> results(thetas.size(), 0);
    const double factor1 = ( 2*l() + 1 )/( 4*M_PI );
    const double factor2 = factorial( l() - m() )/factorial( l() + m() );
    const double preFactor = sqrt( factor1 * factor2 );
    for (decltype(thetas.size()) i = 0; i < thetas.size(); ++i) {
        const double legendre = generalizedLegendre(l(), m(), cos(thetas[i]) );
        std::complex<double> phase(0.0, m()*phis[i] );
        const std::complex<double> expFactor = exp( phase );
        results[i] = preFactor * legendre * expFactor;
    }
    return results;
}

// Generalized Legendre polynomial function P_l^m (x)
double Orbital::generalizedLegendre(const int l, const int m, const double x) 
{
    // handle case if m is outside of range -l to l
    if ( m > l || m< -1*l )
        return 0;
    // calculate the part of preFactor depenting of the sign of m
    int mEff = m;
    double preFactor;
    if ( m < 0 ) {
        mEff = -1*m;
        preFactor = factorial(l - mEff) / factorial(l + mEff);
    } else {
        auto parityVal = [](int m) { if ( m%2 == 0 ) return 1; else return -1; };
        preFactor = parityVal(mEff);
    }
    // calculate the rest of the preFactor
    double twoFactor = pow(2, l);
    double xFactor = pow( (1 - x*x), mEff/2 );
    double factor = preFactor * twoFactor * xFactor;
    // calculate the sum
    double result = 0.0;
    for (int i = mEff; i <= l; ++i) {
        double factorialFactor = factorial(i) / factorial(i - mEff);
        double binomialCoefficient = binomial(l, i);;
        double binomialFactor = binomial((l + i - 1)/2, l);;
        double value = pow(x, i - m);
        result += factorialFactor * binomialCoefficient * binomialFactor * value;
    }
    return factor * result;
}

#endif // ORBITAL_H