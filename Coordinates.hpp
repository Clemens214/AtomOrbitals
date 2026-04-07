#ifndef COORDINATES_H
#define COORDINATES_H

#include <cmath>

struct cartesian {};
struct spherical {};

class Coordinates
{
    public:
        Coordinates () = default;
        Coordinates (double x, double y, double z): mX(x), mY(y), mZ(z) { updateSpherical(); };
        Coordinates (double x, double y, double z, cartesian): mX(x), mY(y), mZ(z) { updateSpherical(); };
        Coordinates (double r, double theta, double phi, spherical): mR(r), mTheta(theta), mPhi(phi) { updateCartesian(); };

        // return cartesian coordiantes
        double x () const { return mX; }
        double y () const { return mY; }
        double z () const { return mZ; }
        // return spherical coordinates
        double r () const { return mR; }
        double theta () const { return mTheta; }
        double phi () const { return mPhi; }

        // set cartesian coordiantes
        double& setX (const double X) { mX = X; updateSpherical (); return mX; }
        double& setY (const double Y) { mY = Y; updateSpherical (); return mY; }
        double& setZ (const double Z) { mZ = Z; updateSpherical (); return mZ; }
        // set spherical coordinates
        double& setR (const double R) { mR = R; updateSpherical (); return mR; }
        double& setTheta (const double Theta) { mTheta = Theta; updateSpherical (); return mTheta; }
        double& setPhi (const double Phi) { mPhi = Phi; updateSpherical (); return mPhi; }

        // convert the coordinates
        void updateSpherical ();
        void updateCartesian ();

    private:
        // cartesian coordinates
        double mX = 0;
        double mY = 0;
        double mZ = 0;
        // spherical coordinates
        double mR = 0;
        double mTheta = 0;
        double mPhi = 0;
};

void Coordinates::updateSpherical () {
    double r = sqrt( x()*x() + y()*y() + z()*z() );
    setR(r);
    double theta = acos( z() / r );
    setTheta(theta);
    auto sign = [](double i) { if ( i>= 0 ) return 1; else return -1; };
    double rho = sqrt( x()*x() + y()*y() );
    double phi = sign( y() ) * acos( x() / rho );
    setPhi(phi);
    return;
}

void Coordinates::updateCartesian () {
    double x = r() * sin( theta() ) * cos( phi() );
    setX(x);
    double y = r() * sin( theta() ) * sin( phi() );
    setY(y);
    double z = r() * cos( theta() );
    setZ(z);
    return;
}

#endif // COORDINATES_H