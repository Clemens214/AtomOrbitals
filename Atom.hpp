#ifndef ATOM_H
#define ATOM_H

#include <vector>

#include "Coordinates.hpp"
#include "Orbital.hpp"

class Atom
{
    struct Core {
        Coordinates location = Coordinates(0, 0, 0);
        int charge = 1;
    };

    public:
        Atom () = default;
        Atom (int c): charge(c) { core.charge = c; };

        
        
    private:
        std::vector<Orbital> orbitals;
        Core core;
        int charge = core.charge;
};

#endif // ATOM_H