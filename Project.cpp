#include <iostream>
#include <cctype>

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

int main()
{
    return 0;
}