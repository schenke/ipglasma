#include "Group.h"

#include <complex>
#include <cstdlib>
#include <iostream>

using std::complex;
using std::cout;
using std::endl;

Group::Group(int N) {
    if (N != 3) {
        std::cerr << "Error: IP-Glasma Group is SU(3)-only; received SU(" << N
                  << "). Exiting." << endl;
        std::exit(1);
    }

    cout << "Initializing group SU(3) ... ";
    // fundamental rep.
    t[0].set(0, 0, 0.);
    t[0].set(0, 1, 0.5);
    t[0].set(0, 2, 0.);
    t[0].set(1, 0, 0.5);
    t[0].set(1, 1, 0.);
    t[0].set(1, 2, 0.);
    t[0].set(2, 0, 0.);
    t[0].set(2, 1, 0.);
    t[0].set(2, 2, 0.);

    t[1].set(0, 0, 0.);
    t[1].set(0, 1, complex<double>(0., -0.5));
    t[1].set(0, 2, 0.);
    t[1].set(1, 0, complex<double>(0., 0.5));
    t[1].set(1, 1, 0.);
    t[1].set(1, 2, 0.);
    t[1].set(2, 0, 0.);
    t[1].set(2, 1, 0.);
    t[1].set(2, 2, 0.);

    t[2].set(0, 0, 0.5);
    t[2].set(0, 1, 0.);
    t[2].set(0, 2, 0.);
    t[2].set(1, 0, 0.);
    t[2].set(1, 1, -0.5);
    t[2].set(1, 2, 0.);
    t[2].set(2, 0, 0.);
    t[2].set(2, 1, 0.);
    t[2].set(2, 2, 0.);

    t[3].set(0, 0, 0.);
    t[3].set(0, 1, 0.);
    t[3].set(0, 2, 0.5);
    t[3].set(1, 0, 0.);
    t[3].set(1, 1, 0.);
    t[3].set(1, 2, 0.);
    t[3].set(2, 0, 0.5);
    t[3].set(2, 1, 0.);
    t[3].set(2, 2, 0.);

    t[4].set(0, 0, 0.);
    t[4].set(0, 1, 0.);
    t[4].set(0, 2, complex<double>(0., -0.5));
    t[4].set(1, 0, 0.);
    t[4].set(1, 1, 0.);
    t[4].set(1, 2, 0.);
    t[4].set(2, 0, complex<double>(0., 0.5));
    t[4].set(2, 1, 0.);
    t[4].set(2, 2, 0.);

    t[5].set(0, 0, 0.);
    t[5].set(0, 1, 0.);
    t[5].set(0, 2, 0.);
    t[5].set(1, 0, 0.);
    t[5].set(1, 1, 0.);
    t[5].set(1, 2, 0.5);
    t[5].set(2, 0, 0.);
    t[5].set(2, 1, 0.5);
    t[5].set(2, 2, 0.);

    t[6].set(0, 0, 0.);
    t[6].set(0, 1, 0.);
    t[6].set(0, 2, 0.);
    t[6].set(1, 0, 0.);
    t[6].set(1, 1, 0.);
    t[6].set(1, 2, complex<double>(0., -0.5));
    t[6].set(2, 0, 0.);
    t[6].set(2, 1, complex<double>(0., 0.5));
    t[6].set(2, 2, 0.);

    t[7].set(0, 0, 1. / (2. * sqrt(3.)));
    t[7].set(0, 1, 0.);
    t[7].set(0, 2, 0.);
    t[7].set(1, 0, 0.);
    t[7].set(1, 1, 1. / (2. * sqrt(3.)));
    t[7].set(1, 2, 0.);
    t[7].set(2, 0, 0.);
    t[7].set(2, 1, 0.);
    t[7].set(2, 2, -1. / (sqrt(3.)));
    cout << "done." << endl;
}
