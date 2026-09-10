#include "Group.h"

#include <complex>
#include <iostream>

using std::complex;
using std::cout;
using std::endl;

Group::Group() {
    cout << "Initializing group SU(3) ... ";
    // fundamental rep.
    t_[0].set(0, 0, 0.);
    t_[0].set(0, 1, 0.5);
    t_[0].set(0, 2, 0.);
    t_[0].set(1, 0, 0.5);
    t_[0].set(1, 1, 0.);
    t_[0].set(1, 2, 0.);
    t_[0].set(2, 0, 0.);
    t_[0].set(2, 1, 0.);
    t_[0].set(2, 2, 0.);

    t_[1].set(0, 0, 0.);
    t_[1].set(0, 1, complex<double>(0., -0.5));
    t_[1].set(0, 2, 0.);
    t_[1].set(1, 0, complex<double>(0., 0.5));
    t_[1].set(1, 1, 0.);
    t_[1].set(1, 2, 0.);
    t_[1].set(2, 0, 0.);
    t_[1].set(2, 1, 0.);
    t_[1].set(2, 2, 0.);

    t_[2].set(0, 0, 0.5);
    t_[2].set(0, 1, 0.);
    t_[2].set(0, 2, 0.);
    t_[2].set(1, 0, 0.);
    t_[2].set(1, 1, -0.5);
    t_[2].set(1, 2, 0.);
    t_[2].set(2, 0, 0.);
    t_[2].set(2, 1, 0.);
    t_[2].set(2, 2, 0.);

    t_[3].set(0, 0, 0.);
    t_[3].set(0, 1, 0.);
    t_[3].set(0, 2, 0.5);
    t_[3].set(1, 0, 0.);
    t_[3].set(1, 1, 0.);
    t_[3].set(1, 2, 0.);
    t_[3].set(2, 0, 0.5);
    t_[3].set(2, 1, 0.);
    t_[3].set(2, 2, 0.);

    t_[4].set(0, 0, 0.);
    t_[4].set(0, 1, 0.);
    t_[4].set(0, 2, complex<double>(0., -0.5));
    t_[4].set(1, 0, 0.);
    t_[4].set(1, 1, 0.);
    t_[4].set(1, 2, 0.);
    t_[4].set(2, 0, complex<double>(0., 0.5));
    t_[4].set(2, 1, 0.);
    t_[4].set(2, 2, 0.);

    t_[5].set(0, 0, 0.);
    t_[5].set(0, 1, 0.);
    t_[5].set(0, 2, 0.);
    t_[5].set(1, 0, 0.);
    t_[5].set(1, 1, 0.);
    t_[5].set(1, 2, 0.5);
    t_[5].set(2, 0, 0.);
    t_[5].set(2, 1, 0.5);
    t_[5].set(2, 2, 0.);

    t_[6].set(0, 0, 0.);
    t_[6].set(0, 1, 0.);
    t_[6].set(0, 2, 0.);
    t_[6].set(1, 0, 0.);
    t_[6].set(1, 1, 0.);
    t_[6].set(1, 2, complex<double>(0., -0.5));
    t_[6].set(2, 0, 0.);
    t_[6].set(2, 1, complex<double>(0., 0.5));
    t_[6].set(2, 2, 0.);

    t_[7].set(0, 0, 1. / (2. * sqrt(3.)));
    t_[7].set(0, 1, 0.);
    t_[7].set(0, 2, 0.);
    t_[7].set(1, 0, 0.);
    t_[7].set(1, 1, 1. / (2. * sqrt(3.)));
    t_[7].set(1, 2, 0.);
    t_[7].set(2, 0, 0.);
    t_[7].set(2, 1, 0.);
    t_[7].set(2, 2, -1. / (sqrt(3.)));
    cout << "done." << endl;
}
