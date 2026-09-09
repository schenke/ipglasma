#include "Cell.h"

#include <cstdlib>
#include <iostream>

Cell::Cell(int Nc)
    : epsilon(0.),
      g2mu2A(0.),
      TpA(0.),
      g2mu2B(0.),
      TpB(0.),
      Ttautau(0.),
      Txx(0.),
      Tyy(0.),
      Txy(0.),
      Tetaeta(0.),
      Ttaux(0.),
      Ttauy(0.),
      Ttaueta(0.),
      Txeta(0.),
      Tyeta(0.),
      pitautau(0.),
      pixx(0.),
      piyy(0.),
      pixy(0.),
      pietaeta(0.),
      pitaux(0.),
      pitauy(0.),
      pitaueta(0.),
      pixeta(0.),
      piyeta(0.),
      utau(0.),
      ux(0.),
      uy(0.),
      ueta(0.) {
    if (Nc != 3) {
        std::cerr << "Error: Cell storage is SU(3)-only; received Nc=" << Nc
                  << ". Exiting." << std::endl;
        std::exit(1);
    }
}
