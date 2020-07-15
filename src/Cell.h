#ifndef Cell_h
#define Cell_h

#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "Matrix.h"

using namespace std;

class Cell {
private:
  int Nc;
  double epsilon; // energy density after collision

  // nucleus A
  double g2mu2A; // color charge density of nucleus A
  double TpA;    // sum over the proton T(b) in this cell for nucleus A
  Matrix *U;     // U is in the fundamental rep. (Nc*Nc matrix) // duobles as x
                 // component of electric field

  // nucleus B
  double g2mu2B; // color charge density of nucleus B
  double TpB;    // sum over the proton T(b) in this cell for nucleus B
  Matrix *U2; // Ui is the initial U in the fundamental rep. (Nc*Nc matrix) //
              // doubles as y component of electric field

  Matrix *Ux; // U is in the fundamental rep. (Nc*Nc matrix)
  Matrix *Uy; // U is in the fundamental rep. (Nc*Nc matrix)

  Matrix *Ux1; // U is in the fundamental rep. (Nc*Nc matrix) nucleus 1 (also
               // room to save g, the gauge fixing matrix)
  Matrix *Uy1; // U is in the fundamental rep. (Nc*Nc matrix) nucleus 1 (also
               // room to save Uplaq, the plaquette)

  Matrix *Ux2; // U is in the fundamental rep. (Nc*Nc matrix) nucleus 2 (doubles
               // as longitudinal electric field pi)
  Matrix *Uy2; // U is in the fundamental rep. (Nc*Nc matrix) nucleus 2 (doubles
               // as scalar field (longitudinal) )

  //  bool parity; // Parity of the cell (needed for Gauge fixing)

  double Ttautau; // energy momentum tensor
  double Txx;     // energy momentum tensor
  double Tyy;     // energy momentum tensor
  double Txy;     // energy momentum tensor
  double Tetaeta; // energy momentum tensor
  double Ttaux;   // energy momentum tensor
  double Ttauy;   // energy momentum tensor
  double Ttaueta; // energy momentum tensor
  double Txeta;   // energy momentum tensor
  double Tyeta;   // energy momentum tensor

  double pitautau; // energy momentum tensor
  double pixx;     // energy momentum tensor
  double piyy;     // energy momentum tensor
  double pixy;     // energy momentum tensor
  double pietaeta; // energy momentum tensor
  double pitaux;   // energy momentum tensor
  double pitauy;   // energy momentum tensor
  double pitaueta; // energy momentum tensor
  double pixeta;   // energy momentum tensor
  double piyeta;   // energy momentum tensor

  double utau; // flow velocity
  double ux;   // flow velocity
  double uy;   // flow velocity
  double ueta; // flow velocity

public:
  Cell(int N);
  ~Cell();

  //  void setParity(bool in) { parity = in; };
  //  bool getParity() { return parity; };

  void setg2mu2A(double in) { g2mu2A = in; };
  void setg2mu2B(double in) { g2mu2B = in; };

  double getg2mu2A() { return g2mu2A; };
  double getg2mu2B() { return g2mu2B; };

  void setTpA(double in) { TpA = in; };
  void setTpB(double in) { TpB = in; };

  double getTpA() { return TpA; };
  double getTpB() { return TpB; };

  void setU(const Matrix &x) { *U = x; };
  void setU2(const Matrix &x) { *U2 = x; };
  void setUplaq(const Matrix &x) {
    *Uy1 = x;
  }; // using unused Uy1 to store Uplaq

  void setUx(const Matrix &x) { *Ux = x; };
  void setUy(const Matrix &x) { *Uy = x; };
  void setUx1(const Matrix &x) { *Ux1 = x; };
  void setUy1(const Matrix &x) { *Uy1 = x; };
  void setUx2(const Matrix &x) { *Ux2 = x; };
  void setUy2(const Matrix &x) { *Uy2 = x; };

  void setEpsilon(const double in) { epsilon = in; };
  double getEpsilon() { return epsilon; };

  void setTtautau(const double in) { Ttautau = in; };
  double getTtautau() { return Ttautau; };
  void setTxx(const double in) { Txx = in; };
  double getTxx() { return Txx; };
  void setTyy(const double in) { Tyy = in; };
  double getTyy() { return Tyy; };
  void setTxy(const double in) { Txy = in; };
  double getTxy() { return Txy; };
  void setTetaeta(const double in) { Tetaeta = in; };
  double getTetaeta() { return Tetaeta; };
  void setTtaux(const double in) { Ttaux = in; };
  double getTtaux() { return Ttaux; };
  void setTtauy(const double in) { Ttauy = in; };
  double getTtauy() { return Ttauy; };
  void setTtaueta(const double in) { Ttaueta = in; };
  double getTtaueta() { return Ttaueta; };
  void setTxeta(const double in) { Txeta = in; };
  double getTxeta() { return Txeta; };
  void setTyeta(const double in) { Tyeta = in; };
  double getTyeta() { return Tyeta; };

  void setpitautau(const double in) { pitautau = in; };
  double getpitautau() { return pitautau; };
  void setpixx(const double in) { pixx = in; };
  double getpixx() { return pixx; };
  void setpiyy(const double in) { piyy = in; };
  double getpiyy() { return piyy; };
  void setpixy(const double in) { pixy = in; };
  double getpixy() { return pixy; };
  void setpietaeta(const double in) { pietaeta = in; };
  double getpietaeta() { return pietaeta; };
  void setpitaux(const double in) { pitaux = in; };
  double getpitaux() { return pitaux; };
  void setpitauy(const double in) { pitauy = in; };
  double getpitauy() { return pitauy; };
  void setpitaueta(const double in) { pitaueta = in; };
  double getpitaueta() { return pitaueta; };
  void setpixeta(const double in) { pixeta = in; };
  double getpixeta() { return pixeta; };
  void setpiyeta(const double in) { piyeta = in; };
  double getpiyeta() { return piyeta; };

  void setutau(const double in) { utau = in; };
  double getutau() { return utau; };
  void setux(const double in) { ux = in; };
  double getux() { return ux; };
  void setuy(const double in) { uy = in; };
  double getuy() { return uy; };
  void setueta(const double in) { ueta = in; };
  double getueta() { return ueta; };

  Matrix &getg() const { return *Ux1; }; // use unused Ux1 to store g
  Matrix &getU() const { return *U; };
  Matrix &getUx() const { return *Ux; };
  Matrix &getUy() const { return *Uy; };
  Matrix &getU2() const { return *U2; };
  Matrix &getUx1() const { return *Ux1; };
  Matrix &getUy1() const { return *Uy1; };
  Matrix &getUx2() const { return *Ux2; };
  Matrix &getUy2() const { return *Uy2; };
  Matrix &getUplaq() const { return *Uy1; }; // use unused Uy1 to store Uplaq

  void setE1(const Matrix &x) { *U = x; }; // use unused U to store E1
  Matrix &getE1() const { return *U; };
  void setE2(const Matrix &x) { *U2 = x; }; // use unused U2 to store E2
  Matrix &getE2() const { return *U2; };
  void setphi(const Matrix &x) { *Uy2 = x; }; // use unused Uy2 to store phi
  Matrix &getphi() const { return *Uy2; };
  void setpi(const Matrix &x) { *Ux2 = x; }; // use unused Ux2 to store pi
  Matrix &getpi() const { return *Ux2; };
  void setg(const Matrix &x) { *Ux1 = x; }; // using unused Ux1 to store g

  //  void computeAdjointU();
};

class SmallCell {
private:
  int Nc;
  Matrix *buffer1;
  Matrix *buffer2;

public:
  SmallCell(int N);
  ~SmallCell();

  Matrix &getbuffer1() const { return *buffer1; };
  void setbuffer1(const Matrix &x) { *buffer1 = x; };
  Matrix &getbuffer2() const { return *buffer2; };
  void setbuffer2(const Matrix &x) { *buffer2 = x; };
};

#endif
