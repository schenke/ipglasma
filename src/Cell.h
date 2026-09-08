#ifndef Cell_h
#define Cell_h

// Scalar per-site state. Fundamental SU(3) matrices are intentionally not
// stored here; they live in contiguous structure-of-arrays fields on Lattice.
class Cell {
  private:
    double epsilon;

    double g2mu2A;
    double TpA;
    double g2mu2B;
    double TpB;

    double Ttautau;
    double Txx;
    double Tyy;
    double Txy;
    double Tetaeta;
    double Ttaux;
    double Ttauy;
    double Ttaueta;
    double Txeta;
    double Tyeta;

    double pitautau;
    double pixx;
    double piyy;
    double pixy;
    double pietaeta;
    double pitaux;
    double pitauy;
    double pitaueta;
    double pixeta;
    double piyeta;

    double utau;
    double ux;
    double uy;
    double ueta;

  public:
    explicit Cell(int Nc);
    ~Cell() = default;

    void setg2mu2A(double in) { g2mu2A = in; }
    void setg2mu2B(double in) { g2mu2B = in; }
    double getg2mu2A() const { return g2mu2A; }
    double getg2mu2B() const { return g2mu2B; }

    void setTpA(double in) { TpA = in; }
    void setTpB(double in) { TpB = in; }
    double getTpA() const { return TpA; }
    double getTpB() const { return TpB; }

    void setEpsilon(double in) { epsilon = in; }
    double getEpsilon() const { return epsilon; }

    void setTtautau(double in) { Ttautau = in; }
    double getTtautau() const { return Ttautau; }
    void setTxx(double in) { Txx = in; }
    double getTxx() const { return Txx; }
    void setTyy(double in) { Tyy = in; }
    double getTyy() const { return Tyy; }
    void setTxy(double in) { Txy = in; }
    double getTxy() const { return Txy; }
    void setTetaeta(double in) { Tetaeta = in; }
    double getTetaeta() const { return Tetaeta; }
    void setTtaux(double in) { Ttaux = in; }
    double getTtaux() const { return Ttaux; }
    void setTtauy(double in) { Ttauy = in; }
    double getTtauy() const { return Ttauy; }
    void setTtaueta(double in) { Ttaueta = in; }
    double getTtaueta() const { return Ttaueta; }
    void setTxeta(double in) { Txeta = in; }
    double getTxeta() const { return Txeta; }
    void setTyeta(double in) { Tyeta = in; }
    double getTyeta() const { return Tyeta; }

    void setpitautau(double in) { pitautau = in; }
    double getpitautau() const { return pitautau; }
    void setpixx(double in) { pixx = in; }
    double getpixx() const { return pixx; }
    void setpiyy(double in) { piyy = in; }
    double getpiyy() const { return piyy; }
    void setpixy(double in) { pixy = in; }
    double getpixy() const { return pixy; }
    void setpietaeta(double in) { pietaeta = in; }
    double getpietaeta() const { return pietaeta; }
    void setpitaux(double in) { pitaux = in; }
    double getpitaux() const { return pitaux; }
    void setpitauy(double in) { pitauy = in; }
    double getpitauy() const { return pitauy; }
    void setpitaueta(double in) { pitaueta = in; }
    double getpitaueta() const { return pitaueta; }
    void setpixeta(double in) { pixeta = in; }
    double getpixeta() const { return pixeta; }
    void setpiyeta(double in) { piyeta = in; }
    double getpiyeta() const { return piyeta; }

    void setutau(double in) { utau = in; }
    double getutau() const { return utau; }
    void setux(double in) { ux = in; }
    double getux() const { return ux; }
    void setuy(double in) { uy = in; }
    double getuy() const { return uy; }
    void setueta(double in) { ueta = in; }
    double getueta() const { return ueta; }
};

#endif
