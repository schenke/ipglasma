#ifndef SRC_CELL_H_
#define SRC_CELL_H_

// Scalar per-site state. Fundamental SU(3) matrices are intentionally not
// stored here; they live in contiguous structure-of-arrays fields on Lattice.
class Cell {
  private:
    double epsilon_;

    double g2mu2A_;
    double TpA_;
    double g2mu2B_;
    double TpB_;

    double Ttautau_;
    double Txx_;
    double Tyy_;
    double Txy_;
    double Tetaeta_;
    double Ttaux_;
    double Ttauy_;
    double Ttaueta_;
    double Txeta_;
    double Tyeta_;

    double pitautau_;
    double pixx_;
    double piyy_;
    double pixy_;
    double pietaeta_;
    double pitaux_;
    double pitauy_;
    double pitaueta_;
    double pixeta_;
    double piyeta_;

    double utau_;
    double ux_;
    double uy_;
    double ueta_;

  public:
    Cell();
    ~Cell() = default;

    void setg2mu2A(double x) { g2mu2A_ = x; }
    void setg2mu2B(double x) { g2mu2B_ = x; }
    double getg2mu2A() const { return g2mu2A_; }
    double getg2mu2B() const { return g2mu2B_; }

    void setTpA(double x) { TpA_ = x; }
    void setTpB(double x) { TpB_ = x; }
    double getTpA() const { return TpA_; }
    double getTpB() const { return TpB_; }

    void setEpsilon(double x) { epsilon_ = x; }
    double getEpsilon() const { return epsilon_; }

    void setTtautau(double x) { Ttautau_ = x; }
    double getTtautau() const { return Ttautau_; }
    void setTxx(double x) { Txx_ = x; }
    double getTxx() const { return Txx_; }
    void setTyy(double x) { Tyy_ = x; }
    double getTyy() const { return Tyy_; }
    void setTxy(double x) { Txy_ = x; }
    double getTxy() const { return Txy_; }
    void setTetaeta(double x) { Tetaeta_ = x; }
    double getTetaeta() const { return Tetaeta_; }
    void setTtaux(double x) { Ttaux_ = x; }
    double getTtaux() const { return Ttaux_; }
    void setTtauy(double x) { Ttauy_ = x; }
    double getTtauy() const { return Ttauy_; }
    void setTtaueta(double x) { Ttaueta_ = x; }
    double getTtaueta() const { return Ttaueta_; }
    void setTxeta(double x) { Txeta_ = x; }
    double getTxeta() const { return Txeta_; }
    void setTyeta(double x) { Tyeta_ = x; }
    double getTyeta() const { return Tyeta_; }

    void setpitautau(double x) { pitautau_ = x; }
    double getpitautau() const { return pitautau_; }
    void setpixx(double x) { pixx_ = x; }
    double getpixx() const { return pixx_; }
    void setpiyy(double x) { piyy_ = x; }
    double getpiyy() const { return piyy_; }
    void setpixy(double x) { pixy_ = x; }
    double getpixy() const { return pixy_; }
    void setpietaeta(double x) { pietaeta_ = x; }
    double getpietaeta() const { return pietaeta_; }
    void setpitaux(double x) { pitaux_ = x; }
    double getpitaux() const { return pitaux_; }
    void setpitauy(double x) { pitauy_ = x; }
    double getpitauy() const { return pitauy_; }
    void setpitaueta(double x) { pitaueta_ = x; }
    double getpitaueta() const { return pitaueta_; }
    void setpixeta(double x) { pixeta_ = x; }
    double getpixeta() const { return pixeta_; }
    void setpiyeta(double x) { piyeta_ = x; }
    double getpiyeta() const { return piyeta_; }

    void setutau(double x) { utau_ = x; }
    double getutau() const { return utau_; }
    void setux(double x) { ux_ = x; }
    double getux() const { return ux_; }
    void setuy(double x) { uy_ = x; }
    double getuy() const { return uy_; }
    void setueta(double x) { ueta_ = x; }
    double getueta() const { return ueta_; }
};

#endif  // SRC_CELL_H_
