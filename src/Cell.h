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

    void setg2mu2A(double in) { g2mu2A_ = in; }
    void setg2mu2B(double in) { g2mu2B_ = in; }
    double getg2mu2A() const { return g2mu2A_; }
    double getg2mu2B() const { return g2mu2B_; }

    void setTpA(double in) { TpA_ = in; }
    void setTpB(double in) { TpB_ = in; }
    double getTpA() const { return TpA_; }
    double getTpB() const { return TpB_; }

    void setEpsilon(double in) { epsilon_ = in; }
    double getEpsilon() const { return epsilon_; }

    void setTtautau(double in) { Ttautau_ = in; }
    double getTtautau() const { return Ttautau_; }
    void setTxx(double in) { Txx_ = in; }
    double getTxx() const { return Txx_; }
    void setTyy(double in) { Tyy_ = in; }
    double getTyy() const { return Tyy_; }
    void setTxy(double in) { Txy_ = in; }
    double getTxy() const { return Txy_; }
    void setTetaeta(double in) { Tetaeta_ = in; }
    double getTetaeta() const { return Tetaeta_; }
    void setTtaux(double in) { Ttaux_ = in; }
    double getTtaux() const { return Ttaux_; }
    void setTtauy(double in) { Ttauy_ = in; }
    double getTtauy() const { return Ttauy_; }
    void setTtaueta(double in) { Ttaueta_ = in; }
    double getTtaueta() const { return Ttaueta_; }
    void setTxeta(double in) { Txeta_ = in; }
    double getTxeta() const { return Txeta_; }
    void setTyeta(double in) { Tyeta_ = in; }
    double getTyeta() const { return Tyeta_; }

    void setpitautau(double in) { pitautau_ = in; }
    double getpitautau() const { return pitautau_; }
    void setpixx(double in) { pixx_ = in; }
    double getpixx() const { return pixx_; }
    void setpiyy(double in) { piyy_ = in; }
    double getpiyy() const { return piyy_; }
    void setpixy(double in) { pixy_ = in; }
    double getpixy() const { return pixy_; }
    void setpietaeta(double in) { pietaeta_ = in; }
    double getpietaeta() const { return pietaeta_; }
    void setpitaux(double in) { pitaux_ = in; }
    double getpitaux() const { return pitaux_; }
    void setpitauy(double in) { pitauy_ = in; }
    double getpitauy() const { return pitauy_; }
    void setpitaueta(double in) { pitaueta_ = in; }
    double getpitaueta() const { return pitaueta_; }
    void setpixeta(double in) { pixeta_ = in; }
    double getpixeta() const { return pixeta_; }
    void setpiyeta(double in) { piyeta_ = in; }
    double getpiyeta() const { return piyeta_; }

    void setutau(double in) { utau_ = in; }
    double getutau() const { return utau_; }
    void setux(double in) { ux_ = in; }
    double getux() const { return ux_; }
    void setuy(double in) { uy_ = in; }
    double getuy() const { return uy_; }
    void setueta(double in) { ueta_ = in; }
    double getueta() const { return ueta_; }
};

#endif  // SRC_CELL_H_
