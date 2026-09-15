#ifndef SRC_CELL_H_
#define SRC_CELL_H_

/**
 * Scalar per-site hydrodynamic and Yang-Mills observables.
 *
 * One Cell exists per transverse lattice site. The fundamental SU(3)
 * matrix fields are intentionally not stored here; they live in
 * Lattice's structure-of-arrays members instead. Cell's own fields are
 * populated in three stages of a run:
 * - \f$g^2\mu_A^2\f$/\f$g^2\mu_B^2\f$ and \f$T_p^A\f$/\f$T_p^B\f$: the
 *   initial-state color-charge density and nuclear thickness for the
 *   projectile/target, set once during initialization (see Init.cpp)
 *   and used to sample the classical fields and the local running
 *   coupling/saturation scale.
 * - The \f$T^{\mu\nu}\f$ components: the local energy-momentum tensor
 *   in Milne (\f$\tau, x, y, \eta\f$) coordinates, computed from the
 *   classical Yang-Mills fields by Evolution::tmunu().
 * - \f$\epsilon\f$, the \f$u^\mu\f$ components, and the
 *   \f$\pi^{\mu\nu}\f$ components: the local rest-frame energy density,
 *   fluid four-velocity, and (traceless) shear-stress tensor, obtained
 *   by diagonalizing \f$T^{\mu\nu}\f$ in
 *   MyEigen::solveFlowVelocityAtCell() and written out as the
 *   hydrodynamic initial condition.
 */
class Cell {
  private:
    /// Local rest-frame energy density \f$\epsilon\f$ [1/fm^4].
    double epsilon_;

    /// Local \f$g^2\mu_A^2\f$ color-charge-density-squared value for the
    /// projectile (nucleus A) [lattice units].
    double g2mu2A_;
    /// Local nuclear thickness function \f$T_p^A\f$ for the projectile
    /// (nucleus A) [1/fm^2].
    double TpA_;
    /// Local \f$g^2\mu_B^2\f$ color-charge-density-squared value for the
    /// target (nucleus B) [lattice units].
    double g2mu2B_;
    /// Local nuclear thickness function \f$T_p^B\f$ for the target
    /// (nucleus B) [1/fm^2].
    double TpB_;

    /// \f$T^{\tau\tau}\f$ component of the energy-momentum tensor
    /// [1/fm^4].
    double Ttautau_;
    /// \f$T^{xx}\f$ component of the energy-momentum tensor [1/fm^4].
    double Txx_;
    /// \f$T^{yy}\f$ component of the energy-momentum tensor [1/fm^4].
    double Tyy_;
    /// \f$T^{xy}\f$ component of the energy-momentum tensor [1/fm^4].
    double Txy_;
    /// \f$T^{\eta\eta}\f$ component of the energy-momentum tensor
    /// [1/fm^4].
    double Tetaeta_;
    /// \f$T^{\tau x}\f$ component of the energy-momentum tensor
    /// [1/fm^4].
    double Ttaux_;
    /// \f$T^{\tau y}\f$ component of the energy-momentum tensor
    /// [1/fm^4].
    double Ttauy_;
    /// \f$T^{\tau\eta}\f$ component of the energy-momentum tensor
    /// [1/fm^4].
    double Ttaueta_;
    /// \f$T^{x\eta}\f$ component of the energy-momentum tensor
    /// [1/fm^4].
    double Txeta_;
    /// \f$T^{y\eta}\f$ component of the energy-momentum tensor
    /// [1/fm^4].
    double Tyeta_;

    /// \f$\pi^{\tau\tau}\f$ component of the local-rest-frame shear-
    /// stress tensor [1/fm^4].
    double pitautau_;
    /// \f$\pi^{xx}\f$ component of the local-rest-frame shear-stress
    /// tensor [1/fm^4].
    double pixx_;
    /// \f$\pi^{yy}\f$ component of the local-rest-frame shear-stress
    /// tensor [1/fm^4].
    double piyy_;
    /// \f$\pi^{xy}\f$ component of the local-rest-frame shear-stress
    /// tensor [1/fm^4].
    double pixy_;
    /// \f$\pi^{\eta\eta}\f$ component of the local-rest-frame shear-
    /// stress tensor [1/fm^4].
    double pietaeta_;
    /// \f$\pi^{\tau x}\f$ component of the local-rest-frame shear-stress
    /// tensor [1/fm^4].
    double pitaux_;
    /// \f$\pi^{\tau y}\f$ component of the local-rest-frame shear-stress
    /// tensor [1/fm^4].
    double pitauy_;
    /// \f$\pi^{\tau\eta}\f$ component of the local-rest-frame shear-
    /// stress tensor [1/fm^4].
    double pitaueta_;
    /// \f$\pi^{x\eta}\f$ component of the local-rest-frame shear-stress
    /// tensor [1/fm^4].
    double pixeta_;
    /// \f$\pi^{y\eta}\f$ component of the local-rest-frame shear-stress
    /// tensor [1/fm^4].
    double piyeta_;

    /// \f$u^\tau\f$ component of the local fluid four-velocity
    /// [dimensionless].
    double utau_;
    /// \f$u^x\f$ component of the local fluid four-velocity
    /// [dimensionless].
    double ux_;
    /// \f$u^y\f$ component of the local fluid four-velocity
    /// [dimensionless].
    double uy_;
    /// \f$u^\eta\f$ component of the local fluid four-velocity
    /// [dimensionless].
    double ueta_;

  public:
    /**
     * Constructs a Cell with every field zero-initialized.
     */
    Cell();
    ~Cell() = default;

    /**
     * Sets the projectile's (nucleus A) local \f$g^2\mu_A^2\f$ value.
     * \param[in] x The new \f$g^2\mu_A^2\f$ value [lattice units].
     */
    void setg2mu2A(double x) { g2mu2A_ = x; }
    /**
     * Sets the target's (nucleus B) local \f$g^2\mu_B^2\f$ value.
     * \param[in] x The new \f$g^2\mu_B^2\f$ value [lattice units].
     */
    void setg2mu2B(double x) { g2mu2B_ = x; }
    /**
     * Returns the projectile's (nucleus A) local \f$g^2\mu_A^2\f$ value,
     * used to sample the classical color source and to set the local
     * running coupling/saturation scale.
     * \return The stored \f$g^2\mu_A^2\f$ value [lattice units].
     */
    double getg2mu2A() const { return g2mu2A_; }
    /**
     * Returns the target's (nucleus B) local \f$g^2\mu_B^2\f$ value,
     * used to sample the classical color source and to set the local
     * running coupling/saturation scale.
     * \return The stored \f$g^2\mu_B^2\f$ value [lattice units].
     */
    double getg2mu2B() const { return g2mu2B_; }

    /**
     * Sets the projectile's (nucleus A) local nuclear thickness
     * function \f$T_p^A\f$.
     * \param[in] x The new \f$T_p^A\f$ value [1/fm^2].
     */
    void setTpA(double x) { TpA_ = x; }
    /**
     * Sets the target's (nucleus B) local nuclear thickness function
     * \f$T_p^B\f$.
     * \param[in] x The new \f$T_p^B\f$ value [1/fm^2].
     */
    void setTpB(double x) { TpB_ = x; }
    /**
     * Returns the projectile's (nucleus A) local nuclear thickness
     * function, used to look up the local saturation scale
     * \f$Q_s^2\f$.
     * \return The stored \f$T_p^A\f$ value [1/fm^2].
     */
    double getTpA() const { return TpA_; }
    /**
     * Returns the target's (nucleus B) local nuclear thickness
     * function, used to look up the local saturation scale
     * \f$Q_s^2\f$.
     * \return The stored \f$T_p^B\f$ value [1/fm^2].
     */
    double getTpB() const { return TpB_; }

    /**
     * Sets the local rest-frame energy density \f$\epsilon\f$.
     * \param[in] x The new energy density [1/fm^4].
     */
    void setEpsilon(double x) { epsilon_ = x; }
    /**
     * Returns the local rest-frame energy density, obtained either
     * directly from \f$T^{\tau\tau}\f$ (before the flow-velocity solve)
     * or from diagonalizing \f$T^{\mu\nu}\f$ in
     * MyEigen::solveFlowVelocityAtCell(). Multiply by \f$\hbar c\f$ to
     * convert to GeV/fm^3.
     * \return The stored energy density [1/fm^4].
     */
    double getEpsilon() const { return epsilon_; }

    /**
     * Sets the \f$T^{\tau\tau}\f$ component of the energy-momentum
     * tensor.
     * \param[in] x The new \f$T^{\tau\tau}\f$ value [1/fm^4].
     */
    void setTtautau(double x) { Ttautau_ = x; }
    /**
     * Returns the \f$T^{\tau\tau}\f$ component of the energy-momentum
     * tensor, set by Evolution::tmunu() from the classical Yang-Mills
     * fields.
     * \return The stored \f$T^{\tau\tau}\f$ value [1/fm^4].
     */
    double getTtautau() const { return Ttautau_; }
    /**
     * Sets the \f$T^{xx}\f$ component of the energy-momentum tensor.
     * \param[in] x The new \f$T^{xx}\f$ value [1/fm^4].
     */
    void setTxx(double x) { Txx_ = x; }
    /**
     * Returns the \f$T^{xx}\f$ component of the energy-momentum tensor,
     * set by Evolution::tmunu() from the classical Yang-Mills fields.
     * \return The stored \f$T^{xx}\f$ value [1/fm^4].
     */
    double getTxx() const { return Txx_; }
    /**
     * Sets the \f$T^{yy}\f$ component of the energy-momentum tensor.
     * \param[in] x The new \f$T^{yy}\f$ value [1/fm^4].
     */
    void setTyy(double x) { Tyy_ = x; }
    /**
     * Returns the \f$T^{yy}\f$ component of the energy-momentum tensor,
     * set by Evolution::tmunu() from the classical Yang-Mills fields.
     * \return The stored \f$T^{yy}\f$ value [1/fm^4].
     */
    double getTyy() const { return Tyy_; }
    /**
     * Sets the \f$T^{xy}\f$ (transverse shear) component of the
     * energy-momentum tensor.
     * \param[in] x The new \f$T^{xy}\f$ value [1/fm^4].
     */
    void setTxy(double x) { Txy_ = x; }
    /**
     * Returns the \f$T^{xy}\f$ (transverse shear) component of the
     * energy-momentum tensor, set by Evolution::tmunu() from the
     * classical Yang-Mills fields.
     * \return The stored \f$T^{xy}\f$ value [1/fm^4].
     */
    double getTxy() const { return Txy_; }
    /**
     * Sets the \f$T^{\eta\eta}\f$ (longitudinal) component of the
     * energy-momentum tensor.
     * \param[in] x The new \f$T^{\eta\eta}\f$ value [1/fm^4].
     */
    void setTetaeta(double x) { Tetaeta_ = x; }
    /**
     * Returns the \f$T^{\eta\eta}\f$ (longitudinal) component of the
     * energy-momentum tensor, set by Evolution::tmunu() from the
     * classical Yang-Mills fields.
     * \return The stored \f$T^{\eta\eta}\f$ value [1/fm^4].
     */
    double getTetaeta() const { return Tetaeta_; }
    /**
     * Sets the \f$T^{\tau x}\f$ (energy flux in x) component of the
     * energy-momentum tensor.
     * \param[in] x The new \f$T^{\tau x}\f$ value [1/fm^4].
     */
    void setTtaux(double x) { Ttaux_ = x; }
    /**
     * Returns the \f$T^{\tau x}\f$ (energy flux in x) component of the
     * energy-momentum tensor, set by Evolution::tmunu() from the
     * classical Yang-Mills fields.
     * \return The stored \f$T^{\tau x}\f$ value [1/fm^4].
     */
    double getTtaux() const { return Ttaux_; }
    /**
     * Sets the \f$T^{\tau y}\f$ (energy flux in y) component of the
     * energy-momentum tensor.
     * \param[in] x The new \f$T^{\tau y}\f$ value [1/fm^4].
     */
    void setTtauy(double x) { Ttauy_ = x; }
    /**
     * Returns the \f$T^{\tau y}\f$ (energy flux in y) component of the
     * energy-momentum tensor, set by Evolution::tmunu() from the
     * classical Yang-Mills fields.
     * \return The stored \f$T^{\tau y}\f$ value [1/fm^4].
     */
    double getTtauy() const { return Ttauy_; }
    /**
     * Sets the \f$T^{\tau\eta}\f$ (longitudinal energy flux) component
     * of the energy-momentum tensor.
     * \param[in] x The new \f$T^{\tau\eta}\f$ value [1/fm^4].
     */
    void setTtaueta(double x) { Ttaueta_ = x; }
    /**
     * Returns the \f$T^{\tau\eta}\f$ (longitudinal energy flux)
     * component of the energy-momentum tensor, set by
     * Evolution::tmunu() from the classical Yang-Mills fields.
     * \return The stored \f$T^{\tau\eta}\f$ value [1/fm^4].
     */
    double getTtaueta() const { return Ttaueta_; }
    /**
     * Sets the \f$T^{x\eta}\f$ (transverse-longitudinal shear)
     * component of the energy-momentum tensor.
     * \param[in] x The new \f$T^{x\eta}\f$ value [1/fm^4].
     */
    void setTxeta(double x) { Txeta_ = x; }
    /**
     * Returns the \f$T^{x\eta}\f$ (transverse-longitudinal shear)
     * component of the energy-momentum tensor, set by
     * Evolution::tmunu() from the classical Yang-Mills fields.
     * \return The stored \f$T^{x\eta}\f$ value [1/fm^4].
     */
    double getTxeta() const { return Txeta_; }
    /**
     * Sets the \f$T^{y\eta}\f$ (transverse-longitudinal shear)
     * component of the energy-momentum tensor.
     * \param[in] x The new \f$T^{y\eta}\f$ value [1/fm^4].
     */
    void setTyeta(double x) { Tyeta_ = x; }
    /**
     * Returns the \f$T^{y\eta}\f$ (transverse-longitudinal shear)
     * component of the energy-momentum tensor, set by
     * Evolution::tmunu() from the classical Yang-Mills fields.
     * \return The stored \f$T^{y\eta}\f$ value [1/fm^4].
     */
    double getTyeta() const { return Tyeta_; }

    /**
     * Sets the \f$\pi^{\tau\tau}\f$ component of the local-rest-frame
     * shear-stress tensor.
     * \param[in] x The new \f$\pi^{\tau\tau}\f$ value [1/fm^4].
     */
    void setpitautau(double x) { pitautau_ = x; }
    /**
     * Returns the \f$\pi^{\tau\tau}\f$ component of the local-rest-
     * frame shear-stress tensor: the traceless part of
     * \f$T^{\mu\nu}\f$ once boosted to the local rest frame via
     * \f$u^\mu\f$, computed by MyEigen::solveFlowVelocityAtCell() and
     * written out as the viscous-hydrodynamics initial condition.
     * \return The stored \f$\pi^{\tau\tau}\f$ value [1/fm^4].
     */
    double getpitautau() const { return pitautau_; }
    /**
     * Sets the \f$\pi^{xx}\f$ component of the local-rest-frame
     * shear-stress tensor.
     * \param[in] x The new \f$\pi^{xx}\f$ value [1/fm^4].
     */
    void setpixx(double x) { pixx_ = x; }
    /**
     * Returns the \f$\pi^{xx}\f$ component of the local-rest-frame
     * shear-stress tensor, computed by
     * MyEigen::solveFlowVelocityAtCell().
     * \return The stored \f$\pi^{xx}\f$ value [1/fm^4].
     */
    double getpixx() const { return pixx_; }
    /**
     * Sets the \f$\pi^{yy}\f$ component of the local-rest-frame
     * shear-stress tensor.
     * \param[in] x The new \f$\pi^{yy}\f$ value [1/fm^4].
     */
    void setpiyy(double x) { piyy_ = x; }
    /**
     * Returns the \f$\pi^{yy}\f$ component of the local-rest-frame
     * shear-stress tensor, computed by
     * MyEigen::solveFlowVelocityAtCell().
     * \return The stored \f$\pi^{yy}\f$ value [1/fm^4].
     */
    double getpiyy() const { return piyy_; }
    /**
     * Sets the \f$\pi^{xy}\f$ component of the local-rest-frame
     * shear-stress tensor.
     * \param[in] x The new \f$\pi^{xy}\f$ value [1/fm^4].
     */
    void setpixy(double x) { pixy_ = x; }
    /**
     * Returns the \f$\pi^{xy}\f$ component of the local-rest-frame
     * shear-stress tensor, computed by
     * MyEigen::solveFlowVelocityAtCell().
     * \return The stored \f$\pi^{xy}\f$ value [1/fm^4].
     */
    double getpixy() const { return pixy_; }
    /**
     * Sets the \f$\pi^{\eta\eta}\f$ component of the local-rest-frame
     * shear-stress tensor.
     * \param[in] x The new \f$\pi^{\eta\eta}\f$ value [1/fm^4].
     */
    void setpietaeta(double x) { pietaeta_ = x; }
    /**
     * Returns the \f$\pi^{\eta\eta}\f$ component of the local-rest-
     * frame shear-stress tensor, computed by
     * MyEigen::solveFlowVelocityAtCell().
     * \return The stored \f$\pi^{\eta\eta}\f$ value [1/fm^4].
     */
    double getpietaeta() const { return pietaeta_; }
    /**
     * Sets the \f$\pi^{\tau x}\f$ component of the local-rest-frame
     * shear-stress tensor.
     * \param[in] x The new \f$\pi^{\tau x}\f$ value [1/fm^4].
     */
    void setpitaux(double x) { pitaux_ = x; }
    /**
     * Returns the \f$\pi^{\tau x}\f$ component of the local-rest-frame
     * shear-stress tensor, computed by
     * MyEigen::solveFlowVelocityAtCell().
     * \return The stored \f$\pi^{\tau x}\f$ value [1/fm^4].
     */
    double getpitaux() const { return pitaux_; }
    /**
     * Sets the \f$\pi^{\tau y}\f$ component of the local-rest-frame
     * shear-stress tensor.
     * \param[in] x The new \f$\pi^{\tau y}\f$ value [1/fm^4].
     */
    void setpitauy(double x) { pitauy_ = x; }
    /**
     * Returns the \f$\pi^{\tau y}\f$ component of the local-rest-frame
     * shear-stress tensor, computed by
     * MyEigen::solveFlowVelocityAtCell().
     * \return The stored \f$\pi^{\tau y}\f$ value [1/fm^4].
     */
    double getpitauy() const { return pitauy_; }
    /**
     * Sets the \f$\pi^{\tau\eta}\f$ component of the local-rest-frame
     * shear-stress tensor.
     * \param[in] x The new \f$\pi^{\tau\eta}\f$ value [1/fm^4].
     */
    void setpitaueta(double x) { pitaueta_ = x; }
    /**
     * Returns the \f$\pi^{\tau\eta}\f$ component of the local-rest-
     * frame shear-stress tensor, computed by
     * MyEigen::solveFlowVelocityAtCell().
     * \return The stored \f$\pi^{\tau\eta}\f$ value [1/fm^4].
     */
    double getpitaueta() const { return pitaueta_; }
    /**
     * Sets the \f$\pi^{x\eta}\f$ component of the local-rest-frame
     * shear-stress tensor.
     * \param[in] x The new \f$\pi^{x\eta}\f$ value [1/fm^4].
     */
    void setpixeta(double x) { pixeta_ = x; }
    /**
     * Returns the \f$\pi^{x\eta}\f$ component of the local-rest-frame
     * shear-stress tensor, computed by
     * MyEigen::solveFlowVelocityAtCell().
     * \return The stored \f$\pi^{x\eta}\f$ value [1/fm^4].
     */
    double getpixeta() const { return pixeta_; }
    /**
     * Sets the \f$\pi^{y\eta}\f$ component of the local-rest-frame
     * shear-stress tensor.
     * \param[in] x The new \f$\pi^{y\eta}\f$ value [1/fm^4].
     */
    void setpiyeta(double x) { piyeta_ = x; }
    /**
     * Returns the \f$\pi^{y\eta}\f$ component of the local-rest-frame
     * shear-stress tensor, computed by
     * MyEigen::solveFlowVelocityAtCell().
     * \return The stored \f$\pi^{y\eta}\f$ value [1/fm^4].
     */
    double getpiyeta() const { return piyeta_; }

    /**
     * Sets the \f$u^\tau\f$ component of the local fluid four-velocity.
     * \param[in] x The new \f$u^\tau\f$ value [dimensionless].
     */
    void setutau(double x) { utau_ = x; }
    /**
     * Returns the \f$u^\tau\f$ component of the local fluid
     * four-velocity, obtained by diagonalizing \f$T^{\mu\nu}\f$ in
     * MyEigen::solveFlowVelocityAtCell().
     * \return The stored \f$u^\tau\f$ value [dimensionless].
     */
    double getutau() const { return utau_; }
    /**
     * Sets the \f$u^x\f$ component of the local fluid four-velocity.
     * \param[in] x The new \f$u^x\f$ value [dimensionless].
     */
    void setux(double x) { ux_ = x; }
    /**
     * Returns the \f$u^x\f$ component of the local fluid four-velocity,
     * obtained by diagonalizing \f$T^{\mu\nu}\f$ in
     * MyEigen::solveFlowVelocityAtCell().
     * \return The stored \f$u^x\f$ value [dimensionless].
     */
    double getux() const { return ux_; }
    /**
     * Sets the \f$u^y\f$ component of the local fluid four-velocity.
     * \param[in] x The new \f$u^y\f$ value [dimensionless].
     */
    void setuy(double x) { uy_ = x; }
    /**
     * Returns the \f$u^y\f$ component of the local fluid four-velocity,
     * obtained by diagonalizing \f$T^{\mu\nu}\f$ in
     * MyEigen::solveFlowVelocityAtCell().
     * \return The stored \f$u^y\f$ value [dimensionless].
     */
    double getuy() const { return uy_; }
    /**
     * Sets the \f$u^\eta\f$ component of the local fluid four-velocity.
     * \param[in] x The new \f$u^\eta\f$ value [dimensionless].
     */
    void setueta(double x) { ueta_ = x; }
    /**
     * Returns the \f$u^\eta\f$ component of the local fluid
     * four-velocity, obtained by diagonalizing \f$T^{\mu\nu}\f$ in
     * MyEigen::solveFlowVelocityAtCell().
     * \return The stored \f$u^\eta\f$ value [dimensionless].
     */
    double getueta() const { return ueta_; }
};

#endif  // SRC_CELL_H_
