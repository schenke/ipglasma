#ifndef DATA_H
#define DATA_H

struct ReturnValue {
  double x;
  double y;
  int collided;
  int acceptances;
};

typedef struct init_data {
  double **gmunu; /* metric */

  int nx;
  int ny;
  int neta;
  int nt;

  double x_size;   /* in fermi -x_size/2 < x < x_size/2 */
  double y_size;   /* in fermi, ditto */
  double eta_size; /* ditto */
  double tau_size; /* tau_0 < tau < tau0+tau_size  */
  double tau0;

  double delta_x;
  double delta_y;
  double delta_eta;
  double delta_tau;

  double epsilon0;
  double rhoB0;
  double eta_fall_off;
  double eta_flat;
  double R_A;
  double a_A;
  double a_short;
  double a_long;
  int rk_order;
  int taufac;

  double minmod_theta;

  double SigmaNN;
  double b;
  double bmin;
  double bmax;
  char *Target;
  char *Projectile;
  char *initName;

  int LexusImax;
  int facTau; // maximally 10

  double deltaY;
  double ymax;

  int mode; // 1: do everything;
            // 2: do hydro evolution only;
            // 3: do calculation of thermal spectra only;
            // 4: do resonance decays only

  // eccentricities and axes with respect to which we have to compute v_2 or v_3
  double ecc2;
  double ecc3;
  double ecc3r3;
  double Psi2;
  double Psi3;
  double Psi3r3;

  int seed;
  double sigma0;
  double TFO;   // freeze-out temperature. Used if useEpsFO=0
  int useEpsFO; // if 1, use energy density value to define freeze out
                // condition, if 0 use temperature in TFO
  int size;
  int rank;
} InitData;

#endif
