/**
 * @file main.cpp
 *
 * @brief
 *  PROGRAM TITLE: MD_bulk_activity
 *  Molecular dynamics code in 2D for a single active flexible polymer
 * (LJ+H.P.+activity) Aitor Martin 2017. 18/10/2017
 *
 */

// TODO check the problen of voroni

#include <float.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
#include "./prng.h"
#include "H5Cpp.h"
using namespace std;

#ifdef _OPENMP
#define N_THREADS 4

#else
#define N_THREADS 1
#endif

// *****************Physical*parameters*and*geometry**************************
uint64_t MD_total_step;
/**< total MPC time*/
// @param MD_total_step_over100 MD_total_step/1000
uint64_t MD_total_step_over100;
// @param t total steps
uint64_t t_step_relax;

uint64_t t_step_relax_over100;
uint64_t MD_step;

bool is_active = true;

// ****************** HDF5 File *********************
// const H5File CfgFile;
const H5std_string FILE_NAME;
const H5std_string DATASET_NAME;

// ***********************output file define*********
// uint64_t FREQ_CFG;
int FREQ_CFG_DETA = 1;
uint64_t FREQ_EVOL;
// frequency of calculate active pressure
// uint64_t FREQ_CFG_vcorr; //freq evol and pressure
// uint64_t limit_vcorr; //freq evol and pressure
int ic;
// ic = 0 random pos and ic = 1 read from file
int aV;
int aV_data;
int cfgs_per_file = 500;
int cfg_conter = 0;
double dt;
// time step Molecular dynamics
uint step_t_1;
uint step_t_10;
uint step_t_100;
uint step_t_1000;
double cfg_dt = 1.;
// uint64_t t_step_relax_vcorr;

// @param Lboxx Size of the box X
const double Lboxx = 200;
// @param Lboxy Size of the box Y
const double Lboxy = 200;
// @param Lboxz Size of the box Z
const double Lboxz = 1.5;
// @param half_Lboxx Half of the Box X
const double half_Lboxx = 0.5 * Lboxx;
// @param half_Lboxy Half of the Box Y
const double half_Lboxy = 0.5 * Lboxy;
// @param half_Lboxz Half of the Box Z
const double half_Lboxz = 0.5 * Lboxz;
// @param inv_Lboxx inerse of box size 1/X
const double inv_Lboxx = 1. / Lboxx;
// @param inv_Lboxy inerse of box size 1/Y
const double inv_Lboxy = 1. / Lboxy;
// @param inv_Lboxz inerse of box size 1/Z
const double inv_Lboxz = 1. / Lboxz;

const double PI = 4.0 * atan(1.0);

int press_count;
int cfg_count;
// @param b_dt beta * dt
double b_dt;
// @param half_dt double dt;
double half_dt;
double b_half_dt;
double b_half_dt2;

const double k_BT = 1.0;  // energy units K_B * Temperature
double Lm = 1.0;          // Kuhn length
double Pe;                // Pe = v_0/(Lm * Dr)
double Pe_over_Delta;     // Pe = v_0/(Lm * Dr)
double Kp;                // Gaussian chain constant
double Kb;                // Bending constant

double density;
double density_prev;
// @param Nm number of monomers
int Nm;
// @param Np number of polymers
int Np;
int Nm1p;  // Number of monomers per polymer
bool is_polymer;
int bond_counter;
// @param mmass monomer mass
const double mmass = 1.;
double Diff_rot = 0.03;
double Delta = 1. / 3.;
// @param alpha 1. / (Delta*Diff_rot);
double alpha = 1. / (Delta * Diff_rot);
double eps;
const double Rh = 0.5;
// @param Rh Hydrodynamic Radius
// TODO change the sigma there
const double sigma = 1.;  // 2*Rh
// const double sigma = pow(32., 1. / 6.);  // 2*Rh
const double sigma2 = sigma * sigma;
// @param sigmawall JL_potential sigma of the wall
const double sigmawall = 0.7;
const double sigmawall2 = sigmawall * sigmawall;
double beta;  ///< position factor
double bta;  // thermal noise
double Ito;   // amplitudse of the gaussian distribution for the active vector
const double rcutLJ = pow(2., 1. / 6.) * sigma;  // cutoff radius in LJ
                                                 // potential
const double rcutLJ2 = rcutLJ * rcutLJ;
const double rcutLJWall = pow(2., 1. / 6.) * sigmawall;
const double rcutLJWall2 = rcutLJWall * rcutLJWall;
const int Ncx = static_cast<int>(Lboxx / rcutLJ);
const int Ncy = static_cast<int>(Lboxy / rcutLJ);
const int Nc = Ncx * Ncy;
const double Nlx = static_cast<double>(Ncx) / Lboxx;
const double Nly = static_cast<double>(Ncy) / Lboxy;
// *****************Computational*Quantities*Energy*&*Transport*Coef***********
const int map_neigh = 4;
const int Nneigh = Nc * map_neigh;
int *map;
int *head;
int *list;
// double endend2, rgyr2;  // end-to-end distance sum
// mean r of Nm
double derr_mol = sigma / 4.;
double derr_mol2 = derr_mol * derr_mol;
// limit distance to compare the errors
double E_kinx_pol, E_kiny_pol, E_kinz_pol;
// Kinetic energy in x,y,z coordinates
double E_pot;
// Potential energy for the polymer
double virial_act;
double virial_pass, vir_pass;

// ***************************************************************************
// *************Positions*velocities*and*displacements************************
// ***************************************************************************

double *pos_x, *pos_y, *pos_z;
// x y & z-coordinate of the fluid+polymer at the present step in time t + dt
double *vel_x, *vel_y, *vel_z;
// x y & z-coordinate of the fluid+polymer at the present step in time t + dt
double *pos_x_prev, *pos_y_prev, *pos_z_prev;
// x y & z-coordinate of the fluid+polymer at the present step in time t + dt
double *f_x, *f_y, *f_z;
// forces in time t + dt
double *fx_prev, *fy_prev, *fz_prev;
// forces in time t + dt
double **fx_priv, **fy_priv, **fz_priv;
// forces in time t + dt
double **px_priv, **py_priv;
// pressure in time t + dt
// double *ex;  ///< active vector in x direction
// double *ey;  ///< active vector in y direction
// double *ez;  ///< actice vector in z direction
// // active vector
// double *ex_prev, *ey_prev, *ez_prev;
// active vector
double *bta_x, *bta_y, *bta_z;

// ***************************************************************************
// **************************RANDOM*GENERATOR*********************************
// ***************************************************************************
cPRNG<MT19937, NORMAL::ZIGGURAT> prng[N_THREADS];
void InitRNG() {
  for (int thid = 0; thid < N_THREADS; ++thid) {
    prng[thid].init(time(NULL) + thid * 2013);
    for (int k = 0; k < 1048756; ++k) prng[thid].gen_int32();
  }
}

// *****************************POINTERS**************************************

void initialization(int Nm, int Nneigh, int Nc) {
  /**
   * @brief the initialize of position, velocity, force,
   */
  pos_x = new double[Nm];
  pos_y = new double[Nm];
  pos_z = new double[Nm];
  vel_x = new double[Nm];
  vel_y = new double[Nm];
  vel_z = new double[Nm];
  pos_x_prev = new double[Nm];
  pos_y_prev = new double[Nm];
  pos_z_prev = new double[Nm];
  f_x = new double[Nm];
  f_y = new double[Nm];
  f_z = new double[Nm];
  fx_prev = new double[Nm];
  fy_prev = new double[Nm];
  fz_prev = new double[Nm];
  fx_priv = new double *[N_THREADS];
  fy_priv = new double *[N_THREADS];
  fz_priv = new double *[N_THREADS];
  for (int i = 0; i < N_THREADS; ++i) {
    fx_priv[i] = new double[Nm];
    fy_priv[i] = new double[Nm];
    fz_priv[i] = new double[Nm];
  }
  px_priv = new double *[N_THREADS];
  py_priv = new double *[N_THREADS];
  for (int i = 0; i < N_THREADS; ++i) {
    px_priv[i] = new double[Nm];
    py_priv[i] = new double[Nm];
  }
  bta_x = new double[Nm];
  bta_y = new double[Nm];
  bta_z = new double[Nm];
  // ex = new double[Nm];
  // ey = new double[Nm];
  // ez = new double[Nm];
  // ex_prev = new double[Nm];
  // ey_prev = new double[Nm];
  // ez_prev = new double[Nm];
  map = new int[Nneigh + 1];
  head = new int[Nc + 1];
  list = new int[Nm + 1];
}

// @brief close the pointers of the whole code
void close_pointers(void) {
  delete[] pos_x;
  delete[] pos_y;
  delete[] pos_z;
  delete[] vel_x;
  delete[] vel_y;
  delete[] vel_z;
  delete[] pos_x_prev;
  delete[] pos_y_prev;
  delete[] pos_z_prev;
  delete[] f_x;
  delete[] f_y;
  delete[] f_z;
  delete[] fx_prev;
  delete[] fy_prev;
  delete[] fz_prev;
  for (int i = 0; i < N_THREADS; ++i) {
    delete[] fx_priv[i];
    delete[] fy_priv[i];
    delete[] fz_priv[i];
  }
  delete[] fx_priv;
  delete[] fy_priv;
  delete[] fz_priv;
  for (int i = 0; i < N_THREADS; ++i) {
    delete[] px_priv[i];
    delete[] py_priv[i];
  }
  delete[] px_priv;
  delete[] py_priv;
  // delete[] ex;
  // delete[] ey;
  // delete[] ez;
  // delete[] ex_prev;
  // delete[] ey_prev;
  // delete[] ez_prev;
  delete[] bta_x;
  delete[] bta_y;
  delete[] bta_z;
  delete[] map;
  delete[] head;
  delete[] list;
}
// ***************************************************************************
// *********************PERIODIC*BOUND*CONDITIONS*****************************
// ***************************************************************************
// computes the x,y-coord of the displacement vector between two disks in a
// torus (PBC)
double vec_PBCx(double x1, double x2) {
  double rx, xij, xi, xj;
  xi = x1 - floor(x1 * inv_Lboxx) * Lboxx;
  xj = x2 - floor(x2 * inv_Lboxx) * Lboxx;
  xij = (xi - xj);
  if (xij < (-1. * half_Lboxx))
    rx = xij + Lboxx;
  else if (xij > half_Lboxx)
    rx = xij - Lboxx;
  else
    rx = xij;

  return rx;
}

double vec_PBCy(double y1, double y2) {
  double ry, yij;
  yij = (y1 - y2) - floor((y1 - y2) * inv_Lboxy) * Lboxy;

  if (yij < (-1. * half_Lboxy))
    ry = yij + Lboxy;
  else if (yij > half_Lboxy)
    ry = yij - Lboxy;
  else
    ry = yij;

  return ry;
}

// ***************************************************************************
// ****************rand*particles*initialization******************************
// ****************RAND*POSITONS*MONOMERS*STRAIGHT********************

void init_pol_straight(void) {
  int row_x, row_y, kk;
  double dr[3];

  row_x = static_cast<int>(Lboxx / static_cast<double>(Nm1p));
  row_y = 1 + static_cast<int>(static_cast<double>(Np) /
                               static_cast<double>(row_x));
  // Nm1p * *3 > Nm

  for (int j = 0; j < row_y; ++j) {             // pick up a disk
    for (int i = 0; i < (row_x * Nm1p); ++i) {  // pick up a disk
      kk = j * (row_x * Nm1p) + i;
      if (kk >= Nm) {
        goto sal;
      }
      pos_x[kk] = static_cast<double>(i) * Lm + 0.5 * sigma +
                  0.25 * (1 - (j % 2)) * Nm1p;
      // save t = 0 values in x,y & z coord
      pos_y[kk] =
          static_cast<double>(j) * Lboxy / static_cast<double>(row_y + 3) +
          0.5 * sigma;
      // save t = 0 values in x,y & z coord
      pos_z[kk] = half_Lboxz;
      // save t = 0 values in x,y & z coord
    }
  }
sal:
  cout << "init config: OK" << endl;
}

void new_init_pol(void) {
  /**
   * @brief initial initialize the position in the slop of arcsin(sigma/Lboxx)
   *
   * like :
   * _________________________
   * |              |
   * |              |
   * |              |
   * |             *|*******....
   * |**************|
   * _______________|__________
   *
   */

  double L_snake_x =
      Lboxy / sigma * sqrt(Lboxx * Lboxx - sigma2) / static_cast<double>(Np);
  double L_snake_y = Lboxy / static_cast<double>(Np);
  double segment_x = Lm * sqrt(Lboxx * Lboxx - sigma2) / Lboxx;
  double segment_y = Lm * sigma / Lboxx;
  int kk;

  for (int i = 0; i < Np; ++i) {
    for (int j = 0; j < Nm1p; ++j) {
      kk = i * Nm1p + j;
      pos_x[kk] = static_cast<double>(i) * L_snake_x +
                  static_cast<double>(j) * segment_x;
      pos_y[kk] = static_cast<double>(i) * L_snake_y +
                  static_cast<double>(j) * segment_y;
      pos_z[kk] = half_Lboxz;
    }
  }

  cout << "init config: OK" << endl;
}
//****************************************************************************
// ****************FILE*FUNCTIONS*READING*&*WRITING***************************
//****************************************************************************
/**
 * @fn read_init_file(void)
 * @brief read the initial data form previous simulation
 */

void read_init_file(void) {
  char file_name[50];
  double x_out, y_out, z_out, vxx, vyy, vzz;
  //  exx, eyy, ezz;

  sprintf(file_name, "../init_cfgs/read_init_cfg_pos_vel_Nm_%d_Pe_%4.2f.xyz",
          Nm, Pe);
  ifstream fifi(file_name, ios::in);

  for (int i = 0; i < Nm; ++i) {
    fifi >> x_out >> y_out >> z_out >> vxx >> vyy >> vzz;
    //  >> exx >> eyy >> ezz;
    pos_x[i] = x_out;
    pos_y[i] = y_out;
    pos_z[i] = z_out;
    vel_x[i] = vxx;
    vel_y[i] = vyy;
    vel_z[i] = vzz;
    // ex[i] = exx;
    // ey[i] = eyy;
    // ez[i] = ezz;
  }
  fifi.close();

  sprintf(file_name, "../init_cfgs/params_Nm_%d_Pe_%4.2f.xyz", Nm, Pe);
  ifstream read_params(file_name, ios::in);
  read_params >> cfg_count >> press_count >> virial_act >> virial_pass;
}

void print_last_cfg(void) {
  char file_name[50];
  sprintf(file_name, "read_init_cfg_pos_vel_Nm_%d_Pe_%4.2f.xyz", Nm, Pe);
  ofstream last_cfg_pos_vel(file_name, ios::trunc);
  last_cfg_pos_vel.precision(7);
  for (int i = 0; i < Nm; ++i) {
    last_cfg_pos_vel << pos_x[i] << " " << pos_y[i] << " " << pos_z[i] << " "
                     << vel_x[i] << " " << vel_y[i] << " " << vel_z[i] << endl;
                    //  << " "
                    //  << ex[i] << " " << ey[i] << " " << ez[i] << endl;
  }
  last_cfg_pos_vel.close();

  sprintf(file_name, "params_Nm_%d_Pe_%4.2f.xyz", Nm, Pe);
  ofstream params(file_name, ios::trunc);
  params << cfg_count << " " << press_count << " " << virial_act << " "
         << virial_pass << " " << aV << endl;
  params.close();
}

void cell_map(void) {
  int imap;
  int cell[Ncx + 2][Ncy + 2];

  /**
   * @brief generate cell index (exp lx=ly=200)
   * cell[0][0] = 31684, cell[1][0] = 31507, cell[2][0] = 31508;
   * cell[0][1] = 178,   cell[1][1] = 1,     cell[2][1] = 2
   * cell[0][2] = 356,   cell[1][2] = 179,   cell[2][2] = 180
   * cell[i][j] > 0
   */
  for (int j = 0; j <= Ncy + 1; ++j) {
    for (int i = 0; i <= Ncx + 1; ++i) {
      cell[i][j] = 1 + ((i - 1 + Ncx) % Ncx) + ((j - 1 + Ncy) % Ncy) * Ncx;
    }
  }
  for (int j = 1; j <= Ncy; ++j) {
    for (int i = 1; i <= Ncx; ++i) {
      imap = (cell[i][j] - 1) * map_neigh;

      map[imap + 1] = cell[i + 1][j];
      map[imap + 2] = cell[i + 1][j + 1];
      map[imap + 3] = cell[i][j + 1];
      map[imap + 4] = cell[i - 1][j + 1];
    }
  }
}

void neigh_cell(void) {
  double xx, yy;
  int icell;

  for (int j = 0; j <= Nc; ++j) {
    head[j] = -1;
  }

  for (int i = 0; i < Nm; ++i) {
    xx = Nlx * (pos_x[i] - floor(pos_x[i] * inv_Lboxx) * Lboxx);
    yy = Nly * (pos_y[i] - floor(pos_y[i] * inv_Lboxy) * Lboxy);
    icell = 1 + static_cast<int>(xx) + static_cast<int>(yy) * Ncx;
    // TODO chang the index start from 0
    list[i + 1] = head[icell];
    head[icell] = i + 1;
  }
}

// ***************************************************************************
// *****************************MOLECULAR*DYNAM*******************************
// ***************************************************************************
// *************************VERLET*LIST*FOR*MONOMER***************************
// *************************************INIT**********************************
void init_velocities(void) {
  double v0_xCM, v0_yCM, v0_zCM;
  // x y z & total velocitie of the CM of the solvent
  double *v_aux_x = new double[Nm], *v_aux_y = new double[Nm],
         *v_aux_z = new double[Nm];
  double initE_kin, mass;
  int i;

  for (i = 0; i < Nm; ++i) {
  redox:
    v_aux_x[i] = prng[0].gen_gauss();
    if (fabs(b_dt * v_aux_x[i]) > 0.1 * derr_mol) {
      goto redox;
    }
  redoy:
    v_aux_y[i] = prng[0].gen_gauss();
    if (fabs(b_dt * v_aux_y[i]) > 0.1 * derr_mol) {
      goto redoy;
    }
  redoz:
    v_aux_z[i] = prng[0].gen_gauss();
    if (fabs(b_dt * v_aux_z[i]) > 0.1 * derr_mol) {
      goto redoz;
    }
  }

  v0_xCM = 0.;
  v0_yCM = 0.;
  v0_zCM = 0.;
  for (i = 0; i < Nm; ++i) {
    v0_xCM += v_aux_x[i];
    v0_yCM += v_aux_y[i];
    v0_zCM += v_aux_z[i];
  }

  initE_kin = 0.;
  for (i = 0; i < Nm; ++i) {
    vel_x[i] = v_aux_x[i] - v0_xCM / static_cast<double>(Nm);
    vel_y[i] = v_aux_y[i] - v0_yCM / static_cast<double>(Nm);
    vel_z[i] = v_aux_z[i] - v0_zCM / static_cast<double>(Nm);
    initE_kin +=
        (vel_x[i] * vel_x[i] + vel_y[i] * vel_y[i] + vel_z[i] * vel_z[i]);
  }

  double rescFAC = sqrt(3.0 * Nm * k_BT / initE_kin);
  for (i = 0; i < Nm; ++i) {
    vel_x[i] *= rescFAC;
    vel_y[i] *= rescFAC;
    vel_z[i] *= rescFAC;
  }

  delete[] v_aux_x;
  delete[] v_aux_y;
  delete[] v_aux_z;
}

/**
 * @brief initialisation activity
 */
// void init_activity(void) {
//   double eax, eay, eaz, enorm2, enorm;
//   for (int i = 0; i < Nm; ++i) {  // for each monomer 'i'
//     do {
//       eax = 1.0 - 2.0 * drand48();
//       eay = 1.0 - 2.0 * drand48();
//       eaz = 1.0 - 2.0 * drand48();

//       enorm2 = eax * eax + eay * eay + eaz * eaz;
//     } while (1.0 < enorm2);

//     enorm = sqrt(enorm2);
//     ex[i] = eax / enorm;
//     ey[i] = eay / enorm;
//     ez[i] = eaz / enorm;
//   }
// }

// ***************************************************************************
// **********************************MD*FORCES********************************
// ***************************************************************************
/**
 * @brief calculate the bond force
 */
void calcBondForces(void) {
  double dr[3], fbond[3];
  double r, coeff;
  int i, j, k, ii;

  for (k = 0; k < Np; ++k) {
    // i-first particle
    for (ii = 0; ii < Nm1p - 1; ++ii) {
      // i-first particle
      i = Nm1p * k + ii;
      j = i + 1;  // second particle
      dr[0] = pos_x[j] - pos_x[i];
      dr[1] = pos_y[j] - pos_y[i];
      dr[2] = pos_z[j] - pos_z[i];
      r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

      coeff = Kp * Lm / r - Kp;
      fbond[0] = coeff * dr[0];
      fbond[1] = coeff * dr[1];
      fbond[2] = coeff * dr[2];

      // Forces
      f_x[i] -= fbond[0];
      f_y[i] -= fbond[1];
      f_z[i] -= fbond[2];

      f_x[j] += fbond[0];
      f_y[j] += fbond[1];
      f_z[j] += fbond[2];

      E_pot += (0.5 * Kp) * (r - Lm) * (r - Lm);
      // sum Bond energy
    }
  }
}

/**
 * @brief calculate bond force between i, j
 */
void calcBondForces_ij(int i, int j) {
  double dr[3], fbond[3];
  double r, coeff;

  dr[0] = pos_x[j] - pos_x[i];
  dr[1] = pos_y[j] - pos_y[i];
  dr[2] = pos_z[j] - pos_z[i];
  r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);

  coeff = Kp * Lm / r - Kp;
  fbond[0] = coeff * dr[0];
  fbond[1] = coeff * dr[1];
  fbond[2] = coeff * dr[2];

  // Forces
  f_x[i] -= fbond[0];
  f_y[i] -= fbond[1];
  f_z[i] -= fbond[2];

  f_x[j] += fbond[0];
  f_y[j] += fbond[1];
  f_z[j] += fbond[2];

  E_pot += (0.5 * Kp) * (r - Lm) * (r - Lm);
}

void calcBendingForces(void) {
  double dr0[2], dr1[2];
  double r0q, r1q, coeff_x, coeff_y;

  for (int k = 0; k < Np; ++k) {  // k-polymer
    int jj = Nm1p * k;

    coeff_x = -Kb * (pos_x[jj + 2] - 2 * pos_x[jj + 1] + pos_x[jj]);
    coeff_y = -Kb * (pos_y[jj + 2] - 2 * pos_y[jj + 1] + pos_y[jj]);
    f_x[jj] += coeff_x;
    f_y[jj] += coeff_y;

    coeff_x = -Kb * (pos_x[jj + 3] - 4 * pos_x[jj + 2] + 5 * pos_x[jj + 1] -
                     2 * pos_x[jj]);

    coeff_y = -Kb * (pos_x[jj + 3] - 4 * pos_x[jj + 2] + 5 * pos_x[jj + 1] -
                     2 * pos_x[jj]);

    f_x[jj + 1] += coeff_x;
    f_y[jj + 1] += coeff_y;

    // middle of the chain
    for (int ii = 2; ii < Nm1p - 2; ++ii) {
      int i = jj + ii;
      coeff_x = -Kb * (pos_x[i + 2] - 4 * pos_x[i + 1] + 6 * pos_x[i] -
                       4 * pos_x[i - 1] + pos_x[i - 2]);

      coeff_y = -Kb * (pos_y[i + 2] - 4 * pos_y[i + 1] + 6 * pos_y[i] -
                       4 * pos_y[i - 1] + pos_y[i - 2]);

      f_x[i] += coeff_x;
      f_y[i] += coeff_y;
    }
    // end of the chain
    coeff_x = -Kb * (pos_x[jj + Nm1p - 4] - 4 * pos_x[jj + Nm1p - 3] +
                     5 * pos_x[jj + Nm1p - 2] - 2 * pos_x[jj + Nm1p - 1]);

    coeff_y = -Kb * (pos_y[jj + Nm1p - 4] - 4 * pos_y[jj + Nm1p - 3] +
                     5 * pos_y[jj + Nm1p - 2] - 2 * pos_y[jj + Nm1p - 1]);

    f_x[jj + Nm1p - 2] += coeff_x;
    f_y[jj + Nm1p - 2] += coeff_y;

    coeff_x = -Kb * (pos_x[jj + Nm1p - 3] - 2 * pos_x[jj + Nm1p - 2] +
                     pos_x[jj + Nm1p - 1]);
    coeff_y = -Kb * (pos_y[jj + Nm1p - 3] - 2 * pos_y[jj + Nm1p - 2] +
                     pos_y[jj + Nm1p - 1]);
    f_x[jj + Nm1p - 1] += coeff_x;
    f_y[jj + Nm1p - 1] += coeff_y;
  }
}

void calcLJForces(void) {
  double fLJ[3], dr[3];
  double r2, r_inv2, sr_inv2, sr_inv6, sr_inv12, coeffLJ;
  int i, j, k, neigh, kneigh, nn, pp, j_nb[2];
  int thread, n_threads;
  double fxj, fyj, fzj, pxj, pyj;

  bond_counter = 0;
  neigh_cell();

#pragma omp parallel private(thread, i, j, k, neigh, kneigh, nn, pp, j_nb,  \
                             fxj, fyj, fzj, pxj, pyj, coeffLJ, fLJ, dr, r2, \
                             r_inv2, sr_inv2, sr_inv6, sr_inv12)            \
    shared(pos_x, pos_y, pos_z)

  n_threads = omp_get_num_threads();
  thread = omp_get_thread_num();
  for (i = 0; i < Nm; ++i) {
    fx_priv[thread][i] = 0.0;
    fy_priv[thread][i] = 0.0;
    fz_priv[thread][i] = 0.0;
    px_priv[thread][i] = 0.0;
    py_priv[thread][i] = 0.0;
  }
#pragma omp for reduction(+ : E_pot, bond_counter)
  for (k = 1; k <= Nc; ++k) {
    j = head[k];
    while (j > 0) {
      // if j == -1 stop the loop (see list)
      // @brief for polymer generate the bonding particle index
      if (is_polymer) {
        pp = (j - 1) % Nm1p;  // @param position in polymer
        j_nb[0] = j - 1;
        j_nb[1] = j + 1;
        if (pp == 0) {
          j_nb[0] = -1;
        } else if (pp == Nm1p - 1) {
          j_nb[1] = -1;
        }
        // cout << "pp " << pp << " j_nb[0] " << j_nb[0] << " j " << j << "
        // j_nb[1] " << j_nb[1] << endl;
      }
      i = list[j];
      fxj = fx_priv[thread][j - 1];
      fyj = fy_priv[thread][j - 1];
      fzj = fz_priv[thread][j - 1];

      pxj = px_priv[thread][j - 1];
      pyj = py_priv[thread][j - 1];

      while (i > 0) {
        // skip the bonding particle

        if (is_polymer && (i == j_nb[0] || i == j_nb[1])) {
          calcBondForces_ij(j - 1, i - 1);
          // cout << "i = " << i - 1 << "j = " << j - 1 << endl;
          i = list[i];
          bond_counter += 1;
          continue;
        }
        dr[0] = vec_PBCx(pos_x[j - 1], pos_x[i - 1]);
        dr[1] = vec_PBCy(pos_y[j - 1], pos_y[i - 1]);
        dr[2] = pos_z[j - 1] - pos_z[i - 1];
        r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
        if (r2 < rcutLJ2) {
          r_inv2 = 1. / r2;
          sr_inv2 = sigma2 * r_inv2;
          sr_inv6 = sr_inv2 * sr_inv2 * sr_inv2;
          sr_inv12 = sr_inv6 * sr_inv6;

          coeffLJ = 24. * eps * r_inv2 * (sr_inv12 + sr_inv12 - sr_inv6);
          fLJ[0] = coeffLJ * dr[0];
          fLJ[1] = coeffLJ * dr[1];
          fLJ[2] = coeffLJ * dr[2];

          px_priv[thread][i - 1] += fLJ[0] * dr[0];
          py_priv[thread][i - 1] += fLJ[1] * dr[1];

          fx_priv[thread][i - 1] -= fLJ[0];
          fy_priv[thread][i - 1] -= fLJ[1];
          fz_priv[thread][i - 1] -= fLJ[2];

          fxj += fLJ[0];
          fyj += fLJ[1];
          fzj += fLJ[2];

          pxj += fLJ[0] * dr[0];
          pyj += fLJ[1] * dr[1];

          E_pot += 4. * eps * (sr_inv12 - sr_inv6) + eps;
        }
        //}
        i = list[i];
      }
      for (neigh = 1; neigh <= map_neigh; ++neigh) {
        nn = (k - 1) * map_neigh + neigh;
        kneigh = map[nn];
        if (kneigh > 0) {
          i = head[kneigh];
          while (i > 0) {
            if (is_polymer && (i == j_nb[0] || i == j_nb[1])) {
              calcBondForces_ij(j - 1, i - 1);
              // cout << "i = " << i - 1 << " j = " << j - 1 << endl;
              i = list[i];
              bond_counter += 1;
              continue;
            }

            dr[0] = vec_PBCx(pos_x[j - 1], pos_x[i - 1]);
            dr[1] = vec_PBCy(pos_y[j - 1], pos_y[i - 1]);
            dr[2] = pos_z[j - 1] - pos_z[i - 1];
            r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

            if (r2 < rcutLJ2) {
              r_inv2 = 1. / r2;
              sr_inv2 = sigma2 * r_inv2;
              sr_inv6 = sr_inv2 * sr_inv2 * sr_inv2;
              sr_inv12 = sr_inv6 * sr_inv6;

              coeffLJ = 24. * eps * r_inv2 * (sr_inv12 + sr_inv12 - sr_inv6);
              fLJ[0] = coeffLJ * dr[0];
              fLJ[1] = coeffLJ * dr[1];
              fLJ[2] = coeffLJ * dr[2];

              px_priv[thread][i - 1] += fLJ[0] * dr[0];
              py_priv[thread][i - 1] += fLJ[1] * dr[1];

              fx_priv[thread][i - 1] -= fLJ[0];
              fy_priv[thread][i - 1] -= fLJ[1];
              fz_priv[thread][i - 1] -= fLJ[2];

              fxj += fLJ[0];
              fyj += fLJ[1];
              fzj += fLJ[2];

              pxj += fLJ[0] * dr[0];
              pyj += fLJ[1] * dr[1];

              E_pot += 4. * eps * (sr_inv12 - sr_inv6) + eps;
            }

            i = list[i];
          }
        }
      }

      fx_priv[thread][j - 1] = fxj;
      fy_priv[thread][j - 1] = fyj;
      fz_priv[thread][j - 1] = fzj;

      px_priv[thread][j - 1] = pxj;
      py_priv[thread][j - 1] = pyj;
      j = list[j];
    }
  }
  vir_pass = 0.0;
  for (i = 0; i < Nm; ++i) {
    for (thread = 0; thread < n_threads; ++thread) {
      f_x[i] += fx_priv[thread][i];
      f_y[i] += fy_priv[thread][i];
      f_z[i] += fz_priv[thread][i];

      vir_pass += px_priv[thread][i] + py_priv[thread][i];
    }
  }
  // cout << "bond counter: " << bond_counter << "/" << Np * (Nm1p - 1) << endl;
}

void calcLJWalls(void) {
  double drz, r2z, r_inv2w, sr_inv2w, sr_inv6w, sr_inv12w, coeffLJWalls;
  int thread, n_threads, j;

#pragma omp parallel private(thread, j, coeffLJWalls, drz, r2z, r_inv2w, \
                             sr_inv2w, sr_inv6w, sr_inv12w) shared(f_z, pos_z)

  n_threads = omp_get_num_threads();
  thread = omp_get_thread_num();
// #pragma omp for
#pragma omp for reduction(+ : E_pot)
  for (j = 0; j < Nm; ++j) {
    if (pos_z[j] > half_Lboxz) {
      drz = Lboxz - pos_z[j];
    } else {
      drz = -pos_z[j];
    }
    r2z = drz * drz;
    if (r2z < rcutLJWall2) {
      r_inv2w = 1. / r2z;
      sr_inv2w = sigmawall2 * r_inv2w;
      sr_inv6w = sr_inv2w * sr_inv2w * sr_inv2w;
      sr_inv12w = sr_inv6w * sr_inv6w;
      coeffLJWalls = 24. * eps * r_inv2w * (sr_inv12w + sr_inv12w - sr_inv6w);
      f_z[j] -= coeffLJWalls * drz;

      E_pot += 4. * eps * (sr_inv12w - sr_inv6w) + eps;
    }
  }
}

// void calcActivity(void) {
//   double phi, theta, s_th, c_th, s_phi, c_phi, e_thx, e_thy, e_thz, e_phix,
//       e_phiy, dnu_th, dnu_phi, ex_aux, ey_aux, ez_aux, norme;

// #ifdef _OPENMP
// #pragma omp parallel for schedule(guided) private(                             \
//     phi, theta, s_th, c_th, s_phi, c_phi, e_thx, e_thy, e_thz, e_phix, e_phiy, \
//     dnu_th, dnu_phi, ex_aux, ey_aux, ez_aux, norme) shared(ex, ey, ez)         \
//     num_threads(N_THREADS)
// #endif
//   for (int i = 0; i < Nm; ++i) {
// #ifdef _OPENMP
//     int thid = omp_get_thread_num();
// #else
//     int thid = 0;
// #endif

//     // int thid=0;
//     phi = atan2(ey[i], ex[i]);
//     theta = acos(ez[i]);
//     s_th = sin(theta);
//     c_th = cos(theta);
//     s_phi = sin(phi);
//     c_phi = cos(phi);
//     /**
//      * @brief get old e and convert to spherical coordinates
//      */

//     e_thx = c_phi * c_th;
//     e_thy = s_phi * c_th;
//     e_thz = -s_th;

//     e_phix = -s_phi;
//     e_phiy = c_phi;
//     /**
//      * @brief normed orthogonal direction of theta and phi
//      */

//     dnu_th = Ito * prng[thid].gen_gauss();
//     dnu_phi = Ito * prng[thid].gen_gauss();

//     ex_aux = ex[i] + e_thx * dnu_th + e_phix * dnu_phi;
//     ey_aux = ey[i] + e_thy * dnu_th + e_phiy * dnu_phi;
//     ez_aux = ez[i] + e_thz * dnu_th;

//     norme = sqrt(ex_aux * ex_aux + ey_aux * ey_aux + ez_aux * ez_aux);
//     ex[i] = ex_aux / norme;
//     ey[i] = ey_aux / norme;
//     ez[i] = ez_aux / norme;
//   }
// }

void forces(void) {
  /**
   * @brief calculate the force
   */
  // update the neighbour list manually
  for (int i = 0; i < Nm; ++i) {
    fx_prev[i] = f_x[i];
    fy_prev[i] = f_y[i];
    fz_prev[i] = f_z[i];
    // save previous forces, verlet

    // ex_prev[i] = ex[i];
    // ey_prev[i] = ey[i];
    // ez_prev[i] = ez[i];
    // save previous forces, verlet

    f_x[i] = 0.0;
    f_y[i] = 0.0;
    f_z[i] = 0.0;
    // restart to 0 the forces
  }
  calcLJWalls();
  // calcBondForces();
  // TODO need to change bond Force
  // calculate harmonic potential force
  if (Nm1p > 2 && Kb > 0.0) {
    calcBendingForces();
  }
  // calculate Bending potential force
  calcLJForces();
  // calculate Lennard-Jones forces
  // calcActivity();
  // calculate active forces
}

void force_no_active(void) {
  for (int i = 0; i < Nm; ++i) {
    fx_prev[i] = f_x[i];
    fy_prev[i] = f_y[i];
    fz_prev[i] = f_z[i];
    // save previous forces, verlet
    f_x[i] = 0.0;
    f_y[i] = 0.0;
    f_z[i] = 0.0;
    // restart to 0 the forces
  }
  calcLJWalls();
  // calcBondForces();
  // TODO need to change bond force
  // calculate harmonic potential force
  if (Nm1p > 2 && Kb > 0.0) {
    calcBendingForces();
  }
  // calculate Bending potential force
  calcLJForces();
  // calculate Lennard-Jones forces
}

void pos_MD() {
  double disx, disy, disz;
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) private(disx, disy, disz) shared( \
    pos_x, pos_y, pos_z, f_x, f_y, f_z) num_threads(N_THREADS)
    // , ex, ey, ez
#endif

  for (int i = 0; i < Nm; ++i) {  // calculate new positions

#ifdef _OPENMP
    int thid = omp_get_thread_num();
#else
    int thid = 0;
#endif

    // beta random force
    bta_x[i] = bta * prng[thid].gen_gauss();
    bta_y[i] = bta * prng[thid].gen_gauss();
    bta_z[i] = bta * prng[thid].gen_gauss();

    disx = vel_x[i] * b_dt + f_x[i] * b_half_dt2 + b_half_dt * bta_x[i];
    disy = vel_y[i] * b_dt + f_y[i] * b_half_dt2 + b_half_dt * bta_y[i];
    disz = vel_z[i] * b_dt + f_z[i] * b_half_dt2 + b_half_dt * bta_z[i];

    // disx += Pe_over_Delta * ex[i] * b_half_dt2;
    // disy += Pe_over_Delta * ey[i] * b_half_dt2;
    // disz += Pe_over_Delta * ez[i] * b_half_dt2;

    pos_x_prev[i] = pos_x[i];
    pos_y_prev[i] = pos_y[i];
    pos_z_prev[i] = pos_z[i];

    pos_x[i] += disx;
    pos_y[i] += disy;
    pos_z[i] += disz;

    if ((disx * disx + disy * disy + disz * disz) > derr_mol2) {
      // too much displacement
      cout << "PUUUUM" << endl;

      cout << "time " << MD_step * dt << ", i = " << i
           << " btax = " << b_half_dt * bta_x[i]
           << " btay = " << b_half_dt * bta_y[i]
           << " btaz = " << b_half_dt * bta_z[i] << endl;

      cout << "time " << MD_step * dt << ", i = " << i
           << " vel_x = " << vel_x[i] * b_dt << " vel_y = " << vel_y[i] * b_dt
           << " vel_z = " << vel_z[i] * b_dt << endl;

      cout << "time " << MD_step * dt << ", i = " << i
           << " fx = " << b_half_dt2 * f_x[i] << " fy = " << b_half_dt2 * f_y[i]
           << " fz = " << b_half_dt2 * f_z[i] << endl;

      cout << "time " << MD_step * dt << ", i = " << i << " x = " << pos_x[i]
           << " y = " << pos_y[i] << " z = " << pos_z[i] << endl;

      cout << "time " << MD_step * dt << ", i = " << i << " dx = " << disx
           << " dy = " << disy << " dz = " << disz << endl;

      exit(1);
    }
  }
}

void pos_MD_no_cative(void) {
  double disx, disy, disz;
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) private(disx, disy, disz) \
    shared(pos_x, pos_y, pos_z, f_x, f_y, f_z) num_threads(N_THREADS)
#endif
  for (int i = 0; i < Nm; ++i) {
  // calculate new positions

#ifdef _OPENMP
    int thid = omp_get_thread_num();
#else
    int thid = 0;
#endif

    disx = vel_x[i] * b_dt + f_x[i] * b_half_dt2;
    disy = vel_y[i] * b_dt + f_y[i] * b_half_dt2;
    disz = vel_z[i] * b_dt + f_z[i] * b_half_dt2;

    pos_x_prev[i] = pos_x[i];
    pos_y_prev[i] = pos_y[i];
    pos_z_prev[i] = pos_z[i];

    pos_x[i] += disx;
    pos_y[i] += disy;
    pos_z[i] += disz;

    if ((disx * disx + disy * disy + disz * disz) > derr_mol2) {
      // too much displacement
      cout << "PUUUUM" << endl;

      cout << "time " << MD_step * dt << ", i = " << i
           << " btax = " << b_half_dt * bta_x[i]
           << " btay = " << b_half_dt * bta_y[i]
           << " btaz = " << b_half_dt * bta_z[i] << endl;

      cout << " vel_x = " << vel_x[i] * b_dt << " vel_y = " << vel_y[i] * b_dt
           << " vel_z = " << vel_z[i] * b_dt << endl;

      cout << " fx = " << b_half_dt2 * f_x[i] << " fy = " << b_half_dt2 * f_y[i]
           << " fz = " << b_half_dt2 * f_z[i] << endl;

      cout << " x = " << pos_x[i] << " y = " << pos_y[i] << " z = " << pos_z[i]
           << endl;

      cout << " dx = " << disx << " dy = " << disy << " dz = " << disz << endl;

      exit(1);
    }
  }
}

void vel_MD() {
  for (int i = 0; i < Nm; ++i) {
    // calculate new velocities
    vel_x[i] +=
        (fx_prev[i] + f_x[i] 
        // + Pe_over_Delta * (ex[i] + ex_prev[i])
        ) * half_dt -
        alpha * (pos_x[i] - pos_x_prev[i]) + bta_x[i];

    vel_y[i] +=
        (fy_prev[i] + f_y[i]
        //  + Pe_over_Delta * (ey[i] + ey_prev[i])
        ) * half_dt -
        alpha * (pos_y[i] - pos_y_prev[i]) + bta_y[i];

    vel_z[i] +=
        (fz_prev[i] + f_z[i] 
        // + Pe_over_Delta * (ez[i] + ez_prev[i])
        ) * half_dt -
        alpha * (pos_z[i] - pos_z_prev[i]) + bta_z[i];
  }
}

void vel_MD_no_active(void) {
  for (int i = 0; i < Nm; ++i) {
    // calculate new velocities
    vel_x[i] += (fx_prev[i] + f_x[i]) * half_dt;
    vel_y[i] += (fy_prev[i] + f_y[i]) * half_dt;
    vel_z[i] += (fz_prev[i] + f_z[i]) * half_dt;
  }
}

void print_ener_pol(void) {
  double E_kin = 0.0;
  for (int i = 0; i < Nm; ++i) {
    E_kin += vel_x[i] * vel_x[i] + vel_y[i] * vel_y[i] + vel_z[i] * vel_z[i];
  }

  /**
   * @brief write energy and protential to ener_pol file with column:
   * time, average potential, average kinetic energy, total energy
   */
  cout.width(8);
  cout.precision(8);

  cout << static_cast<double>(MD_step * dt) << "\t"
       << E_pot / static_cast<double>(Nm) << "\t"
       << 0.5 * E_kin / static_cast<double>(Nm) << "\t"
       << (E_pot + 0.5 * E_kin) / static_cast<double>(Nm) << endl;
}

void mol_dyn(void) {
  /**
   * @brief the main molecular dynamic simulation process
   */
  E_pot = 0.;
  if (is_active) {
    pos_MD();
    // update positions
    forces();
    // calculate force from potential
    vel_MD();
    // update velocities
  } else {
    pos_MD_no_cative();
    force_no_active();
    vel_MD_no_active();
  }
  if (MD_step % static_cast<int>(0.1 / dt) == 0) {
    print_ener_pol();
  }
}

// ***************************************************************************
// ***************************************************************************
// *****************DYNAMIC/STATIC*PROPERTIES*&*CORRELATIONS******************
// ***************************************************************************
void evol(ofstream &ener_pol) {
  /**
   * @brief write to file ener_pol, e_act, endtoend
   */

  // cell-level scaling thermostat in the canonical ensemble
  cfg_count += 1;
  // saving configuration of the polymer
  double E_kin = 0.0;
  for (int i = 0; i < Nm; ++i) {
    E_kin += vel_x[i] * vel_x[i] + vel_y[i] * vel_y[i] + vel_z[i] * vel_z[i];
  }

  /**
   * @brief write energy and protential to ener_pol file with column:
   * time, average potential, average kinetic energy, total energy
   */
  ener_pol << static_cast<double>(MD_step * dt) << " "
           << E_pot / static_cast<double>(Nm) << " "
           << 0.5 * E_kin / static_cast<double>(Nm) << " "
           << (E_pot + 0.5 * E_kin) / static_cast<double>(Nm) << endl;
}

// void cfg_files(int aV_sf, int numm, ofstream &cfg) {
//   for (int i = 0; i < Nm; ++i) {
//     cfg << pos_x[i] << " " << pos_y[i] << " "
//         << pos_z[i]
//         // << " " << vel_x[i] << " " << vel_y[i] << " " << vel_z[i]
//         // << " " << ex[i] << " " << ey[i] << " " << ez[i]
//         << endl;
//   }
//   cfg << setw(0) << endl;

//   if ((aV_sf % cfgs_per_file) == 0) {
//     aV += 1;
//     cfg.close();
//     char file_name[80];
//     sprintf(file_name,
//             "cfg_vel_s%d_aV_%d_phi_%5.4f_Nm_%d_Nm1p_%d_Pe_%4.2f_stf_%4.2f.dat",
//             numm, aV, density, Nm, Nm1p, Pe, Kb);
//     cfg.open(file_name);
//   }
// }

void cfg_data(ofstream &cfg_data_file, ofstream &cfg_time_file) {
  for (int i = 0; i < Nm; ++i) {
    cfg_data_file << pos_x[i] << " " << pos_y[i] << " " << pos_z[i] << " "
                  << vel_x[i] << " " << vel_y[i] << " " << vel_z[i] << endl;
                  // << " "
                  // << ex[i] << " " << ey[i] << " " << ez[i] << endl;
  }
  cfg_data_file << setw(0) << endl;
  cfg_time_file << static_cast<double>(MD_step) * dt << endl;
}

// void pressure(ofstream &act_pres) {
//   press_count += 1;
//   double vir_act = 0.;
//   for (int i = 0; i < Nm; ++i) {
//     // just one/the first polymer
//     vir_act += ex[i] * pos_x[i] + ey[i] * pos_y[i];
//   }

//   virial_act =
//       virial_act * (press_count - 1) / static_cast<double>(press_count) +
//       vir_act / static_cast<double>(press_count);

//   virial_pass =
//       virial_pass * (press_count - 1) / static_cast<double>(press_count) +
//       vir_pass / static_cast<double>(press_count);

//   if ((MD_step % FREQ_EVOL) == 0) {
//     act_pres << static_cast<double>(MD_step * dt) << " "
//              << Pe_over_Delta * virial_act / (2. * Lboxx * Lboxy) << " "
//              << virial_pass / (4. * Lboxx * Lboxy) << endl;
//   }
// }

// *****************************************************************************
// *****************************************************************************
// ***********************************MAIN**************************************
// *****************************************************************************

int main(int argc, char *argv[]) {
  InitRNG();

  // declarations of flie names
  char file_name[80];
  cout << "running thread" << N_THREADS << endl;
  // time counter
  clock_t start = clock();
  time_t sec;
  sec = time(NULL);
  // set the start time
  time_t now = time(0);
  // current date/time based on current system
  char *real_time = ctime(&now);
  // convert now to string form
  // initialization of variables
  sscanf(argv[1], "%lf", &Pe);
  // Pe
  sscanf(argv[2], "%lf", &dt);
  // dt
  MD_total_step = static_cast<int>(atof(argv[3]) / dt);
  // MD_total_step
  sscanf(argv[4], "%lf", &density_prev);
  // desnsity
  Nm1p = atoi(argv[5]);
  is_polymer = (Nm1p > 1);
  // is_polymer = false;
  Kb = atof(argv[6]);
  ic = atoi(argv[7]);
  // @param num series number
  int num = atoi(argv[8]);

  step_t_1 = static_cast<int>(1. / dt);
  step_t_10 = static_cast<int>(10. / dt);
  step_t_100 = static_cast<int>(100. / dt);
  step_t_1000 = static_cast<int>(1000. / dt);

  cout << "ic = 0 (no reset), ic = 1 (reset) //  ic = " << ic << endl;
  cout << "total time = " << MD_total_step << ", Pe = " << Pe << ", "
       << "dt = " << dt << ", Nmp = " << Nm1p << ", Kb = " << Kb << endl;
  cout << "calculate polymere bond skip LJ force " << is_polymer << endl;
  cout << "open active force " << is_active << endl;

  if (ic == 0) {
    t_step_relax = static_cast<uint64_t>((3000. - 2. * Pe) / dt);
    cout << "t_step_relax = " << t_step_relax << endl;
    // TODO change there to change the relax time
    // t_step_relax = 0.;
  } else {
    t_step_relax = 0.1 * MD_total_step;
  }

  Pe_over_Delta = Pe / Delta;
  if (is_active) {
    beta = 1. / (1. + 0.5 * alpha * dt);
  } else {
    beta = 1;
  }

  b_dt = dt * beta;
  half_dt = 0.5 * dt;
  b_half_dt = 0.5 * dt * beta;
  b_half_dt2 = 0.5 * dt * dt * beta;

  // FREQ_CFG = (MD_total_step - t_step_relax) / NUM_CFG;
  // FREQ_CFG = static_cast<int>(cfg_dt / dt);
  FREQ_EVOL = static_cast<int>(0.1 / dt);
  // TODO FREQ_CFG_vcorr = static_cast<int>((0.2 / (2. * Diff_rot)) / dt);
  // TODO limit_vcorr = t_step_relax_vcorr + FREQ_CFG_vcorr * NUM_CFG_vcorr;
  MD_total_step_over100 = MD_total_step / 100;
  t_step_relax_over100 = t_step_relax / 100;

  eps = Pe_over_Delta;
  if (Pe <= 1.0 || !is_active) {
    eps = 1.0;
  }
  bta = sqrt(2. * alpha * dt);
  // thermal noise
  Ito = sqrt(2. * Diff_rot * dt);
  Kp = 1000. + 500. * Pe;
  // Gaussian chain constant (50+2*Pe)*1000/(Pe/Delta)
  Np = static_cast<int>(density_prev * (Lboxx * Lboxy) / (PI * Rh * Rh) /
                        static_cast<double>(Nm1p));
  // fit the Np form input density
  Nm = Np * Nm1p;
  density = Nm * PI * Rh * Rh / (Lboxx * Lboxy);
  initialization(Nm, Nneigh, Nc);

  sprintf(file_name,
          "ener_poly_s_%d_phi_%5.4f_Nm_%d_Nm1p_%d_Pe_%4.2f_stf_%4.2f.dat", num,
          density, Nm, Nm1p, Pe, Kb);
  ofstream ener_pol(file_name, ios::app);

  // sprintf(file_name,
  //         "active_pressure_s_%d_phi_%5.4f_Nm_%d_Nm1p_%d_Pe_%4.2f_stf_%4.2f.dat",
  //         num, density, Nm, Nm1p, Pe, Kb);
  // ofstream act_pres(file_name, ios::app);

  sprintf(file_name,
          "ReadMe_In_s_%d_phi_%5.4f_Nm_%d_Nm1p_%d_Pe_%4.2f_stf_%4.2f.txt", num,
          density, Nm, Nm1p, Pe, Kb);
  ofstream entrada(file_name, ios::app);

  sprintf(file_name,
          "ReadMe_Out_s_%d_phi_%5.4f_Nm_%d_Nm1p_%d_Pe_%4.2f_stf_%4.2f.txt", num,
          density, Nm, Nm1p, Pe, Kb);
  ofstream salida(file_name, ios::app);

  // sprintf(file_name,
  //         "cfg_vel_s%d_aV_1_phi_%5.4f_Nm_%d_Nm1p_%d_Pe_%4.2f_stf_%4.2f.dat",
  //         num, density, Nm, Nm1p, Pe, Kb);
  // sprintf(FILE_NAME,
  //         "cfg_s%d_phi_%5.4f_Nm_%d_Nm1p_%d_Pe_%4.2f_stf_%4.2f.h5",
  //         num, density, Nm, Nm1p, Pe, Kb);
  // ofstream cfg(file_name, ios::app);

  sprintf(file_name,
          "cfg_data_s%d_aV_%d_phi_%5.4f_Nm_%d_Nm1p_%d_Pe_%4.2f_stf_%4.2f.dat",
          num, aV_data, density, Nm, Nm1p, Pe, Kb);
  ofstream cfg_data_file(file_name, ios::app);
  cfg_data_file.precision(7);

  sprintf(file_name,
          "cfg_time_s%d_aV_%d_phi_%5.4f_Nm_%d_Nm1p_%d_Pe_%4.2f_stf_%4.2f.dat",
          num, aV_data, density, Nm, Nm1p, Pe, Kb);
  ofstream cfg_time_file(file_name, ios::app);
  cfg_time_file.precision(7);

  // endtoend.precision(5);

  entrada << "current version: BD_POLYMERS_2D_active" << endl;

  entrada << "density = " << density << ", eps = " << eps << endl;

  entrada << "Nm = " << Nm << " Nm1p = " << Nm1p << " Np = " << Np
          << ",  Kuhn_length = " << Lm << endl;

  entrada << "Lx = " << Lboxx << " Ly = " << Lboxy << ",  Lz = " << Lboxz
          << endl;

  entrada << "K_bond = " << Kp << " K_bend = " << Kb << ", Pe = " << Pe << endl;

  entrada << "gamma = " << alpha << " b = " << beta << ", DR = " << Diff_rot
          << endl;

  entrada << ", Delta = " << Delta << ", rot-noise Ito = " << Ito
          << ", MD time step, dt = " << dt << ", MD_iter = " << MD_total_step
          << endl;

  entrada << "sigma = " << sigma << ", r_cutLJ = " << rcutLJ
          << " thermal noise = " << bta << endl;

  entrada << "MD_time = " << setprecision(15)
          << static_cast<double>(MD_total_step * dt) << " t_step_relax = "
          << t_step_relax  // <<  " t_step_relax_vcor = " << t_step_relax_vcorr
          << endl;

  entrada
      // << "Cfg_dt" << cfg_dt << ", Freq_cfg = " << FREQ_CFG
      << ", Freq_evol = " << FREQ_EVOL << endl;

  entrada << "Running from: " << real_time << endl;

  for (int i = 0; i < Nm; ++i) {
    f_x[i] = 0.;
    f_y[i] = 0.;
    f_z[i] = 0.;
  }

  if (ic == 0) {
    // new_init_pol();
    /**
     * TODO
     * There is some problem in new_init_pol
     * cause end_to_end 10 peak in 1
     */
    // initial position
    init_pol_straight();
    init_velocities();
    // init_activity();
    // init ex,ey,ez
    press_count = 0.;
    cfg_count = 0.;
    virial_act = 0.;
    virial_pass = 0.;
    aV = 1;
  } else if (ic == 1) {
    // TODO(chen): need to get aV
    read_init_file();  // initialize from a file
  }

  int aV_single_file = 0;
  cell_map();
  forces();

  // @brief the MD loop for relaxation
  for (MD_step = 0; MD_step < t_step_relax; ++MD_step) {
    // START MD-TIME EVOLUTION
    if ((MD_step % t_step_relax_over100) == 0) {
      cout << "relax: "
           << 100. * static_cast<double>(MD_step) /
                  static_cast<double>(t_step_relax)
           << " % , "
           << "MD_step = " << MD_step << endl;
    }
    mol_dyn();
  }
  // @brief the MD loop for data output
  for (MD_step = 0; MD_step <= MD_total_step; ++MD_step) {
    if ((MD_step % MD_total_step_over100) == 0) {
      cout << 100. * static_cast<double>(MD_step) /
                  static_cast<double>(MD_total_step)
           << " % , "
           << "MD_step = " << MD_step << endl;
    }
    mol_dyn();
    // pressure(act_pres);
    if ((MD_step % FREQ_EVOL) == 0) {
      evol(ener_pol);
      print_last_cfg();
    }

    // if ((MD_step % FREQ_CFG) == 0) {
    //   aV_single_file += 1;
    //   cfg_files(aV_single_file, num, cfg);
    // }
    /**
     * manipulation the fequce of cfg data output
     */
    if (MD_step == step_t_1) {
      FREQ_CFG_DETA = 10;
    }
    if (MD_step == step_t_10) {
      FREQ_CFG_DETA = 100;
    }
    if (MD_step == step_t_100) {
      FREQ_CFG_DETA = 1000;
    }
    if (MD_step == step_t_1000) {
      FREQ_CFG_DETA = 10000;
    }

    if (MD_step % FREQ_CFG_DETA == 0) {
      if (cfg_conter >= cfgs_per_file) {
        cfg_conter = 0;

        cfg_data_file.close();
        cfg_time_file.close();

        aV_data += 1;
        sprintf(file_name,
                "cfg_data_s%d_aV_%d_phi_%5.4f_Nm_%d_Nm1p_%d_Pe_%4.2f_stf_%4."
                "2f.dat",
                num, aV_data, density, Nm, Nm1p, Pe, Kb);
        cfg_data_file.open(file_name, ios::app);

        sprintf(file_name,
                "cfg_time_s%d_aV_%d_phi_%5.4f_Nm_%d_Nm1p_%d_Pe_%4.2f_stf_%4."
                "2f.dat",
                num, aV_data, density, Nm, Nm1p, Pe, Kb);
        cfg_time_file.open(file_name, ios::app);
      }

      cfg_data(cfg_data_file, cfg_time_file);
      cfg_conter += 1;
    }
  }

  salida << "time= " << time(NULL) - sec << endl;

  // close files
  entrada.close();
  // e_act.close();
  // vcorr.close();
  ener_pol.close();
  // endtoend.close();
  salida.close();
  // act_pres.close();
  close_pointers();
  // cfg.close();
  cfg_data_file.close();
  cfg_time_file.close();
  // kill memory space of the pointers

  return 0;
}
