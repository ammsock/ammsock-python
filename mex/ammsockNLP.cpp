/*
  Copyright (C) 2017 Julian Späth
  ----------------------------------------------------------------------------
  This file is part of AMMSoCK.

  AMMSoCK is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version. AMMSoCK is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
  details. You should have received a copy of the GNU General Public License
  along with AMMSoCK. If not, see <https://www.gnu.org/licenses/>.
*/

#include "ammsockNLP.hpp"

#include <cassert>

// include mech depending files
#include "auxiliaries.hpp"
#include "eval_f.hpp"
#include "eval_g.hpp"
#include "eval_grad_f.hpp"
#include "eval_h.hpp"
#include "eval_jac_g.hpp"
#include "get_bounds_info.hpp"
#include "get_nlp_info.hpp"
#include "rates.hpp"
#include "resolution.hpp"
#include <iostream>

using namespace Ipopt;

/* Constructor. */
ammsockNLP::ammsockNLP() // Methodenaufruf der Klasse ammsockNLP.hpp
{}

ammsockNLP::~ammsockNLP() {}
/*bool ammsockNLP::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // initialize to the given starting point
  x[0] = 1.0;
  x[1] = 5.0;
  x[2] = 5.0;
  x[3] = 1.0;

  return true;
}*/

bool ammsockNLP::get_starting_point(Index n, bool init_x, Number *x, //
                                    bool init_z, Number *z_L, Number *z_U,
                                    Index m, bool init_lambda, Number *lambda) {
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // set initial values for each stage
  for (int i = 0; i < (NSTAGES + 1) * (NOP + 1) + NINT * NSTAGES * (NOP + 1);
       ++i) {
    x[i] = initialGuess[i];
  }

  return true;
}

bool ammsockNLP::setParameter(Number *rpvfix, Number *init, Number _hcoll,
                              Number _hFixed, Number _scale,
                              Number _consAtom[NATOM]) {
  // set constants
  h = _hcoll;
  hFixed = _hFixed;
  scale = _scale;
  for (int i = 0; i < NATOM; ++i) {
    consAtom[i] = _consAtom[i];
  }

  // set rpv values
  for (int i = 0; i < NRPV; ++i) {
    rpv_fixed[i] = rpvfix[i];
  }

  // set constant inital guess for each stage
  for (int i = 0; i <= NSTAGES; ++i) {
    for (int s = 0; s <= NOP; ++s) {
      initialGuess[i * (NOP + 1) + s] = init[s];
    }
  }
  for (int k = 1; k <= NINT; ++k) {
    for (int i = 1; i <= NSTAGES; ++i) {
      for (int s = 0; s <= NOP; ++s) {
        initialGuess[k * (NSTAGES * (NOP + 1)) + i * (NOP + 1) + s] = init[s];
      }
    }
  }

  return true;
}

void ammsockNLP::getSolution(Number *x) {
  // set solution
  for (int i = 0; i < (NOP + 1); ++i) {
    x[i] = sol[i];
  }
}

void ammsockNLP::finalize_solution(SolverReturn status, Index n,
                                   const Number *x, const Number *z_L,
                                   const Number *z_U, Index m, const Number *g,
                                   const Number *lambda, Number obj_value,
                                   const IpoptData *ip_data,
                                   IpoptCalculatedQuantities *ip_cq) {
  // set solution
  //	std::cout << "Writing Solution: " << std::endl;
  for (int i = 0; i < (NOP + 1); ++i) {
    //		std::cout << "sol["<<i<<"] = x["<<
    //(NOP+1)+NINT*NSTAGES*(NOP+1)+(NSTAGES-1)*(NOP+1)+i<<"] = " <<
    //x[(NOP+1)+NINT*NSTAGES*(NOP+1)+(NSTAGES-1)*(NOP+1)+i] << std::endl;
    sol[i] = x[(NOP + 1) + NINT * NSTAGES * (NOP + 1) +
               (NSTAGES - 1) * (NOP + 1) + i];
  }
  /*	std::cout << std::endl << std::endl << "Solution of the primal
    variables, x" << std::endl; for (Index i=0; i<n; i++) { std::cout << "x[" <<
    i << "] = " << x[i] << std::endl;
    }
    std::cout << std::endl << std::endl << "Solution of the bound multipliers,
    z_L and z_U" << std::endl; for (Index i=0; i<n; i++) { std::cout << "z_L["
    << i << "] = " << z_L[i] << std::endl;
    }
    for (Index i=0; i<n; i++) {
      std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
    }

    std::cout << std::endl << std::endl << "Objective value" << std::endl;
    std::cout << "f(x*) = " << obj_value << std::endl;

    std::cout << std::endl << "Final value of the constraints:" << std::endl;
    for (Index i=0; i<m ;i++) {
      std::cout << "g(" << i << ") = " << g[i] << std::endl;
    }
  */
}
/*
void ammsockNLP::finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L,
const Number* z_U, Index m, const Number* g, const Number* lambda, Number
obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file,
etc
  // so we could use the solution.

  // For this example, we write the solution to the console
  std::cout << std::endl << std::endl << "Solution of the primal variables, x"
<< std::endl; for (Index i=0; i<n; i++) { std::cout << "x[" << i << "] = " <<
x[i] << std::endl;
  }

  std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L
and z_U" << std::endl; for (Index i=0; i<n; i++) { std::cout << "z_L[" << i <<
"] = " << z_L[i] << std::endl;
  }
  for (Index i=0; i<n; i++) {
    std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
  }

  std::cout << std::endl << std::endl << "Objective value" << std::endl;
  std::cout << "f(x*) = " << obj_value << std::endl;

  std::cout << std::endl << "Final value of the constraints:" << std::endl;
  for (Index i=0; i<m ;i++) {
    std::cout << "g(" << i << ") = " << g[i] << std::endl;
  }
}*/

void ammsockNLP::checkMethods() {
  std::cout << std::setprecision(15) << "Check Methods" << std::endl;

  int n = ((NINT + 1) * NSTAGES + 1) * (NOP + 1);
  int m = 0, nnz_jac_g = 0;
  if (NNOP > 0) {
    m = NRPV + (NINT + 1) * NSTAGES * (NOP + 1) + 1 +
        NNOP * ((NINT + 1) * NSTAGES + 1);
    nnz_jac_g = NSTAGES * NOP + (NOP + 1) + NINT * NSTAGES * NOP +
                (NINT + 1) * (NSTAGES * NSTAGES * NOP * (NOP + 1) +
                              (NOP + 1) * NSTAGES) +
                NRPV + (NSTAGES * (NINT + 1) + 1) * NOP * NNOP;
  } else {
    m = NRPV + (NINT + 1) * NSTAGES * (NOP + 1) + 1 + NATOM;
    nnz_jac_g = NSTAGES * NOP + (NOP + 1) + NINT * NSTAGES * NOP +
                (NINT + 1) * (NSTAGES * NSTAGES * NOP * (NOP + 1) +
                              (NOP + 1) * NSTAGES) +
                NRPV + NATOM * NOP;
  }
  int nnz_h_lag =
      (NOP + 1) * (NOP + 1) + (NINT + 1) * NSTAGES * (NOP + 1) * (NOP + 1);

  // allocate memory
  Number f[NOP], dYf[NOP][NOP], dTf[NOP], dYYf[NOP][NOP][NOP], dYTf[NOP][NOP],
      dTTf[NOP], dYYYf[NOP][NOP][NOP][NOP], dYYTf[NOP][NOP][NOP],
      dYTTf[NOP][NOP];
  Number rho, dYrho[NOP], dTrho, dYYrho[NOP][NOP], dYTrho[NOP], dTTrho,
      dYYYrho[NOP][NOP][NOP];
  Number Mbar, dYMbar[NOP], dYYMbar[NOP][NOP];
  Number q[NREAC], dYq[NREAC][NOP], dTq[NREAC], dYYq[NREAC][NOP][NOP],
      dYTq[NREAC][NOP], dTTq[NREAC], dYYYq[NREAC][NOP][NOP][NOP],
      dYYTq[NREAC][NOP][NOP], dYTTq[NREAC][NOP];
  Number qtilde[NREAC], dYqtilde[NREAC][NOP], dTqtilde[NREAC],
      dYYqtilde[NREAC][NOP][NOP], dYTqtilde[NREAC][NOP], dTTqtilde[NREAC],
      dYYYqtilde[NREAC][NOP][NOP][NOP], dYYTqtilde[NREAC][NOP][NOP],
      dYTTqtilde[NREAC][NOP];
  Number omega[NOP], dYomega[NOP][NOP], dTomega[NOP], dYYomega[NOP][NOP][NOP],
      dYTomega[NOP][NOP], dTTomega[NOP], dYYYomega[NOP][NOP][NOP][NOP],
      dYYTomega[NOP][NOP][NOP], dYTTomega[NOP][NOP];
  Number M[NTB], dYM[NTB][NOP], dTM[NTB], dYYM[NTB][NOP][NOP], dYTM[NTB][NOP],
      dTTM[NTB], dYYYM[NTB][NOP][NOP][NOP], dYYTM[NTB][NOP][NOP],
      dYTTM[NTB][NOP];
  Number Gp[NREAC], dYGp[NREAC][NOP], dYYGp[NREAC][NOP][NOP],
      dYYYGp[NREAC][NOP][NOP][NOP];
  Number Gm[NREAC], dYGm[NREAC][NOP], dYYGm[NREAC][NOP][NOP],
      dYYYGm[NREAC][NOP][NOP][NOP];
  Number Rp[NREAC], dYRp[NREAC][NOP], dYYRp[NREAC][NOP][NOP],
      dYYYRp[NREAC][NOP][NOP][NOP];
  Number Rm[NREAC], dYRm[NREAC][NOP], dYYRm[NREAC][NOP][NOP],
      dYYYRm[NREAC][NOP][NOP][NOP];
  Number kf[NREAC], dTkf[NREAC], dTTkf[NREAC];
  Number kr[NREAC], dTkr[NREAC], dTTkr[NREAC];
  Number H[NSPEC], dTH[NSPEC], dTTH[NSPEC];
  double Y[NOP], T;

  /* SMOOKE
          if (NNOP > 0) {
             Y[0] = 0.000436994141148;
             Y[1] = 0.000112506816684;
             Y[2] = 0.000443278880829;
             Y[3] = 0.139018933204699;
             Y[4] = 0.000045534053435;
             Y[5] = 0.002511084155766;
             Y[6] = 0.055881327051107;
             Y[7] = 0.001284231169616;
             Y[8] = 0.037605284785079;
             Y[9] = 0.015441699147586;
             Y[10] = 0.000000358318762;
             Y[11] = 0.000004064738122;
          } else {
                  Y[0] = 0.000436994141148;
                     Y[1] = 0.000112506816684;
                     Y[2] = 0.000443278880829;
                     Y[3] = 0.000079028348710;
                     Y[4] = 0.139018933204699;
                     Y[5] = 0.001821330761267;
                     Y[6] = 0.000045534053435;
                     Y[7] = 0.002511084155766;
                     Y[8] = 0.055881327051107;
                     Y[9] = 0.001284231169616;
                     Y[10] = 0.000043121715668;
                     Y[11] = 0.037605284785079;
                     Y[12] = 0.015441699147586;
                     Y[13] = 0.000000358318762;
                     Y[14] = 0.000004064738122;
                     Y[15] = 0.745271222711521;
          }
          T = 1.660206066877847e+03;
  */

  if (NNOP > 0) {
    Y[0] = 0.065160244954992;
    Y[1] = 0.010852820557299;
    Y[2] = 0.006976335346012;
    Y[3] = 0.000057664436232;
    Y[4] = 0.000000083461521;
    Y[5] = 0.154921223995219;
  } else {
    Y[0] = 0.013133867263872;
    Y[1] = 0.065160244954992;
    Y[2] = 0.004118491065228;
    Y[3] = 0.010852820557299;
    Y[4] = 0.006976335346012;
    Y[5] = 0.000057664436232;
    Y[6] = 0.000000083461521;
    Y[7] = 0.154921223995219;
    Y[8] = 0.744779415313207;
  }
  T = 1.620023072534618e+03;

  Number Yres[NNOP];

  get_f(T, Y, f, H, Yres);
  std::cout << "get_f:" << std::endl;
  std::cout << "\tf = [";
  for (int i = 0; i < NOP; ++i) {
    std::cout << std::setprecision(15) << f[i] << ", ";
  }
  std::cout << "]" << std::endl;
  std::cout << "\tH = [";
  for (int i = 0; i < 9; ++i) {
    std::cout << std::setprecision(15) << H[i] << ", ";
  }
  std::cout << "]" << std::endl;

  get_f(T, Y, f, rho, Mbar, q, qtilde, omega, M, Gp, Gm, Rp, Rm, kf, dTkf,
        dTTkf, kr, dTkr, dTTkr, H, dTH, dTTH);
  std::cout << "get_f(param):" << std::endl;
  std::cout << "\tf = [";
  for (int i = 0; i < NOP; ++i) {
    std::cout << std::setprecision(15) << f[i] << ", ";
  }
  std::cout << "]" << std::endl;

  std::cout << "\tkf = [";
  for (int i = 0; i < NREAC; ++i) {
    std::cout << std::setprecision(15) << kf[i] << ", ";
  }
  std::cout << "]" << std::endl;
  std::cout << "\tkr = [";
  for (int i = 0; i < NREAC; ++i) {
    std::cout << std::setprecision(15) << kr[i] << ", ";
  }
  std::cout << "]" << std::endl;

  getFirstDerivative_f(T, Y, f, dYf, dTf, H, dTH);

  std::cout << "getFirstDerivative_f:" << std::endl;
  std::cout << "\tf = [";
  for (int i = 0; i < NOP; ++i) {
    std::cout << std::setprecision(15) << f[i] << ", ";
  }
  std::cout << "]" << std::endl;

  std::cout << "\tdTf = [";
  for (int i = 0; i < NOP; ++i) {
    std::cout << std::setprecision(15) << dTf[i] << ", ";
  }
  std::cout << "]" << std::endl;

  std::cout << "\tdYf = [" << std::endl;
  for (int i = 0; i < NOP; ++i) {
    std::cout << "\t       [";
    for (int s = 0; s < NOP; ++s) {
      std::cout << std::setprecision(15) << dYf[i][s] << ", ";
    }
    std::cout << "]," << std::endl;
  }
  std::cout << "\t       ]" << std::endl;

  /*

          getFirstDerivative_f(T,Y,f,dYf,dTf,rho,dYrho,dTrho,Mbar,dYMbar,q,dYq,dTq,qtilde,dYqtilde,dTqtilde,omega,dYomega,dTomega,M,dYM,dTM,Gp,dYGp,Gm,dYGm,Rp,dYRp,Rm,dYRm,kf,dTkf,dTTkf,kr,dTkr,dTTkr,H,dTH,dTTH);

          std::cout << "getFirstDerivative_f(param):" << std::endl;

          std::cout << "\tf = [";
          for (int i = 0; i < NOP; ++i) {
                  std::cout << std::setprecision(15) << f[i] << ", ";
          }
          std::cout << "]" << std::endl;

          std::cout << "\tdTf = [";
          for (int i = 0; i < NOP; ++i) {
                  std::cout << std::setprecision(15) << dTf[i] << ", ";
          }
          std::cout << "]" << std::endl;

          std::cout << "\tdYGp = [" << std::endl;
          for (int i = 0; i < 28; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << std::setprecision(15) << dYGp[i][s] << ",
     ";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYGm = [" << std::endl;
          for (int i = 0; i < 28; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << std::setprecision(15) << dYGm[i][s] << ",
     ";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYf = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << std::setprecision(15) << dYf[i][s] << ",
     ";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;
  */

  getSecondDerivative_f(
      T, Y, f, dYf, dTf, dYYf, dYTf, dTTf, rho, dYrho, dTrho, dYYrho, dYTrho,
      dTTrho, Mbar, dYMbar, dYYMbar, q, dYq, dTq, dYYq, dYTq, dTTq, qtilde,
      dYqtilde, dTqtilde, dYYqtilde, dYTqtilde, dTTqtilde, omega, dYomega,
      dTomega, dYYomega, dYTomega, dTTomega, M, dYM, dTM, dYYM, dYTM, dTTM, Gp,
      dYGp, dYYGp, Gm, dYGm, dYYGm, Rp, dYRp, dYYRp, Rm, dYRm, dYYRm, kf, dTkf,
      dTTkf, kr, dTkr, dTTkr, H, dTH, dTTH);
  std::cout << "getSecondDerivative_f(param):" << std::endl;

  std::cout << "\tf = [";
  for (int i = 0; i < NOP; ++i) {
    std::cout << std::setprecision(15) << f[i] << ", ";
  }
  std::cout << "]" << std::endl;
  /*
          std::cout << "\tdTf = [";
          for (int i = 0; i < NOP; ++i) {
                  std::cout << std::setprecision(15) << dTf[i] << ", ";
          }
          std::cout << "]" << std::endl;

          std::cout << "\tdYYf = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << "[";
                                          for (int m = 0; m < NOP; ++m) {
                                                  std::cout <<
     std::setprecision(15) << dYYf[i][s][m] << ", ";
                                          }
                          std::cout << "],";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYTf = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << std::setprecision(15) << dYTf[i][s] << ",
     ";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdTTf = [";
          for (int i = 0; i < NOP; ++i) {
                          std::cout << std::setprecision(15) << dTTf[i] << ", ";
          }
          std::cout << " ]:" << std::endl;


          getSecondDerivative_f(T,Y,f,dYf,dTf,dYYf,dYTf,H,dTH,dTTH);
          std::cout << "getSecondDerivative_f:" << std::endl;

          std::cout << "\tf = [";
          for (int i = 0; i < NOP; ++i) {
                  std::cout << std::setprecision(15) << f[i] << ", ";
          }
          std::cout << "]" << std::endl;

          std::cout << "\tdTf = [";
          for (int i = 0; i < NOP; ++i) {
                  std::cout << std::setprecision(15) << dTf[i] << ", ";
          }
          std::cout << "]" << std::endl;

          std::cout << "\tdYf = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << std::setprecision(15) << dYf[i][s] << ",
     ";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYTf = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << std::setprecision(15) << dYTf[i][s] << ",
     ";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;


          std::cout << "\tdYYf = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << "[";
                                          for (int m = 0; m < NOP; ++m) {
                                                  std::cout <<
     std::setprecision(15) << dYYf[i][s][m] << ", ";
                                          }
                          std::cout << "],";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYTf = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << std::setprecision(20) << dYTf[i][s] << ",
     ";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;
  */
  /*
          getThirdDerivative_f(T,Y,f,dYf,dTf,dYYf,dYTf,dTTf,dYYYf,dYYTf,dYTTf,rho,dYrho,dTrho,dYYrho,dYTrho,dTTrho,dYYYrho,Mbar,dYMbar,dYYMbar,q,dYq,dTq,dYYq,dYTq,dTTq,dYYYq,dYYTq,dYTTq,qtilde,dYqtilde,dTqtilde,dYYqtilde,dYTqtilde,dTTqtilde,dYYYqtilde,dYYTqtilde,dYTTqtilde,omega,dYomega,dTomega,dYYomega,dYTomega,dTTomega,dYYYomega,dYYTomega,dYTTomega,M,dYM,dTM,dYYM,dYTM,dTTM,dYYYM,dYYTM,dYTTM,Gp,dYGp,dYYGp,dYYYGp,Gm,dYGm,dYYGm,dYYYGm,Rp,dYRp,dYYRp,dYYYRp,Rm,dYRm,dYYRm,dYYYRm,kf,dTkf,dTTkf,kr,dTkr,dTTkr,H,dTH,dTTH);
          std::cout << "getThirdDerivative_f:" << std::endl;

          std::cout << "\tf = [";
          for (int i = 0; i < NOP; ++i) {
                  std::cout << std::setprecision(15) << f[i] << ", ";
          }
          std::cout << "]" << std::endl;
          std::cout << "\tdYYf = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << "[";
                                          for (int m = 0; m < NOP; ++m) {
                                                  std::cout <<
     std::setprecision(15) << dYYf[i][s][m] << ", ";
                                          }
                          std::cout << "],";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;
          std::cout << "\tdYTf = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << std::setprecision(20) << dYTf[i][s] << ",
     ";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYYYrho = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << "[";
                          for (int k = 0; k < NOP; ++k) {
                                  std::cout << std::setprecision(20) <<
     dYYYrho[i][s][k] << ", ";
                          }
                          std::cout << "],";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYTTM = [" << std::endl;
          for (int i = 0; i < NTB; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << std::setprecision(20) << dYTTM[i][s] <<
     ", ";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYYTM = [" << std::endl;
          for (int i = 0; i < NTB; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << "[";
                          for (int k = 0; k < NOP; ++k) {
                                  std::cout << std::setprecision(20) <<
     dYYTM[i][s][k] << ", ";
                          }
                          std::cout << "],";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYYYM = [" << std::endl;
          for (int i = 0; i < NTB; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << "[";
                          for (int k = 0; k < NOP; ++k) {
                                  std::cout << "[";
                                  for (int m = 0; m < NOP; ++m) {
                                          std::cout << std::setprecision(20) <<
     dYYYM[i][s][k][m] << ", ";
                                  }
                                  std::cout << "],";
                          }
                          std::cout << "],";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYYYq = [" << std::endl;
          for (int i = 0; i < NREAC; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << "[";
                          for (int k = 0; k < NOP; ++k) {
                                  std::cout << "[";
                                  for (int m = 0; m < NOP; ++m) {
                                          std::cout << std::setprecision(20) <<
     dYYYq[i][s][k][m] << ", ";
                                  }
                                  std::cout << "],";
                          }
                          std::cout << "],";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;


          std::cout << "\tdYYTqtilde = [" << std::endl;
          for (int i = 0; i < NREAC; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << "[";
                          for (int k = 0; k < NOP; ++k) {
                                  std::cout << std::setprecision(20) <<
     dYYTqtilde[i][s][k] << ", ";
                          }
                          std::cout << "],";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYYTq = [" << std::endl;
          for (int i = 0; i < NREAC; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << "[";
                          for (int k = 0; k < NOP; ++k) {
                                  std::cout << std::setprecision(20) <<
     dYYTq[i][s][k] << ", ";
                          }
                          std::cout << "],";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYTTqtilde = [" << std::endl;
          for (int i = 0; i < NREAC; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << std::setprecision(20) << dYTTqtilde[i][s]
     << ", ";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYTTq = [" << std::endl;
          for (int i = 0; i < NREAC; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << std::setprecision(20) << dYTTq[i][s] <<
     ", ";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYYYomega = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << "[";
                          for (int k = 0; k < NOP; ++k) {
                                  std::cout << "[";
                                  for (int m = 0; m < NOP; ++m) {
                                          std::cout << std::setprecision(20) <<
     dYYYomega[i][s][k][m] << ", ";
                                  }
                                  std::cout << "],";
                          }
                          std::cout << "],";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYYTomega = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << "[";
                          for (int k = 0; k < NOP; ++k) {
                                  std::cout << std::setprecision(20) <<
     dYYTomega[i][s][k] << ", ";
                          }
                          std::cout << "],";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYTTomega = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << std::setprecision(20) << dYTTomega[i][s]
     << ", ";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYYYf = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << "[";
                          for (int k = 0; k < NOP; ++k) {
                                  std::cout << "[";
                                  for (int m = 0; m < NOP; ++m) {
                                          std::cout << std::setprecision(20) <<
     dYYYf[i][s][k][m] << ", ";
                                  }
                                  std::cout << "],";
                          }
                          std::cout << "],";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYYTf = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << "[";
                          for (int k = 0; k < NOP; ++k) {
                                  std::cout << std::setprecision(20) <<
     dYYTf[i][s][k] << ", ";
                          }
                          std::cout << "],";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

          std::cout << "\tdYTTf = [" << std::endl;
          for (int i = 0; i < NOP; ++i) {
                  std::cout << "\t       [";
                  for (int s = 0; s < NOP; ++s) {
                          std::cout << std::setprecision(20) << dYTTf[i][s] <<
     ", ";
                  }
                  std::cout << "]," << std::endl;
          }
          std::cout << "\t       ]" << std::endl;

  */

  /*
          Number x[NOP+1], grad_f[NOP+1], objvalue;
          for (int i = 0; i < NOP; ++i) {
                  x[i] = Y[i];
          }
          x[NOP] = T;
  */

  // define distinct points
  Number x[n], grad_f[n], objvalue;
  x[0] = 0.06516036588982107580;
  x[1] = 0.01085284247054557269;
  x[2] = 0.00697634938291600420;
  x[3] = 0.00005766459509550730;
  x[4] = 0.00000008346158013730;
  x[5] = 0.15492092560675799207;
  x[6] = 1620.02298516142309381394;
  x[7] = 0.04532604475926410581;
  x[8] = 0.01480552101610078247;
  x[9] = 0.00540273399912655997;
  x[10] = 0.00003142121870769591;
  x[11] = 0.00000009574191704357;
  x[12] = 0.17176129294499814515;
  x[13] = 1815.97265504894630794297;
  x[14] = 0.04254313771312739911;
  x[15] = 0.01599351787506801681;
  x[16] = 0.00540133550526537744;
  x[17] = 0.00002603777546671426;
  x[18] = 0.00000011085095364815;
  x[19] = 0.17468109656927380269;
  x[20] = 1899.43227417249704558344;
  x[21] = 0.03909048175968590777;
  x[22] = 0.01726836695536585614;
  x[23] = 0.00534646274375098155;
  x[24] = 0.00002022372991004300;
  x[25] = 0.00000012916388470504;
  x[26] = 0.17872295152508163585;
  x[27] = 2005.45913449510567261314;
  x[28] = 0.02555983740805415835;
  x[29] = 0.01868100153044827474;
  x[30] = 0.00450167727988897998;
  x[31] = 0.00000684772032255288;
  x[32] = 0.00000018577511168378;
  x[33] = 0.19906508789617133326;
  x[34] = 2418.15420995643762580585;
  x[35] = 0.02457541228125392427;
  x[36] = 0.01851006364284429537;
  x[37] = 0.00440268847115802848;
  x[38] = 0.00000646687510630283;
  x[39] = 0.00000019381770471409;
  x[40] = 0.20081310513042407995;
  x[41] = 2447.12082114346139860572;
  x[42] = 0.02191069135995103029;
  x[43] = 0.01784495364881045001;
  x[44] = 0.00411143116296858787;
  x[45] = 0.00000582030238549453;
  x[46] = 0.00000022828387810440;
  x[47] = 0.20570704296647226994;
  x[48] = 2524.27209875927246685023;
  x[49] = 0.02110864963521254592;
  x[50] = 0.01758688358972280985;
  x[51] = 0.00401733080476635637;
  x[52] = 0.00000573043993633865;
  x[53] = 0.00000024375265659530;
  x[54] = 0.20722420125621648923;
  x[55] = 2547.09738519242955590016;
  x[56] = 0.01997993609471254359;
  x[57] = 0.01717866957162951513;
  x[58] = 0.00388005836546617298;
  x[59] = 0.00000568046479168769;
  x[60] = 0.00000027096151863868;
  x[61] = 0.20939186604388113078;
  x[62] = 2578.87968521139146105270;
  x[63] = 0.01768425225683526586;
  x[64] = 0.01618870067869662699;
  x[65] = 0.00358411723730947497;
  x[66] = 0.00000582880715911408;
  x[67] = 0.00000035192278633877;
  x[68] = 0.21391016468008355877;
  x[69] = 2642.19575649198031896958;

  eval_f(n, x, true, objvalue);
  eval_grad_f(n, x, true, grad_f);

  std::cout << "eval_f:" << std::endl;
  std::cout << "\tobjValue = " << objvalue << std::endl;
  std::cout << std::endl;

  std::cout << "eval_grad_f:" << std::endl;
  std::cout << "\tgrad_f = [";
  for (int i = 0; i <= NOP; ++i) {
    std::cout << grad_f[i] << ", ";
  }
  std::cout << "]:" << std::endl;

  Number g[m], values[nnz_jac_g];
  Index iRow[nnz_jac_g], jCol[nnz_jac_g];

  /*for (int i = 0; i <= NSTAGES; ++i) {
          for (int s = 0; s <= NOP; ++s) {
                  x_g[i*(NOP+1)+s] = x[s];
  //			std::cout << "x_g["<<i*(NOP+1)+s<<"] = " <<
  x_g[i*(NOP+1)+s] << std::endl;
          }
  }

  for (int k = 1; k <= NINT; ++k) {
          for (int i = 1; i <= NSTAGES; ++i) {
                  for (int s = 0; s <= NOP; ++s) {
                          x_g[k*(NSTAGES*(NOP+1))+i*(NOP+1)+s] = x[s];
  //			std::cout << "x_g["<<k*(NSTAGES*(NOP+1))+i*(NOP+1)+s<<"] =
  " << x_g[k*(NSTAGES*(NOP+1))+i*(NOP+1)+s] << std::endl;
                  }
          }
  }
  */
  eval_g(n, x, true, m, g);
  eval_jac_g(n, x, true, m, nnz_jac_g, iRow, jCol, values);

  std::cout << "eval_g := [ ";
  for (int i = 0; i < m; ++i) {
    std::cout << std::setprecision(25) << g[i] << ", ";
  }
  std::cout << "]:" << std::endl;

  std::cout << "eval_jac_g := " << std::endl;
  for (int i = 0; i < nnz_jac_g; ++i) {
    std::cout << "values(" << i + 1 << ") = " << values[i] << ";" << std::endl;
    //		std::cout << "\tvalues["<<i<<"] = " << values[i] << std::endl;
  }
  std::cout << std::endl;

  Number valuesH[nnz_h_lag];
  Index iRowH[nnz_h_lag], jColH[nnz_h_lag];

  std::cout << "nnz_h_lag = " << nnz_h_lag << std::endl;

  Number lambda[m];
  for (int i = 0; i < m; ++i) {
    lambda[i] = i + 1;
  }

  std::cout << "lambda = [";
  for (int i = 0; i < m; ++i) {
    std::cout << "," << lambda[i];
  }
  std::cout << "]:" << std::endl;

  eval_h(n, x, true, 1.0, m, lambda, true, nnz_h_lag, iRowH, jColH, NULL);
  eval_h(n, x, true, 1.0, m, lambda, true, nnz_h_lag, iRowH, jColH, valuesH);

  /*
          std::cout << "iRowH, jColh = ";
          for (int i = 0; i < nnz_h_lag; ++i) {
                  std::cout << "iRow("<<i+1<<") = " << iRowH[i] << "; " <<
     "jCol("<<i+1<<") = " << jColH[i] << ";" << std::endl;
          }
          std::cout << std::endl;
  */

  std::cout << "eval_h:" << std::endl;
  for (int i = 0; i < nnz_h_lag; ++i) {
    std::cout << "values(" << i + 1 << ") = " << valuesH[i] << ";" << std::endl;
  }
  std::cout << std::endl;

  //	eval_h(n, x_g, true, 0.0, m, lambda, true, nnz_h_lag, iRowH, jColH,
  //valuesH);
}
