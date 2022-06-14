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

#ifndef __AMMSOCKNLP_HPP__
#define __AMMSOCKNLP_HPP__

#include "IpTNLP.hpp"
#include "generationStats.hpp" //lädt das in AMMSoCK.m erzeugte File
#include <cmath> //lädt mathematische Opterationen,Funktionen,Konstanten etc
#include <iomanip>
#include <iostream> // zuständig für input und output (Bildschirmausgabe etc)

using namespace Ipopt;

class ammsockNLP : public TNLP {
public:
  ammsockNLP();
  virtual ~ammsockNLP(); // destructor, zerstört die Funktion

  virtual bool get_nlp_info(Index &n, Index &m, Index &nnz_jac_g,
                            Index &nnz_h_lag, IndexStyleEnum &index_style);
  virtual bool setParameter(Number *p, Number *init, Number _hcoll,
                            Number _hFixed, Number _scale,
                            Number _consAtom[NATOM]);
  virtual void checkMethods();
  virtual void getSolution(Number *x);
  virtual bool get_bounds_info(Index n, Number *x_l, Number *x_u, Index m,
                               Number *g_l, Number *g_u);
  virtual bool get_starting_point(Index n, bool init_x, Number *x, bool init_z,
                                  Number *z_L, Number *z_U, Index m,
                                  bool init_lambda, Number *lambda);
  virtual bool eval_f(Index n, const Number *x, bool new_x, Number &obj_value);
  virtual bool eval_grad_f(Index n, const Number *x, bool new_x,
                           Number *grad_f);
  virtual bool eval_g(Index n, const Number *x, bool new_x, Index m, Number *g);
  virtual bool eval_jac_g(Index n, const Number *x, bool new_x, Index m,
                          Index nele_jac, Index *iRow, Index *jCol,
                          Number *values);
  virtual bool eval_h(Index n, const Number *x, bool new_x, Number obj_factor,
                      Index m, const Number *lambda, bool new_lambda,
                      Index nele_hess, Index *iRow, Index *jCol,
                      Number *values);

  virtual void rates(Number T, Number kf[NREAC], Number kr[NREAC],
                     Number dTkf[NREAC], Number dTkr[NREAC],
                     Number d2Tkf[NREAC], Number d2Tkr[NREAC],
                     Number Hbar[NSPEC], Number dTHbar[NSPEC],
                     Number d2THbar[NSPEC]);

  virtual void get_f(const Number T, const Number Y[NOP], Number f[NOP]);

  virtual void get_f(const Number T, const Number Y[NOP], Number f[NOP],
                     Number H[NSPEC], Number Yres[NNOP]);

  virtual void get_f(const Number T, const Number Y[NOP], Number f[NOP],
                     Number &rho, Number &Mbar, Number q[NREAC],
                     Number qtilde[NREAC], Number omega[NOP], Number M[NTB],
                     Number Gp[NREAC], Number Gm[NREAC], Number Rp[NREAC],
                     Number Rm[NREAC], Number kf[NREAC], Number dTkf[NREAC],
                     Number dTTkf[NREAC], Number kr[NREAC], Number dTkr[NREAC],
                     Number dTTkr[NREAC], Number H[NSPEC], Number dTH[NSPEC],
                     Number dTTH[NSPEC]);

  virtual void getFirstDerivative_f(const Number T, const Number Y[NOP],
                                    Number f[NOP], Number dYf[NOP][NOP],
                                    Number dTf[NOP], Number H[NSPEC],
                                    Number dTH[NSPEC]);

  virtual void getFirstDerivative_f(const Number T, const Number Y[NOP],
                                    Number f[NOP], Number dYf[NOP][NOP],
                                    Number dTf[NOP], Number H[NSPEC],
                                    Number dTH[NSPEC], Number Yres[NNOP],
                                    Number dYYres[NNOP][NOP]);

  virtual void getFirstDerivative_f(
      const Number T, const Number Y[NOP], Number f[NOP], Number dYf[NOP][NOP],
      Number dTf[NOP], Number &rho, Number dYrho[NOP], Number &dTrho,
      Number &Mbar, Number dYMbar[NOP], Number q[NREAC], Number dYq[NREAC][NOP],
      Number dTq[NREAC], Number qtilde[NREAC], Number dYqtilde[NREAC][NOP],
      Number dTqtilde[NREAC], Number omega[NOP], Number dYomega[NOP][NOP],
      Number dTomega[NOP], Number M[NTB], Number dYM[NTB][NOP], Number dTM[NTB],
      Number Gp[NREAC], Number dYGp[NREAC][NOP], Number Gm[NREAC],
      Number dYGm[NREAC][NOP], Number Rp[NREAC], Number dYRp[NREAC][NOP],
      Number Rm[NREAC], Number dYRm[NREAC][NOP], Number kf[NREAC],
      Number dTkf[NREAC], Number dTTkf[NREAC], Number kr[NREAC],
      Number dTkr[NREAC], Number dTTkr[NREAC], Number H[NSPEC],
      Number dTH[NSPEC], Number dTTH[NSPEC]);

  virtual void getSecondDerivative_f(const Number T, const Number Y[NOP],
                                     Number f[NOP], Number dYf[NOP][NOP],
                                     Number dTf[NOP],
                                     Number dYYf[NOP][NOP][NOP],
                                     Number dYTf[NOP][NOP], Number H[NSPEC],
                                     Number dTH[NSPEC], Number dTTH[NSPEC]);

  virtual void getSecondDerivative_f(
      const Number T, const Number Y[NOP], Number f[NOP], Number dYf[NOP][NOP],
      Number dTf[NOP], Number dYYf[NOP][NOP][NOP], Number dYTf[NOP][NOP],
      Number dTTf[NOP], Number H[NSPEC], Number dTH[NSPEC], Number dTTH[NSPEC],
      Number Yres[NNOP], Number dYYres[NNOP][NOP]);

  virtual void getSecondDerivative_f(
      const Number T, const Number Y[NOP], Number f[NOP], Number dYf[NOP][NOP],
      Number dTf[NOP], Number dYYf[NOP][NOP][NOP], Number dYTf[NOP][NOP],
      Number dTTf[NOP], Number &rho, Number dYrho[NOP], Number &dTrho,
      Number dYYrho[NOP][NOP], Number dYTrho[NOP], Number &dTTrho, Number &Mbar,
      Number dYMbar[NOP], Number dYYMbar[NOP][NOP], Number q[NREAC],
      Number dYq[NREAC][NOP], Number dTq[NREAC], Number dYYq[NREAC][NOP][NOP],
      Number dYTq[NREAC][NOP], Number dTTq[NREAC], Number qtilde[NREAC],
      Number dYqtilde[NREAC][NOP], Number dTqtilde[NREAC],
      Number dYYqtilde[NREAC][NOP][NOP], Number dYTqtilde[NREAC][NOP],
      Number dTTqtilde[NREAC], Number omega[NOP], Number dYomega[NOP][NOP],
      Number dTomega[NOP], Number dYYomega[NOP][NOP][NOP],
      Number dYTomega[NOP][NOP], Number dTTomega[NOP], Number M[NTB],
      Number dYM[NTB][NOP], Number dTM[NTB], Number dYYM[NTB][NOP][NOP],
      Number dYTM[NTB][NOP], Number dTTM[NTB], Number Gp[NREAC],
      Number dYGp[NREAC][NOP], Number dYYGp[NREAC][NOP][NOP], Number Gm[NREAC],
      Number dYGm[NREAC][NOP], Number dYYGm[NREAC][NOP][NOP], Number Rp[NREAC],
      Number dYRp[NREAC][NOP], Number dYYRp[NREAC][NOP][NOP], Number Rm[NREAC],
      Number dYRm[NREAC][NOP], Number dYYRm[NREAC][NOP][NOP], Number kf[NREAC],
      Number dTkf[NREAC], Number dTTkf[NREAC], Number kr[NREAC],
      Number dTkr[NREAC], Number dTTkr[NREAC], Number H[NSPEC],
      Number dTH[NSPEC], Number dTTH[NSPEC]);

  virtual void
  getThirdDerivative_f(const Number T, const Number Y[NOP], Number f[NOP],
                       Number dYf[NOP][NOP], Number dTf[NOP],
                       Number dYYf[NOP][NOP][NOP], Number dYTf[NOP][NOP],
                       Number dTTf[NOP], Number dYYYf[NOP][NOP][NOP][NOP],
                       Number dYYTf[NOP][NOP][NOP], Number dYTTf[NOP][NOP],
                       Number H[NSPEC], Number dTH[NSPEC], Number dTTH[NSPEC],
                       Number Yres[NNOP], Number dYYres[NNOP][NOP]);

  virtual void getThirdDerivative_f(
      const Number T, const Number Y[NOP], Number f[NOP], Number dYf[NOP][NOP],
      Number dTf[NOP], Number dYYf[NOP][NOP][NOP], Number dYTf[NOP][NOP],
      Number dTTf[NOP], Number dYYYf[NOP][NOP][NOP][NOP],
      Number dYYTf[NOP][NOP][NOP], Number dYTTf[NOP][NOP], Number &rho,
      Number dYrho[NOP], Number &dTrho, Number dYYrho[NOP][NOP],
      Number dYTrho[NOP], Number &dTTrho, Number dYYYrho[NOP][NOP][NOP],
      Number &Mbar, Number dYMbar[NOP], Number dYYMbar[NOP][NOP],
      Number q[NREAC], Number dYq[NREAC][NOP], Number dTq[NREAC],
      Number dYYq[NREAC][NOP][NOP], Number dYTq[NREAC][NOP], Number dTTq[NREAC],
      Number dYYYq[NREAC][NOP][NOP][NOP], Number dYYTq[NREAC][NOP][NOP],
      Number dYTTq[NREAC][NOP], Number qtilde[NREAC],
      Number dYqtilde[NREAC][NOP], Number dTqtilde[NREAC],
      Number dYYqtilde[NREAC][NOP][NOP], Number dYTqtilde[NREAC][NOP],
      Number dTTqtilde[NREAC], Number dYYYqtilde[NREAC][NOP][NOP][NOP],
      Number dYYTqtilde[NREAC][NOP][NOP], Number dYTTqtilde[NREAC][NOP],
      Number omega[NOP], Number dYomega[NOP][NOP], Number dTomega[NOP],
      Number dYYomega[NOP][NOP][NOP], Number dYTomega[NOP][NOP],
      Number dTTomega[NOP], Number dYYYomega[NOP][NOP][NOP][NOP],
      Number dYYTomega[NOP][NOP][NOP], Number dYTTomega[NOP][NOP],
      Number M[NTB], Number dYM[NTB][NOP], Number dTM[NTB],
      Number dYYM[NTB][NOP][NOP], Number dYTM[NTB][NOP], Number dTTM[NTB],
      Number dYYYM[NTB][NOP][NOP][NOP], Number dYYTM[NTB][NOP][NOP],
      Number dYTTM[NTB][NOP], Number Gp[NREAC], Number dYGp[NREAC][NOP],
      Number dYYGp[NREAC][NOP][NOP], Number dYYYGp[NREAC][NOP][NOP][NOP],
      Number Gm[NREAC], Number dYGm[NREAC][NOP], Number dYYGm[NREAC][NOP][NOP],
      Number dYYYGm[NREAC][NOP][NOP][NOP], Number Rp[NREAC],
      Number dYRp[NREAC][NOP], Number dYYRp[NREAC][NOP][NOP],
      Number dYYYRp[NREAC][NOP][NOP][NOP], Number Rm[NREAC],
      Number dYRm[NREAC][NOP], Number dYYRm[NREAC][NOP][NOP],
      Number dYYYRm[NREAC][NOP][NOP][NOP], Number kf[NREAC], Number dTkf[NREAC],
      Number dTTkf[NREAC], Number kr[NREAC], Number dTkr[NREAC],
      Number dTTkr[NREAC], Number H[NSPEC], Number dTH[NSPEC],
      Number dTTH[NSPEC]);

  virtual void finalize_solution(SolverReturn status, Index n, const Number *x,
                                 const Number *z_L, const Number *z_U, Index m,
                                 const Number *g, const Number *lambda,
                                 Number obj_value, const IpoptData *ip_data,
                                 IpoptCalculatedQuantities *ip_cq);

private:
  ammsockNLP(const ammsockNLP &);
  ammsockNLP &operator=(const ammsockNLP &);

  Number sol[NOP + 1];
  Number initialGuess[(NSTAGES + 1) * (NOP + 1) + NSTAGES * NINT * (NOP + 1)];

  Number Iop[NOP] = IOP;
  Number Inop[NNOP] = INOP;

  Number rpv_fixed[NRPV];
  Index rpv_index[NRPV] = INDRPV;

  Number h;
  Number hFixed;
  Number consAtom[NATOM];
  Number scale;
};

#endif
