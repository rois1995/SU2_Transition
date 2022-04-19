//
// Created by marcocera on 25/02/22.
//

#pragma once

#include "CTurbVariable.hpp"

class CTransLMVariable final : public CTurbVariable {
 protected:  // elencati in base all'ordine del sito NASA
  VectorType F_length, F_onset;
  VectorType F_turb;
  VectorType F_onset1, reV, F_onset2, R_T, F_onset3, F_length1, F_sublayer, rew, rethetac;
  VectorType T, rethetat_eq, F_thetat;
  VectorType Velocity_Mag, delta_param, F_wake, lambda_theta, Turb_Intens, du_dx, du_dy, du_dz, dU_ds, F_lambda, thetat;
  VectorType gamma_sep, gamma_eff;


 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] intermittency - Intermittency (gamma) (initialization value).
   * \param[in] Re_theta - theta Reynolds number (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] constants   // for now it is off
   * \param[in] config - Definition of the particular problem.
   */
  CTransLMVariable(su2double intermittency, su2double Re_theta, const su2double* constants, unsigned long npoint,
                   unsigned long ndim, unsigned long nvar, /*const su2double* constants*/ CConfig* config);

  /*!
     * \brief Destructor of the class.
   */
  ~CTransLMVariable() override = default;

  /*!
   * \brief Set the quantities for the terms calculation.
   * \param[in] val_viscosity - Value of the viscosity.
   * \param[in] val_dist - Value of the distance to the wall.
   * \param[in] val_density - Value of the density.
   * \param[in] val_desiredquantity - Value of the desired quantity.
   */
  void SetQuantities(unsigned long iPoint, su2double* constants, su2double val_viscosity, su2double val_dist, su2double val_density,
                     su2double val_vort, su2double StrainMag, CMatrixView<const su2double> Velocity_Gradient,
                     su2double* Velocity, su2double* TurbVars) override;

  // come ordine tengo quello in cui le ho definite nel .cpp -> probabilmente non tutti sono necessari
  /*!
   * \brief Get the various quantities.
   */
  inline su2double Getrew(unsigned long iPoint) const override { return rew(iPoint); }
  inline su2double GetF_length1(unsigned long iPoint) const override { return F_length1(iPoint); }
  inline su2double GetF_sublayer(unsigned long iPoint) const override { return F_sublayer(iPoint); }
  inline su2double GetF_length(unsigned long iPoint) const override { return F_length(iPoint); }
  inline su2double GetreV(unsigned long iPoint) const override { return reV(iPoint); }
  inline su2double GetF_rethetac(unsigned long iPoint) const override { return rethetac(iPoint); }
  inline su2double GetF_onset1(unsigned long iPoint) const override { return F_onset1(iPoint); }
  inline su2double GetR_T(unsigned long iPoint) const override { return R_T(iPoint); }
  inline su2double GetF_onset2(unsigned long iPoint) const override { return F_onset2(iPoint); }
  inline su2double GetF_onset3(unsigned long iPoint) const override { return F_onset3(iPoint); }
  inline su2double GetF_onset(unsigned long iPoint) const override { return F_onset(iPoint); }
  inline su2double GetF_turb(unsigned long iPoint) const override { return F_turb(iPoint); }
  inline su2double GetVelocity_Mag(unsigned long iPoint) const override { return Velocity_Mag(iPoint); }
  inline su2double GetT(unsigned long iPoint) const override { return T(iPoint); }
  inline su2double Getdu_dx(unsigned long iPoint) const override { return du_dx(iPoint); }
  inline su2double Getdu_dy(unsigned long iPoint) const override { return du_dy(iPoint); }
  inline su2double Getdu_dz(unsigned long iPoint) const override { return du_dz(iPoint); }
  inline su2double GetdU_ds(unsigned long iPoint) const override { return dU_ds(iPoint); }
  inline su2double GetTurb_Intens(unsigned long iPoint) const override { return Turb_Intens(iPoint); }
  inline su2double GetF_lambda(unsigned long iPoint) const override { return F_lambda(iPoint); }
  inline su2double Getrethetat_eq(unsigned long iPoint) const override { return rethetat_eq(iPoint); }
  inline su2double Getrethetac(unsigned long iPoint) const override { return rethetac(iPoint); }
  inline su2double Getthetat(unsigned long iPoint) const override { return thetat(iPoint); }
  inline su2double Getlambda_theta(unsigned long iPoint) const override { return lambda_theta(iPoint); }
  inline su2double Getdelta_param(unsigned long iPoint) const override { return delta_param(iPoint); }
  inline su2double GetF_wake(unsigned long iPoint) const override { return F_wake(iPoint); }
  inline su2double GetF_thetat(unsigned long iPoint) const override { return F_thetat(iPoint); }



  /*!
   * \brief Compute the correction for separation-induced transition.
   */
  inline void SetGammaSep(unsigned long iPoint, su2double val_gamma_sep) override {
    gamma_sep(iPoint) = val_gamma_sep;
  }

  /*!
   * \brief Correction for separation-induced transition.
   */
   // Non sono sicuro che sia il modo giusto. Forse conviene lasciare la soluzione come quella che è
  // ed usare un vettore di gammaEff?
  inline void SetGammaEff(unsigned long iPoint) override {
    gamma_eff(iPoint) = max(Solution(iPoint, 0), gamma_sep(iPoint));
  }

  inline void CorrectGamma(unsigned long iPoint) override {
    Solution(iPoint, 0) = gamma_eff(iPoint);
  }

  /*!
   * \brief Get intermittency value
   */
//  inline su2double GetIntermittency(unsigned long iPoint) const override { return Solution(iPoint, 0); }
  inline su2double GetIntermittency(unsigned long iPoint) const override { return gamma_eff(iPoint); }


  inline su2double GetGammaSep(unsigned long iPoint) const override { return gamma_sep(iPoint);}
  inline su2double GetGammaEff(unsigned long iPoint) const override { return gamma_eff(iPoint);}
};