//
// Created by marcocera on 24/02/22.
//


#ifndef TRANSITION_MY_TRANS_SOURCES_HPP
#define TRANSITION_MY_TRANS_SOURCES_HPP

#endif //TRANSITION_MY_TRANS_SOURCES_HPP   // questi 3 li ha aggiunti automaticamente Clion, non so cosa siano


#pragma once

#include "../scalar/scalar_sources.hpp"

template <class FlowIndices>
class CSourcePiecewise_TransLM final : public CNumerics {
private:

    const FlowIndices idx;

 // costanti sono nel variable

    su2double gammaAmb,
              rethetaAmb;

    su2double Residual[2],
              *Jacobian_i[2] = {nullptr},
              Jacobian_Buffer[4] = {0.0};

    su2double F_length_i,
              F_onset_i,
              F_thetat_i,
              F_turb_i,
              rethetat_eq_i,
              T_param_i,
              F_wake_i,
              delta_param_i,
              U_mag_i;

    su2double F_length_j,
              F_onset_j,
              F_thetat_j,
              F_turb_j,
              rethetat_eq_j,
              T_param_j,
              F_wake_j,
              delta_param_j,
              U_mag_j;

    su2double c_a1,
              c_a2,
              c_e1,
              c_e2,
              c_thetat;



    bool incompressible;

    // eventuali parti aggiuntive

public:
    /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   CSourcePiecewise_TransLM(unsigned short val_nDim, unsigned short val_nVar,
                               const su2double *constants, su2double val_gamma_Inf,
                               su2double val_rethetat_Inf, const CConfig* config):
                               CNumerics(val_nDim, val_nVar, config),
                               idx(val_nDim, config->GetnSpecies()),
                               incompressible(config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE),
                                                                               gammaAmb(val_gamma_Inf),
                                                                               rethetaAmb(val_rethetat_Inf),
                                                                               c_a1(constants[0]),
                                                                               c_a2(constants[1]),
                                                                               c_e1(constants[2]),
                                                                               c_e2(constants[3]),
                                                                               c_thetat(constants[4]) {
     /*--- "Allocate" the Jacobian using the static buffer. ---*/
     Jacobian_i[0] = Jacobian_Buffer;
     Jacobian_i[1] = Jacobian_Buffer + 2;

   }


    /*!
    * \brief Residual for source term integration.
    * \param[in] config - Definition of the particular problem.
    * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
    */
    ResidualType<> ComputeResidual(const CConfig* config) override{

//      if(config->dummyVar == 0)
//        cout << "trans_sources::ComputeResidual" << endl;

      unsigned short iDim;
      // quelli in comune con SST non li ridefinisco => ci sarà da collegarsi a quello quindi
      su2double P_gamma, P_rethetat, D_gamma;
      su2double VorticityMag = sqrt(Vorticity_i[0]*Vorticity_i[0] +
                                    Vorticity_i[1]*Vorticity_i[1] +
                                    Vorticity_i[2]*Vorticity_i[2]);


      Density_i = V_i[idx.Density()];
      Laminar_Viscosity_i = V_i[idx.LaminarViscosity()];
      Eddy_Viscosity_i = V_i[idx.EddyViscosity()];



      Residual[0] = 0.0;            Residual[1] = 0.0;
      Jacobian_i[0][0] = 0.0;       Jacobian_i[0][1] = 0.0;
      Jacobian_i[1][0] = 0.0;       Jacobian_i[1][1] = 0.0;


    if (dist_i > 1e-10) {
      // Questo deve essere preso da nodes. è già una variabile dentro questa classe
      //    su2double U_mag;
      //    su2double T_param = 500 * Laminar_Viscosity_i / (Density_i * pow(U_mag, 2));

      /*--- Production ---*/
      P_gamma = F_length_i * c_a1 * Density_i * StrainMag_i * pow(ScalarVar_i[0] * F_onset_i, 0.5) *
                (1 - c_e1 * ScalarVar_i[0]);
      P_rethetat = c_thetat * Density_i / T_param_i * (rethetat_eq_i - ScalarVar_i[1]) * (1.0 - F_thetat_i);
      D_gamma = c_a2 * Density_i * VorticityMag * ScalarVar_i[0] * F_turb_i * (c_e2 * ScalarVar_i[0] - 1);

      int iPoint2print = 22;
//      if(config->dummyVar == iPoint2print) {
//        cout << "Dentro source residuals" << endl;
//        cout << "c_a2 = " << c_a2 << endl;
//        cout << "Density_i = " << Density_i << endl;
//        cout << "VorticityMag = " << VorticityMag << endl;
//        cout << "F_turb_i = " << F_turb_i << endl;
//        cout << "P_gamma = " << P_gamma << endl;
//        cout << "D_gamma = " << D_gamma << endl;
//        cout << "P_rethetat = " << P_rethetat << endl;
//      }
//      cout << "F_turb_i = " << F_turb_i << endl;


      /*--- Add the production terms to the residuals ---*/
      Residual[0] += P_gamma * Volume;
      Residual[0] -= D_gamma * Volume;
      Residual[1] += P_rethetat * Volume;

      /*--- Crossflow extension ---*/
      // DA QUI

      /*--- Implicit part ---*/

      // Allora, ScalarVar sono direttamente gamma e Retheta.
      // lo jacobiano va fatto rispetto alle variabili conservative rho*gamma e rho*Retheta!
      // lo jacobiano è approssimato, quindi se il termine di produzione ha derivate brutte lo lasciamo stare.
      // l'importante è che il residuo sia giusto!
      Jacobian_i[0][0] = (F_length_i*c_a1*StrainMag_i*sqrt(F_onset_i)*(0.5*pow(ScalarVar_i[0], -0.5) -1.5*c_e1*pow(ScalarVar_i[0], 0.5))
                          - c_a2 * VorticityMag*F_turb_i*(2.0*c_e2*ScalarVar_i[0]-1.0) )*Volume;
      Jacobian_i[0][1] = 0.0;
      Jacobian_i[1][0] = 0.0;
      Jacobian_i[1][1] = -c_thetat/T_param_i*(1.0-F_thetat_i)*Volume;


      //Termine [1][0] è la derivata del termine di produzione (non ha distruzione) di Retheta rispetto a rho*gamma.
      //Ha dei livelli di approssimazione
      int approxLevel = 2;
//      if (approxLevel == 1)
//        Jacobian_i[1][0] = 0.0;
//      else if (approxLevel == 2) {
//        su2double DF_theta_t_D_rhoGamma =
//            -2 * (c_e2 / Density_i) * (c_e2 * ScalarVar_i[0] - 1) / pow((c_e2 - 1.0), 2.0);
//        Jacobian_i[1][0] =
//            -DF_theta_t_D_rhoGamma * c_thetat * (Density_i / T_param_i) * (rethetat_eq_i - ScalarVar_i[1]);
//      }

      // Termine [1][1] è la derivata del termine di produzione (non ha distruzione) di Retheta rispetto a rho*Retheta. Ha 3 livelli di approssimazione
//      approxLevel = 3;
//      Jacobian_i[1][1] = 0.0;
//      if (approxLevel == 1)
//        Jacobian_i[1][1] = 0.0;
//      else if (approxLevel == 2) {
//        Jacobian_i[1][1] = -c_thetat * (1.0 - F_thetat_i) * Volume/ T_param_i;
//      } else {
//        Jacobian_i[1][1] = -c_thetat * (1.0 - F_thetat_i) * Volume / T_param_i;
//        su2double FirstMaxTerm = F_wake_i * exp(-pow(dist_i / delta_param_i, 4));
//        su2double SecondMaxTerm = 1.0 - pow((c_e2 * ScalarVar_i[0] - 1) / (c_e2 - 1), 2);
//        // Only if true then the derivative of F_theta_t wrt rho*Re_theta_t is different from 0
//        if (FirstMaxTerm > SecondMaxTerm && FirstMaxTerm < 1.0) {
//          su2double additiveTerm = c_thetat * Density_i * ScalarVar_i[1] / T_param_i;
//          su2double tmp = dist_i / delta_param_i;
//          su2double F_theta_t_derivative =
//              4.0 * pow(tmp, 3) * F_wake_i * exp(-pow(tmp, 4.0)) * tmp / (Density_i * ScalarVar_i[1]);
//          additiveTerm *= F_theta_t_derivative;
//          Jacobian_i[1][1] += additiveTerm * Volume;
//        }
//      }
    }

      return ResidualType<>(Residual, Jacobian_i, nullptr);

   }

    inline void SetF_length(su2double val_F_length_i, su2double val_F_length_j) override {
      F_length_i = val_F_length_i;
      F_length_j = val_F_length_j;
    }

    inline void SetF_onset(su2double val_F_onset_i, su2double val_F_onset_j) override {
      F_onset_i = val_F_onset_i;
      F_onset_j = val_F_onset_j;
    }

    inline void SetF_thetat(su2double val_F_thetat_i, su2double val_F_thetat_j) override {
      F_thetat_i = val_F_thetat_i;
      F_thetat_j = val_F_thetat_j;
    }

    inline void SetF_turb(su2double val_F_turb_i, su2double val_F_turb_j) override {
      F_turb_i = val_F_turb_i;
      F_turb_j = val_F_turb_j;
    }

    inline void SetF_wake(su2double val_F_wake_i, su2double val_F_wake_j) override {
      F_wake_i = val_F_wake_i;
      F_wake_j = val_F_wake_j;
    }

    inline void Setrethetat_eq(su2double val_rethetat_eq_i, su2double val_rethetat_eq_j) override {
      rethetat_eq_i = val_rethetat_eq_i;
      rethetat_eq_j = val_rethetat_eq_j;
    }

    inline void SetT_param(su2double val_T_param_i, su2double val_T_param_j) override {
      T_param_i = val_T_param_i;
      T_param_j = val_T_param_j;
    }

    inline void Setdelta_param(su2double val_delta_param_i, su2double val_delta_param_j) override {
      delta_param_i = val_delta_param_i;
      delta_param_j = val_delta_param_j;
    }

    inline void SetU_mag(su2double val_U_mag_i, su2double val_U_mag_j) override {
      U_mag_i = val_U_mag_i;
      U_mag_j = val_U_mag_j;
    }

};
