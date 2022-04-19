//
// Created by marcocera on 04/03/22.
//

#ifndef TRANSITION_MY_TRANS_DIFFUSION_HPP
#define TRANSITION_MY_TRANS_DIFFUSION_HPP

#endif //TRANSITION_MY_TRANS_DIFFUSION_HPP   // questi 3 li ha aggiunti automaticamente Clion, non so cosa siano

// nei relativi turb dice che tutto ciò che riguarda i viscous residual si può fare estendendo al caso interessato quelle classi!

#pragma once

#include "../CNumerics.hpp"
#include "../scalar/scalar_diffusion.hpp"

/*!
 * \class CAvgGrad_Scalar
 * \brief Template class for computing viscous residual of scalar values
 * \details This class serves as a template for the scalar viscous residual
 *   classes.  The general structure of a viscous residual calculation is the
 *   same for many different  models, which leads to a lot of repeated code.
 *   By using the template design pattern, these sections of repeated code are
 *   moved to a shared base class, and the specifics of each model
 *   are implemented by derived classes.  In order to add a new residual
 *   calculation for a viscous residual, extend this class and implement
 *   the pure virtual functions with model-specific behavior.
 * \ingroup ViscDiscr
 * \author C. Pederson, A. Bueno, and F. Palacios
 */
 // Aggiunto da me
 // Perchè la crei di nuovo? c'è già in scalar_diffusion.hpp
//class CAvgGrad_Scalar : public CNumerics {
//protected:
//    su2double
//    *Proj_Mean_GradTransVar_Normal = nullptr,  /*!< \brief Mean_gradTransVar DOT normal. */
//    *Proj_Mean_GradTransVar = nullptr,         /*!< \brief Mean_gradTransVar DOT normal, corrected if required. */
//    proj_vector_ij = 0.0,                     /*!< \brief (Edge_Vector DOT normal)/|Edge_Vector|^2 */
//    *Flux = nullptr,                          /*!< \brief Final result, diffusive flux/residual. */
//    **Jacobian_i = nullptr,                   /*!< \brief Flux Jacobian w.r.t. node i. */
//    **Jacobian_j = nullptr;                   /*!< \brief Flux Jacobian w.r.t. node j. */
//
//    const bool correct_gradient = false, implicit = false, incompressible = false;
//
//    // saltate le parti di AD
//
//    /*!
//   * \brief Model-specific steps in the ComputeResidual method, derived classes
//   *        should compute the Flux and Jacobians (i/j) inside this method.
//   * \param[in] config - Definition of the particular problem.
//   */
//    virtual void FinishResidualCalc(const CConfig* config) = 0;
//
//public:
///*!
// * \brief Constructor of the class.
// * \param[in] val_nDim - Number of dimensions of the problem.
// * \param[in] val_nVar - Number of variables of the problem.
// * \param[in] correct_gradient - Whether to correct gradient for skewness.
// * \param[in] config - Definition of the particular problem.
// */
//    CAvgGrad_Scalar(unsigned short val_nDim, unsigned short val_nVar, bool correct_gradient, const CConfig* config);
//
//    /*!
//     * \brief Destructor of the class.
//     */
//    ~CAvgGrad_Scalar(void) override;
//
//    /*!
//   * \brief Compute the viscous residual using an average of gradients without correction.
//   * \param[in] config - Definition of the particular problem.
//   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
//   */
//    ResidualType<> ComputeResidual(const CConfig* config) override;
//
//};

// saltato SA

// saltato SA_NEG

/*!
 * \class my_CAvgGrad_TransLM
 * \brief Class for computing viscous term using average of gradient with correction (Langtry-Menter gamma-Re_theta transition model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 */
template <class FlowIndices>
class CAvgGrad_TransLM final : public CAvgGrad_Scalar<FlowIndices> {
private:
    // saltate le costanti -> dovrebbe bastare averle nel variable

    // saltate le parti di AD

   using Base = CAvgGrad_Scalar<FlowIndices>;
   using Base::Laminar_Viscosity_i;
   using Base::Laminar_Viscosity_j;
   using Base::Eddy_Viscosity_i;
   using Base::Eddy_Viscosity_j;
   using Base::Density_i;
   using Base::Density_j;
   using Base::ScalarVar_i;
   using Base::ScalarVar_j;
   using Base::Proj_Mean_GradScalarVar;
   using Base::proj_vector_ij;
   using Base::implicit;
   using Base::Flux;
   using Base::Jacobian_i;
   using Base::Jacobian_j;

   su2double sigma_f, sigma_thetat;

   /*!
   * \brief Adds any extra variables to AD
    */
   void ExtraADPreaccIn() override {
//     AD::SetPreaccIn(F1_i, F1_j);
   }

    /*!
   * \brief LM specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
    void FinishResidualCalc(const CConfig *config) override {

//      if(config->dummyVar == 0)
//        cout << "CAvgGrad_TransLM::FinishResidualCalc" << endl;


      /*--- Compute mean effective viscosity ---*/
      const su2double diff_i_intermittency = Laminar_Viscosity_i + Eddy_Viscosity_i / sigma_f;
      const su2double diff_j_intermittency = Laminar_Viscosity_j + Eddy_Viscosity_j / sigma_f;
      const su2double diff_i_rethetat = sigma_thetat * (Laminar_Viscosity_i + Eddy_Viscosity_i);
      const su2double diff_j_rethetat = sigma_thetat * (Laminar_Viscosity_j + Eddy_Viscosity_j);

      const su2double diff_intermittency = 0.5*(diff_i_intermittency + diff_j_intermittency);
      const su2double diff_rethetat = 0.5*(diff_i_rethetat + diff_j_rethetat);

      Flux[0] = diff_intermittency*Proj_Mean_GradScalarVar[0];
      Flux[1] = diff_rethetat*Proj_Mean_GradScalarVar[1];

      /*--- For Jacobians -> Use of TSL (Thin Shear Layer) approx. to compute derivatives of the gradients ---*/
      if (implicit) {
        const su2double proj_on_rho_i = proj_vector_ij/Density_i;

        Jacobian_i[0][0] = -diff_intermittency*proj_on_rho_i;       Jacobian_i[0][1] = 0.0;
        Jacobian_i[1][0] = 0.0;                                     Jacobian_i[1][1] = -diff_rethetat*proj_on_rho_i;

        const su2double proj_on_rho_j = proj_vector_ij/Density_j;

        Jacobian_j[0][0] = diff_intermittency*proj_on_rho_j;        Jacobian_j[0][1] = 0.0;
        Jacobian_j[1][0] = 0.0;                                     Jacobian_j[1][1] = diff_rethetat*proj_on_rho_j;
      }


    }

public:
/*!
 * \brief Constructor of the class.
 * \param[in] val_nDim - Number of dimensions of the problem.
 * \param[in] val_nVar - Number of variables of the problem.
 * \param[in] constants - Constants of the model.
 * \param[in] correct_grad - Whether to correct gradient for skewness.
 * \param[in] config - Definition of the particular problem.
 */
    CAvgGrad_TransLM(unsigned short val_nDim, unsigned short val_nVar, const su2double* constants, bool correct_grad,
                        const CConfig *config) : CAvgGrad_Scalar<FlowIndices>(val_nDim, val_nVar, correct_grad, config),
    sigma_f(constants[6]),
    sigma_thetat(constants[7]) {
 }

    // saltata la parte di F1, che non esiste nel trans

    // è molto copia e incolla rispetto al turb, se è corretto come approccio allora bisognerà ripensare anche il sources!

};


