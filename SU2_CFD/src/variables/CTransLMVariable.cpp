//
// Created by marcocera on 25/02/22.
//

#include "../../include/variables/CTransLMVariable.hpp"

// bisogna capire se va fatto come "if(LM)" in SST o se va fatto a parte (come questo file nuovo)

CTransLMVariable::CTransLMVariable(su2double intermittency, su2double Re_theta, const su2double* constants, unsigned long npoint, unsigned long ndim,
                                         unsigned long nvar, CConfig *config) :
                                         CTurbVariable(npoint, ndim, nvar, config) {

    for (unsigned long iPoint=0; iPoint<nPoint; ++iPoint) {
        Solution(iPoint,0) = intermittency;
        Solution(iPoint,1) = Re_theta;
    }

    Solution_Old = Solution;   // queste assegnazioni a Solution e Solution_Old (Runge-Kutta problem) sono copiate da SST =>
                               // bisogna capire se cambiare nome o se il vector è pensato per contenere un numero qualsiasi
                               // di elementi (forse contiene nvar elementi)


    F_length.resize(nPoint) = su2double(0.0);
    F_onset.resize(nPoint) = su2double(0.0);
    F_turb.resize(nPoint) = su2double(0.0);
    F_onset1.resize(nPoint) = su2double(0.0);
    reV.resize(nPoint) = su2double(0.0);
    F_onset2.resize(nPoint) = su2double(0.0);
    R_T.resize(nPoint) = su2double(0.0);
    F_onset3.resize(nPoint) = su2double(0.0);
    F_length1.resize(nPoint) = su2double(0.0);
    F_sublayer.resize(nPoint) = su2double(0.0);
    rew.resize(nPoint) = su2double(0.0);
    rethetac.resize(nPoint) = su2double(0.0);
    T.resize(nPoint) = su2double(0.0);
    rethetat_eq.resize(nPoint) = su2double(0.0);
    F_thetat.resize(nPoint) = su2double(0.0);
    Velocity_Mag.resize(nPoint) = su2double(0.0);
    delta_param.resize(nPoint) = su2double(0.0);
    F_wake.resize(nPoint) = su2double(0.0);
    lambda_theta.resize(nPoint) = su2double(0.0);
    Turb_Intens.resize(nPoint) = su2double(0.0);
    du_dx.resize(nPoint) = su2double(0.0);
    du_dy.resize(nPoint) = su2double(0.0);
    du_dz.resize(nPoint) = su2double(0.0);
    dU_ds.resize(nPoint) = su2double(0.0);
    F_lambda.resize(nPoint) = su2double(0.0);
    thetat.resize(nPoint) = su2double(0.0);
    gamma_sep.resize(nPoint) = su2double(1.0);
    gamma_eff.resize(nPoint) = su2double(1.0);


}

void CTransLMVariable::SetQuantities(unsigned long iPoint, su2double* constants,  su2double val_viscosity,
                                        su2double val_dist, su2double val_density, su2double val_vort, su2double StrainMag, CMatrixView<const su2double> Velocity_Gradient, su2double *Velocity, su2double *TurbVars) {



    /*--- Source and Sink terms ---*/
//    cout << "entro" << endl;

    /*--- F_length ---*/
    rew(iPoint) = val_density*TurbVars[1]*val_dist*val_dist/val_viscosity;   // assunto "val_viscosity" come la laminar
//    cout << "Dopo rew" << endl;
    if (Solution(iPoint, 1)<400) {
        F_length1(iPoint) = 39.8189+(-119.270*pow(10,-4))*Solution(iPoint, 1)+(-132.567*pow(10,-6))*Solution(iPoint, 1)*Solution(iPoint, 1);
    }
    else if (Solution(iPoint, 1) < 596) {
        F_length1(iPoint) = 263.404+(-123.939*pow(10,-2))*Solution(iPoint, 1)+(194.548*pow(10,-5))*pow(Solution(iPoint, 1),2)+
                    (-101.695*pow(10,-8))*pow(Solution(iPoint, 1),3);
    }
    else if (Solution(iPoint, 1) < 1200) {
        F_length1(iPoint) = 0.5-(3.0*pow(10,-4))*(Solution(iPoint, 1)-596.0);
    }
    else {
        F_length1(iPoint) = 0.3188;
    }
    F_sublayer(iPoint) = exp(-pow(0.005*rew(iPoint),2));
//    cout << "Dopo F_sublayer" << endl;
    F_length(iPoint) = F_length1(iPoint)*(1-F_sublayer(iPoint))+40.0*F_sublayer(iPoint);

//    cout << "Dopo F_Length" << endl;

    /*--- F_onset ---*/
    reV(iPoint) = val_density*StrainMag*val_dist*val_dist/val_viscosity;

    if (Solution(iPoint, 1)<=1870){
      su2double FirstTerm = (-396.035*pow(10,-2));
      su2double SecondTerm = (10120.656*pow(10,-4))*Solution(iPoint, 1);
      su2double ThirdTerm = (-868.230*pow(10,-6))*pow(Solution(iPoint, 1),2);
      su2double ForthTerm = (696.506*pow(10,-9))*pow(Solution(iPoint, 1),3);
      su2double FifthTerm = (-174.105*pow(10,-12))*pow(Solution(iPoint, 1),4);
        rethetac(iPoint) = FirstTerm + SecondTerm + ThirdTerm + ForthTerm + FifthTerm;
    }
    else {
        rethetac(iPoint) = Solution(iPoint, 1)-(593.11+0.482*(Solution(iPoint, 1)-1870.0));
    }

//    cout << "Dopo rethetac" << endl;

    /*--- F_onset ---*/
    F_onset1(iPoint) = reV(iPoint)/(2.193*rethetac(iPoint));
    R_T(iPoint) = val_density*TurbVars[0]/(val_viscosity*TurbVars[1]);   // assunto "val_viscosity" come la laminar
    F_onset2(iPoint) = min(max(F_onset1(iPoint),pow(F_onset1(iPoint),4)),2.0);
    F_onset3(iPoint) = max(1.0-pow(0.4*R_T(iPoint),3),0.0);
    F_onset(iPoint) = max(F_onset2(iPoint)-F_onset3(iPoint),0.0);
//    cout << "Dopo F_onset" << endl;

    /*--- F_turb ---*/
    F_turb(iPoint) = exp(-pow(0.25*R_T(iPoint),4));

    /*--- T ---*/
    if (nDim==2) {
        Velocity_Mag(iPoint) = sqrt(Velocity[0]*Velocity[0]+Velocity[1]*Velocity[1]);
    }
    else if (nDim==3) {
      Velocity_Mag(iPoint) = sqrt(Velocity[0]*Velocity[0]+Velocity[1]*Velocity[1]+Velocity[2]*Velocity[2]);
    }
    T(iPoint) = 500*val_viscosity/(val_density*Velocity_Mag(iPoint)*Velocity_Mag(iPoint));   // assunto "val_viscosity" come la laminar
//    cout << "Dopo T" << endl;

    /*--- rethetat_eq ---*/

    // Aggiunto da me
    // Calcolate come nel paper

    for(auto iDim = 0u; iDim < nDim; iDim++){
      du_dx(iPoint) += 2 * Velocity[iDim] * Velocity_Gradient[iDim][0];
    }
    du_dx(iPoint) = du_dx(iPoint) * 0.5 / Velocity_Mag(iPoint);

    for(auto iDim = 0u; iDim < nDim; iDim++){
      du_dy(iPoint) += 2 * Velocity[iDim] * Velocity_Gradient[iDim][1];
    }
    du_dy(iPoint) = du_dy(iPoint) * 0.5 / Velocity_Mag(iPoint);

    if(nDim == 3) {
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        du_dz(iPoint) += 2 * Velocity[iDim] * Velocity_Gradient[iDim][2];
      }
      du_dz(iPoint) = du_dz(iPoint) * 0.5 / Velocity_Mag(iPoint);
    }
//    cout << "Dopo du_dz" << endl;



    dU_ds(iPoint) = Velocity[0] * du_dx(iPoint) / Velocity_Mag(iPoint) + Velocity[1] * du_dy(iPoint) / Velocity_Mag(iPoint);
    if (nDim==3) {
        dU_ds(iPoint) += Velocity[2] * du_dz(iPoint) / Velocity_Mag(iPoint);
    }
//    cout << "Dopo dU_ds" << endl;

    Turb_Intens(iPoint) = 100*sqrt(2*TurbVars[0]/3)/Velocity_Mag(iPoint);
    Turb_Intens(iPoint) = max(Turb_Intens(iPoint), 0.027);
//    cout << "Dopo F_Length" << endl;

    // Guess iniziale fatta con lambda = 0.
    F_lambda(iPoint) = 1.0;
    su2double toll = 1e-5;
    su2double rethetat_eq_old = 1.0;
    int nMax = 100;
//    cout << "Dopo F_lambda" << endl;

    // Aggiustare assolutamente questa soluzione
    for (int iter=0;iter<nMax;iter++) {   // quante iterazioni? non è meglio un while che abbia una tolleranza sensata? nel caso sì, quale valore?

      thetat(iPoint) = rethetat_eq(iPoint)*val_viscosity/(val_density*Velocity_Mag(iPoint));
      lambda_theta(iPoint) = val_density*thetat(iPoint)*thetat(iPoint)*dU_ds(iPoint)/val_viscosity;
      // lambda_theta dovrebbe essere limitata tra -0.1 e 0.1
      lambda_theta(iPoint) = max(-0.1, lambda_theta(iPoint));
      lambda_theta(iPoint) = min(0.1, lambda_theta(iPoint));

      if (lambda_theta(iPoint)<=0.0) {
        su2double FirstTerm = 12.986*lambda_theta(iPoint);
        su2double SecondTerm = 123.66*pow(lambda_theta(iPoint),2);
        su2double ThirdTerm = 405.689*pow(lambda_theta(iPoint),3);
        F_lambda(iPoint) = 1 + ( FirstTerm + SecondTerm + ThirdTerm ) * exp(- pow(Turb_Intens(iPoint)/1.5,1.5));
      }
      else {
        su2double FirstTerm = exp(-35.0*lambda_theta(iPoint));
        F_lambda(iPoint) = 1.0 + 0.275*(1.0 - FirstTerm) * exp(-2.0*Turb_Intens(iPoint));
      }   // F_lambda andrà a essere usato per ricalcolare rethetat_eq (a inizio for)

      if (Turb_Intens(iPoint)<=1.3) {
            su2double FirstTerm = 589.428*Turb_Intens(iPoint);
            su2double SecondTerm = 0.2196/(Turb_Intens(iPoint)*Turb_Intens(iPoint));
            rethetat_eq(iPoint) = (1173.51 - FirstTerm + SecondTerm) * F_lambda(iPoint);
        }
        else {
            rethetat_eq(iPoint) = 331.5*pow(Turb_Intens(iPoint)-0.5658,-0.671)*F_lambda(iPoint);
        }

        rethetat_eq(iPoint) = max(20.0, rethetat_eq(iPoint));   // limite inferiore preso da pag 4 del sito NASA


        if(abs(rethetat_eq(iPoint) - rethetat_eq_old)/rethetat_eq_old < toll){
          iter = nMax+1;
        } else
          rethetat_eq_old = rethetat_eq(iPoint);

    }

//    cout << "Dopo rethetat_eq_old" << endl;

    /*--- F_thetat ---*/
    su2double theta_BL = Solution(iPoint, 1) * val_viscosity / (val_density * Velocity_Mag(iPoint));
    su2double delta_BL = 15.0 * theta_BL / 2.0;
    delta_param(iPoint) = 50.0 * val_vort * val_dist * delta_BL / Velocity_Mag(iPoint);
    F_wake(iPoint) = exp(-pow(0.00001*rew(iPoint),2));
    su2double c_e2 = constants[3];
    su2double FirstMaxTerm = F_wake(iPoint) * exp(-pow(val_dist/delta_param(iPoint), 4));
    su2double SecondMaxTerm = 1 - pow((c_e2 * Solution(iPoint, 0) -1)/(c_e2-1), 2);
//    cout << "Dopo SecondMaxTerm" << endl;

    F_thetat(iPoint) = min( max(FirstMaxTerm, SecondMaxTerm), 1.0);



    /* ------ Computation of separation induced intermittency ------ */

    su2double FReattach = exp(-pow(R_T(iPoint)/20.0, 4));

    su2double maxInnerTerm = max(0.0, (reV(iPoint)/(3.235*rethetac(iPoint)))-1);
    gamma_sep(iPoint) = min(constants[5] * maxInnerTerm*FReattach, 2.0);
    gamma_sep(iPoint) = gamma_sep(iPoint) * F_thetat(iPoint);


    bool print = false;
    if(print) {
      int iPoint2print = 200;
      if (iPoint == iPoint2print) {
        cout << "F_length(0) = " << F_length(iPoint2print) << endl;
        cout << "F_onset(0) = " << F_onset(iPoint2print) << endl;
        cout << "F_turb(0) = " << F_turb(iPoint2print) << endl;
        cout << "val_density = " << val_density << endl;
        cout << "val_viscosity = " << val_viscosity << endl;
        cout << "TurbVars(0) = " << TurbVars[0] << endl;
        cout << "TurbVars(1) = " << TurbVars[1] << endl;
        cout << "F_onset1(0) = " << F_onset1(iPoint2print) << endl;
        cout << "reV(0) = " << reV(iPoint2print) << endl;
        cout << "F_onset2(0) = " << F_onset2(iPoint2print) << endl;
        cout << "R_T(0) = " << R_T(iPoint2print) << endl;
        cout << "F_onset3(0) = " << F_onset3(iPoint2print) << endl;
        cout << "F_length1(0) = " << F_length1(iPoint2print) << endl;
        cout << "F_sublayer(0) = " << F_sublayer(iPoint2print) << endl;
        cout << "rew(0) = " << rew(iPoint2print) << endl;
        cout << "rethetac(0) = " << rethetac(iPoint2print) << endl;
        cout << "T(0) = " << T(iPoint2print) << endl;
        cout << "rethetat_eq(0) = " << rethetat_eq(iPoint2print) << endl;
        cout << "F_thetat(0) = " << F_thetat(iPoint2print) << endl;
        cout << "Velocity_Mag(0) = " << Velocity_Mag(iPoint2print) << endl;
        cout << "delta_param(0) = " << delta_param(iPoint2print) << endl;
        cout << "F_wake(0) = " << F_wake(0) << endl;
        cout << "lambda_theta(0) = " << lambda_theta(iPoint2print) << endl;
        cout << "Turb_Intens(0) = " << Turb_Intens(iPoint2print) << endl;
        cout << "du_dx(0) = " << du_dx(iPoint2print) << endl;
        cout << "du_dy(0) = " << du_dy(iPoint2print) << endl;
        cout << "du_dz(0) = " << du_dz(iPoint2print) << endl;
        cout << "dU_ds(0) = " << dU_ds(iPoint2print) << endl;
        cout << "F_lambda(0) = " << F_lambda(iPoint2print) << endl;
        cout << "thetat(0) = " << thetat(iPoint2print) << endl;
        cout << "gamma_sep(0) = " << gamma_sep(iPoint2print) << endl;
        cout << "gamma_eff(0) = " << gamma_eff(iPoint2print) << endl;
      }
    }


    // saltato la parte di AD

}
