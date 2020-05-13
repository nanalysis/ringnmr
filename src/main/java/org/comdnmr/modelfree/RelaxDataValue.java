/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.modelfree;

import static org.comdnmr.modelfree.RelaxEquations.I;
import static org.comdnmr.modelfree.RelaxEquations.ImS;
import static org.comdnmr.modelfree.RelaxEquations.IpS;

/**
 *
 * @author brucejohnson
 */
public class RelaxDataValue {

    final int resIndex;
    final String resSpecifier;
    final double R1;
    final double R1err;
    final double R2;
    final double R2err;
    final double NOE;
    final double NOEerr;
    final RelaxEquations relaxObj;
    double rhoError = 0.0;

    public RelaxDataValue(int resIndex, String resSpecifier, double r1,
            double r1Error, double r2, double r2Error, double noe, double noeError,
            RelaxEquations relaxObj) {
        this.resIndex = resIndex;
        this.resSpecifier = resSpecifier;
        this.R1 = r1;
        this.R1err = r1Error;
        this.R2 = r2;
        this.R2err = r2Error;
        this.NOE = noe;
        this.NOEerr = noeError;
        this.relaxObj = relaxObj;
    }

    public void updateRhoError(double[] J) {
        double gammaS = relaxObj.getGammaS();
        double gammaI = relaxObj.getGammaI();
        double w1 = (J[ImS] + 6 * J[IpS]) / (6 * J[IpS] - J[ImS]);
        double w2 = (6 * J[I]) / (6 * J[IpS] - J[ImS]);
        double f = (gammaS / gammaI) * R1 * (NOE - 1);
        double rhoExp1 = 2 * R2 - R1 - w2 * f;
        double rhoExp2 = R1 - w1 * f;
        double rhoExpErr1 = Math.sqrt(2 * 2 * R2err * R2err + R1err * R1err
                + w2 * w2 * (f * f * ((R1err / R1) * (R1err / R1) + (NOEerr / NOE) * (NOEerr / NOE)) + R1err * R1err));
        double rhoExpErr2 = Math.sqrt(R1err * R1err
                + w1 * w1 * (f * f * ((R1err / R1) * (R1err / R1) + (NOEerr / NOE) * (NOEerr / NOE)) + R1err * R1err));
        double rhoExp = rhoExp1 / rhoExp2;
        rhoError = Math.sqrt(rhoExp * rhoExp * ((rhoExpErr1 / rhoExp1) * (rhoExpErr1 / rhoExp1) + (rhoExpErr2 / rhoExp2) * (rhoExpErr2 / rhoExp2)));

    }

    public String getResSpecifier() {
        return resSpecifier;
    }

    public double calcExpRho(double[] J) {
        double rhoExp = relaxObj.calcRhoExp(R1, R2, NOE, J);
        return rhoExp;
    }

    public double calcPredRho(double[] J) {
        double rhoPred = relaxObj.calcRhoPred(J);
        return rhoPred;
    }

    public double score2(double r1P, double r2P, double noeP) {
        double r1D = (r1P - R1) / R1err;
        double r2D = (r2P - R2) / R2err;
        double noeD = (noeP - NOE) / NOEerr;
        double score = r1D * r1D + r2D * r2D + noeD * noeD;
        return score;
    }

}
