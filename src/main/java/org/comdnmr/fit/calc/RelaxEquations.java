package org.comdnmr.fit.calc;

/**
 *
 * @author brucejohnson
 */
public class RelaxEquations {

    double mu0 = 4.0e-7 * Math.PI;
    double gammaI = -2.71e7;
    double gammaS = 2.68e8;
    double planck = 1.054e-34;
    double r = 1.01e-10;
    // don't use this yet.  Various inconsistencies with various different presentations of equations
    //   consider using scaled versions (smaller exponents)

    public RelaxEquations() {

    }

    public double JModelFree(double w, double tau, double taui, double s2) {
        double value1 = s2 / (1.0 + w * w * taui * taui);
        double value2 = ((1.0 - s2) * (tau + taui) * tau) / ((tau + taui) * (tau + taui) + w * w * taui * taui * tau * tau);
        double value = 0.4 * taui * (value1 + value2);
        return value;
    }

    public double J(double w, double tau) {
        double value = tau / (1.0 + w * w * tau * tau);
        return value;
    }

    public double T1(double tau) {
        double wH = 800.0e6 * 2.0 * Math.PI;
        double wN = Math.abs(wH * gammaI / gammaS);
        double A = (J(wH - wN, tau) + 3.0 * J(wN, tau) + 6.0 * J(wH + wN, tau));
        double d = mu0 * (gammaS * gammaI * planck) / (4.0 * Math.PI * r * r * r);
        double d2 = 0.1 * d * d;

        double invT1 = d2 * A;
        return invT1;

    }

    public double T2(double tau) {
        double wH = 800.0e6 * 2.0 * Math.PI;
        double wN = Math.abs(wH * gammaI / gammaS);
        double A = (4.0 * J(0.0, tau) + J(wH - wN, tau) + 3.0 * J(wN, tau) + 6.0 * J(wH, tau) + 6.0 * J(wH + wN, tau));
        double d = mu0 * (gammaS * gammaI * planck) / (4.0 * Math.PI * r * r * r);
        double d2 = 0.1 * d * d;

        double invT1 = 0.5 * d2 * A;
        return invT1;

    }

    public double T2s(double tau) {
        double A = (1.0 / Math.pow((2.0 * Math.PI), 2)) * Math.pow(mu0 / (4.0 * Math.PI), 2);
        double x = Math.pow(gammaI * gammaS * planck, 2) / Math.pow(r, 6);
        double invT2 = A * x * tau;
        return invT2;
    }
}
