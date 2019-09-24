/*
 * CoMD/NMR Software : A Program for Analyzing NMR Dynamics Data
 * Copyright (C) 2018-2019 Bruce A Johnson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
package org.comdnmr.modelfree;

import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author brucejohnson
 */
public class RelaxEquations {

    static final int S = 1;
    static final int ImS = 2;
    static final int I = 3;
    static final int IpS = 4;

    public final static double MU0 = 4.0e-7 * Math.PI;
    public final static double GAMMA_N = -2.71e7;
    public final static double GAMMA_H = 2.68e8;
    public final static double PLANCK = 1.054e-34;
    public final static double R_HN = 1.01e-10;
    public final static double SIGMA = 165.0e-6;
    public final static Map<String, Double> GAMMA_MAP = new HashMap<>();
    public final static Map<String, Double> R_MAP = new HashMap<>();

    static {
        GAMMA_MAP.put("H", GAMMA_H);
        GAMMA_MAP.put("N", GAMMA_N);
        R_MAP.put("HN", R_HN);
        R_MAP.put("NH", R_HN);
    }

    private final double r;
    private final double d;
    private final double d2;
    private final double c;
    private final double c2;
    private final double gammaS;
    private final double gammaI;
    private final double sf;
    private final double wI;
    private final double wS;

    //   consider using scaled versions (smaller exponents)
    public RelaxEquations(double sf, String elem1, String elem2) {
        gammaI = GAMMA_MAP.get(elem1);
        gammaS = GAMMA_MAP.get(elem2);
        wI = sf * 2.0 * Math.PI;
        wS = wI * gammaS / gammaI;
        r = R_MAP.get(elem1 + elem2);
        d = MU0 * (gammaI * gammaS * PLANCK) / (4.0 * Math.PI * r * r * r);
        d2 = d * d;
        c = wS * SIGMA / Math.sqrt(3.0);
        c2 = c * c;

        this.sf = sf;
    }

    // Note: taui = tm in Art Palmer's code, and taui in Relax. 
    // tau = ts in Art Palmer's code (taue in the paper: Phys Chem Chem Phys, 2016, 18, 5839-5849), and taue in Relax.
    public double JModelFree(double w, double tau, double taui, double s2) {
        double value1 = s2 / (1.0 + w * w * taui * taui);
        double value2 = ((1.0 - s2) * (tau + taui) * tau) / ((tau + taui) * (tau + taui) + w * w * taui * taui * tau * tau);
        double value = 0.4 * taui * (value1 + value2);
        return value;
    }
    
    // Note: taui = tm in Art Palmer's code. tau = ts in Art Palmer's code.
    public double JModelFree(double w, double tau, double taui, double s2, double sf2) {
        double value1 = s2 / (1.0 + w * w * taui * taui);
        double value2 = ((sf2 - s2) * (tau + taui) * tau) / ((tau + taui) * (tau + taui) + w * w * taui * taui * tau * tau);
        double value = 0.4 * taui * (value1 + value2);
        return value;
    }
    
    // Note: taui = tm in Art Palmer's code. tau = tf in Art Palmer's code. tauj = ts in Art Palmer's code.
    public double JModelFree(double w, double tau, double taui, double tauj, double s2, double sf2) {
        double value1 = s2 / (1.0 + w * w * taui * taui);
        double value2 = ((1 - sf2) * (tau + taui) * tau) / ((tau + taui) * (tau + taui) + w * w * taui * taui * tau * tau);
        double value3 = ((sf2 - s2) * (tauj + taui) * tauj) / ((tauj + taui) * (tauj + taui) + w * w * taui * taui * tauj * tauj);
        double value = 0.4 * taui * (value1 + value2 + value3);
        return value;
    }

    public double J(double w, double tau) {
        double value = 0.4 * tau / (1.0 + w * w * tau * tau);
        return value;
    }

    public double J(double w, double tau, double S2) {
        double value = 0.4 * S2 * tau / (1.0 + w * w * tau * tau);
        return value;
    }

    public double r2r1Ratio(double tau) {
        double num = 4.0 * J(0, tau) + J(wS - wI, tau) + 3.0 * J(wS, tau)
                + 6.0 * J(wI, tau) + 6.0 * J(wS + wI, tau)
                + (c2 / (3.0 * d2) * (4.0 * J(0, tau) + 3.0 * J(wS, tau)));

        double denom = 2.0 * J(wS - wI, tau) + 6.0 * J(wS, tau) + 12.0 * J(wS + wI, tau) + 2.0 * (c2 / (3.0 * d2) * J(wS, tau));
        return num / denom;
    }

    public double R1(double tau) {
        double dipolarContrib = d2 / 4.0 * (J(wI - wS, tau) + 3.0 * J(wS, tau) + 6.0 * J(wI + wS, tau));
        double csaContrib = c2 * J(wS, tau);
        return dipolarContrib + csaContrib;
    }

    public double R2(double tau) {
        double dipolarContrib = d2 / 8.0 * (4.0 * J(0.0, tau) + J(wI - wS, tau) + 3.0 * J(wS, tau)
                + 6.0 * J(wI, tau) + 6.0 * J(wI + wS, tau));
        double csaContrib = c2 / 6 * (4.0 * J(0.0, tau) + 3.0 * J(wS, tau));
        return dipolarContrib + csaContrib;
    }

    public double NOE(double tau) {
        double R1 = R1(tau);
        return 1.0 + (d2 / (4.0 * R1)) * (gammaI / gammaS)
                * (6.0 * J(wI + wS, tau) - J(wI - wS, tau));
    }

    public double R1(double[] J) {
        double dipolarContrib = d2 / 4.0 * (J[ImS] + 3.0 * J[S]
                + 6.0 * J[IpS]);
        double csaContrib = c2 * J[S];
        return dipolarContrib + csaContrib;
    }

    public double R2(double[] J, double Rex) {
        double dipolarContrib = d2 / 8.0 * (4.0 * J[0] + J[ImS]
                + 3.0 * J[S]
                + 6.0 * J[I] + 6.0 * J[IpS]);
        double csaContrib = c2 / 6 * (4.0 * J[0] + 3.0 * J[S]);
        return dipolarContrib + csaContrib + Rex;
    }

    public double NOE(double[] J) {
        double R1 = R1(J);
        return 1.0 + (d2 / (4.0 * R1)) * (gammaI / gammaS)
                * (6.0 * J[IpS] - J[ImS]);
    }
    
    public double sigmaSI(double[] J) {
        double R1 = R1(J);
        double NOE = NOE(J);
        return R1*(NOE - 1)*(gammaS / gammaI);
    }
    
    public double gamma(double[] J, double Rex) {
        double R1 = R1(J);
        double R2 = R2(J, Rex);
        double sigmaSI = sigmaSI(J);
        return R2 - 0.5*R1 - 0.454*sigmaSI;
    }
    
    public double[] getJ(double tau) {
        double J0 = J(0.0, tau); // R2
        double JIplusS = J(wI + wS, tau); //R1, R2, NOE
        double JIminusS = J(wI - wS, tau); // R1, R2, NOE
        double JI = J(wI, tau); // R2
        double JS = J(wS, tau); //R1, R2
        double[] result = {J0, JS, JIminusS, JI, JIplusS};
        return result;
    }
    
    // Note: taui = tm in Art Palmer's code, and taui in Relax. tau = ts in Art Palmer's code (taue in the paper), and taue in Relax.
    public double[] getJModelFree(double tau, double taui, double s2) {
        double J0 = JModelFree(0.0, tau, taui, s2);
        double JIplusS = JModelFree(wI + wS, tau, taui, s2); //R1, R2, NOE
        double JIminusS = JModelFree(wI - wS, tau, taui, s2); // R1, R2, NOE
        double JI = JModelFree(wI, tau, taui, s2); // R2
        double JS = JModelFree(wS, tau, taui, s2); //R1, R2
        double[] result = {J0, JS, JIminusS, JI, JIplusS};
        return result;
    }
    
    // Note: taui = tm in Art Palmer's code. tau = ts in Art Palmer's code.
    public double[] getJModelFree(double tau, double taui, double s2, double sf2) {
        double J0 = JModelFree(0.0, tau, taui, s2, sf2);
        double JIplusS = JModelFree(wI + wS, tau, taui, s2, sf2); //R1, R2, NOE
        double JIminusS = JModelFree(wI - wS, tau, taui, s2, sf2); // R1, R2, NOE
        double JI = JModelFree(wI, tau, taui, s2, sf2); // R2
        double JS = JModelFree(wS, tau, taui, s2, sf2); //R1, R2
        double[] result = {J0, JS, JIminusS, JI, JIplusS};
        return result;
    }
    
    // Note: taui = tm in Art Palmer's code. tau = tf in Art Palmer's code. tauj = ts in Art Palmer's code.
    public double[] getJModelFree(double tau, double taui, double tauj, double s2, double sf2) {
        double J0 = JModelFree(0.0, tau, taui, tauj, s2, sf2);
        double JIplusS = JModelFree(wI + wS, tau, taui, tauj, s2, sf2); //R1, R2, NOE
        double JIminusS = JModelFree(wI - wS, tau, taui, tauj, s2, sf2); // R1, R2, NOE
        double JI = JModelFree(wI, tau, taui, tauj, s2, sf2); // R2
        double JS = JModelFree(wS, tau, taui, tauj, s2, sf2); //R1, R2
        double[] result = {J0, JS, JIminusS, JI, JIplusS};
        return result;
    }

    public double[] getJ(double tau, double S) {
        double[] result = getJ(tau);
        for (int i = 0; i < result.length; i++) {
            result[i] *= S;
        }
        return result;
    }

}
