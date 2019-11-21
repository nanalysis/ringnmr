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
package org.comdnmr.eqnfit;

import org.comdnmr.util.MtxExp;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.ejml.data.Complex_F64;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.decomposition.eig.WatchedDoubleStepQRDecomposition_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;

/**
 *
 * @author brucejohnson
 */
public class R1RhoEquations {

    public static double[] r1rhoLaguerre(double[] omega, double pb, double kex, double[] deltaA, double[] deltaB, double R1A, double R1B, double R2A, double R2B) {
        // Calculates the Miloushev-Palmer Laguerre second-order approximation to the eigenvalue and returns R1rho.
        // Assumes the intrinsic relaxation rate constants for sites A and B are nearly identical, so that only
        // the average rates are calculated for projection from the laboratory frame to the tilted frame
        // (and are not incorporated into the exchange formalism).
        //
        // omega: B1 field strength (1/s)
        // pb: population of minor state
        // kex: k12+k21 (1/s)
        // deltaA: offset of A state (angular units, 1/s)
        // deltaB: offset of B state (angular units, 1/s)
        // R1A, R1B: R10 relaxation rate constants of A and B states
        // R2A, R2B: R20 relaxation rate constants of A and B states
        double pa = 1.0 - pb;
        int size = omega.length;
        double R1Bar = pa * R1A + pb * R1B;
        double R2Bar = pa * R2A + pb * R2B;
        double[] r1rho = new double[size];
        for (int i = 0; i < size; i++) {
            double dw = deltaB[i] - deltaA[i];
            double omegaBar = pa * deltaA[i] + pb * deltaB[i];
            double weA = Math.sqrt(omega[i] * omega[i] + deltaA[i] * deltaA[i]);
            double weB = Math.sqrt(omega[i] * omega[i] + deltaB[i] * deltaB[i]);
            double we = Math.sqrt(omega[i] * omega[i] + omegaBar * omegaBar);
            double sin2t = (omega[i] / we) * (omega[i] / we);
            double x = (pa * pb * (dw * dw)) * sin2t;
            double y = (weA * weA) * (weB * weB) / (we * we) + kex * kex;
            double z = x * (1 + 2 * (kex * kex) * (pa * (weA * weA) + pb * (weB * weB)) / ((weA * weA) * (weB * weB) + (we * we) * (kex * kex)));
            double REx = kex * x / (y - z);
            r1rho[i] = (1 - sin2t) * R1Bar + sin2t * R2Bar + REx;
        }
        return r1rho;
    }

    public static double[] r1rhoBaldwinKay(double[] omega, double pb, double kex, double[] deltaA, double[] deltaB, double R1A, double R1B, double R2A, double R2B) {
        // Calculates the Baldwin-Kay first-order approximation to the eigenvalue and returns R1rho.
        // Allows the intrinsic relaxation rate constants for sites A and B to differ.
        //
        // omega: B1 field strength (1/s)
        // pb: population of minor state
        // kex: k12+k21 (1/s)
        // deltaA: offset of A state (angular units, 1/s)
        // deltaB: offset of B state (angular units, 1/s)
        // R1A, R1B: R10 relaxation rate constants of A and B states
        // R2A, R2B: R20 relaxation rate constants of A and B states
        int size = omega.length;
        double pa = 1.0 - pb;
        double dR = R2B - R2A;
        double[] r1rho = new double[size];
        for (int i = 0; i < size; i++) {
            double dw = deltaB[i] - deltaA[i];
            double omegaBar = pa * deltaA[i] + pb * deltaB[i];
            double weA = Math.sqrt(omega[i] * omega[i] + deltaA[i] * deltaA[i]);
            double weB = Math.sqrt(omega[i] * omega[i] + deltaB[i] * deltaB[i]);
            double we = Math.sqrt(omega[i] * omega[i] + omegaBar * omegaBar);
            double sin2t = (omega[i] / we) * (omega[i] / we);
            double cos2t = 1 - sin2t;
            double tan2t = sin2t / cos2t;
            double f1p = (pa * pb * (dw * dw));
            double f2p = kex * kex + omega[i] * omega[i] + (deltaA[i] * deltaA[i]) * (deltaB[i] * deltaB[i]) / (omegaBar * omegaBar);
            double dp = kex * kex + (weA * weA) * (weB * weB) / (we * we);
            double f1 = pb * (weA * weA + kex * kex + dR * pa * kex);
            double f2 = 2 * kex + omega[i] * omega[i] / kex + dR * pa;
            double f3 = 3 * pb * kex + (2 * pa * kex + omega[i] * omega[i] / kex + dR + dR * (pb * pb) * (kex * kex) / (weA * weA)) * ((weA * weA) / (omega[i] * omega[i]));
            double c1 = (f2p + (f1p + dR * (f3 - f2)) * tan2t) / (dp + dR * f3 * sin2t);
            double c2 = (dp / sin2t - f2p / tan2t - f1p + dR * f2) / (dp + dR * f3 * sin2t);
            double rex = (f1p * kex + dR * f1) / (dp + dR * f3 * sin2t);
            r1rho[i] = c1 * R1A * cos2t + sin2t * (c2 * R2A + rex);
        }
        return r1rho;
    }

    public static double[] r1rhoPerturbationNoEx(double[] omega, double[] deltaA, double R1A, double R2A) {
        int size = omega.length;
        double[] r1rho = new double[size];
        for (int i = 0; i < size; i++) {
            double weA = Math.sqrt(omega[i] * omega[i] + deltaA[i] * deltaA[i]);
            double sin2t = (omega[i] / weA) * (omega[i] / weA);
            r1rho[i] = (1 - sin2t) * R1A + sin2t * R2A;
        }
        return r1rho;
    }

    public static double r1rhoExact0(DMatrixRMaj Z, double[] m0, double[] m1, double tdelay) {//(double tdelay, double omega, double pB, double kex, double deltaA, double deltaB, double R1A, double R1B, double R2A, double R2B) {
        // Performs an exact numerical calculation and returns CEST intensity ratio.
        //
        // X: array containing two arrays:
        //  omegarf: CEST irradiation frequency (ppm)
        //  omega1: B1 field strength (1/s)
        // 
        // pb: population of minor state
        // kex: k12+k21 (1/s)
        // deltaA: offset of A state (angular units, 1/s)
        // deltaB: offset of B state (angular units, 1/s)
        // R1A, R1B: R10 relaxation rate constants of A and B states
        // R2A, R2B: R20 relaxation rate constants of A and B states

        // time delay is hard-coded below

        DMatrixRMaj at = MtxExp.matrixExp(Z);
        DMatrixRMaj v = new DMatrixRMaj(m0);
        CommonOps_DDRM.mult(at, new DMatrixRMaj(m0), v);
        double magA = 0.0;
        double magA0 = 0.0;
        for (int i = 0; i < v.getNumRows(); i++) {
            magA += m1[i] * v.get(i, 0);
            magA0 += m1[i] * m0[i];
        }

//        double magA0 = m0[2] * m1[2] + m0[5] * m1[5];
//        double magA = at[2][2] * m0[2] * m1[2] + at[2][5] * m0[5] * m1[5];
//        magA = Math.abs(magA);
        double r1rho = -Math.log(magA / magA0) / tdelay;
//        System.out.println(omega + " " + at[2][2] + " " + at[2][5] + " " + magA0 + " " + magA + " " + r1rho);

        return r1rho;
    }

    public static double r1rhoExact(DMatrixRMaj Z) {//(double omega, double pb, double kex, double deltaA, double deltaB, double R1A, double R1B, double R2A, double R2B) {
        // Performs an exact numerical calculation of the eigenvalue and returns R1rho.
        //
        // omega: B1 field strength (1/s)
        // pb: population of minor state
        // kex: k12+k21 (1/s)
        // deltaA: offset of A state (angular units, 1/s)
        // deltaB: offset of B state (angular units, 1/s)
        // R1A, R1B: R10 relaxation rate constants of A and B states
        // R2A, R2B: R20 relaxation rate constants of A and B states
        List<Double> r1rho1 = new ArrayList<>();
        WatchedDoubleStepQRDecomposition_DDRM eig = (WatchedDoubleStepQRDecomposition_DDRM) DecompositionFactory_DDRM.eig(Z.getNumCols(), true, false);
        eig.decompose(Z);
        int numEigVec = eig.getNumberOfEigenvalues();
        for (int v=0; v<numEigVec; v++) {
            Complex_F64 ZEigVal = eig.getEigenvalue(v);
            if (ZEigVal.isReal()) {
                r1rho1.add(ZEigVal.getReal());
            }
        }
        double r1rho = Math.abs(Collections.max(r1rho1));
        // original Python code for the above, starting from Z = ...:
        // Z = La+Lb+K
        // lam = np.linalg.eigvals(Z)
        // condition = np.imag(lam) == 0
        // r1rho = np.extract(condition, lam)
        // r1rho = np.abs(np.max(r1rho))
        return r1rho;
    }

    public static double[] r1rhoPerturbation(double[] omega, double pb, double kex, double[] deltaA, double[] deltaB, double R1A, double R1B, double R2A, double R2B) {
        // Calculates the Trott-Palmer perturbation approximation to the eigenvalue and returns R1rho.
        // Allows the intrinsic relaxation rate constants for sites A and B to differ.
        // This result is not as accurate as the first-order Baldwin-Kay result, but simpler.
        //
        // omega: B1 field strength (1/s)
        // pb: population of minor state
        // kex: k12+k21 (1/s)
        // deltaA: offset of A state (angular units, 1/s)
        // deltaB: offset of B state (angular units, 1/s)
        // R1A, R1B: R10 relaxation rate constants of A and B states
        // R2A, R2B: R20 relaxation rate constants of A and B states
        int size = omega.length;
        double pa = 1.0 - pb;
        double k1 = pb * kex;
        double km1 = pa * kex;
        double dR = Math.abs(R2B - R2A);
        double[] r1rho = new double[size];
        for (int i = 0; i < size; i++) {
            double dw = deltaB[i] - deltaA[i];
            double weA = Math.sqrt(omega[i] * omega[i] + deltaA[i] * deltaA[i]);
            double weB = Math.sqrt(omega[i] * omega[i] + deltaB[i] * deltaB[i]);
            double sin2t = (omega[i] / weA) * (omega[i] / weA);
            double x = ((dw * dw) + (dR * dR)) * km1 + dR * ((weA * weA) + (km1 * km1));
            double y = km1 * (weB * weB + (km1 + dR) * (km1 + dR)) + dR * (omega[i] * omega[i]);
            double REx = k1 * x / y;
            r1rho[i] = (1 - sin2t) * R1A + sin2t * R2A + sin2t * REx;
        }
        return r1rho;
    }

}
