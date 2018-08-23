/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;

/**
 *
 * @author Martha Beckwith
 *
 *
 * Original Python script by Arthur Palmer.
 *
 *
 * Function definitions:
 *
 * R1rho functions:
 *
 * r1rhoLaguerre Second-order approximation, assuming R1A=R1B and R2A=R2B
 *
 * r1rhoPerturbation Perturbation analysis, assuming R1A=R1B, but R2A != R2B
 *
 * r1rhoBaldwinKay First-order approximation, assuming R1A=R1B, but R2A != R2B
 *
 * r1rhoExact Numerical determination of least negative eigenvalue of
 * Bloch-McConnell rate matrix. No restriction on R10 or R20
 *
 * CEST functions
 *
 * cestExact0 Numerical integration of thermalized Bloch-McConnell rate matrix.
 * No restriction on R10 or R20
 *
 * cestExact1 Numerical integration of thermalized Bloch-McConnell rate matrix.
 * R1A = R1B
 *
 * cestExact2 Numerical integration of thermalized Bloch-McConnell rate matrix.
 * R2A = R2B
 *
 * cestR1rhoN Uses Laguerre R1rho approximation to calculate CEST Second-order
 * approximation, assuming R1A=R1B and R2A=R2B
 *
 * cestR1rhoPerturbation Uses Trott perturbation R1rho approximation to
 * calculate CEST Assumes R1A=R1B
 *
 * cestR1rhoSD Uses perturbation R1rho approximation to calculate CEST
 * Second-order approximation, assuming R1A=R1B Averages over B1 inhomogeneity
 *
 *
 *
 * cestR1rhoBaldwinKay Uses first-order R1rho approximation to calculate CEST
 * Assumes R1A=R1B
 *
 * cestR1rhoExact1 Uses numerical determination of least negative eigenvalue to
 * calculate CEST Assumes R1A=R1B
 */
public class CESTEquations {

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
        double[] dw = new double[deltaA.length];
        double[] omegaBar = new double[deltaA.length];
        for (int i = 0; i < deltaA.length; i++) {
            dw[i] = deltaB[i] - deltaA[i];
            omegaBar[i] = pa * deltaA[i] + pb * deltaB[i];
        }

        double[] weA = new double[omega.length];
        double[] weB = new double[omega.length];
        double[] we = new double[omega.length];
        double[] sin2t = new double[omega.length];
        for (int i = 0; i < omega.length; i++) {
            weA[i] = Math.sqrt(omega[i] * omega[i] + deltaA[i] * deltaA[i]);
            weB[i] = Math.sqrt(omega[i] * omega[i] + deltaB[i] * deltaB[i]);
            we[i] = Math.sqrt(omega[i] * omega[i] + omegaBar[i] * omegaBar[i]);
            sin2t[i] = (omega[i] / we[i]) * (omega[i] / we[i]);
        }

        double k1 = pb * kex;
        double km1 = pa * kex;

        double R1Bar = pa * R1A + pb * R1B;
        double R2Bar = pa * R2A + pb * R2B;

        double[] x = new double[weA.length];
        double[] y = new double[weB.length];
        double[] z = new double[weB.length];
        double[] r1rho = new double[x.length];
        for (int i = 0; i < weA.length; i++) {
            x[i] = (pa * pb * (dw[i] * dw[i])) * sin2t[i];
            y[i] = (weA[i] * weA[i]) * (weB[i] * weB[i]) / (we[i] * we[i]) + kex * kex;
            z[i] = x[i] * (1 + 2 * (kex * kex) * (pa * (weA[i] * weA[i]) + pb * (weB[i] * weB[i])) / ((weA[i] * weA[i]) * (weB[i] * weB[i]) + (we[i] * we[i]) * (kex * kex)));
            r1rho[i] = kex * x[i] / (y[i] - z[i]);
            r1rho[i] = (1 - sin2t[i]) * R1Bar + sin2t[i] * R2Bar + r1rho[i];
        }

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

        double pa = 1.0 - pb;
        double[] dw = new double[deltaA.length];
        double[] omegaBar = new double[deltaA.length];
        for (int i = 0; i < deltaA.length; i++) {
            dw[i] = deltaB[i] - deltaA[i];
            omegaBar[i] = pa * deltaA[i] + pb * deltaB[i];
        }
        double dR = R2B - R2A;

        double[] weA = new double[omega.length];
        double[] weB = new double[omega.length];
        double[] we = new double[omega.length];
        double[] sin2t = new double[omega.length];
        double[] cos2t = new double[omega.length];
        double[] tan2t = new double[omega.length];
        for (int i = 0; i < omega.length; i++) {
            weA[i] = Math.sqrt(omega[i] * omega[i] + deltaA[i] * deltaA[i]);
            weB[i] = Math.sqrt(omega[i] * omega[i] + deltaB[i] * deltaB[i]);
            we[i] = Math.sqrt(omega[i] * omega[i] + omegaBar[i] * omegaBar[i]);
            sin2t[i] = (omega[i] / we[i]) * (omega[i] / we[i]);
            cos2t[i] = 1 - sin2t[i];
            tan2t[i] = sin2t[i] / cos2t[i];
        }

        double k1 = pb * kex;
        double km1 = pa * kex;

        double[] f1p = new double[dw.length];
        double[] f2p = new double[omega.length];
        double[] dp = new double[weA.length];
        double[] f1 = new double[weA.length];
        double[] f2 = new double[omega.length];
        double[] f3 = new double[omega.length];
        double[] c1 = new double[tan2t.length];
        double[] c2 = new double[tan2t.length];
        double[] rex = new double[f1.length];
        double[] r1rho = new double[sin2t.length];
        for (int i = 0; i < dw.length; i++) {
            f1p[i] = (pa * pb * (dw[i] * dw[i]));
            f2p[i] = kex * kex + omega[i] * omega[i] + (deltaA[i] * deltaA[i]) * (deltaB[i] * deltaB[i]) / (omegaBar[i] * omegaBar[i]);
            dp[i] = kex * kex + (weA[i] * weA[i]) * (weB[i] * weB[i]) / (we[i] * we[i]);
            f1[i] = pb * (weA[i] * weA[i] + kex * kex + dR * pa * kex);
            f2[i] = 2 * kex + omega[i] * omega[i] / kex + dR * pa;
            f3[i] = 3 * pb * kex + (2 * pa * kex + omega[i] * omega[i] / kex + dR + dR * (pb * pb) * (kex * kex) / (weA[i] * weA[i])) * ((weA[i] * weA[i]) / (omega[i] * omega[i]));
            c1[i] = (f2p[i] + (f1p[i] + dR * (f3[i] - f2[i])) * tan2t[i]) / (dp[i] + dR * f3[i] * sin2t[i]);
            c2[i] = (dp[i] / sin2t[i] - f2p[i] / tan2t[i] - f1p[i] + dR * f2[i]) / (dp[i] + dR * f3[i] * sin2t[i]);
            rex[i] = (f1p[i] * kex + dR * f1[i]) / (dp[i] + dR * f3[i] * sin2t[i]);
            r1rho[i] = c1[i] * R1A * cos2t[i] + sin2t[i] * (c2[i] * R2A + rex[i]);
        }

        return r1rho;
    }

    public static double r1rhoExact(double omega, double pb, double kex, double deltaA, double deltaB, double R1A, double R1B, double R2A, double R2B) {
        // Performs an exact numerical calculation of the eigenvalue and returns R1rho.
        //
        // omega: B1 field strength (1/s)
        // pb: population of minor state
        // kex: k12+k21 (1/s)
        // deltaA: offset of A state (angular units, 1/s)
        // deltaB: offset of B state (angular units, 1/s)
        // R1A, R1B: R10 relaxation rate constants of A and B states
        // R2A, R2B: R20 relaxation rate constants of A and B states

        double k1 = pb * kex;
        double km1 = (1 - pb) * kex;

        double[][] La = {{-R2A, -deltaA, 0, 0, 0, 0},
        {deltaA, -R2A, -omega, 0, 0, 0},
        {0, omega, -R1A, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0}};

        double[][] Lb = {{0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, -R2B, -deltaB, 0},
        {0, 0, 0, deltaB, -R2B, -omega},
        {0, 0, 0, 0, omega, -R1B}};

        double[][] K = {{-k1, 0, 0, km1, 0, 0},
        {0, -k1, 0, 0, km1, 0},
        {0, 0, -k1, 0, 0, km1},
        {k1, 0, 0, -km1, 0, 0},
        {0, k1, 0, 0, -km1, 0},
        {0, 0, k1, 0, 0, -km1}};

        double[][] Z = new double[La.length][La[0].length];
        for (int i = 0; i < La.length; i++) {
            for (int j = 0; j < La[i].length; j++) {
                Z[i][j] = La[i][j] + Lb[i][j] + K[i][j];
            }
        }

        RealMatrix zR = new Array2DRowRealMatrix(Z);
        EigenDecomposition Zeig = new EigenDecomposition(zR);
        double[] Zeigreal = Zeig.getRealEigenvalues();
        double[] Zeigimag = Zeig.getImagEigenvalues();

        List<Double> r1rho1 = new ArrayList<>();
        for (int i = 0; i < Zeigimag.length; i++) {
            if (Zeigimag[i] == 0) {
                r1rho1.add(Zeigreal[i]);
                // i++;
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

    public static double[] cestExact0(double[][] X, double pb, double kex, double deltaA0, double deltaB0, double R1A, double R1B, double R2A, double R2B) {
        // Performs an exact numerical calculation and returns CEST intensity ratio.
        //
        // X: array containing two arrays:
        //  omegarf: CEST irradiation frequency (1/s)
        //  omega1: B1 field strength (1/s)
        // 
        // pb: population of minor state
        // kex: k12+k21 (1/s)
        // deltaA: offset of A state (angular units, 1/s)
        // deltaB: offset of B state (angular units, 1/s)
        // R1A, R1B: R10 relaxation rate constants of A and B states
        // R2A, R2B: R20 relaxation rate constants of A and B states

        double[] omegarf = X[0];
        double[] omega1 = X[1];
        double[] Tex = X[2];

        // time delay is hard-coded below
        double tdelay = Tex[0];//0.3;

        double[] cest = new double[omegarf.length];

        double k1 = pb * kex;
        double km1 = (1 - pb) * kex;

        double[] m0 = {0.5, 0, 0, 1 - pb, 0, 0, pb};
        double[] m1 = {0.5, 0, 0, -(1 - pb), 0, 0, -pb};

        double[][] K = {{0, 0, 0, 0, 0, 0, 0},
        {0, -k1, 0, 0, km1, 0, 0},
        {0, 0, -k1, 0, 0, km1, 0},
        {0, 0, 0, -k1, 0, 0, km1},
        {0, k1, 0, 0, -km1, 0, 0},
        {0, 0, k1, 0, 0, -km1, 0},
        {0, 0, 0, k1, 0, 0, -km1}};

        for (int i = 0; i < omegarf.length; i++) {
            double deltaA = deltaA0 - omegarf[i];
            double deltaB = deltaB0 - omegarf[i];
            double omegaBar = (1 - pb) * deltaA + pb * deltaB;
            double we = Math.sqrt(omega1[i] * omega1[i] + omegaBar * omegaBar);

            double sint = omega1[i] / we;
            double cost = omegaBar / we;

            double[][] La = {{0, 0, 0, 0, 0, 0, 0},
            {0, -R2A, -deltaA, 0, 0, 0, 0},
            {0, deltaA, -R2A, -omega1[i], 0, 0, 0},
            {2 * R1A * (1 - pb), 0, omega1[i], -R1A, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0}};

            double[][] Lb = {{0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, -R2B, -deltaB, 0},
            {0, 0, 0, 0, deltaB, -R2B, -omega1[i]},
            {2 * R1B * pb, 0, 0, 0, 0, omega1[i], -R1B}};

            double[][] Z = new double[La.length][La[0].length];
            for (int k = 0; k < La.length; k++) {
                for (int j = 0; j < La[k].length; j++) {
                    Z[k][j] = La[k][j] + Lb[k][j] + K[k][j];
                }
            }

            double[][] at = new double[Z.length][Z[0].length];
            for (int k = 0; k < Z.length; k++) {
                for (int j = 0; j < Z[k].length; j++) {
                    at[k][j] = tdelay * Z[k][j];
                }
            }

            at = MtxExp.matrixExp(at);

            double magA = at[3][0] * m0[0] + at[3][3] * m0[3] + at[3][6] * m0[6];
            magA = magA - (at[3][0] * m1[0] + at[3][3] * m1[3] + at[3][6] * m1[6]);
            magA = magA / 2;

            cest[i] = magA;
        }

        return cest;
    }

    public static double[] cestR1rhoExact1(double[][] X, double pb, double kex, double deltaA0, double deltaB0, double R1A, double R1B, double R2A, double R2B) {
        // Performs an exact numerical calculation and returns CEST intensity ratio.
        // Assumes R1A = R1B.
        //
        // X: array containing two arrays:
        //  omegarf: CEST irradiation frequency (1/s)
        //  omega1: B1 field strength (1/s)
        // 
        // pb: population of minor state
        // kex: k12+k21 (1/s)
        // deltaA: offset of A state (angular units, 1/s)
        // deltaB: offset of B state (angular units, 1/s)
        // R1A, R1B: R10 relaxation rate constants of A and B states
        // R2A, R2B: R20 relaxation rate constants of A and B states

        double[] omegarf = X[0];
        double[] omega1 = X[1];
        double[] Tex = X[2];

        //double R1B = R1A;
        // time delay is hard-coded below
        double tdelay = Tex[0];//0.3;

        double[] cest = new double[omegarf.length];

        double k1 = pb * kex;
        double km1 = (1 - pb) * kex;

        double[][] K = {{-k1, 0, 0, km1, 0, 0},
        {0, -k1, 0, 0, km1, 0},
        {0, 0, -k1, 0, 0, km1},
        {k1, 0, 0, -km1, 0, 0},
        {0, k1, 0, 0, -km1, 0},
        {0, 0, k1, 0, 0, -km1}};

        for (int i = 0; i < omegarf.length; i++) {
            double deltaA = deltaA0 - omegarf[i];
            double deltaB = deltaB0 - omegarf[i];
            double omegaBar = (1 - pb) * deltaA + pb * deltaB;
            double we = Math.sqrt(omega1[i] * omega1[i] + omegaBar * omegaBar);

            double sint = omega1[i] / we;
            double cost = omegaBar / we;

            double[][] La = {{-R2A, -deltaA, 0, 0, 0, 0},
            {deltaA, -R2A, -omega1[i], 0, 0, 0},
            {0, omega1[i], -R1A, 0, 0, 0},
            {0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0}};

            double[][] Lb = {{0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0},
            {0, 0, 0, -R2B, -deltaB, 0},
            {0, 0, 0, deltaB, -R2B, -omega1[i]},
            {0, 0, 0, 0, omega1[i], -R1B}};

            double[][] Z = new double[La.length][La[0].length];
            for (int k = 0; k < La.length; k++) {
                for (int j = 0; j < La[k].length; j++) {
                    Z[k][j] = La[k][j] + Lb[k][j] + K[k][j];
                }
            }

            RealMatrix zR = new Array2DRowRealMatrix(Z);
            EigenDecomposition Zeig = new EigenDecomposition(zR);
            double[] Zeigreal = Zeig.getRealEigenvalues();
            double[] Zeigimag = Zeig.getImagEigenvalues();

            List<Double> r1rho1 = new ArrayList<>();
            for (int j = 0; j < Zeigimag.length; j++) {
                if (Zeigimag[j] == 0) {
                    r1rho1.add(Zeigreal[j]);
                    // i++;
                }
            }

            double r1rho = Math.abs(Collections.max(r1rho1));

            cest[i] = cost * cost * Math.exp(-tdelay * r1rho);
        }

        return cest;
    }

    public static double[] cestR1rhoApprox(String approx, double[][] X, double pb, double kex, double deltaA0, double deltaB0, double R1A, double R1B, double R2A, double R2B) {

        // X: array containing two arrays:
        //  omegarf: CEST irradiation frequency (1/s)
        //  omega1: B1 field strength (1/s)
        //
        // pb: population of minor state
        // kex: k12+k21 (1/s)
        // deltaA: offset of A state (angular units, 1/s)
        // deltaB: offset of B state (angular units, 1/s)
        // R1A, R1B: R10 relaxation rate constants of A and B states
        // R2A, R2B: R20 relaxation rate constants of A and B states
        // In the present implementation, the irradiation time is hard-coded below
        double[] omegarf = X[0];
        double[] omega1 = X[1];
        double[] Tex = X[2];

        double trad = Tex[0];//0.3;

        double pa = 1.0 - pb;
        double[] deltaA = new double[omegarf.length];
        double[] deltaB = new double[omegarf.length];
        double[] omegaBar = new double[omegarf.length];
        for (int i = 0; i < omegarf.length; i++) {
            deltaA[i] = deltaA0 - omegarf[i];
            deltaB[i] = deltaB0 - omegarf[i];
            omegaBar[i] = pa * deltaA[i] + pb * deltaB[i];
        }

        double[] we = new double[omegaBar.length];
        double[] cos2t = new double[omegaBar.length];
        for (int i = 0; i < omegaBar.length; i++) {
            we[i] = Math.sqrt(omega1[i] * omega1[i] + omegaBar[i] * omegaBar[i]);
            cos2t[i] = (omegaBar[i] / we[i]) * (omegaBar[i] / we[i]);
        }

        double[] cest = new double[cos2t.length];

        if (approx == "laguerre") {
            // if(R2.length == 1){
            double[] r1rho = r1rhoLaguerre(omega1, pb, kex, deltaA, deltaB, R1A, R1B, R2A, R2B);
            //double[] cest = new double[cos2t.length];
            for (int i = 0; i < cos2t.length; i++) {
                cest[i] = cos2t[i] * Math.exp(-trad * r1rho[i]);
            }
//            } else {
//                System.out.print("Error: Incorrect number of values for R2. Input one value for R2.");
//            } 
        } else if (approx == "trott") {
            double[] r1rho = r1rhoPerturbation(omega1, pb, kex, deltaA, deltaB, R1A, R1B, R2A, R2B);
            //double[] cest = new double[cos2t.length];
            for (int i = 0; i < cos2t.length; i++) {
                cest[i] = cos2t[i] * Math.exp(-trad * r1rho[i]);
            }

        } else if (approx == "baldwinkay") {
            double[] r1rho = r1rhoBaldwinKay(omega1, pb, kex, deltaA, deltaB, R1A, R1B, R2A, R2B);
            //double[] cest = new double[cos2t.length];
            for (int i = 0; i < cos2t.length; i++) {
                cest[i] = cos2t[i] * Math.exp(-trad * r1rho[i]);
            }

        } else if (approx == "sd") {
            int wlen = 11;

            double omegaSD = 0.2;   // fractional variation in B1
            double[] omwt = {0.022, 0.0444, 0.0777, 0.1159, 0.1473, 0.1596, 0.1473, 0.1159, 0.0777, 0.0444, 0.0216};
            double omwtsum = 0;
            for (int i = 0; i < omwt.length; i++) {
                omwtsum += omwt[i];
            }
            for (int i = 0; i < omwt.length; i++) {
                omwt[i] = omwt[i] / omwtsum;
            }

            double[] omegagauss = new double[wlen];
            for (int i = 0; i < wlen; i++) {
                omegagauss[i] = -2 * omegaSD + i * (2 * omegaSD - (-2 * omegaSD)) / (wlen - 1);
            }

            double[] magA = new double[omega1.length];
            double[] omegatmp = new double[omega1.length];
            double[] r1rho = new double[omegatmp.length];
            //double[] cest = new double[r1rho.length];
            for (int i = 0; i < wlen; i++) {
                for (int j = 0; j < omegatmp.length; j++) {
                    omegatmp[j] = omega1[j] * (1 + omegagauss[i]);
                    we[j] = Math.sqrt(omegatmp[j] * omegatmp[j] + omegaBar[j] * omegaBar[j]);
                    cos2t[j] = (omegaBar[j] / we[j]) * (omegaBar[j] / we[j]);
                }

                r1rho = r1rhoPerturbation(omegatmp, pb, kex, deltaA, deltaB, R1A, R1B, R2A, R2B);

                for (int j = 0; j < cest.length; j++) {
                    cest[j] = cos2t[j] * Math.exp(-trad * r1rho[j]);
                    magA[j] = magA[j] + omwt[i] * cest[j];
                }
            }

            cest = magA;
        }

        return cest;
    }

    public static double[][] cestPeakGuess(double[][] xvals, double[] yvals) {
        // Estimates CEST peak positions for initial guesses for before fitting.

        List<Double> xmin = new ArrayList<>();
        List<Double> ymin = new ArrayList<>();
        List<Double> fwhm = new ArrayList<>();

        double previousYdiff = 0;

        double B1val = xvals[1][0];

        for (int i = 3; i < yvals.length - 2; i++) {
            double ydiff = yvals[i] - yvals[i - 1];
            double ydiff2l = yvals[i - 3] - yvals[i - 2];
            double ydiff2r = yvals[i + 2] - yvals[i + 1];
            if (Math.abs(yvals[1] - yvals[i - 1]) / yvals[1] > 0.02) { //peak intensity threshold
                if (ydiff > 0 & previousYdiff < 0) { //look for sign changes going from - to +
                    if (ydiff2l > 0 & ydiff2r > 0) {
                        ymin.add(yvals[i - 1]);
                        xmin.add(xvals[0][i - 1]);
                    } else if (ydiff2l > 0) {
                        ymin.add(yvals[i - 1]);
                        xmin.add(xvals[0][i - 1]);
                    } else if (ydiff2r > 0) {
                        ymin.add(yvals[i - 1]);
                        xmin.add(xvals[0][i - 1]);
                    }
                } else if (ydiff == 0) {
                    if (ydiff2l > 0 & ydiff2r > 0) {
                        double yavg = (yvals[i] + yvals[i - 1]) / 2;
                        double xavg = (xvals[0][i] + xvals[0][i - 1]) / 2;
                        ymin.add(yavg);
                        xmin.add(xavg);
                    } //else if (ydiff2l > 0) {
//                        double yavg = (yvals[i] + yvals[i - 1]) / 2;
//                        double xavg = (xvals[i] + xvals[i - 1]) / 2;
//                        ymin.add(yavg);
//                        xmin.add(xavg);
//                    } else if (ydiff2r > 0) {
//                        double yavg = (yvals[i] + yvals[i - 1]) / 2;
//                        double xavg = (xvals[i] + xvals[i - 1]) / 2;
//                        ymin.add(yavg);
//                        xmin.add(xavg);
//                    }
                }
                if (fwhm.size() < ymin.size()) {
                    // FWHM calculation
                    double halfinten = (yvals[1] - yvals[i - 1]) / 2 + yvals[i - 1];
                    double distance = Math.abs(yvals[1] - halfinten);
                    int idx = 0;
                    int idx1 = 0;
                    for (int c = i - 10; c < i; c++) {
                        double cdistance = Math.abs(yvals[c] - halfinten);
                        if (cdistance < distance) {
                            idx = c;
                            distance = cdistance;
                        }
                    }
                    double halfleft = xvals[0][idx];// /(2*Math.PI);
                    double distance1 = Math.abs(yvals[1] - halfinten);
                    for (int c = i + 1; c < i + 10; c++) {
                        double cdistance1 = Math.abs(yvals[c] - halfinten);
                        if (cdistance1 < distance1) {
                            idx1 = c;
                            distance1 = cdistance1;
                        }
                    }
                    double halfright = xvals[0][idx1];// /(2*Math.PI);
                    fwhm.add(Math.abs(halfright - halfleft));
                }
            }
            previousYdiff = ydiff;
            if (xvals[1][i] != B1val) {
                break;
            }
        }

//        System.out.println("ymin = " + ymin);
        double[][] peaks = new double[2][3];
        double ymin0 = ymin.get(0);

        for (int i = 0; i < xmin.size(); i++) {
            if (ymin.get(i) < ymin0) {
                peaks[0][0] = xmin.get(ymin.indexOf(ymin0));
                peaks[0][1] = ymin.get(ymin.indexOf(ymin0));
                peaks[0][2] = fwhm.get(ymin.indexOf(ymin0));
                peaks[1][0] = xmin.get(i);
                peaks[1][1] = ymin.get(i);
                peaks[1][2] = fwhm.get(i);
                ymin0 = ymin.get(i);
            }
        }
//        for (int i=0; i<peaks.length; i++) {
//           for (int j=0; j<peaks[i].length; j++) {
//                System.out.println("peaks guess " + i + " " + j + " " + peaks[i][j]);
//            } 
//        }
        return peaks;
    }

    public static double cestPbGuess(double[][] peaks) {
        // Estimates CEST pb values from peak intensities for initial guesses for before fitting.
        // Uses the output from cestPeakGuess as the input.

        if (peaks.length > 1) {
            double[] pb = new double[peaks.length / 2];

            if (peaks.length == 2) {
                pb[0] = peaks[0][1] - peaks[1][1];
            } else {
                for (int i = 0; i < pb.length; i++) {
                    pb[i] = peaks[2 * i][1] - peaks[2 * i + 1][1];
                }
            }
//            System.out.println("pb guess = " + pb[0]);
            return pb[0];
        } else {
            return 0.0;
        }

    }

    public static double[][] cestR2Guess(double[][] peaks) {
        // Estimates CEST R2A and R2B values from peak widths for initial guesses for before fitting.
        // Uses the output from cestPeakGuess as the input.

        double pb = cestPbGuess(peaks);
        for (int i = 0; i < peaks.length; i++) {
            double p = peaks[i][2];
            double ppi = p / (2.0 * Math.PI);
            System.out.println(i + " " + peaks[i][0] + " " + peaks[i][1] + " " + p + " " + ppi);
        }

        if (peaks.length > 1) {
            double[][] r2 = new double[2][peaks.length / 2];

            for (int i = 0; i < r2[0].length; i++) {
                double awidth = peaks[2 * i + 1][2] / (2 * Math.PI);
                double bwidth = peaks[2 * i][2] / (2 * Math.PI);
                double kex = (awidth + bwidth) / 2;
                double kb = pb * kex;
                double ka = (1 - pb) * kex;
                r2[0][i] = Math.abs(awidth - ka); //R2A
                r2[1][i] = Math.abs(bwidth - kb); //R2B
            }
//            for (int i=0; i<r2.length; i++) {
//                for (int j=0; j<r2[i].length; j++) {
//                    System.out.println("R2 guess " + i + " " + j + " " + r2[i][j]);
//                }
//            }
            return r2;
        } else {
            double[][] r2 = new double[2][1];

            for (int i = 0; i < r2[0].length; i++) {
                double awidth = peaks[0][2] / (2 * Math.PI);
                r2[0][0] = awidth; //R2A
                r2[1][0] = 0.0; //R2B
            }
            return r2;
        }
    }

    public static double cestKexGuess(double[][] peaks) {
        // Estimates CEST kex values from peak widths for initial guesses for before fitting.
        // Uses the output from cestPeakGuess as the input.

        if (peaks.length > 1) {
            double[] kex = new double[peaks.length / 2];

            for (int i = 0; i < kex.length; i++) {
                double awidth = peaks[2 * i + 1][2] / (2 * Math.PI);
                double bwidth = peaks[2 * i][2] / (2 * Math.PI);
                kex[i] = (awidth + bwidth) / 2; //peaks[2 * i][2]/(2*Math.PI); //Kex
            }
//            for (int i=0; i<kex.length; i++) {
//                System.out.println("Kex guess " + i + " " + kex[i]);
//            }
            return kex[0];
        } else {
            return peaks[0][2] / (2 * Math.PI); //Kex;
        }

    }

    public static double[] cestR1Guess(double[] yvals, Double Tex) {
        // Estimates CEST R1 values from data baseline intensity and Tex for initial guesses for before fitting.
        // Reference: Palmer, A. G. "Chemical exchange in biomacromolecules: Past, present, and future." J. Mag. Res. 241 (2014) 3-17.

        double[] r1 = {-Math.log(yvals[1]) / 0.3, -Math.log(yvals[1]) / 0.3}; //{R1A, R1B}
        if (Tex != null) {
            r1[0] = -Math.log(yvals[1]) / Tex; //R1A
            r1[1] = -Math.log(yvals[1]) / Tex; //R1B
            if (Tex.equals(0.0)) {
                r1[0] = 0.0;
                r1[1] = 0.0;
            }
        }
//        for (int i=0; i<r1.length; i++) {
//           System.out.println("R1 guess " + i + " " + r1[i]); 
//        }
        return r1;
    }
}
