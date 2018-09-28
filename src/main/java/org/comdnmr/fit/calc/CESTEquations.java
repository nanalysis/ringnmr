/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

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

    public static double[] r1rhoPerturbationNoEx(double[] omega, double pb, double kex, double[] deltaA, double R1A, double R2A) {
        int size = omega.length;
        double pa = 1.0 - pb;

        double k1 = pb * kex;
        double km1 = pa * kex;

        double dR = Math.abs(R2A);

        double[] r1rho = new double[size];
        for (int i = 0; i < size; i++) {
            double dw = deltaA[i];
            double weA = Math.sqrt(omega[i] * omega[i] + deltaA[i] * deltaA[i]);
//            double weB = Math.sqrt(omega[i] * omega[i] + deltaB[i] * deltaB[i]);
            double sin2t = (omega[i] / weA) * (omega[i] / weA);
            double x = ((dw * dw) + (dR * dR)) * km1 + dR * ((weA * weA) + (km1 * km1));
            double y = dR * (omega[i] * omega[i]);
            double REx = k1 * x / y;
            r1rho[i] = (1 - sin2t) * R1A + sin2t * R2A + sin2t * REx;
        }
        return r1rho;
    }

    public static double[] r1rhoPerturbationNoEx2(double[] omega, double[] deltaA, double R1A, double R2A) {
        int size = omega.length;

        double dR = Math.abs(R2A);

        double[] r1rho = new double[size];
        for (int i = 0; i < size; i++) {
            double dw = deltaA[i];
            double weA = Math.sqrt(omega[i] * omega[i] + deltaA[i] * deltaA[i]);
//            double weB = Math.sqrt(omega[i] * omega[i] + deltaB[i] * deltaB[i]);
            double sin2t = (omega[i] / weA) * (omega[i] / weA);
            double x = dR * ((weA * weA));
            double y = dR * (omega[i] * omega[i]);
            r1rho[i] = (1 - sin2t) * R1A + sin2t * R2A;
        }
        return r1rho;
    }

    public static double[] cestR1rhoPerturbationNoEx2(double[][] X, double deltaA0, double R1A, double R2A) {
        double[] omegarf = X[0];
        double[] omega1 = X[1];
        double[] Tex = X[2];

        double trad = Tex[0];//0.3;

        double pa = 1.0;
        double[] deltaA = new double[omegarf.length];
        double[] deltaB = new double[omegarf.length];
        double[] omegaBar = new double[omegarf.length];
        for (int i = 0; i < omegarf.length; i++) {
            deltaA[i] = deltaA0 - omegarf[i];
            omegaBar[i] = pa * deltaA[i];
        }

        double[] we = new double[omegaBar.length];
        double[] cos2t = new double[omegaBar.length];
        for (int i = 0; i < omegaBar.length; i++) {
            we[i] = Math.sqrt(omega1[i] * omega1[i] + omegaBar[i] * omegaBar[i]);
            cos2t[i] = (omegaBar[i] / we[i]) * (omegaBar[i] / we[i]);
        }

        double[] cest = new double[cos2t.length];

        double[] r1rho = r1rhoPerturbationNoEx2(omega1, deltaA, R1A, R2A);
        //double[] cest = new double[cos2t.length];
        for (int i = 0; i < cos2t.length; i++) {
            cest[i] = cos2t[i] * Math.exp(-trad * r1rho[i]);
        }
        return cest;
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

        } else if (approx == "trottnoex") {
            double[] r1rho = r1rhoPerturbationNoEx(omega1, pb, kex, deltaA, R1A, R2A);
            //double[] cest = new double[cos2t.length];
            for (int i = 0; i < cos2t.length; i++) {
                cest[i] = cos2t[i] * Math.exp(-trad * r1rho[i]);
            }

        } else if (approx == "trottnoex2") {
            double[] r1rho = r1rhoPerturbationNoEx2(omega1, deltaA, R1A, R2A);
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

    static class Peak {

        final double position;
        final double depth;
        final double width;
        final double widthLB;
        final double widthUB;
        final int pkInd;

        Double getDepth() {
            return depth;
        }

        Double getPosition() {
            return depth;
        }

        Peak(double position, double depth, double width, double widthLB, double widthUB, int pkInd) {
            this.position = position;
            this.depth = depth;
            this.width = width;
            this.widthLB = widthLB;
            this.widthUB = widthUB;
            this.pkInd = pkInd;
        }
    }

    public static double[] smoothCEST(double[] vec, int j1, int j2, int order, int smoothSize, double[] X) {
        // Applies a Savitzky-Golay filter to the data for smoothing.
        SavitzkyGolay sg = new SavitzkyGolay();
        double[] smoothed = new double[X.length];
        try {
            smoothed = sg.runningSavitzkyGolaySmoothing(vec, j1, j2, order, smoothSize, X);
        } catch (Exception e) {
            System.out.println("Smooth-size should be one of 5,7,9,...,25");
        }
        return smoothed;
    }

    static double[] getBaseline(double[] vec) {
        int winSize = 8;
        double maxValue = Double.NEGATIVE_INFINITY;
        double sDev = 0.0;
        DescriptiveStatistics stat = new DescriptiveStatistics(winSize);
        for (int i = 0; i < vec.length; i++) {
            stat.addValue(vec[i]);
            if (i >= (winSize - 1)) {
                double mean = stat.getMean();
                if (mean > maxValue) {
                    maxValue = mean;
                    sDev = stat.getStandardDeviation();
                }
            }
        }
        double[] result = {maxValue, sDev};
        return result;
    }

    public static double[][] cestPeakGuess(double[][] xvals, double[] yvals) {
        // Estimates CEST peak positions for initial guesses for before fitting.

        List<Peak> peaks = new ArrayList<>();
        double B1val = xvals[1][0];

        double[] syvals = new double[yvals.length];
        double[] baseValues = getBaseline(yvals);
        double baseline = baseValues[0];

        yvals = smoothCEST(yvals, 0, yvals.length, 3, 11, syvals);

        // A point must have a lower value than this number of points on each
        // side of the point in order to be a peak.
        int nP = 2;
        double baseRatio = 3.0;

        // A threshold to use when deciding if a point is deep enough
        // calculated as baseline - a multiple of the standard deviation estimate
        double threshold = baseline - baseValues[1] * baseRatio;
        double yMin = Double.MAX_VALUE;
        double xAtYMin = 0.0;
        for (int i = nP; i < yvals.length - nP; i++) {
            if (yvals[i] < yMin) {
                yMin = yvals[i];
                xAtYMin = xvals[0][i];
            }
            if (yvals[i] < threshold) {
                boolean ok = true;
                for (int j = i - nP; j <= (i + nP); j++) {
                    if (yvals[i] > yvals[j]) {
                        ok = false;
                        break;
                    }
                }
                int iCenter = i;
                if (ok) {
                    double halfinten = (baseline - yvals[iCenter]) / 2 + yvals[iCenter];
                    double[] halfPos = new double[2];

                    // search from peak center in both directions to find
                    // the peak width.  Find a value above and below the 
                    // half-height and interpolate to get the width on each
                    // side.
                    for (int k = 0; k < 2; k++) {
                        int iDir = k * 2 - 1; // make iDir -1, 1
                        int j = iCenter + iDir;
                        double dUp = Double.MAX_VALUE;
                        double dLow = Double.MAX_VALUE;
                        int iUp = 0;
                        int iLow = 0;
                        while ((j >= 0) && (j < yvals.length)) {
                            double delta = yvals[j] - halfinten;
                            if (delta < 0.0) {
                                if (Math.abs(delta) < dUp) {
                                    dUp = Math.abs(delta);
                                    iUp = j;
                                }
                            } else {
                                if (Math.abs(delta) < dLow) {
                                    dLow = Math.abs(delta);
                                    iLow = j;
                                }
                                break;
                            }
                            j += iDir;
                        }
                        if ((dLow == Double.MAX_VALUE) || (dUp == Double.MAX_VALUE)) {
                            ok = false;
                            break;
                        }
                        double delta = dLow + dUp;
                        halfPos[k] = xvals[0][iLow] * dUp / delta + xvals[0][iUp] * dLow / delta;
                    }
                    if (ok) {
                        double xCenter = xvals[0][iCenter];
                        double yCenter = yvals[iCenter];
                        double width = Math.abs(halfPos[0] - halfPos[1]);
                        double widthL = Math.abs(halfPos[0] - xCenter);
                        double widthR = Math.abs(halfPos[1] - xCenter);
                        Peak peak = new Peak(xCenter, yCenter, width, widthL, widthR, iCenter);
                        peaks.add(peak);
                    }
                }
            }

            if (xvals[1][i] != B1val) {
                break;
            }
        }

        peaks.sort(Comparator.comparingDouble(Peak::getDepth));
        System.out.println("min at " + xAtYMin + " " + yMin);
        for (int i = 0; i < peaks.size(); i++) {
            System.out.println("peaks guess " + i + " x = " + peaks.get(i).position);
            System.out.println("peaks guess " + i + " y = " + peaks.get(i).depth);
            System.out.println("peaks guess " + i + " fwhm = " + peaks.get(i).width);
        }
        List<Peak> peaks2 = peaks;
        if (peaks.size() >= 2) {
            peaks2 = peaks.subList(0, 2);
        } else if (peaks.size() == 1) {
            // If there is only one peak found add another peak on the side
            // with the largest width
            peaks2 = peaks.subList(0, 1);
            Peak peak = peaks2.get(0);
            double newCenter;
            if (peak.widthLB > peak.widthUB) {
                newCenter = peak.position - peak.widthLB / 2.0;
            } else {
                newCenter = peak.position + peak.widthUB / 2.0;
            }
            double newDepth = (baseline + peak.depth) / 2.0;
            Peak newPeak = new Peak(newCenter, newDepth, peak.width, peak.widthLB, peak.widthUB, peak.pkInd);
            peaks2.add(newPeak);
        }

        peaks2.sort(Comparator.comparingDouble(Peak::getDepth));
        for (int i = 0; i < peaks2.size(); i++) {
            System.out.println("peaks guess " + i + " x = " + peaks2.get(i).position);
            System.out.println("peaks guess " + i + " y = " + peaks2.get(i).depth);
            System.out.println("peaks guess " + i + " fwhm = " + peaks2.get(i).width);
        }

        double peak1diff = Math.abs(peaks2.get(0).depth - baseline);
        double peak2diff = peak1diff;
        if (peaks2.size() == 2) {
            peak2diff = Math.abs(peaks2.get(1).depth - baseline);
        }
        if (peak1diff < 0.05 || peak2diff < 0.05) {
            double[][] peaks1 = new double[1][3];
            int index = 0;
            if (peaks2.get(0).depth > peaks2.get(1).depth && peaks2.get(1).depth != 0) {
                index = 1;
            }
            peaks1[0][0] = peaks2.get(index).position;
            peaks1[0][1] = peaks2.get(index).depth;
            peaks1[0][2] = peaks2.get(index).width;
            return peaks1;
        } else {
            double[][] peaks1 = new double[peaks2.size()][3];
            for (int i = 0; i < peaks1.length; i++) {
                peaks1[i][0] = peaks2.get(i).position;
                peaks1[i][1] = peaks2.get(i).depth;
                peaks1[i][2] = peaks2.get(i).width;
            }
            return peaks1;
        }
    }

    public static double cestPbGuess(double[][] peaks, double[] yvals) {
        // Estimates CEST pb values from peak intensities for initial guesses for before fitting.
        // Uses the output from cestPeakGuess as the input.

//        System.out.println(peaks.length + " peaks found.");
//        for (int i = 0; i < peaks.length; i++) {
//            for (int j = 0; j < peaks[i].length; j++) {
//                System.out.println(i + " " + j + " " + peaks[i][j]);
//            }
//        }
        double[] baseValues = getBaseline(yvals);
        double baseline = baseValues[0];

        if (peaks.length > 1) {
            double[] pb = new double[peaks.length / 2];

            if (peaks.length == 2) {
                if (peaks[0][1] > peaks[1][1]) {
                    pb[0] = ((baseline - peaks[0][1]) / (baseline - peaks[1][1])) / 4;
                } else {
                    pb[0] = ((baseline - peaks[1][1]) / (baseline - peaks[0][1])) / 4;
                }
            } else {
                for (int i = 0; i < pb.length; i++) {
                    pb[i] = (baseline - peaks[2 * i][1]) / (baseline - peaks[2 * i + 1][1]) / 4;
                }
            }
//            System.out.println("pb guess = " + pb[0]);
            return pb[0];
        } else {
            return 0.1;
        }

    }

    public static double[][] cestR2Guess(double[][] peaks, double[] yvals) {
        // Estimates CEST R2A and R2B values from peak widths for initial guesses for before fitting.
        // Uses the output from cestPeakGuess as the input.

        double pb = cestPbGuess(peaks, yvals);

        if (peaks.length > 1) {
            double[][] r2 = new double[2][peaks.length / 2];
            double awidth = 0;
            double bwidth = 0;

            if (peaks.length == 2) {
                if (peaks[0][1] > peaks[1][1]) {
                    awidth = peaks[1][2] / (2 * Math.PI);
                    bwidth = peaks[0][2] / (2 * Math.PI);
                } else {
                    awidth = peaks[0][2] / (2 * Math.PI);
                    bwidth = peaks[1][2] / (2 * Math.PI);
                }
                double kex = (awidth + bwidth) / 2;
                double kb = pb * kex;
                double ka = (1 - pb) * kex;
                r2[0][0] = Math.abs(awidth - ka); //R2A
                r2[1][0] = Math.abs(bwidth - kb); //R2B
            } else {
                for (int i = 0; i < r2[0].length; i++) {
                    awidth = peaks[2 * i + 1][2] / (2 * Math.PI);
                    bwidth = peaks[2 * i][2] / (2 * Math.PI);
                    double kex = (awidth + bwidth) / 2;
                    double kb = pb * kex;
                    double ka = (1 - pb) * kex;
                    r2[0][i] = Math.abs(awidth - ka); //R2A
                    r2[1][i] = Math.abs(bwidth - kb); //R2B
                }
            }
//            for (int i = 0; i < r2.length; i++) {
//                for (int j = 0; j < r2[i].length; j++) {
//                    System.out.println("R2 guess " + i + " " + j + " " + r2[i][j]);
//                }
//            }
            return r2;
        } else {
            double[][] r2 = new double[2][1];

            for (int i = 0; i < r2[0].length; i++) {
                double awidth = peaks[0][2] / (2 * Math.PI);
                r2[0][0] = awidth; //R2A
                r2[1][0] = awidth; //R2B
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

        double[] baseValues = getBaseline(yvals);
        double baseline = baseValues[0];
        double[] r1 = {-Math.log(baseline) / 0.3, -Math.log(baseline) / 0.3}; //{R1A, R1B}
        if (Tex != null) {
            r1[0] = -Math.log(baseline) / Tex; //R1A
            r1[1] = -Math.log(baseline) / Tex; //R1B
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
