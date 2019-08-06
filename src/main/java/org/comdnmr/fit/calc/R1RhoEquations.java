package org.comdnmr.fit.calc;

import org.comdnmr.util.MtxExp;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

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

    public static double r1rhoExact0(double tdelay, double omega, double pB, double kex, double deltaA, double deltaB, double R1A, double R1B, double R2A, double R2B) {
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
        double pA = 1.0 - pB;
        double kAB = pB * kex;
        double kBA = pA * kex;

        double theta = Math.atan2(omega, deltaA);
        double cosA = Math.cos(theta);
        double sinA = Math.sin(theta);

        double[] m0 = {pA * sinA, 0.0, pA * cosA, 0.0, 0.0, 0.0};
        double[] m1 = {sinA, 0.0, cosA, 0.0, 0.0, 0.0};

        double[][] K = {
            {-kAB, 0, 0, kBA, 0, 0},
            {0, -kAB, 0, 0, kBA, 0},
            {0, 0, -kAB, 0, 0, kBA},
            {kAB, 0, 0, -kBA, 0, 0},
            {0, kAB, 0, 0, -kBA, 0},
            {0, 0, kAB, 0, 0, -kBA}};

        double[][] La = {
            {-R2A, -deltaA, 0, 0, 0, 0},
            {deltaA, -R2A, -omega, 0, 0, 0},
            {0, omega, -R1A, 0, 0, 0},
            {0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0}};

        double[][] Lb = {
            {0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0},
            {0, 0, 0, -R2B, -deltaB, 0},
            {0, 0, 0, deltaB, -R2B, -omega},
            {0, 0, 0, 0, omega, -R1B}};

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
        Array2DRowRealMatrix aM = new Array2DRowRealMatrix(at);
        double[] v = aM.operate(m0);
        double magA = 0.0;
        double magA0 = 0.0;
        for (int i = 0; i < v.length; i++) {
            magA += m1[i] * v[i];
            magA0 += m1[i] * m0[i];
        }

//        double magA0 = m0[2] * m1[2] + m0[5] * m1[5];
//        double magA = at[2][2] * m0[2] * m1[2] + at[2][5] * m0[5] * m1[5];
//        magA = Math.abs(magA);
        double r1rho = -Math.log(magA / magA0) / tdelay;
//        System.out.println(omega + " " + at[2][2] + " " + at[2][5] + " " + magA0 + " " + magA + " " + r1rho);

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
        double[][] La = {{-R2A, -deltaA, 0, 0, 0, 0}, {deltaA, -R2A, -omega, 0, 0, 0}, {0, omega, -R1A, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
        double[][] Lb = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, -R2B, -deltaB, 0}, {0, 0, 0, deltaB, -R2B, -omega}, {0, 0, 0, 0, omega, -R1B}};
        double[][] K = {{-k1, 0, 0, km1, 0, 0}, {0, -k1, 0, 0, km1, 0}, {0, 0, -k1, 0, 0, km1}, {k1, 0, 0, -km1, 0, 0}, {0, k1, 0, 0, -km1, 0}, {0, 0, k1, 0, 0, -km1}};
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
