package org.comdnmr.fit.calc;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

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
            r1rho[i] = kex * x / (y - z);
            r1rho[i] = (1 - sin2t) * R1Bar + sin2t * R2Bar + r1rho[i];
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
