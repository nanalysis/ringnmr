package org.comdnmr.modelfree;

import org.apache.commons.math3.linear.*;
import org.comdnmr.modelfree.models.MFModelIso2sf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.comdnmr.modelfree.models.MFModelIso2sf.TAU_PRIME;

public class CRLBCalc {
    private static final double RCOND = 1.0e-12;  // eigenvalue floor, relative to the largest
    private static final double VIF_MAX = 1.0e8;    // variance-inflation cutoff

    private boolean atLowerBound(MFModelIso2sf.ORDERPARS op, MFModelIso2sf model, double tol) {

        double v = op.getParam(model);
        double lo = op.getLowerBound();          // wherever your optimiser bounds live
        return (v - lo) <= tol * Math.max(Math.abs(lo), 1e-3);
    }

    /**
     * Lower Cholesky factor of a 3x3 correlation matrix; null if not positive definite.
     */
    private static double[][] chol3(double[][] r) {
        double[][] l = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j <= i; j++) {
                double s = r[i][j];
                for (int k = 0; k < j; k++) {
                    s -= l[i][k] * l[j][k];
                }
                if (i == j) {
                    if (!(s > 0.0)) {
                        return null;
                    }
                    l[i][i] = Math.sqrt(s);
                } else {
                    l[i][j] = s / l[j][j];
                }
            }
        }
        return l;
    }

    /**
     * CRLB per parameter. Dimensionless for sf2/ss2, Nanoseconds for the taus.
     */
    public double[] cramerRao(MFModelIso2sf model, double[] omegas,
                               double[][][] cov3, boolean fitTauM) {
        final double EPS = 1e-9;
        final double BOUND_TOL = 1e-6;

        boolean hasFast = (model.getTauF() > 0.0)
                && (model.getSf2() / model.getSN() < 1.0 - EPS)
                && !atLowerBound(MFModelIso2sf.ORDERPARS.TAUF, model, BOUND_TOL);

        boolean inModelSlow = model.getTauS() > 0.0;
        boolean hasSlowAmp = inModelSlow
                && ((model.getSf2() / model.getSN()) * (1.0 - model.getSs2()) > 1.0e-4)
                && !atLowerBound(MFModelIso2sf.ORDERPARS.TAUS, model, BOUND_TOL);

        List<MFModelIso2sf.ORDERPARS> active = new ArrayList<>();
        active.add(MFModelIso2sf.ORDERPARS.SF2);
        if (inModelSlow) active.add(MFModelIso2sf.ORDERPARS.SS2);
        if (hasSlowAmp) active.add(MFModelIso2sf.ORDERPARS.TAUS);
        if (hasFast) active.add(MFModelIso2sf.ORDERPARS.TAUF);
        if (fitTauM) active.add(MFModelIso2sf.ORDERPARS.TAUM);

        int p = active.size();
        int m = omegas.length;
        double[] sigmaJ = new double[m];
        for (int i = 0; i < m; i++) {
            sigmaJ[i] = Math.sqrt(cov3[i / 3][i % 3][i % 3]);
        }        double[][] jac = new double[m][p];
        double[][] w = new double[m][p];

        for (int c = 0; c < p; c++) {
            MFModelIso2sf.ORDERPARS op = active.get(c);
            boolean isOrderPar = op.name().startsWith("S");
            double v = op.getParam(model);
            double h = isOrderPar ? 1.0e-6 : 1.0e-4 * Math.abs(v);
            if (!isOrderPar && (h < 1.0e-4 * TAU_PRIME)) {
                h = 1.0e-4 * TAU_PRIME;
            }

            double up = v + h;
            double dn = v - h;
            if (isOrderPar) {
                double hi = op == MFModelIso2sf.ORDERPARS.SF2 ? model.getSN() : 1.0;
                up = Math.min(up, hi);
                dn = Math.max(dn, 1e-6);
            }
            op.setParam(model, up);
            double[] ju = model.calc(omegas);
            op.setParam(model, dn);
            double[] jd = model.calc(omegas);
            op.setParam(model, v);

            double denom = up - dn;
            for (int i = 0; i < m; i++) {
                jac[i][c] = (ju[i] - jd[i]) / denom;
                w[i][c] = jac[i][c] / sigmaJ[i];
            }
        }
        if (cov3.length * 3 != m) {
            throw new IllegalArgumentException(
                    "cov3 has " + cov3.length + " fields, omegas has " + m + " entries");
        }

        double[][] F = new double[p][p];
        for (int f = 0; f < cov3.length; f++) {
            double[][] corrF = new double[3][3];
            for (int i = 0; i < 3; i++) {
                corrF[i][i] = 1.0;                       // exactly 1, for Cholesky
                for (int j = i + 1; j < 3; j++) {
                    double r = cov3[f][i][j]
                            / Math.sqrt(cov3[f][i][i] * cov3[f][j][j]);
                    corrF[i][j] = r;
                    corrF[j][i] = r;
                }
            }
            double[][] lF = chol3(corrF);      // null -> fall back to diagonal weighting
            int base = f * 3;

            double[][] u = new double[3][p];   // u = L^-1 * (this field's rows of w)
            for (int c = 0; c < p; c++) {
                for (int i = 0; i < 3; i++) {
                    double s = w[base + i][c];
                    if (lF == null) {
                        u[i][c] = s;
                    } else {
                        for (int k = 0; k < i; k++) {
                            s -= lF[i][k] * u[k][c];
                        }
                        u[i][c] = s / lF[i][i];
                    }
                }
            }
            for (int a = 0; a < p; a++) {
                for (int b = a; b < p; b++) {
                    double s = 0.0;
                    for (int i = 0; i < 3; i++) {
                        s += u[i][a] * u[i][b];
                    }
                    F[a][b] += s;
                    if (b != a) {
                        F[b][a] += s;
                    }
                }
            }
        }

        double[] nrm = new double[p];
        for (int c = 0; c < p; c++) {
            nrm[c] = Math.sqrt(F[c][c]);
            if (!(nrm[c] > 0.0)) {
                nrm[c] = 1.0;
            }
        }
        for (int a = 0; a < p; a++) {
            for (int b = 0; b < p; b++) {
                F[a][b] /= (nrm[a] * nrm[b]);
            }
        }
        int nPar = 5;
        double[] crlb = new double[nPar];
        Arrays.fill(crlb, Double.POSITIVE_INFINITY);

        EigenDecomposition eig =
                new EigenDecomposition(new Array2DRowRealMatrix(F, false));
        double[] lambda = eig.getRealEigenvalues();
        double lamMax = 0.0;
        for (double l : lambda) {
            lamMax = Math.max(lamMax, l);
        }
        if (!(lamMax > 0.0)) {
            return crlb;                       // nothing is determined
        }
        double floor = lamMax * RCOND;

        for (int c = 0; c < p; c++) {
            double variance = 0.0;
            for (int i = 0; i < p; i++) {
                double v = eig.getEigenvector(i).getEntry(c);
                variance += v * v / Math.max(lambda[i], floor);
            }
            if (variance <= VIF_MAX) {
                crlb[active.get(c).index()] = Math.sqrt(variance) / nrm[c];
            }
        }
        return crlb;
    }
}
