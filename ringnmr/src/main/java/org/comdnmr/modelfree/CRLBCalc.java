package org.comdnmr.modelfree;

import org.apache.commons.math3.linear.*;
import org.comdnmr.modelfree.models.MFModelIso2sf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.comdnmr.modelfree.models.MFModelIso2sf.TAU_PRIME;

public class CRLBCalc {

    private boolean atLowerBound(MFModelIso2sf.ORDERPARS op, MFModelIso2sf model, double tol) {

        double v  = op.getParam(model);
        double lo = op.getLowerBound();          // wherever your optimiser bounds live
        return (v - lo) <= tol * Math.max(Math.abs(lo), 1e-3);
    }
    /**
     * CRLB per parameter. Dimensionless for sf2/ss2, SECONDS for the taus.
     */
    public double[] cramerRao(MFModelIso2sf model, double[] omegas,
                              double[] jValues, double[] sigmaJ, boolean fitTauM) {
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
        if (hasSlowAmp)  active.add(MFModelIso2sf.ORDERPARS.TAUS);
        if (hasFast)     active.add(MFModelIso2sf.ORDERPARS.TAUF);
        if (fitTauM)     active.add(MFModelIso2sf.ORDERPARS.TAUM);

        int p = active.size();
        int m = omegas.length;
        double[][] jac = new double[m][p];
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
                up = Math.min(up, 1.0);
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
        double[] nrm = new double[p];
        for (int c = 0; c < p; c++) {
            double s = 0.0;
            for (int i = 0; i < m; i++) {
                s += w[i][c] * w[i][c];
            }
            nrm[c] = Math.sqrt(s);
            if (nrm[c] <= 0.0) {
                nrm[c] = 1.0;
            }
            for (int i = 0; i < m; i++) {
                w[i][c] /= nrm[c];
            }
        }

        double[][] F = new double[p][p];
        for (int a = 0; a < p; a++) {
            for (int b = a; b < p; b++) {
                double s = 0.0;
                for (int i = 0; i < m; i++) {
                    s += w[i][a] * w[i][b];
                }
                F[a][b] = s;
                F[b][a] = s;
            }
        }
        RealMatrix fisher = new Array2DRowRealMatrix(F, false);
        double cond = new SingularValueDecomposition(fisher).getConditionNumber();
        if (cond > 1e10) {
            return null;
        }
        RealMatrix cov;
        try {
            cov = new CholeskyDecomposition(fisher).getSolver().getInverse();
        } catch (NonPositiveDefiniteMatrixException | SingularMatrixException _) {
            return null;
        }

        int nPar = 5;
        double[] crlb = new double[nPar];
        Arrays.fill(crlb, Double.POSITIVE_INFINITY);
        for (int c = 0; c < p; c++) {
            int k = active.get(c).index(model.fitTau());
            crlb[k] = Math.sqrt(cov.getEntry(c, c)) / nrm[c];
        }
        return crlb;
    }
}
