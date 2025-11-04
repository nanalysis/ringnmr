package org.comdnmr.modelfree.models;

public class MFModelIso2sfx extends MFModelIso2sf {

    @Override
    public double[] calc(double[] omegas) {
        double tauMx = tauM * 1.0e-9;
        double tauFx = tauF * 1.0e-9;
        double tauSx = tauS * 1.0e-9;
        double[] J = new double[omegas.length];
        int j = 0;
        double ss2 = this.ss2;
        double sf2 = this.sf2 / sN;
        double s2 = ss2 * sf2;
        for (double omega : omegas) {
            double omega2 = omega * omega;
            double vM = s2 * tauMx / (1.0 + omega2 * tauMx * tauMx);
            double vMS = sf2 * (1.0 - ss2) * (tauMx * tauSx * (tauMx + tauSx))
                    / (tauMx * tauMx * tauSx * tauSx * omega2 + (tauMx + tauSx) * (tauMx + tauSx));
            double vMF = (1.0 - sf2) * ss2 * (tauMx * tauFx * (tauMx + tauFx))
                    / (tauMx * tauMx * tauFx * tauFx * omega2 + (tauMx + tauFx) * (tauMx + tauFx));
            double tauMFS = tauFx * (tauMx + tauSx) + tauMx * tauSx;
            double vMFS = (1.0 - sf2) * (1.0 - ss2) * (tauFx * tauMx * tauSx * tauMFS)
                    / (tauFx * tauFx * tauMx * tauMx * tauSx * tauSx * omega2
                    + tauMFS * tauMFS);
            J[j++] = 0.4 * (vM + vMS + vMF + vMFS);
        }
        complexityS = Math.abs(1.0 - sf2) + Math.abs(1.0 - ss2);
        complexityTau =
                Math.log10((tauSx + tauPrime) / tauPrime) +
                        Math.log10((tauFx + tauPrime) / tauPrime);
        return J;

    }
}
