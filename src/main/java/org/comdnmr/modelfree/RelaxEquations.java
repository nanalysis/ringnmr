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

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import static org.comdnmr.modelfree.RelaxFit.DiffusionType;
import static org.comdnmr.modelfree.RelaxFit.DiffusionType.ANISOTROPIC;
import static org.comdnmr.modelfree.RelaxFit.DiffusionType.OBLATE;
import static org.comdnmr.modelfree.RelaxFit.DiffusionType.PROLATE;
import static org.comdnmr.modelfree.RelaxFit.DiffusionType.ISOTROPIC;

/**
 *
 * @author brucejohnson
 */
public class RelaxEquations {

    static final int S = 1;
    static final int ImS = 2;
    static final int I = 3;
    static final int IpS = 4;
    public final static double SQRT2 = Math.sqrt(2.0);

    public final static double MU0 = 4.0e-7 * Math.PI;
    public final static double GAMMA_N = -2.7116e7;
    public final static double GAMMA_C = 6.72828e7;
    public final static double GAMMA_H = 2.6752218744e8;

    public final static double PLANCK = 1.054e-34;
    public final static double R_HN = 1.02e-10;
    public final static double R_HC = 1.09e-10;
    public final static double SIGMA = 160.0e-6;
    public final static Map<String, Double> GAMMA_MAP = new HashMap<>();
    public final static Map<String, Double> R_MAP = new HashMap<>();

    static {
        GAMMA_MAP.put("H", GAMMA_H);
        GAMMA_MAP.put("N", GAMMA_N);
        GAMMA_MAP.put("C", GAMMA_C);
        R_MAP.put("HN", R_HN);
        R_MAP.put("NH", R_HN);
        R_MAP.put("HC", R_HC);
        R_MAP.put("CH", R_HC);
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
    /**
     *
     * @param sf double. 1H NMR Spectrometer frequency.
     * @param elem1 String. First element ("H" for 1H NMR).
     * @param elem2 String. Second element (C, N, etc.)
     */
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

    public static double getSF(double sf, String elemX) {
        return Math.abs(sf * GAMMA_MAP.get(elemX) / GAMMA_H);
    }

    // Note: tauM = tm in Art Palmer's code, and taui in Relax. 
    /**
     * Model Free spectral density function, J(omega), calculation using Model
     * 1.
     *
     * @param w double. The frequency, omega.
     * @param tauM double. The overall correlation time.
     * @param s2 double. The order parameter S^2.
     * @return J(w) value.
     */
    public double JModelFree(double w, double tauM, double s2) {
        double value1 = s2 / (1.0 + w * w * tauM * tauM);
        double value = 0.4 * tauM * (value1);
        return value;
    }

    // Note: tauM = tm in Art Palmer's code, and taui in Relax. 
    // tau = ts in Art Palmer's code (taue in the paper: Phys Chem Chem Phys, 2016, 18, 5839-5849), and taue in Relax.
    /**
     * Model Free spectral density function, J(omega), calculation using Model
     * 2.
     *
     * @param w double. The frequency, omega.
     * @param tau double. The internal correlation time.
     * @param tauM double. The overall correlation time.
     * @param s2 double. The order parameter S^2.
     * @return J(w) value.
     */
    public double JModelFree(double w, double tau, double tauM, double s2) {
        double value1 = s2 / (1.0 + w * w * tauM * tauM);
        double value2 = ((1.0 - s2) * (tau + tauM) * tau) / ((tau + tauM) * (tau + tauM) + w * w * tauM * tauM * tau * tau);
        double value = 0.4 * tauM * (value1 + value2);
        return value;
    }

    // Note: tauM = tm in Art Palmer's code. tau = ts in Art Palmer's code.
    /**
     * Model Free spectral density function, J(omega), calculation using Model
     * 5.
     *
     * @param w double. The frequency, omega.
     * @param tau double. The internal correlation time.
     * @param tauM double. The overall correlation time.
     * @param s2 double. The order parameter S^2.
     * @param sf2 double. The order parameter for intramolecular motions with
     * fast correlation times.
     * @return J(w) value.
     */
    public double JModelFree(double w, double tau, double tauM, double s2, double sf2) {
        double value1 = s2 / (1.0 + w * w * tauM * tauM);
        double value2 = ((sf2 - s2) * (tau + tauM) * tau) / ((tau + tauM) * (tau + tauM) + w * w * tauM * tauM * tau * tau);
        double value = 0.4 * tauM * (value1 + value2);
        return value;
    }

    /**
     * Model Free spectral density function, J(omega), calculation using Model
     * 6.
     *
     * @param w double. The frequency, omega.
     * @param tauF double. The internal fast correlation time.
     * @param tauM double. The overall correlation time.
     * @param tauS double. The internal slow correlation time.
     * @param s2 double. The order parameter S^2.
     * @param sf2 double. The order parameter for intramolecular motions with
     * fast correlation times.
     * @return J(w) value.
     */
    public double JModelFree(double w, double tauF, double tauM, double tauS, double s2, double sf2) {
        double value1 = s2 / (1.0 + w * w * tauM * tauM);
        double value2 = ((1 - sf2) * (tauF + tauM) * tauF) / ((tauF + tauM) * (tauF + tauM) + w * w * tauM * tauM * tauF * tauF);
        double value3 = ((sf2 - s2) * (tauS + tauM) * tauS) / ((tauS + tauM) * (tauS + tauM) + w * w * tauM * tauM * tauS * tauS);
        double value = 0.4 * tauM * (value1 + value2 + value3);
        return value;
    }

    // Note: tauM = tm in Art Palmer's code, and taui in Relax. 
    /**
     * Model Free spectral density function, J(omega), calculations using Model
     * 1.
     *
     * @param tauM double. The overall correlation time.
     * @param s2 double. The order parameter S^2.
     * @return double[]. Array of J(w) values.
     */
    public double[] getJModelFree(double tauM, double s2) {
        double J0 = JModelFree(0.0, tauM, s2);
        double JIplusS = JModelFree(wI + wS, tauM, s2); //R1, R2, NOE
        double JIminusS = JModelFree(wI - wS, tauM, s2); // R1, R2, NOE
        double JI = JModelFree(wI, tauM, s2); // R2
        double JS = JModelFree(wS, tauM, s2); //R1, R2
        double[] result = {J0, JS, JIminusS, JI, JIplusS};
        return result;
    }

    // Note: tauM = tm in Art Palmer's code, and taui in Relax. tau = ts in Art Palmer's code (taue in the paper), and taue in Relax.
    /**
     * Model Free spectral density function, J(omega), calculations using Model
     * 2.
     *
     * @param tau double. The internal correlation time.
     * @param tauM double. The overall correlation time.
     * @param s2 double. The order parameter S^2.
     * @return double[]. Array of J(w) values.
     */
    public double[] getJModelFree(double tau, double tauM, double s2) {
        double J0 = JModelFree(0.0, tau, tauM, s2);
        double JIplusS = JModelFree(wI + wS, tau, tauM, s2); //R1, R2, NOE
        double JIminusS = JModelFree(wI - wS, tau, tauM, s2); // R1, R2, NOE
        double JI = JModelFree(wI, tau, tauM, s2); // R2
        double JS = JModelFree(wS, tau, tauM, s2); //R1, R2
        double[] result = {J0, JS, JIminusS, JI, JIplusS};
        return result;
    }

    // Note: tauM = tm in Art Palmer's code. tau = ts in Art Palmer's code.
    /**
     * Model Free spectral density function, J(omega), calculations using Model
     * 5.
     *
     * @param tau double. The internal correlation time.
     * @param tauM double. The overall correlation time.
     * @param s2 double. The order parameter S^2.
     * @param sf2 double. The order parameter for intramolecular motions with
     * fast correlation times.
     * @return double[]. Array of J(w) values.
     */
    public double[] getJModelFree(double tau, double tauM, double s2, double sf2) {
        double J0 = JModelFree(0.0, tau, tauM, s2, sf2);
        double JIplusS = JModelFree(wI + wS, tau, tauM, s2, sf2); //R1, R2, NOE
        double JIminusS = JModelFree(wI - wS, tau, tauM, s2, sf2); // R1, R2, NOE
        double JI = JModelFree(wI, tau, tauM, s2, sf2); // R2
        double JS = JModelFree(wS, tau, tauM, s2, sf2); //R1, R2
        double[] result = {J0, JS, JIminusS, JI, JIplusS};
        return result;
    }

    /**
     * Model Free spectral density function, J(omega), calculations using Model
     * 6.
     *
     * @param tauF double. The internal fast correlation time.
     * @param tauM double. The overall correlation time.
     * @param tauS double. The internal slow correlation time.
     * @param s2 double. The order parameter S^2.
     * @param sf2 double. The order parameter for intramolecular motions with
     * fast correlation times.
     * @return double[]. Array of J(w) values.
     */
    public double[] getJModelFree(double tauF, double tauM, double tauS, double s2, double sf2) {
        double J0 = JModelFree(0.0, tauF, tauM, tauS, s2, sf2);
        double JIplusS = JModelFree(wI + wS, tauF, tauM, tauS, s2, sf2); //R1, R2, NOE
        double JIminusS = JModelFree(wI - wS, tauF, tauM, tauS, s2, sf2); // R1, R2, NOE
        double JI = JModelFree(wI, tauF, tauM, tauS, s2, sf2); // R2
        double JS = JModelFree(wS, tauF, tauM, tauS, s2, sf2); //R1, R2
        double[] result = {J0, JS, JIminusS, JI, JIplusS};
        return result;
    }

    /**
     * Calculate the d array for the diffusion J(w) calculations. Equations from
     * the SI of Berlin K.; Longhini, A.; Dayie, T. K. and Fushman, D., J.
     * Biomol NMR, 2013.
     *
     * @param diffType DiffusionType. The type of diffusion: anisotropic,
     * prolate, or oblate.
     * @param D double[][]. The diffusion matrix.
     * @return double[]. The d array.
     */
    public double[] calcDiffusiond(DiffusionType diffType, double[][] D) {
        double Dx = D[0][0];
        double Dy = D[1][1];
        double Dz = D[2][2];
//        System.out.println("Dx = " + Dx + " Dy = " + Dy + " Dz = " + Dz);
        double[] dDiff = new double[3];
        switch (diffType) {
            case ANISOTROPIC:
                dDiff = new double[5];
                double[] k = {Dy - Dx, Dz - Dx, (Dx + Dy + Dz) / 3.0, 0};
                k[3] = Math.sqrt(k[0] * k[0] - k[0] * k[1] + k[1] * k[1]);
                dDiff[0] = 4.0 * Dx + Dy + Dz;
                dDiff[1] = Dx + 4.0 * Dy + Dz;
                dDiff[2] = Dx + Dy + 4.0 * Dz;
                dDiff[3] = 6.0 * k[2] + 2.0 * k[3];
                dDiff[4] = 6.0 * k[2] - 2.0 * k[3];
                break;
            case PROLATE: {
                //Dxx = Dyy
                double Dpar = Dz;
                double Dperp = Dx;
                dDiff[0] = 5.0 * Dperp + Dpar;
                dDiff[1] = 2.0 * Dperp + 4.0 * Dpar;
                dDiff[2] = 6.0 * Dperp;
                break;
            }
            case OBLATE: {
                //Dyy = Dzz
                double Dpar = Dx;
                double Dperp = Dz;
                dDiff[0] = 5.0 * Dperp + Dpar;
                dDiff[1] = 2.0 * Dperp + 4.0 * Dpar;
                dDiff[2] = 6.0 * Dperp;
                break;
            }
            default:
                break;
        }
//        for (int i=0; i<dDiff.length; i++) {
//            System.out.println("dDiff " + i + ": " + dDiff[i]);
//        }
        return dDiff;
    }

    /**
     * Calculate the a array for the diffusion J(w) calculations. Equations from
     * the SI of Berlin K.; Longhini, A.; Dayie, T. K. and Fushman, D., J.
     * Biomol NMR, 2013.
     *
     * @param diffType DiffusionType. The type of diffusion: anisotropic,
     * prolate, or oblate.
     * @param D double[][]. The diffusion matrix.
     * @param vec double[]. The unit vector for the SI bond.
     * @return double[]. The a array.
     */
    public double[] calcDiffusiona(DiffusionType diffType, double[][] D, double[][] VT, double[] vec) {
        double Dx = D[0][0];
        double Dy = D[1][1];
        double Dz = D[2][2];
        DMatrixRMaj vec1 = new DMatrixRMaj(vec);
        DMatrixRMaj VT1 = new DMatrixRMaj(VT);
        DMatrixRMaj v = new DMatrixRMaj(vec1.numRows, vec1.numCols);
        CommonOps_DDRM.mult(VT1, vec1, v);
        double vx = v.get(0, 0); //vec[0]; //vec.getX();
        double vy = v.get(1, 0); //vec[1]; //vec.getY();
        double vz = v.get(2, 0); //vec[2]; //vec.getZ();
        double vx2 = vx * vx;
        double vy2 = vy * vy;
        double vz2 = vz * vz;
        double[] a = new double[3];
        switch (diffType) {
            case ANISOTROPIC:
                a = new double[5];
                double[] k = {Dy - Dx, Dz - Dx, (Dx + Dy + Dz) / 3.0, 0.0};
                k[3] = Math.sqrt(k[0] * k[0] - k[0] * k[1] + k[1] * k[1]);
                double[] delta = {(-k[0] - k[1]) / k[3], (2.0 * k[0] - k[1]) / k[3], (2.0 * k[1] - k[0]) / k[3]};
                if (k[0] <= 1e-12 && k[1] <= 1e-12 && k[3] <= 1e-12) {
                    delta = new double[3];
                }
//                for (int i=0; i<delta.length; i++) {
//                    System.out.print("del " + i + " = " + delta[i] + " ");
//                }
//                System.out.println();
                a[0] = 3.0 * (vy2) * (vz2);
                a[1] = 3.0 * (vx2) * (vz2);
                a[2] = 3.0 * (vx2) * (vy2);
                double p1 = 0.25 * (3.0 * ((vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2)) - 1.0);
                double val1 = delta[0] * (3 * vx2 * vx2 + 2 * a[0] - 1.0);
                double val2 = delta[1] * (3 * vy2 * vy2 + 2 * a[1] - 1.0);
                double val3 = delta[2] * (3 * vz2 * vz2 + 2 * a[2] - 1.0);
                double p2 = (1.0 / 12.0) * (val1 + val2 + val3);
                a[3] = p1 - p2;
                a[4] = p1 + p2;
                break;
            case PROLATE:
                //Dxx = Dyy
                a[0] = 3.0 * (vz2) * (1.0 - vz2);
                a[1] = 0.75 * (1.0 - vz2) * (1.0 - vz2);
                a[2] = 0.25 * (3.0 * vz2 - 1.0) * (3.0 * vz2 - 1.0);
                break;
            case OBLATE:
                //Dyy = Dzz
                a[0] = 3.0 * (vx2) * (1 - vx2);
                a[1] = 0.75 * (1 - vx2) * (1 - vx2);
                a[2] = 0.25 * (3 * vx2 - 1) * (3 * vx2 - 1);
                break;
            default:
                break;
        }
//        for (int i=0; i<a.length; i++) {
//            System.out.println("a " + i + ": " + a[i]);
//        }
        return a;
    }

    /**
     * Calculate the e array for the diffusion J(w) calculations. Equations from
     * the SI of Berlin K.; Longhini, A.; Dayie, T. K. and Fushman, D., J.
     * Biomol NMR, 2013.
     *
     * @param dDiff double[]. The d array.
     * @param tauLoc double. The internal correlation time.
     * @return double[]. The e array.
     */
    public double[] calcDiffusione(double[] dDiff, double tauLoc) {
        double[] e = new double[dDiff.length];
        for (int i = 0; i < e.length; i++) {
            e[i] = tauLoc / (dDiff[i] * tauLoc + 1);
        }
        return e;
    }

    /**
     * Model Free spectral density function, J(omega), diffusion calculation
     * using ModelFree Model 1, 2, 5, or 6. Equations from the SI of Berlin K.;
     * Longhini, A.; Dayie, T. K. and Fushman, D., J. Biomol NMR, 2013.
     *
     * @param diffType DiffusionType. The type of diffusion: anisotropic,
     * prolate, or oblate.
     * @param w double. The frequency, omega.
     * @param D double[][]. The diagonalized diffusion matrix.
     * @param VT double[][]. The transposed orthonormal matrix of the
     * eigenvectors of D.
     * @param v double[]. The unit vector for the SI bond.
     * @param s2 double. The order parameter S^2.
     * @param tauF Double. The internal fast correlation time.
     * @param sf2 Double. The order parameter for intramolecular motions with
     * fast correlation times.
     * @param tauS Double. The internal slow correlation time.
     * @return J(w) value.
     */
    public double JDiffusion(DiffusionType diffType, double w, double[][] D, double[][] VT, double[] v, double s2, Double tauF, Double sf2, Double tauS) {
        double w2 = w * w;
        double[] dDiff = calcDiffusiond(diffType, D);
        double[] a = calcDiffusiona(diffType, D, VT, v);
        double[] eF = new double[dDiff.length];
        double[] eS = new double[dDiff.length];
        if (tauF != null) {
            eF = calcDiffusione(dDiff, tauF);
            if (tauS != null) {
                eS = calcDiffusione(dDiff, tauS);
            }
        }
//        double value1 = 0.0;
//        for (int i=0; i<dDiff.length; i++) {
//            value1 += s2*(dDiff[i]*a[i])/(dDiff[i]*dDiff[i] + w*w);
//        }
//        double value = 0.4 * (value1);
        double sum = 0.0;
        double[] Df = new double[dDiff.length];
        for (int d = 0; d < Df.length; d++) {
            if (w2 > 0.0) {
                Df[d] = dDiff[d] / (dDiff[d] * dDiff[d] + w2);
            } else {
                Df[d] = 1.0 / dDiff[d];
            }
        }
        if (diffType == ISOTROPIC) {
            double Diso = (D[0][0] + D[1][1] + D[2][2]) / 3.0;
            double tauC = 1.0 / (6 * Diso);
            double value1 = s2 * tauC / (1.0 + w2 * tauC * tauC);
            double value2 = 0.0;
            double value3 = 0.0;
            if (tauF != null && sf2 == null && tauS == null) {
                double tauf = tauC * tauF / (tauC + tauF);
                value2 = (1.0 - s2) * (tauf) / (1.0 + w2 * tauf * tauf);
            } else if (tauF != null && sf2 != null && tauS == null) {
                double tauf = tauC * tauF / (tauC + tauF);
                value2 = (sf2 - s2) * (tauf) / (1.0 + w2 * tauf * tauf);
            } else if (tauF != null && sf2 != null && tauS != null) {
                double taue = tauC * tauS / (tauC + tauS);
                double tauf = tauC * tauF / (tauC + tauF);
                value2 = (sf2 - s2) * (taue) / (1.0 + w2 * taue * taue);//(eS[i]*eS[i] + w2);
                value3 = (1.0 - sf2) * (tauf) / (1.0 + w2 * tauf * tauf);//(eF[i]*eF[i] + w2);
            }
            sum = value1 + value2 + value3;
        } else {
            for (int i = 0; i < dDiff.length; i++) {
                double value1 = s2 * (Df[i] * a[i]);
                double value2 = 0.0;
                double value3 = 0.0;
                if (tauF != null && sf2 == null && tauS == null) {
                    value2 = (1.0 - s2) * (eF[i] * a[i]) / (1.0 + w2 * eF[i] * eF[i]);//(eF[i]*eF[i] + w2);
                } else if (tauF != null && sf2 != null && tauS == null) {
                    value2 = (sf2 - s2) * (eF[i] * a[i]) / (1.0 + w2 * eF[i] * eF[i]);//(eF[i]*eF[i] + w2);
                } else if (tauF != null && sf2 != null && tauS != null) {
                    value2 = (sf2 - s2) * (eS[i] * a[i]) / (1.0 + w2 * eS[i] * eS[i]);//(eS[i]*eS[i] + w2);
                    value3 = (1.0 - sf2) * (eF[i] * a[i]) / (1.0 + w2 * eF[i] * eF[i]);//(eF[i]*eF[i] + w2);
                }
                sum += value1 + value2 + value3;
            }
        }
        double value = 0.4 * sum;
        return value;
    }

    // Note: tauM = tm in Art Palmer's code, and taui in Relax. 
    /**
     * Model Free spectral density function, J(omega), calculations using
     * ModelFree Model 1, 2, 5, or 6.
     *
     * @param diffType DiffusionType. The type of diffusion: anisotropic,
     * prolate, or oblate.
     * @param D double[][]. The diagonalized diffusion matrix.
     * @param VT double[][]. The transposed orthonormal matrix of the
     * eigenvectors of D.
     * @param v double[]. The unit vector for the SI bond.
     * @param s2 double. The order parameter S^2.
     * @param tauF Double. The internal fast correlation time.
     * @param sf2 Double. The order parameter for intramolecular motions with
     * fast correlation times.
     * @param tauS Double. The internal slow correlation time.
     * @return double[]. Array of J(w) values.
     */
    public double[] getJDiffusion(DiffusionType diffType, double[][] D, double[][] VT, double[] v, double s2, Double tauF, Double sf2, Double tauS) {
        double J0 = JDiffusion(diffType, 0.0, D, VT, v, s2, tauF, sf2, tauS);
        double JIplusS = JDiffusion(diffType, wI + wS, D, VT, v, s2, tauF, sf2, tauS); //R1, R2, NOE
        double JIminusS = JDiffusion(diffType, wI - wS, D, VT, v, s2, tauF, sf2, tauS); // R1, R2, NOE
        double JI = JDiffusion(diffType, wI, D, VT, v, s2, tauF, sf2, tauS); // R2
        double JS = JDiffusion(diffType, wS, D, VT, v, s2, tauF, sf2, tauS); //R1, R2
        double[] result = {J0, JS, JIminusS, JI, JIplusS};
        return result;
    }

    /**
     * Spectral density function calculation.
     *
     * @param w double. The frequency, omega.
     * @param tau double. The correlation time.
     * @return double. The spectral density, J(w).
     */
    public double J(double w, double tau) {
        double value = 0.4 * tau / (1.0 + w * w * tau * tau);
        return value;
    }

    /**
     * Spectral density function calculation.
     *
     * @param w double. The frequency, omega.
     * @param tau double. The correlation time.
     * @param S2 double. The order parameter.
     * @return double. The spectral density, J(w).
     */
    public double J(double w, double tau, double S2) {
        double value = 0.4 * S2 * tau / (1.0 + w * w * tau * tau);
        return value;
    }

    public double[] getDiffusionConstants(String type) {
        String[] types = {"sphere", "spheroid", "ellipsoid"};
        double[][] constants = {{1}};//, 
//            {0.25*(3.0*dz2 - 1)*(3.0*dz2 - 1), 3*dz2*(1 - dz2), 0.75*(3.0*dz2 - 1)*(3.0*dz2 - 1)},
//            {0.25*(dtot - e), 3*dy2*dz2, 3*dx2*dz2, 3*dx2*dy2, 0.25*(dtot + e)}};
        return constants[Arrays.asList(types).indexOf(type)];
    }

    public double[] getCorrelationTimes(String type, double Diso, double Da, double Dr) {
        String[] types = {"sphere", "spheroid", "ellipsoid"};
        double[][] tauInv = {{6 * Diso}};//, 
//            {6*Diso - 2*Da, 6*Diso - Da, 6*Diso + 2*Da},
//            {6*Diso - 2*Da*R, 6*Diso - Da*(1 + 3*Dr), 6*Diso - Da*(1 - 3*Dr), 6*Diso + 2*Da, 6*Diso + 2*Da*R}};
        int index = Arrays.asList(types).indexOf(type);
        double[] taus = new double[tauInv[index].length];
        for (int i = 0; i < taus.length; i++) {
            taus[i] = 1 / tauInv[index][i];
        }
        return taus;
    }

    public double calcS2(double J0, double bN, double mN) {
        return (5 / 2) * Math.sqrt((J0 - bN) * mN);
    }

    /**
     * Spectral density function calculations.
     *
     * @param tau double. The correlation time.
     * @return double[]. Array of J(w) values.
     */
    public double[] getJ(double tau) {
        double J0 = J(0.0, tau); // R2
        double JIplusS = J(wI + wS, tau); //R1, R2, NOE
        double JIminusS = J(wI - wS, tau); // R1, R2, NOE
        double JI = J(wI, tau); // R2
        double JS = J(wS, tau); //R1, R2
        double[] result = {J0, JS, JIminusS, JI, JIplusS};
        return result;
    }

    /**
     * Spectral density function calculations.
     *
     * @param tau double. The correlation time.
     * @param S double. The order parameter.
     * @return double[]. Array of J(w) values.
     */
    public double[] getJ(double tau, double S) {
        double[] result = getJ(tau);
        for (int i = 0; i < result.length; i++) {
            result[i] *= S;
        }
        return result;
    }

    /**
     * R2:R1 ratio calculation.
     *
     * @param tau double. The correlation time.
     * @return double. R2/R1 value.
     */
    public double r2r1Ratio(double tau) {
        double num = 4.0 * J(0, tau) + J(wS - wI, tau) + 3.0 * J(wS, tau)
                + 6.0 * J(wI, tau) + 6.0 * J(wS + wI, tau)
                + (c2 / (3.0 * d2) * (4.0 * J(0, tau) + 3.0 * J(wS, tau)));

        double denom = 2.0 * J(wS - wI, tau) + 6.0 * J(wS, tau) + 12.0 * J(wS + wI, tau) + 2.0 * (c2 / (3.0 * d2) * J(wS, tau));
        return num / denom;
    }

    /**
     * R1 calculation
     *
     * @param tau double. The correlation time.
     * @return double. R1 value.
     */
    public double R1(double tau) {
        double dipolarContrib = d2 / 4.0 * (J(wI - wS, tau) + 3.0 * J(wS, tau) + 6.0 * J(wI + wS, tau));
        double csaContrib = c2 * J(wS, tau);
        return dipolarContrib + csaContrib;
    }

    /**
     * R2 calculation
     *
     * @param tau double. The correlation time.
     * @return double. R2 value.
     */
    public double R2(double tau) {
        double dipolarContrib = d2 / 8.0 * (4.0 * J(0.0, tau) + J(wI - wS, tau) + 3.0 * J(wS, tau)
                + 6.0 * J(wI, tau) + 6.0 * J(wI + wS, tau));
        double csaContrib = c2 / 6 * (4.0 * J(0.0, tau) + 3.0 * J(wS, tau));
        return dipolarContrib + csaContrib;
    }

    /**
     * NOE calculation
     *
     * @param tau double. The correlation time.
     * @return double. NOE value.
     */
    public double NOE(double tau) {
        double R1 = R1(tau);
        return 1.0 + (d2 / (4.0 * R1)) * (gammaI / gammaS)
                * (6.0 * J(wI + wS, tau) - J(wI - wS, tau));
    }

    /**
     * R1 calculation
     *
     * @param J double[]. Array of spectral density function values, J(w).
     * @return double. R1 value.
     */
    public double R1(double[] J) {
        double dipolarContrib = d2 / 4.0 * (J[ImS] + 3.0 * J[S]
                + 6.0 * J[IpS]);
        double csaContrib = c2 * J[S];
        return dipolarContrib + csaContrib;
    }

    /**
     * R2 calculation
     *
     * @param J double[]. Array of spectral density function values, J(w).
     * @param Rex double. Rate of exchange, Rex, value.
     * @return double. R2 value.
     */
    public double R2(double[] J, double Rex) {
        double dipolarContrib = d2 / 8.0 * (4.0 * J[0] + J[ImS]
                + 3.0 * J[S]
                + 6.0 * J[I] + 6.0 * J[IpS]);
        double csaContrib = c2 / 6 * (4.0 * J[0] + 3.0 * J[S]);
        return dipolarContrib + csaContrib + Rex;
    }

    /**
     * NOE calculation
     *
     * @param J double[]. Array of spectral density function values.
     * @return double. NOE value.
     */
    public double NOE(double[] J) {
        double R1 = R1(J);
        return 1.0 + (d2 / (4.0 * R1)) * (gammaI / gammaS)
                * (6.0 * J[IpS] - J[ImS]);
    }

    /**
     * Sigma calculation
     *
     * @param J double[]. Array of spectral density function values.
     * @return double. Sigma value.
     */
    public double sigmaSI(double[] J) {
        double R1 = R1(J);
        double NOE = NOE(J);
        return R1 * (NOE - 1) * (gammaS / gammaI);
    }

    /**
     * Gamma calculation
     *
     * @param J double[]. Array of spectral density function values, J(w).
     * @param Rex double. Rate of exchange, Rex, value.
     * @return double. Gamma value.
     */
    public double gamma(double[] J, double Rex) {
        double R1 = R1(J);
        double R2 = R2(J, Rex);
        double sigmaSI = sigmaSI(J);
        return R2 - 0.5 * R1 - 0.454 * sigmaSI;
    }

    public double TRACTdeltaAlphaBeta(double tauC) {
        double B0 = wI / GAMMA_H;

        double ddN = 160.0e-6;
        double theta = 17.0 * Math.PI / 180.0;
        double p = MU0 * gammaI * gammaS * PLANCK / (8.0 * Math.PI * SQRT2 * r * r * r);
        double dN = gammaS * B0 * ddN / (3.0 * SQRT2);

        double cosTheta = Math.cos(theta);
        double[] J = getJ(tauC);
        double nuxy2 = 2.0 * p * dN * (4.0 * J[0] + 3.0 * J[S]) * (3.0 * cosTheta * cosTheta - 1.0);
        return nuxy2;
    }

    /**
     * Calculate rhoExp from Eqn 6 in Berlin K.; Longhini, A.; Dayie, T. K. and
     * Fushman, D., J. Biomol NMR, 2013.
     *
     * @param R1 double. Experimentally measured R1 value.
     * @param R2 double. Experimentally measured R2 value.
     * @param NOE double. Experimentally measured NOE value.
     * @param J double[]. Array of spectral density function values, J(w).
     * @return double. RhoExp.
     */
    public double calcRhoExp(double R1, double R2, double NOE, double[] J) {
        double w1 = (J[ImS] + 6 * J[IpS]) / (6 * J[IpS] - J[ImS]);
        double w2 = (6 * J[I]) / (6 * J[IpS] - J[ImS]);
        double f = (gammaS / gammaI) * R1 * (NOE - 1);
        double rhoExp = (2 * R2 - R1 - w2 * f) / (R1 - w1 * f);
        return rhoExp;
    }

    /**
     * Calculate the standard deviation error in rhoExp.
     *
     * @param R1 double. Experimentally measured R1 value.
     * @param R2 double. Experimentally measured R2 value.
     * @param NOE double. Experimentally measured NOE value.
     * @param J double[]. Array of spectral density function values, J(w).
     * @param R1err double. Error value for the experimentally measured R1
     * value.
     * @param R2err double. Error value for the experimentally measured R2
     * value.
     * @param NOEerr double. Error value for the experimentally measured NOE
     * value.
     * @param rhoExp double. Calculated rhoExp value.
     * @return double. RhoExp sigma squared error.
     */
    public double calcRhoExpError(double R1, double R2, double NOE, double[] J, double R1err, double R2err, double NOEerr, double rhoExp) {
        double w1 = (J[ImS] + 6 * J[IpS]) / (6 * J[IpS] - J[ImS]);
        double w2 = (6 * J[I]) / (6 * J[IpS] - J[ImS]);
        double f = (gammaS / gammaI) * R1 * (NOE - 1);
        double rhoExp1 = 2 * R2 - R1 - w2 * f;
        double rhoExp2 = R1 - w1 * f;
        double rhoExpErr1 = Math.sqrt(2 * 2 * R2err * R2err + R1err * R1err
                + w2 * w2 * (f * f * ((R1err / R1) * (R1err / R1) + (NOEerr / NOE) * (NOEerr / NOE)) + R1err * R1err));
        double rhoExpErr2 = Math.sqrt(R1err * R1err
                + w1 * w1 * (f * f * ((R1err / R1) * (R1err / R1) + (NOEerr / NOE) * (NOEerr / NOE)) + R1err * R1err));
        double rhoExpErr = Math.sqrt(rhoExp * rhoExp * ((rhoExpErr1 / rhoExp1) * (rhoExpErr1 / rhoExp1) + (rhoExpErr2 / rhoExp2) * (rhoExpErr2 / rhoExp2)));
        return rhoExpErr;
    }

    /**
     * Calculate rhoExp from function computeRho() in RotDif code
     * (..../relax/RelaxationDatum.java).
     *
     * @param R1 double. Experimentally measured R1 value.
     * @param R2 double. Experimentally measured R2 value.
     * @param NOE double. Experimentally measured NOE value.
     * @param J double[]. Array of spectral density function values, J(w).
     * @return double. RhoExp.
     */
    public double calcRhoExpCode(double R1, double R2, double NOE, double[] J) {
        double w1 = (J[ImS] + 6 * J[IpS]) / (6 * J[IpS] - J[ImS]);
        double w2 = 0.5 * (6 * J[I] + 6 * J[IpS] + J[ImS]) / (6 * J[IpS] - J[ImS]);
        double f = (gammaS / gammaI) * R1 * (NOE - 1);
        double R1c = R1 - w1 * f;
        double R2c = R2 - w2 * f;
        double rhoExp = (2 * R2c - R1c) / (R1c);
        return rhoExp;
    }

    /**
     * Calculate rhoPred from Eqn 7 in Berlin K.; Longhini, A.; Dayie, T. K. and
     * Fushman, D., J. Biomol NMR, 2013.
     *
     * @param J double[]. Array of spectral density function values, J(w).
     * @return double. RhoExp.
     */
    public double calcRhoPred(double[] J) {
        double rhoPred = (4.0 / 3.0) * (J[0] / J[S]);
        return rhoPred;
    }

}
