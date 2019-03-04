package org.comdnmr.fit.calc;

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

    static double[] getBaseline(double[] vec) {
        int winSize = 8;
        double maxValue = Double.POSITIVE_INFINITY;
        double sDev = 0.0;
        DescriptiveStatistics stat = new DescriptiveStatistics(winSize);
        for (int i = 0; i < vec.length; i++) {
            stat.addValue(vec[i]);
            if (i >= (winSize - 1)) {
                double mean = stat.getMean();
                if (mean < maxValue) {
                    maxValue = mean;
                    sDev = stat.getStandardDeviation();
                }
            }
        }
        double[] result = {maxValue, sDev};
        return result;
    }

    public static double[][] r1rhoPeakGuess(double[] xvals, double[] yvals, double field) {
        // Estimates CEST peak positions for initial guesses for before fitting.

        List<CESTEquations.Peak> peaks = new ArrayList<>();

        double[] syvals = new double[yvals.length];
        double[] baseValues = getBaseline(yvals);
        double baseline = baseValues[0];
        int smoothSize;
        if (yvals.length < 20) {
            smoothSize = 0;
        } else if (yvals.length < 30) {
            smoothSize = 5;
        } else if (yvals.length < 40) {
            smoothSize = 7;
        } else if (yvals.length < 50) {
            smoothSize = 9;
        } else {
            smoothSize = 11;
        }

//        if (smoothSize != 0) {
//            yvals = smoothCEST(yvals, 0, yvals.length, 3, smoothSize, syvals);
//        }
        // A point must have a higher value than this number of points on each
        // side of the point in order to be a peak.
        int nP = 2;
        double baseRatio = 3.0;

        // A threshold to use when deciding if a point is deep enough
        // calculated as baseline + a multiple of the standard deviation estimate
        double threshold = baseline + baseValues[1] * baseRatio;
//        System.out.println("baseline = " + baseline);
//        System.out.println("threshold = " + threshold);
        double yMax = Double.MIN_VALUE;
        double xAtYMax = 0.0;
        for (int i = nP; i < yvals.length - nP; i++) {
            if (yvals[i] > yMax) {
                yMax = yvals[i];
                xAtYMax = xvals[i];
            }
            if (yvals[i] > threshold) {
                boolean ok = true;
                for (int j = i - nP; j <= (i + nP); j++) {
//                    System.out.println("i, xi, yi; j, xj, yj = " + i + ", " + xvals[0][i] + ", " + yvals[i] + "; " + j + ", " + xvals[0][j] + ", " + yvals[j]);
                    if (yvals[i] < yvals[j]) {
                        ok = false;
                        break;
                    }
                }
                int iCenter = i;
                if (ok) {
//                    System.out.println("iCenter, xCenter = " + iCenter + ", " + xvals[0][iCenter]);
                    double halfinten = (yvals[iCenter] - baseline) / 2 + baseline;
//                    System.out.println("halfinten = " + halfinten);
                    double[] halfPos = new double[2];

                    // search from peak center in both directions to find
                    // the peak width.  Find a value above and below the 
                    // half-height and interpolate to get the width on each
                    // side.
                    for (int k = 0; k < 2; k++) {
                        int iDir = k * 2 - 1; // make iDir -1, 1
                        int j = iCenter + iDir;
                        double dUp = Double.MIN_VALUE;
                        double dLow = Double.MIN_VALUE;
                        int iUp = 0;
                        int iLow = 0;
                        while ((j >= 0) && (j < yvals.length)) {
                            double delta = yvals[j] - halfinten;
//                            System.out.println("delta = " + delta);
                            if (delta > 0.0) {
                                if (Math.abs(delta) > dUp) {
                                    dUp = Math.abs(delta);
                                    iUp = j;
                                }
                            } else {
                                if (Math.abs(delta) > dLow) {
                                    dLow = Math.abs(delta);
                                    iLow = j;
                                }
                                break;
                            }
                            j += iDir;
                        }
                        if ((dLow == Double.MIN_VALUE) || (dUp == Double.MIN_VALUE)) {
                            ok = false;
                            break;
                        }
                        double delta = dLow + dUp;
                        halfPos[k] = xvals[iLow] * dUp / delta + xvals[iUp] * dLow / delta;
                    }
                    if (ok) {
                        double xCenter = xvals[iCenter];
                        double yCenter = yvals[iCenter];
                        double width = Math.abs(halfPos[0] - halfPos[1]) * field;
                        double widthL = Math.abs(halfPos[0] - xCenter) * field;
                        double widthR = Math.abs(halfPos[1] - xCenter) * field;
                        CESTEquations.Peak peak = new CESTEquations.Peak(xCenter, yCenter, width, widthL, widthR, iCenter);
                        peaks.add(peak);
                    }
                }
            }
        }

        peaks.sort(Comparator.comparingDouble(CESTEquations.Peak::getDepth));
        System.out.println("max at " + xAtYMax + " " + yMax);
        for (int i = 0; i < peaks.size(); i++) {
            System.out.println("orig peaks guess " + i + " x = " + peaks.get(i).position);
            System.out.println("orig peaks guess " + i + " y = " + peaks.get(i).depth);
            System.out.println("orig peaks guess " + i + " fwhm = " + peaks.get(i).width);
        }
        List<CESTEquations.Peak> peaks2 = peaks;
        if (peaks.size() >= 2) {
            peaks2 = peaks.subList(0, 2);
        } else if (peaks.size() == 1) {
            // If there is only one peak found add another peak on the side
            // with the largest width
            peaks2 = peaks.subList(0, 1);
            CESTEquations.Peak peak = peaks2.get(0);
            double newCenter;
            if (peak.widthLB > peak.widthUB) {
                newCenter = peak.position - peak.widthLB / field / 2.0;
            } else {
                newCenter = peak.position + peak.widthUB / field / 2.0;
            }
            double newDepth = (baseline + peak.depth) / 2.0;
            CESTEquations.Peak newPeak = new CESTEquations.Peak(newCenter, newDepth, peak.width, peak.widthLB, peak.widthUB, peak.pkInd);
            peaks2.add(newPeak);
        }

        peaks2.sort(Comparator.comparingDouble(CESTEquations.Peak::getDepth));
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

    public static double r1rhoPbGuess(double[][] peaks, double[] yvals) {
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
                    pb[0] = ((peaks[0][1] - baseline) / (peaks[1][1] - baseline)) / 40;
                } else {
                    pb[0] = ((peaks[1][1] - baseline) / (peaks[0][1] - baseline)) / 40;
                }
            } else {
                for (int i = 0; i < pb.length; i++) {
                    pb[i] = (peaks[2 * i][1] - baseline) / (peaks[2 * i + 1][1] - baseline) / 40;
                }
            }
//            System.out.println("pb guess = " + pb[0]);
            if (pb[0] > 0.25) {
                pb[0] = 0.25;
            }
            return pb[0];
        } else {
            return 0.1;
        }

    }

    public static double[] r1rhoR1Guess(double[] yvals, Double Tex) {
        // Estimates CEST R1 values from data baseline intensity and Tex for initial guesses for before fitting.
        // Reference: Palmer, A. G. "Chemical exchange in biomacromolecules: Past, present, and future." J. Mag. Res. 241 (2014) 3-17.

        double[] baseValues = getBaseline(yvals);
        double baseline = baseValues[0];
        double[] r1 = {baseline, baseline}; //{R1A, R1B}
//        if (Tex != null) {
//            r1[0] = -Math.log(baseline) / Tex; //R1A
//            r1[1] = -Math.log(baseline) / Tex; //R1B
//            if (Tex.equals(0.0)) {
//                r1[0] = 0.0;
//                r1[1] = 0.0;
//            }
//        }
//        for (int i=0; i<r1.length; i++) {
//           System.out.println("R1 guess " + i + " " + r1[i]); 
//        }
        return r1;
    }

    public static double[][] r1rhoR2Guess(double[][] peaks, double[] yvals) {
        // Estimates CEST R2A and R2B values from peak widths for initial guesses for before fitting.
        // Uses the output from cestPeakGuess as the input.

        double pb = r1rhoPbGuess(peaks, yvals);

        if (peaks.length > 1) {
            double[][] r2 = new double[2][peaks.length / 2];
            double awidth;
            double bwidth;

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
                r2[0][0] = Math.abs(awidth - ka) / 12; //R2A
                r2[1][0] = Math.abs(bwidth - kb) / 6; //R2B
            } else {
                for (int i = 0; i < r2[0].length; i++) {
                    awidth = peaks[2 * i + 1][2] / (2 * Math.PI);
                    bwidth = peaks[2 * i][2] / (2 * Math.PI);
                    double kex = (awidth + bwidth) / 2;
                    double kb = pb * kex;
                    double ka = (1 - pb) * kex;
                    r2[0][i] = Math.abs(awidth - ka) / 12; //R2A
                    r2[1][i] = Math.abs(bwidth - kb) / 6; //R2B
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
                r2[0][0] = awidth / 3; //R2A
                r2[1][0] = awidth / 3; //R2B
            }
            return r2;
        }
    }

    public static double r1rhoKexGuess(double[][] peaks) {
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
            return kex[0] / 3;
        } else {
            return peaks[0][2] / (2 * Math.PI) / 3; //Kex;
        }

    }

    public static double[] r1Boundaries(double r1, double tex, double delta) {
        double baseline = Math.exp(-r1 * tex);
        double r1Low = -Math.log(baseline + 0.1) / tex;
        double r1Up = -Math.log(baseline - delta) / tex;
        double[] result = {r1Low, r1Up};
        return result;
    }

}
