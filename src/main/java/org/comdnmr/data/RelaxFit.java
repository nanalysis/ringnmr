package org.comdnmr.data;

import java.util.ArrayList;
import java.util.stream.IntStream;
import org.apache.commons.math3.optim.PointValuePair;
import org.comdnmr.modelfree.RelaxEquations;

/**
 *
 * @author brucejohnson
 */
public class RelaxFit {

    boolean reportFitness = true;
    int reportAt = 10;
    long startTime = 0;
    double B0;
    double[][] xValues;
    double[] yValues;
    double[] errValues;
    double[] bestPars;
    double[] parErrs;
    double bestAIC;
    double bestChiSq;
    ArrayList<RelaxEquations> relaxObjs = new ArrayList<>();

    public void setXYE(double[][] xValues, double[] yValues, double[] errValues) {
        this.xValues = xValues;
        this.yValues = yValues;
        this.errValues = errValues;
    }

    public enum ExptType {
        R1, R2, NOE;
    }

    public void makeRelaxObjs(double[] fields, String elem1, String elem2) {
        for (int i = 0; i < fields.length; i++) {
            relaxObjs.add(new RelaxEquations(fields[i], elem1, elem2));
        }
    }

    public double[] getJ(double[] pars, RelaxEquations relaxObj) {
        double tauM = pars[0];//4.5e-9;
        double s2 = pars[1];
        double[] J = new double[5];
        switch (pars.length) {
            case 2:
                J = relaxObj.getJModelFree(tauM, s2);
                break;
            case 3:
                double tau = pars[2];
                J = relaxObj.getJModelFree(tau, tauM, s2);
                break;
            case 4:
                tau = pars[2];
                double sf2 = pars[3];
                J = relaxObj.getJModelFree(tau, tauM, s2, sf2);
                break;
            case 5:
                tau = pars[2];
                sf2 = pars[3];
                double tauS = pars[4];
//                System.out.println("tau, sf2, tauS = " + tau + " " + sf2 + " " + tauS);
                J = relaxObj.getJModelFree(tau, tauM, tauS, s2, sf2);
                break;
            default:
                break;
        }
//        for (double Jval : J) {
//            System.out.println("J: " + Jval);
//        }
        return J;
    }

    public double getYVal(double[] pars, RelaxEquations relaxObj, double[] J, ExptType type) {
        double y = 0.0;
        switch (type) {
            case R1:
                y = relaxObj.R1(J);
//                System.out.print("R1 = " + y + " ");
                break;
            case R2:
//                double Rex = pars[0];
                y = relaxObj.R2(J, 0.0);
//                System.out.print("R2 = " + y + " ");
                break;
            case NOE:
                y = relaxObj.NOE(J);
//                System.out.print("NOE = " + y + " ");
                break;
            default:
                break;
        }
        return y;
    }

    public double[][] getSimValues(double first, double last, int n, boolean adjust) {
        double[][] result = new double[2][n];
        double delta = (last - first) / (n - 1);
        int i = 0;
        for (RelaxEquations relaxObj : relaxObjs) {
            double[] J = getJ(bestPars, relaxObj);
            for (ExptType type : ExptType.values()) {
                double x = first + delta * i;
                double y = getYVal(bestPars, relaxObj, J, type);
                result[0][i] = x;
                result[1][i] = y;
                i++;
            }
        }
        return result;
    }

    public double value(double[] pars, double[][] values) {
        int n = values[0].length;
        double sum = 0.0;
        int i = 0;
//        for (int j=0;j<pars.length;j++) {
//            System.out.print(" " + pars[j]);
//        }
//        System.out.println("");
        for (RelaxEquations relaxObj : relaxObjs) {
            double[] J = getJ(pars, relaxObj);
            for (ExptType type : ExptType.values()) {
                double y = getYVal(pars, relaxObj, J, type);
                double delta = y - values[2][i];
                sum += (delta * delta) / (errValues[i] * errValues[i]);
//                System.out.print(" " + y + " " + values[2][i]);
                i++;
            }
//            System.out.println(" " + Math.sqrt(sum/n));
        }
        double chisq = Math.sqrt(sum / n);

        return chisq;

    }

    public double valueMultiResidue(double[] pars, double[][] values) {
        int n = values[0].length;
        double sum = 0.0;
        int i = 0;
//        for (int j=0;j<pars.length;j++) {
//            System.out.print(" " + pars[j]);
//        }
//        System.out.println("");
        int nParsPerResidue = 1;
        int nParsPerModel = 2;// model 1
        int nTypes = 3;  // r1 + r2 +  noe
        int nRes = n / nTypes;
        double[] resPars = new double[nParsPerModel];
        resPars[0] = pars[0];
        for (RelaxEquations relaxObj : relaxObjs) {
            for (int iRes = 0; iRes < nRes; iRes++) {
                System.arraycopy(pars, iRes * nParsPerResidue + 1, resPars, 1, nParsPerResidue);
//                for (double par : resPars) {
//                    System.out.print("resPars = " + par + " ");
//                }
//                System.out.println("");
                double[] J = getJ(resPars, relaxObj);
                for (ExptType type : ExptType.values()) {
                    double y = getYVal(resPars, relaxObj, J, type);
                    double delta = y - values[2][i];
                    sum += (delta * delta) / (errValues[i] * errValues[i]);
//                    System.out.print(" " + y + " " + values[2][i]);
                    i++;
                }
            }
//            System.out.println(" " + Math.sqrt(sum/n));
        }
        double chisq = Math.sqrt(sum / n);

        return chisq;

    }

    public double[] calcYVals(double[] pars) {
        double[] calcY = new double[3];
        for (RelaxEquations relaxObj : relaxObjs) {
            double[] J = getJ(pars, relaxObj);
            for (ExptType type : ExptType.values()) {
                double y = getYVal(pars, relaxObj, J, type);
                calcY[type.ordinal()] = y;
            }
        }
        return calcY;
    }

    public double[] getPars() {
        return bestPars;
    }

    public double getAIC() {
        return bestAIC;
    }

    public double[] getParErrs() {
        return parErrs;
    }

    public double[] getYVals() {
        return yValues;
    }

    public double[] getErrVals() {
        return errValues;
    }

    public PointValuePair fit(double sf, double[] guesses, boolean globalFit) {
        double max0 = IntStream.range(0, yValues.length).filter(i -> (i % 2) == 0).mapToDouble(i -> yValues[i]).max().getAsDouble();
        double max1 = IntStream.range(0, yValues.length).filter(i -> (i % 2) == 1).mapToDouble(i -> yValues[i]).max().getAsDouble();
        double halfMax = max0 / 2.0;
        double midDelta = Double.MAX_VALUE;
        double midX = 0.0;
        for (int i = 0; i < yValues.length; i += 2) {
            double y = yValues[i];
            double delta = Math.abs(y - halfMax);
            if (delta < midDelta) {
                midDelta = delta;
                midX = xValues[0][i];
            }
        }
        Fitter fitter;
        this.B0 = sf * 2.0 * Math.PI / RelaxEquations.GAMMA_H;
        if (globalFit) {
            fitter = Fitter.getArrayFitter(this::valueMultiResidue);
        } else {
            fitter = Fitter.getArrayFitter(this::value);
        }
        double[][] xValues2 = {xValues[0], xValues[1]};
        fitter.setXYE(xValues2, yValues, errValues);
        double[] start = guesses; //{max0, r0, max1, 10.0};
        double[] lower = new double[guesses.length]; //{max0 / 2.0, r0 / 2.0, max1 / 2.0, 1.0};
        double[] upper = new double[guesses.length]; //{max0 * 2.0, r0 * 2.0, max1 * 2.0, 100.0};
        for (int i = 0; i < lower.length; i++) {
            lower[i] = guesses[i] / 10.0;
            upper[i] = guesses[i] * 10.0;
        }
        lower[0] = 3.8e-9;
        upper[0] = 4.5e-9;
        lower[1] = 0.0;
        upper[1] = 1.0;
        for (int i = 1; i < lower.length; i += (lower.length - 1) / (yValues.length / 3)) {
            lower[i] = 0.0;
            upper[i] = 1.0;
        }
        if (lower.length > 3) {
            for (int i = 3; i < lower.length; i += (lower.length - 1) / (yValues.length / 3)) {
                lower[i] = 0.0;
                upper[i] = 1.0;
            }
        }
//        for (int i=0; i<guesses.length; i++) {
//            System.out.println("guess lower upper " + i + " " + guesses[i] + " " + lower[i] + " " + upper[i]);
//        }
        try {
            PointValuePair result = fitter.fit(start, lower, upper, 10.0);
            bestPars = result.getPoint();
            bestChiSq = result.getValue();
            int k = bestPars.length;
            bestAIC = bestChiSq + 2 * k;
            double[] yCalc = calcYVals(bestPars);
            parErrs = fitter.bootstrap(result.getPoint(), 300, false, yCalc);
            return result;
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }

    }
}
