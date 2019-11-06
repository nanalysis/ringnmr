package org.comdnmr.data;

import java.util.ArrayList;
import org.apache.commons.math3.optim.PointValuePair;
import static org.comdnmr.data.RelaxFit.ExptType.NOE;
import static org.comdnmr.data.RelaxFit.ExptType.R1;
import static org.comdnmr.data.RelaxFit.ExptType.R2;
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
    Fitter bestFitter;
    int[] residueModels;
    int[] resParStart;
    int[] nParsPerModel = {0, 1, 2, 0, 0, 3, 4};
    static ExptType[] expTypeArray = {R1, R2, NOE};
    double[] tauMBounds = {3.8e-9, 4.5e-9};
    
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

    public double[] getJ(double[] pars, RelaxEquations relaxObj, int modelNum) {
        double tauM = pars[0];//4.5e-9;
        double s2 = pars[1];
        double[] J = new double[5];
        switch (modelNum) {
            case 1:
                J = relaxObj.getJModelFree(tauM, s2);
                break;
            case 2:
                double tau = pars[2];
                J = relaxObj.getJModelFree(tau, tauM, s2);
                break;
            case 5:
                tau = pars[2];
                double sf2 = pars[3];
                J = relaxObj.getJModelFree(tau, tauM, s2, sf2);
                break;
            case 6:
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
//        int i = 0;
        for (int i=0; i<n; i++) {
            int iField = (int) xValues[i][0];
            RelaxEquations relaxObj = relaxObjs.get(iField);
            int iRes = (int) xValues[i][1];
            int iExpType = (int) xValues[i][2];
            int modelNum = residueModels[iRes];
            double[] J = getJ(bestPars, relaxObj, modelNum);
            ExptType type = expTypeArray[iExpType];
            double x = first + delta * i;
            double y = getYVal(bestPars, relaxObj, J, type);
            result[0][i] = x;
            result[1][i] = y;
        }
        return result;
    }

    public double value(double[] pars, double[][] values) {
        int n = values[0].length;
        double sum = 0.0;
        for (int i=0; i<n; i++) {
            int iField = (int) xValues[i][0];
            RelaxEquations relaxObj = relaxObjs.get(iField);
            int iRes = (int) xValues[i][1];
            int iExpType = (int) xValues[i][2];
            int modelNum = residueModels[iRes];
            double[] J = getJ(pars, relaxObj, modelNum);
            ExptType type = expTypeArray[iExpType];
            double y = getYVal(pars, relaxObj, J, type);
            double delta = y - values[2][i];
            sum += (delta * delta) / (errValues[i] * errValues[i]);
        }
        double chisq = Math.sqrt(sum / n);

        return chisq;

    }

    public double valueMultiResidue(double[] pars, double[][] values) {
        int n = values[0].length;
        double sum = 0.0;
        for (int i=0; i<values[2].length; i++) {
            int iField = (int) xValues[i][0];
            RelaxEquations relaxObj = relaxObjs.get(iField);
            int iRes = (int) xValues[i][1];
            int iExpType = (int) xValues[i][2];
            int modelNum = residueModels[iRes];
            double[] resPars = new double[nParsPerModel[modelNum] + 1];
            resPars[0] = pars[0];
            int parStart = resParStart[iRes] + 1;
            System.arraycopy(pars, parStart, resPars, 1, nParsPerModel[modelNum]);
//            parStart += resParStart[i];
            double[] J = getJ(resPars, relaxObj, modelNum);
            ExptType type = expTypeArray[iExpType];
            double y = getYVal(resPars, relaxObj, J, type);
            double delta = y - values[2][i];
            sum += (delta * delta) / (errValues[i] * errValues[i]);
        }

        double chisq = Math.sqrt(sum / n);

        return chisq;

    }

    public double[] calcYVals(double[] pars, boolean multiRes) {
        int n = yValues.length;
        double[] calcY = new double[n];
        for (int i=0; i<n; i++) {
            int iField = (int) xValues[i][0];
            RelaxEquations relaxObj = relaxObjs.get(iField);
            int iRes = (int) xValues[i][1];
            int iExpType = (int) xValues[i][2];
            int modelNum = residueModels[iRes];
            double[] J = getJ(pars, relaxObj, modelNum);
            ExptType type = expTypeArray[iExpType];
            double y = getYVal(pars, relaxObj, J, type);
            if (multiRes) {
                double[] resPars = new double[nParsPerModel[modelNum] + 1];
                resPars[0] = pars[0];
                int parStart = resParStart[iRes] + 1;
                System.arraycopy(pars, parStart, resPars, 1, nParsPerModel[modelNum]);
    //            parStart += resParStart[i];
                J = getJ(resPars, relaxObj, modelNum);
                y = getYVal(resPars, relaxObj, J, type);
            }
            calcY[i] = y;
        }

        return calcY;
    }

    public double[] getPars() {
        return bestPars;
    }

    public double getAIC() {
        return bestAIC;
    }

    public double getChiSq() {
        return bestChiSq;
    }

    public double[] getParErrs() {
        return parErrs;
    }
    
    public double[][] getXVals() {
        return xValues;
    }

    public double[] getYVals() {
        return yValues;
    }

    public double[] getErrVals() {
        return errValues;
    }

    public Fitter getFitter() {
        return bestFitter;
    }
    
    public int[] getResParStart() {
        return resParStart;
    }
    
    public int[] getResidueModels() {
        return residueModels;
    }
    
    public void setResidueModels(int[] bestModels) {
        residueModels = bestModels;
        resParStart = new int[residueModels.length];
        int nParsPrev = 0;
        resParStart[0] = 0;
        for (int i=1; i<residueModels.length; i++) {
            resParStart[i] = nParsPerModel[residueModels[i-1]] + nParsPrev;
            nParsPrev = resParStart[i];
//            System.out.println(i + " resModel " + residueModels[i] + " resParStart " + resParStart[i]);
        }
    }
    
    public void setTauMBounds(double[] bounds) {
        tauMBounds[0] = bounds[0];
        tauMBounds[1] = bounds[1];
    }

    public PointValuePair fit(double sf, double[] guesses, boolean globalFit) {
        Fitter fitter;
        this.B0 = sf * 2.0 * Math.PI / RelaxEquations.GAMMA_H;
        if (globalFit) {
            fitter = Fitter.getArrayFitter(this::valueMultiResidue);
        } else {
            fitter = Fitter.getArrayFitter(this::value);
        }
        double[] xVals0 = new double[yValues.length];
        double[][] xValues2 = {xVals0, xVals0};//{xValues[0], xValues[1]};
        fitter.setXYE(xValues2, yValues, errValues);
        double[] start = guesses; //{max0, r0, max1, 10.0};
        double[] lower = new double[guesses.length]; //{max0 / 2.0, r0 / 2.0, max1 / 2.0, 1.0};
        double[] upper = new double[guesses.length]; //{max0 * 2.0, r0 * 2.0, max1 * 2.0, 100.0};
        for (int i = 0; i < lower.length; i++) {
            lower[i] = guesses[i] / 10.0;
            upper[i] = guesses[i] * 10.0;
        }
        lower[0] = tauMBounds[0];
        upper[0] = tauMBounds[1];
        lower[1] = 0.0;
        upper[1] = 1.0;
        if (globalFit) {
            for (int i=0; i<resParStart.length; i++) {
                int modelNum = residueModels[i];
                int iS2 = resParStart[i]+1;
                lower[iS2] = 0.0;
                upper[iS2] = 1.0;
                if (modelNum >= 5) {
                    lower[iS2+2] = 0.0;
                    upper[iS2+2] = 1.0;
                }
            }
        } else {
            if (lower.length > 3) {
                lower[3] = 0.0;
                upper[3] = 1.0;
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
            int n = yValues.length;
            bestAIC = n * Math.log(bestChiSq) + 2 * k;
            bestFitter = fitter;
            return result;
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }

    }

    public void calcParErrs(PointValuePair result, Fitter fitter, double[] fitPars, boolean globalFit) {
        double[] yCalc;
        yCalc = calcYVals(fitPars, globalFit);
        parErrs = fitter.bootstrap(result.getPoint(), 300, true, yCalc);
    }
}
