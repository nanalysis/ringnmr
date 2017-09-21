package org.comdnmr.cpmgfit2.calc;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Bruce Johnson
 */
public class CPMGFitResult {

    /**
     * @return the nGroupPars
     */
    public int getnGroupPars() {
        return nGroupPars;
    }

    /**
     * @return the nCurves
     */
    public int getNCurves() {
        return curveFits.size();
    }

    /**
     * @return the equationName
     */
    public String getEquationName() {
        return equationName;
    }

    /**
     * @return the fitParNames
     */
    public String[] getFitParNames() {
        return fitParNames;
    }

    /**
     * @return the aicc
     */
    public double getAicc() {
        return aicc;
    }

    /**
     * @return the rms
     */
    public double getRms() {
        return rms;
    }

    public CurveFit getCurveFit(int iCurve) {
        return curveFits.get(iCurve);
    }

    public double[] getPars(int iCurve) {
        return curveFits.get(iCurve).plotEquation.pars;
    }

    public double[] getErrs(int iCurve) {
        return curveFits.get(iCurve).plotEquation.errs;
    }

    public String getResidueNumber(int iCurve) {
        return curveFits.get(iCurve).resNum;
    }

    public boolean exchangeValid() {
        boolean valid = true;
        double[] pars = getPars(0);
        double[] errs = getErrs(0);
        switch (equationName) {
            case "NOEX":
                break;
            case "CPMGFAST":
            case "CPMGSLOW":
                if (errs[0] > pars[0]) {
                    valid = false;
                }
            default: {

            }
        }
        return valid;
    }

    private final String[] fitParNames;
    private final List<CurveFit> curveFits = new ArrayList<>();
    private final double aicc;
    private final double rms;
    private final String equationName;
    private final int nGroupPars;

    public CPMGFitResult(String[] fitParNames, List<CurveFit> curveFits, String equationName, int nGroupPars, double aicc, double rms) {
        this.curveFits.addAll(curveFits);
        this.fitParNames = fitParNames.clone();
        this.equationName = equationName;
        this.nGroupPars = nGroupPars;
        this.aicc = aicc;
        this.rms = rms;
    }
}
