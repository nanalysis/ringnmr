/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 *
 * @author Bruce Johnson
 */
public interface EquationFitter {

    public List<String> getEquationNameList();

    public FitModel getFitModel();

    public CPMGFitResult doFit(String eqn, boolean absMode, boolean nonParBootStrap, double[] sliderGuesses);

    public void setupFit(String eqn, boolean absMode);
    
    public List<ParValueInterface> guessPars(String eqn, boolean absMode);
    
    public double rms(double[] pars);

    public void setData(ResidueProperties resProps, String[] resNums);
    
    public int[] getStateCount();
    
    public int[][] getStates();

    public default CPMGFitResult getResults(EquationFitter fitter, String eqn, String[] parNames, String[] resNums, int[][] map, int[][] states, double[] usedFields, int nGroupPars, double[] pars, double[] errEstimates, double aic, double rms, double[][] simPars) {
        int nNonGroup = parNames.length - nGroupPars;
        List<CurveFit> curveFits = new ArrayList<>();
//        System.out.println("ning " + nCurves);
        int nCurves = states.length;
        for (int iCurve = 0; iCurve < nCurves; iCurve++) {
            String stateString = ResidueProperties.getStateString(states[iCurve]);
            double[] parArray = new double[parNames.length];
            double[] errArray = new double[parNames.length];
            List<ParValueInterface> parValues = new ArrayList<>();
            for (int i = 0; i < nGroupPars; i++) {
                ParValue parValue = new ParValue(resNums[states[iCurve][0]], stateString, parNames[i], pars[i], errEstimates[i]);
                parValues.add(parValue);
                parArray[i] = pars[i];
                errArray[i] = errEstimates[i];
            }
            for (int j = 0; j < nNonGroup; j++) {
                int k = map[iCurve][nGroupPars + j];
                ParValue parValue = new ParValue(resNums[states[iCurve][0]], stateString, parNames[nGroupPars + j], pars[k], errEstimates[k]);
                parValues.add(parValue);
                parArray[nGroupPars + j] = pars[k];
                errArray[nGroupPars + j] = errEstimates[k];
            }
//            System.out.println("res " + resNums[states[iCurve][0]] + " " + parValues.toString());

            HashMap<String, Double> parMap = new HashMap<>();
            for (ParValueInterface parValue : parValues) {
                parMap.put(parValue.getName(), parValue.getValue());
                parMap.put(parValue.getName() + ".sd", parValue.getError());
            }
            parMap.put("AIC", aic);
            parMap.put("RMS", rms);
            FitModel model = getFitModel();

            parMap.put("Equation", 1.0 + fitter.getEquationNameList().indexOf(eqn));
            // fixme
            double[] extras = new double[2];
            extras[0] = usedFields[0];
            extras[1] = 17.0 * 2 * Math.PI;
            //System.out.println("getResults got called with extras length = " + extras.length);
            PlotEquation plotEquation = new PlotEquation(eqn, parArray, errArray, extras);
            CurveFit curveFit = new CurveFit(stateString, resNums[states[iCurve][0]], parMap, plotEquation);
            curveFits.add(curveFit);
        }
        CPMGFitResult fitResult = new CPMGFitResult(parNames, curveFits, eqn, nGroupPars, aic, rms, simPars);
        return fitResult;
    }
}
