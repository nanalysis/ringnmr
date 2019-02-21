/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

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

    public void setData(List<Double>[] allXValues, List<Double> yValues, List<Double> errValues);

    public int[] getStateCount();

    public int[][] getStates();

    public double[] getSimX();

    public default CPMGFitResult getResults(EquationFitter fitter, String eqn, String[] parNames, String[] resNums, int[][] map, int[][] states, double[] usedFields, int nGroupPars, double[] pars, double[] errEstimates, double aic, double rms, double rChiSq, double[][] simPars, boolean hasExchange) {
        int nNonGroup = parNames.length - nGroupPars;
        List<CurveFit> curveFits = new ArrayList<>();
//        System.out.println("ning " + nCurves);
        Map<String, double[]> simsMap = new HashMap<>();
        for (int i = 0; i < map.length; i++) {
            for (int j = 0; j < map[i].length; j++) {
                String key = parNames[j] + " " + (map[i][j] - map[0][j]);
//                System.out.println(i + " " + j + " " + key + " " + map[i][j] + " " + pars[map[i][j]]);
//                DescriptiveStatistics dStat = new DescriptiveStatistics(simPars[map[i][j]]);
                simsMap.put(key, simPars[map[i][j]]);
//                System.out.println(dStat);
            }
        }
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
            parMap.put("rChiSq", rChiSq);
            FitModel model = getFitModel();

            parMap.put("Equation", 1.0 + fitter.getEquationNameList().indexOf(eqn));
            // fixme
            double[] extras = new double[1];
            extras[0] = usedFields[states[iCurve][1]];
            //System.out.println("getResults got called with extras length = " + extras.length);
            PlotEquation plotEquation = new PlotEquation(eqn, parArray, errArray, extras);
            CurveFit curveFit = new CurveFit(stateString, resNums[states[iCurve][0]], parMap, plotEquation);
            curveFits.add(curveFit);
        }
        CPMGFitResult fitResult = new CPMGFitResult(parNames, curveFits, eqn, nGroupPars, aic, rms, rChiSq, simsMap, hasExchange);
        return fitResult;
    }
}
