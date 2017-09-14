package org.comdnmr.cpmgfit2.calc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.optim.PointValuePair;
import static org.comdnmr.cpmgfit2.calc.ResidueFitter.equationNames;

/**
 *
 * @author Bruce Johnson
 */
public class CPMGFit {

    List<Double> xValues = new ArrayList<>();
    List<Double> yValues = new ArrayList<>();
    List<Double> errValues = new ArrayList<>();
    List<Double> fieldValues = new ArrayList<>();
    List<Integer> idValues = new ArrayList<>();
    double[] usedFields = null;
    int nCurves = 1;
    int nResidues = 1;
    Map<Double, Integer> fieldMap = new LinkedHashMap();
    Map<Double, Integer> tempMap = new LinkedHashMap();
    Map<String, Integer> nucMap = new LinkedHashMap();
    int[][] states;
    int[] stateCount;
    String[] resNums;
    List<String> equationNameList = Arrays.asList(equationNames);

    class StateCount {

        int[][] states;
        int[] stateCount;

        StateCount(int[][] states, int[] stateCount) {
            this.states = states;
            this.stateCount = stateCount;
        }

        int getResIndex(int i) {
            return states[i][0];
        }

        int getTempIndex(int i) {
            return states[i][2];
        }
    }

    public static int getMapIndex(int[] state, int[] stateCount, int... mask) {
        int index = 0;
        System.out.println(state.length + " mask " + mask.length);
        for (int i = 0; i < state.length; i++) {
            System.out.print(" " + state[i]);
        }
        System.out.println("");
        double mult = 1.0;
        for (int i = 0; i < mask.length; i++) {
            System.out.println("mask:" + mask[i] + " state[mask]:" + state[mask[i]] + " count:" + stateCount[mask[i]]);
            index += mult * state[mask[i]];
            mult *= stateCount[mask[i]];
        }
        return index;
    }

    int[] getStateIndices(int resIndex, ExperimentData expData) {
        int[] state = new int[4];
        state[0] = resIndex;
        state[1] = fieldMap.get(Math.floor(expData.field));
        state[2] = tempMap.get(Math.floor(expData.temperature));
        state[3] = nucMap.get(expData.nucleus);
        System.out.println(resIndex + " " + expData.field + " " + expData.temperature + " " + expData.nucleus);
        System.out.println("state index residue:" + state[0] + " field:" + state[1] + " temp:" + state[2] + " nuc:" + state[3]);

        return state;
    }

    int[] getStateCount(int nResidues) {
        int[] state = new int[4];
        state[0] = nResidues;
        state[1] = fieldMap.size();
        state[2] = tempMap.size();
        state[3] = nucMap.size();
        System.out.println("state count residues:" + state[0] + " fields:" + state[1] + " temps:" + state[2] + " nucs:" + state[3]);
        return state;
    }

    public void setData(Collection<ExperimentData> expDataList, String[] resNums) {
        this.resNums = resNums.clone();
        nResidues = resNums.length;
        int id = 0;
        fieldMap.clear();
        tempMap.clear();
        nucMap.clear();
        for (ExperimentData expData : expDataList) {
            if (!fieldMap.containsKey(Math.floor(expData.field))) {
                fieldMap.put(Math.floor(expData.field), fieldMap.size());
            }
            if (!tempMap.containsKey(Math.floor(expData.temperature))) {
                tempMap.put(Math.floor(expData.temperature), tempMap.size());
            }
            if (!nucMap.containsKey(expData.nucleus)) {
                nucMap.put(expData.nucleus, nucMap.size());
            }
        }
        stateCount = getStateCount(resNums.length);
        nCurves = resNums.length * expDataList.size();
        states = new int[nCurves][];
        int k = 0;
        int resIndex = 0;
        for (String resNum : resNums) {
            for (ExperimentData expData : expDataList) {
                states[k++] = getStateIndices(resIndex, expData);
                ResidueData resData = expData.getResidueData(resNum);
                //  need peakRefs
                double field = expData.getField();
                double[] x = resData.getXValues();
                double[] y = resData.getYValues();
                double[] err = resData.getErrValues();
                for (int i = 0; i < x.length; i++) {
                    xValues.add(x[i]);
                    yValues.add(y[i]);
                    errValues.add(err[i]);
                    fieldValues.add(field);
                    idValues.add(id);
                }
                id++;

            }
            resIndex++;
        }
        usedFields = new double[expDataList.size()];
        int iExp = 0;
        for (ExperimentData expData : expDataList) {
            usedFields[iExp++] = expData.getField();
        }
    }

    public static String getStateString(int[] state) {
        StringBuilder builder = new StringBuilder();
        builder.append(state[0]);
        for (int i = 1; i < state.length; i++) {
            builder.append(':');
            builder.append(state[i]);
        }
        return builder.toString();
    }

    public CPMGFitResult doFit(String eqn) {
        double[] x = new double[xValues.size()];
        double[] y = new double[xValues.size()];
        double[] err = new double[xValues.size()];
        int[] idNums = new int[xValues.size()];
        double[] fields = new double[xValues.size()];
        for (int i = 0; i < x.length; i++) {
            x[i] = xValues.get(i);
            y[i] = yValues.get(i);
            err[i] = errValues.get(i);
            fields[i] = fieldValues.get(i);
            idNums[i] = idValues.get(i);
        }
        CalcRDisp calcR = new CalcRDisp();
        calcR.setEquation(eqn);

        calcR.setXY(x, y);
        calcR.setIds(idNums);
        calcR.setErr(err);
        calcR.setFieldValues(fields);
        calcR.setFields(usedFields);
        calcR.setMap(stateCount, states);
        int[][] map = calcR.getMap();
        double[] guesses = calcR.guess();
        double[][] boundaries = calcR.boundaries();
        double[] sigma = new double[guesses.length];
        for (int i = 0; i < guesses.length; i++) {
            sigma[i] = (boundaries[1][i] - boundaries[0][i]) / 10.0;
        }
        PointValuePair result = calcR.refine(guesses, boundaries[0], boundaries[1], sigma);

        double[] pars = result.getPoint();
        System.out.print("Fit pars ");
        for (int i = 0; i < pars.length; i++) {
            System.out.printf(" %.3f", pars[i]);
        }
        System.out.println("");
        double aic = calcR.getAICc(pars);
        double rms = calcR.getRMS(pars);
        System.out.println("rms " + rms);
        int nGroupPars = calcR.getNGroupPars();

        String[] parNames = calcR.getParNames();
        double[] errEstimates;
        if (false) {
            errEstimates = calcR.simBoundsBootstrapStream(pars.clone(), boundaries[0], boundaries[1], sigma);
        } else {
            errEstimates = calcR.simBoundsStream(pars.clone(), boundaries[0], boundaries[1], sigma);

        }
        int nNonGroup = parNames.length - nGroupPars;
        List<int[]> stateList = new ArrayList<>();
        List<String> residueNumbers = new ArrayList<>();
        List<CurveFit> curveFits = new ArrayList<>();
        System.out.println("ning " + nCurves);
        for (int iCurve = 0; iCurve < nCurves; iCurve++) {
            double[] parArray = new double[parNames.length];
            List<ParValueInterface> parValues = new ArrayList<>();
            for (int i = 0; i < nGroupPars; i++) {
                ParValue parValue = new ParValue(parNames[i], pars[i], errEstimates[i]);
                parValues.add(parValue);
                parArray[i] = pars[i];
            }
            for (int j = 0; j < nNonGroup; j++) {
                int k = map[iCurve][nGroupPars + j];
                ParValue parValue = new ParValue(parNames[nGroupPars + j], pars[k], errEstimates[k]);
                parValues.add(parValue);
                parArray[nGroupPars + j] = pars[k];
            }
            System.out.println("res " + resNums[states[iCurve][0]] + " " + parValues.toString());
            stateList.add(states[iCurve].clone());
            residueNumbers.add(resNums[states[iCurve][0]]);

            HashMap<String, Double> parMap = new HashMap<>();
            for (ParValueInterface parValue : parValues) {
                parMap.put(parValue.getName(), parValue.getValue());
                parMap.put(parValue.getName() + ".sd", parValue.getError());
            }
            parMap.put("AIC", aic);
            parMap.put("RMS", rms);
            parMap.put("Equation", 1.0 + equationNameList.indexOf(eqn));

            PlotEquation plotEquation = new PlotEquation(eqn, parArray, usedFields);
            CurveFit curveFit = new CurveFit(getStateString(states[iCurve]), resNums[states[iCurve][0]], parMap, plotEquation);
            curveFits.add(curveFit);
        }
        CPMGFitResult fitResult = new CPMGFitResult(parNames, curveFits, eqn, nGroupPars, aic, rms);
        return fitResult;
    }

}
