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
package org.comdnmr.eqnfit;

import org.comdnmr.data.ResidueProperties;
import org.comdnmr.data.ResidueData;
import org.comdnmr.data.ExperimentData;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.inference.TTest;
import org.comdnmr.util.CoMDOptions;

/**
 *
 * @author Bruce Johnson
 */
public class CPMGFitter implements EquationFitter {

    public static final double[] SIMX = {25.0, 50.0, 100.0, 150.0, 250.0, 350.0, 500.0, 750.0, 1000.0};

    CPMGFitFunction calcR;
    CoMDOptions options;
    static List<String> equationNameList = Arrays.asList(CPMGEquation.getEquationNames());
    List<Double> xValues = new ArrayList<>();
    List<Double> yValues = new ArrayList<>();
    List<Double> errValues = new ArrayList<>();
    List<Double> fieldValues = new ArrayList<>();
    List<Integer> idValues = new ArrayList<>();
    int nCurves = 1;
    int nResidues = 1;
    int[][] states;
    int[] stateCount;
    String[] resNums;
    long errTime;
    static final String expType = "cpmg";

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

    public CPMGFitter(CoMDOptions options) {
        calcR = new CPMGFitFunction(options);
        this.options = options;
    }

    public static int getMapIndex(int[] state, int[] stateCount, int... mask) {
        int index = 0;
//        System.out.println(state.length + " mask " + mask.length);
//        for (int i = 0; i < state.length; i++) {
//            System.out.print(" " + state[i]);
//        }
//        System.out.println("");
        double mult = 1.0;
        for (int i = 0; i < mask.length; i++) {
//            System.out.println("mask:" + mask[i] + " state[mask]:" + state[mask[i]] + " count:" + stateCount[mask[i]]);
            index += mult * state[mask[i]];
            mult *= stateCount[mask[i]];
        }
        return index;
    }

    @Override
    public void setData(List<Double>[] allXValues, List<Double> yValues, List<Double> errValues, List<Double> fieldValues) {
        xValues.clear();
        xValues.addAll(allXValues[0]);
        this.yValues.clear();
        this.yValues.addAll(yValues);
        this.errValues.clear();
        this.errValues.addAll(errValues);
        this.fieldValues.clear();
        this.idValues.clear();
        this.fieldValues.addAll(fieldValues);
        for (Double yValue : yValues) {
            this.idValues.add(0);
        }
        resNums = new String[1];
        resNums[0] = "0";
        nCurves = 1;
        stateCount = new int[4];
        stateCount[0] = nResidues;
        stateCount[1] = 1;
        stateCount[2] = 1;
        stateCount[3] = 1;
        states = new int[1][4];
    }

    @Override
    public void setData(ResidueProperties resProps, String[] resNums) {
        xValues.clear();
        this.resNums = resNums.clone();
        nResidues = resNums.length;

        stateCount = resProps.getStateCount(nResidues);
        Collection<ExperimentData> expDataList = resProps.getExperimentData();
        nCurves = resNums.length * expDataList.size();
        states = new int[nCurves][];
        int k = 0;
        int resIndex = 0;
        int id = 0;
        for (String resNum : resNums) {
            for (ExperimentData expData : expDataList) {
                ResidueData resData = expData.getResidueData(resNum);
                if (resData != null) {

                    states[k++] = resProps.getStateIndices(resIndex, expData);
                    //  need peakRefs
                    double field = expData.getNucleusField();
                    double[][] x = resData.getXValues();
                    double[] y = resData.getYValues();
                    double[] err = resData.getErrValues();

                    for (int i = 0; i < y.length; i++) {
                        xValues.add(x[0][i]);
                        yValues.add(y[i]);
                        errValues.add(err[i]);
                        fieldValues.add(field);
                        idValues.add(id);
                    }
                    id++;
                }
            }
            resIndex++;
        }
    }

    @Override
    public FitFunction getFitModel() {
        return calcR;
    }

    @Override
    public List<String> getEquationNameList() {
        return getEquationNames();
    }
    
    @Override
    public String getExpType() {
        return expType;
    }

    public static List<String> getEquationNames() {
        return equationNameList;
    }

    @Override
    public int[] getStateCount() {
        return stateCount;
    }

    @Override
    public int[][] getStates() {
        return states;
    }

    @Override
    public void setupFit(String eqn) {
        double[][] x = new double[1][yValues.size()];
        double[] y = new double[yValues.size()];
        double[] err = new double[yValues.size()];
        int[] idNums = new int[yValues.size()];
        double[] fields = new double[yValues.size()];
        for (int i = 0; i < x[0].length; i++) {
            x[0][i] = xValues.get(i);
            y[i] = yValues.get(i);
            err[i] = errValues.get(i);
            //System.out.println(x[0][i]+", "+x[0][i]+", "+x[0][i]+", "+x[0][i]);
            fields[i] = fieldValues.get(i);
            idNums[i] = idValues.get(i);
        }
        calcR.setEquation(eqn);

        calcR.setXY(x, y);
        calcR.setIds(idNums);
        calcR.setErr(err);
        calcR.setFieldValues(fields);
        calcR.setMap(stateCount, states);
    }

    @Override
    public List<ParValueInterface> guessPars(String eqn) {
        setupFit(eqn);
        double[] guesses = calcR.guess();
        String[] parNames = calcR.getParNames();
        int[][] map = calcR.getMap();
        List<ParValueInterface> parValues = new ArrayList<>();
        for (int i = 0; i < parNames.length; i++) {
            double guess = guesses[map[0][i]];
            ParValueInterface parValue = new ParValue(parNames[i], guess);
            parValues.add(parValue);
        }
        return parValues;
    }

    @Override
    public double rms(double[] pars) {
        double rms = calcR.getRMS(pars);
        return rms;
    }

    @Override
    public FitResult doFit(String eqn, double[] sliderguesses, CoMDOptions options) {
        setupFit(eqn);
        int[][] map = calcR.getMap();
        double[] guesses;
        if (sliderguesses != null) {
            //fixme
            guesses = sliderguesses;
        } else {
            guesses = calcR.guess();
        }
//        System.out.println("dofit guesses = " + guesses);
        double[][] boundaries = calcR.boundaries(guesses);
        double sigma = options.getStartRadius();
        PointValuePair result = calcR.refine(guesses, boundaries[0],
                boundaries[1], sigma, options.getOptimizer());
        double[] pars = result.getPoint();
        System.out.print("Fit pars \n");
        for (int i = 0; i < pars.length; i++) {
            System.out.printf("%d %.3f %.3f %.3f %.3f\n", i, guesses[i], boundaries[0][i], pars[i], boundaries[1][i]);
        }
        System.out.println("");
        int nCurves = states.length;

        /*
        for (int i = 0; i < map.length; i++) {
            for (int j = 0; j < map[i].length; j++) {
                System.out.printf(" %3d", map[i][j]);
            }
            System.out.println("");
        }

      
        System.out.print("Fit " + x.length + " points");
        System.out.print("Fit pars");
        for (int i = 0; i < pars.length; i++) {
            System.out.printf(" %.3f", pars[i]);
        }
        System.out.println("");
         */
        double aic = calcR.getAICc(pars);
        double rms = calcR.getRMS(pars);
        double rChiSq = calcR.getReducedChiSq(pars);
        System.out.printf("%.3f %.3f %.3f\n", aic, rms, rChiSq);
//        System.out.println("rms " + rms);
        int nGroupPars = calcR.getNGroupPars();
        sigma /= 2.0;

        String[] parNames = calcR.getParNames();
        double[] errEstimates;
        double[][] simPars = null;
        double[] rexValues = calcR.getRex(pars);
        boolean okRex = false;
        double rexRatio = options.getRexRatio();
        for (double rexValue : rexValues) {
            if (rexValue > rexRatio * rms) {
                okRex = true;
                break;
            }
        }
        boolean exchangeValid = okRex;

        if (FitFunction.getCalcError()) {
            long startTime = System.currentTimeMillis();
            errEstimates = calcR.simBoundsStream(pars.clone(),
                    boundaries[0], boundaries[1], sigma, options);
            long endTime = System.currentTimeMillis();
            errTime = endTime - startTime;
            simPars = calcR.getSimPars();
            TTest tTest = new TTest();
            for (String parName : parNames) {
                int parIndex = -1;
                if (parName.equals("Kex")) {
                    parIndex = 0;
                    if (pars[parIndex] < errEstimates[parIndex]) {
                        exchangeValid = false;
                    }
                }
                if (parName.equals("dPPMmin")) {
                    parIndex = 2;
                }
                if (parIndex != -1) {
                    boolean valid = tTest.tTest(0.0, simPars[map[0][parIndex]], 0.02);
                    SummaryStatistics sStat = new SummaryStatistics();
                    for (double v : simPars[map[0][parIndex]]) {
                        sStat.addValue(v);
                    }
                    // System.out.println(sStat.toString());
                    double alpha = tTest.tTest(0.0, simPars[map[0][parIndex]]);
                    double mean = sStat.getMean();
                    double sdev = sStat.getStandardDeviation();
                    if (!valid) {
                        exchangeValid = false;
                    }
//                    System.out.println(parName + " " + parIndex + " " + valid + " " + alpha + " " + mean + " " + sdev);
                }
            }
        } else {
            errEstimates = new double[pars.length];
        }
        String refineOpt = options.getOptimizer();
        String bootstrapOpt = options.getBootStrapOptimizer();
        long fitTime = calcR.fitTime;
        long bootTime = errTime;
        int nSamples = options.getSampleSize();
        boolean useAbs = options.getAbsValueFit();
        boolean useNonParametric = options.getNonParametricBootstrap();
        double sRadius = options.getStartRadius();
        double fRadius = options.getFinalRadius();
        double tol = options.getTolerance();
        boolean useWeight = options.getWeightFit();
        CurveFit.CurveFitStats curveStats = new CurveFit.CurveFitStats(refineOpt, bootstrapOpt, fitTime, bootTime, nSamples, useAbs,
                useNonParametric, sRadius, fRadius, tol, useWeight);
        double[] usedFields = getFields(fieldValues, idValues);
        return getResults(this, eqn, parNames, resNums, map, states, usedFields, nGroupPars, pars, errEstimates, aic, rms, rChiSq, simPars, exchangeValid, curveStats);
    }

    @Override
    public double[] getSimX(int nPts, double xLB, double xUB) {
        int nPoints = nPts;
        double[] x = new double[nPoints];
        double firstValue = xLB;
        double lastValue = xUB;
        double delta = (lastValue - firstValue) / (nPoints + 1);
        double value = firstValue;
        for (int i = 0; i < nPoints; i++) {
            x[i] = value;
            value += delta;

        }
        return x;
    }

    @Override
    public double[] getSimXDefaults() {
        return SIMX;
    }
}
