package org.comdnmr.cpmgfit2.calc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import org.apache.commons.math3.optim.PointValuePair;

/**
 *
 * @author Bruce Johnson
 */
public class CESTFit implements EquationFitter {

    FitModel calcCEST = new CalcCEST();
    List<Double>[] xValues;
    List<Double> yValues = new ArrayList<>();
    List<Double> errValues = new ArrayList<>();
    List<Double> fieldValues = new ArrayList<>();
    List<Integer> idValues = new ArrayList<>();
    double[] usedFields = null;
    int nCurves = 1;
    int nResidues = 1;
    int[][] states;
    int[] stateCount;
    String[] resNums;
    static List<String> equationNameList = Arrays.asList(ExpEquation.getEquationNames());

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

    // public void setData(Collection<ExperimentData> expDataList, String[] resNums) {
    @Override
    public void setData(ResidueProperties resProps, String[] resNums) {
        xValues = new ArrayList[2];
        xValues[0] = new ArrayList<>();
        xValues[1] = new ArrayList<>();
        this.resNums = resNums.clone();
        nResidues = resNums.length;
        int id = 0;
        resProps.setupMaps();
        stateCount = resProps.getStateCount(resNums.length);
        Collection<ExperimentData> expDataList = resProps.getExperimentData();
        nCurves = resNums.length * expDataList.size();
        states = new int[nCurves][];
        int k = 0;
        int resIndex = 0;
        for (String resNum : resNums) {
            for (ExperimentData expData : expDataList) {
                states[k++] = resProps.getStateIndices(resIndex, expData);
                ResidueData resData = expData.getResidueData(resNum);
                //  need peakRefs
                double field = expData.getField();
                double[][] x = resData.getXValues();
                double[] y = resData.getYValues();
                double[] err = resData.getErrValues();
                for (int i = 0; i < y.length; i++) {
                    xValues[0].add(x[0][i]);
                    xValues[1].add(x[1][i]);
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

    public void setData(List<Double> xValues0, List<Double> xValues1, List<Double> yValues, List<Double> errValues) {
        xValues = new ArrayList[2];
        xValues[0] = new ArrayList<>();
        xValues[0].addAll(xValues0);
        xValues[1] = new ArrayList<>();
        xValues[1].addAll(xValues1);
        this.yValues.addAll(yValues);
        this.errValues.addAll(errValues);
        for (int i = 0; i < yValues.size(); i++) {
            fieldValues.add(1.0);
            idValues.add(0);
        }
        usedFields = new double[1];
        usedFields[0] = 1.0;
        resNums = new String[1];
        resNums[0] = "0";
        //states = new int[1][];
        //states[0] = new int[7];
        //stateCount = new int[7];
        //for (int i=0;i<states.length;i++) {
        //    states[0][i] = i;
        //}

        stateCount = new int[4];
        stateCount[0] = 1;
        stateCount[1] = 1;
        stateCount[2] = 1;
        stateCount[3] = 1;

        states = new int[1][4];
        states[0][0] = 0;
        states[0][1] = 0;
        states[0][2] = 0;
        states[0][3] = 0;

        // states
        // stateCount
        //System.out.println(xValues[0]);
        //System.out.println(xValues[1]);
        //System.out.println(this.yValues);
        //System.out.println(this.errValues);
        //System.out.print(xValues[0].size() + "\n");
        //System.out.print(xValues[1].size() + "\n");
        //System.out.print(this.yValues.size() + "\n");
    }

    public FitModel getFitModel() {
        return calcCEST;
    }

    public static List<String> getEquationNames() {
        return equationNameList;
    }

    @Override
    public CPMGFitResult doFit(String eqn, boolean absMode, boolean nonParBootStrap) {
        double[][] x = new double[2][yValues.size()];
        double[] y = new double[yValues.size()];
        double[] err = new double[yValues.size()];
        int[] idNums = new int[yValues.size()];
        double[] fields = new double[yValues.size()];
        for (int i = 0; i < x[0].length; i++) {
            x[0][i] = xValues[0].get(i);
            x[1][i] = xValues[1].get(i);
            y[i] = yValues.get(i);
            err[i] = errValues.get(i);
            fields[i] = fieldValues.get(i);
            idNums[i] = idValues.get(i);
        }
        calcCEST.setEquation(eqn);
        calcCEST.setAbsMode(absMode);

        calcCEST.setXY(x, y);
        calcCEST.setIds(idNums);
        calcCEST.setErr(err);
        calcCEST.setFieldValues(fields);
        calcCEST.setFields(usedFields);
        calcCEST.setMap(stateCount, states);
        int[][] map = calcCEST.getMap();
        double[] guesses = calcCEST.guess();
        double[][] boundaries = calcCEST.boundaries();
        double[] sigma = new double[guesses.length];
        for (int i = 0; i < guesses.length; i++) {
            sigma[i] = (boundaries[1][i] - boundaries[0][i]) / 10.0;
            //System.out.println(i + " map " + map[0][i] + " bou0 " + boundaries[0][i] + " bou1 " + boundaries[1][i] + " sig " + sigma[i] + " gue " + guesses[i]);
        }
        PointValuePair result = calcCEST.refine(guesses, boundaries[0], boundaries[1], sigma);
        double[] pars = result.getPoint();
        /*
        for (int i = 0; i < map.length; i++) {
            for (int j = 0; j < map[i].length; j++) {
                System.out.printf(" %3d", map[i][j]);
            }
            System.out.println("");
        }

        System.out.print("Fit pars ");
        for (int i = 0; i < pars.length; i++) {
            System.out.printf(" %.3f", pars[i]);
        }
        System.out.println("");
         */
        double aic = calcCEST.getAICc(pars);
        double rms = calcCEST.getRMS(pars);
//        System.out.println("rms " + rms);
        int nGroupPars = calcCEST.getNGroupPars();
        for (int i = 0; i < guesses.length; i++) {
            sigma[i] /= 2.0;
        }

        String[] parNames = calcCEST.getParNames();
        double[] errEstimates;
        if (nonParBootStrap) {
            errEstimates = calcCEST.simBoundsBootstrapStream(pars.clone(), boundaries[0], boundaries[1], sigma);
        } else {
            errEstimates = calcCEST.simBoundsStream(pars.clone(), boundaries[0], boundaries[1], sigma);

        }
        return getResults(eqn, parNames, resNums, map, states, usedFields, nGroupPars, pars, errEstimates, aic, rms);
    }

}
