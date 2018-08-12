package org.comdnmr.fit.calc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import org.apache.commons.math3.optim.PointValuePair;

/**
 *
 * @author Bruce Johnson
 */
public class ExpFit implements EquationFitter {

    FitModel expModel = new CalcExpDecay();
    List<Double> xValues = new ArrayList<>();
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
    public void setData(ResidueProperties resProps, String[] resNums) {
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
                    xValues.add(x[0][i]);
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

    public FitModel getFitModel() {
        return expModel;
    }

    @Override
    public List<String> getEquationNameList() {
        return getEquationNames();
    }

    public static List<String> getEquationNames() {
        return equationNameList;
    }

    public int[] getStateCount() {
        return stateCount;
    }

    public int[][] getStates() {
        return states;
    }

    public void setupFit(String eqn, boolean absMode) {
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
        expModel.setEquation(eqn);
        expModel.setAbsMode(absMode);

        expModel.setXY(x, y);
        expModel.setIds(idNums);
        expModel.setErr(err);
        expModel.setFieldValues(fields);
        expModel.setFields(usedFields);
        expModel.setMap(stateCount, states);
    }

    public List<ParValueInterface> guessPars(String eqn, boolean absMode) {
        setupFit(eqn, absMode);
        double[] guesses = expModel.guess();
        //System.out.println("guesses length = " + guesses.length);
        String[] parNames = expModel.getParNames();
        //System.out.println("parNames length = " + parNames.length);
        List<ParValueInterface> parValues = new ArrayList<>();
        for (int i = 0; i < parNames.length; i++) {
            ParValueInterface parValue = new ParValue(parNames[i], guesses[i]);
            parValues.add(parValue);
        }
        return parValues;
    }

    public double rms(double[] pars) {
        double rms = expModel.getRMS(pars);
        return rms;
    }

    public CPMGFitResult doFit(String eqn, boolean absMode, boolean nonParBootStrap, double[] sliderguesses) {
        setupFit(eqn, absMode);

        int[][] map = expModel.getMap();
        double[] guesses;
        if (sliderguesses != null) {
            //fixme
            guesses = sliderguesses;
        } else {
            guesses = expModel.guess();
        }
//        System.out.println("dofit guesses = " + guesses);
        double[][] boundaries = expModel.boundaries(guesses);
        double[] sigma = new double[guesses.length];
        for (int i = 0; i < guesses.length; i++) {
            sigma[i] = (boundaries[1][i] - boundaries[0][i]) / 10.0;
//            System.out.println(i + " " + boundaries[0][i] + " " + boundaries[1][i] + " " + sigma[i]);
        }
        PointValuePair result = expModel.refine(guesses, boundaries[0], boundaries[1], sigma);
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
        double aic = expModel.getAICc(pars);
        double rms = expModel.getRMS(pars);
//        System.out.println("rms " + rms);
        int nGroupPars = expModel.getNGroupPars();
        for (int i = 0; i < guesses.length; i++) {
            sigma[i] /= 2.0;
        }

        String[] parNames = expModel.getParNames();
        double[] errEstimates;
        if (nonParBootStrap) {
            errEstimates = expModel.simBoundsBootstrapStream(pars.clone(), boundaries[0], boundaries[1], sigma);
        } else {
            errEstimates = expModel.simBoundsStream(pars.clone(), boundaries[0], boundaries[1], sigma);

        }
        double[][] simPars = expModel.getSimPars();
        return getResults(this, eqn, parNames, resNums, map, states, usedFields, nGroupPars, pars, errEstimates, aic, rms, simPars);
    }

}
