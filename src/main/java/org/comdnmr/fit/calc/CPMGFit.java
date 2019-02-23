package org.comdnmr.fit.calc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.inference.TTest;

/**
 *
 * @author Bruce Johnson
 */
public class CPMGFit implements EquationFitter {

    public static double REF_FIELD = 500.0;
    public static final double[] SIMX = {25.0, 50.0, 100.0, 150.0, 250.0, 350.0, 500.0, 750.0, 1000.0};

    CalcRDisp calcR = new CalcRDisp();
    static List<String> equationNameList = Arrays.asList(CPMGEquation.getEquationNames());
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

    public void setData(List<Double>[] allXValues, List<Double> yValues, List<Double> errValues) {
        xValues.clear();
        xValues.addAll(allXValues[0]);
        this.yValues.clear();
        this.yValues.addAll(yValues);
        this.errValues.clear();
        this.errValues.addAll(errValues);
        this.fieldValues.clear();
        this.idValues.clear();
        for (int i = 0; i < yValues.size(); i++) {
            this.fieldValues.add(500.0);
            this.idValues.add(0);
        }
        resNums = new String[1];
        resNums[0] = "0";
        usedFields = new double[1];
        usedFields[0] = 500.0;
        nCurves = 1;
        stateCount = new int[4];
        stateCount[0] = nResidues;
        stateCount[1] = 1;
        stateCount[2] = 1;
        stateCount[3] = 1;
        states = new int[1][4];
    }

    // public void setData(Collection<ExperimentData> expDataList, String[] resNums) {
    public void setData(ResidueProperties resProps, String[] resNums) {
        this.resNums = resNums.clone();
        nResidues = resNums.length;
        int id = 0;
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
        return calcR;
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
        calcR.setEquation(eqn);
        calcR.setAbsMode(absMode);

        calcR.setXY(x, y);
        calcR.setIds(idNums);
        calcR.setErr(err);
        calcR.setFieldValues(fields);
        calcR.setFields(usedFields);
        calcR.setMap(stateCount, states);
    }

    public List<ParValueInterface> guessPars(String eqn, boolean absMode) {
        setupFit(eqn, absMode);
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

    public double rms(double[] pars) {
        double rms = calcR.getRMS(pars);
        return rms;
    }

    public CPMGFitResult doFit(String eqn, boolean absMode, boolean nonParBootStrap, double[] sliderguesses) {
        setupFit(eqn, absMode);
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
        double sigma = FitModel.SIGMA_DEFAULT;
        PointValuePair result = calcR.refine(guesses, boundaries[0],
                boundaries[1], sigma, CoMDPreferences.getOptimizer());
        double[] pars = result.getPoint();
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
//        System.out.println("rms " + rms);
        int nGroupPars = calcR.getNGroupPars();
        sigma /= 2.0;

        String[] parNames = calcR.getParNames();
        double[] errEstimates;
        double[][] simPars = null;
        boolean exchangeValid = true;
        double[] rexValues = calcR.getRex(pars);
        boolean okRex = false;
        double rexRatio = CoMDPreferences.getRexRatio();
        for (double rexValue : rexValues) {
            if (rexValue > rexRatio * rms) {
                okRex = true;
                break;
            }
        }
        exchangeValid = okRex;

        if (FitModel.getCalcError()) {
            if (nonParBootStrap) {
                errEstimates = calcR.simBoundsBootstrapStream(pars.clone(), boundaries[0], boundaries[1], sigma);
            } else {
                errEstimates = calcR.simBoundsStream(pars.clone(), boundaries[0], boundaries[1], sigma);

            }
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
                if (parName.equals("Rex")) {
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
        return getResults(this, eqn, parNames, resNums, map, states, usedFields, nGroupPars, pars, errEstimates, aic, rms, rChiSq, simPars, exchangeValid);
    }

    public double[] getSimX() {
        return SIMX;
    }

}
