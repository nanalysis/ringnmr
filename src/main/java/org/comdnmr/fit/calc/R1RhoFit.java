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
public class R1RhoFit implements EquationFitter {

    FitModel calcR1Rho = new CalcR1Rho();
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
    static List<String> equationNameList = Arrays.asList(R1RhoEquation.getEquationNames());

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
        xValues = new ArrayList[3];
        xValues[0] = new ArrayList<>();
        xValues[1] = new ArrayList<>();
        xValues[2] = new ArrayList<>();
        this.resNums = resNums.clone();
        nResidues = resNums.length;
        int id = 0;
        resProps.setupMaps();
        stateCount = resProps.getStateCount(resNums.length);
        Collection<ExperimentData> expDataList = resProps.getExperimentData();
        nCurves = resNums.length * expDataList.size();

        // fixme??
        nCurves = resNums.length;
        states = new int[nCurves][];
        int k = 0;
        int resIndex = 0;
        List<Double> fieldList = new ArrayList<>();
        for (String resNum : resNums) {
            for (ExperimentData expData : expDataList) {
                states[k] = resProps.getStateIndices(resIndex, expData);
                ResidueData resData = expData.getResidueData(resNum);
                if (resData != null) {
                    //  need peakRefs
                    double field = expData.getField();
                    fieldList.add(field);
                    double[][] x = resData.getXValues();
//                System.out.println("setData x length = " + x.length);
                    double[] y = resData.getYValues();
                    double[] err = resData.getErrValues();
                    for (int i = 0; i < y.length; i++) {
                        xValues[0].add(x[0][i]);
                        xValues[1].add(x[1][i]);
                        xValues[2].add(x[2][i]);
                        yValues.add(y[i]);
                        errValues.add(err[i]);
                        fieldValues.add(field);
                        idValues.add(id);
                    }
                    // fixme ?? id++;
                }
            }
            k++;
            resIndex++;
        }
        usedFields = new double[fieldList.size()];
        int iExp = 0;
        for (Double field : fieldList) {
            usedFields[iExp++] = field;
        }
    }

    @Override
    public void setData(List<Double>[] allXValues, List<Double> yValues, List<Double> errValues) {
        setData(allXValues[0], allXValues[1], allXValues[2], yValues, errValues);
    }

    public void setData(List<Double> xValues0, List<Double> xValues1, List<Double> xValues2, List<Double> yValues, List<Double> errValues) {
        xValues = new ArrayList[3];
        xValues[0] = new ArrayList<>();
        xValues[0].addAll(xValues0);
        xValues[1] = new ArrayList<>();
        xValues[1].addAll(xValues1);
        xValues[2] = new ArrayList<>();
        xValues[2].addAll(xValues2);
        this.yValues.addAll(yValues);
        this.errValues.addAll(errValues);
        for (Double yValue : yValues) {
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

    @Override
    public FitModel getFitModel() {
        return calcR1Rho;
    }

    @Override
    public List<String> getEquationNameList() {
        return getEquationNames();
    }

    public static List<String> getEquationNames() {
        return equationNameList;
    }

    public static List<String> setEquationNames(List<String> eqnNames) {
        equationNameList = eqnNames;
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
    public void setupFit(String eqn, boolean absMode) {
        double[][] x = new double[3][yValues.size()];
        double[] y = new double[yValues.size()];
        double[] err = new double[yValues.size()];
        int[] idNums = new int[yValues.size()];
        double[] fields = new double[yValues.size()];
        for (int i = 0; i < x[0].length; i++) {
            x[0][i] = xValues[0].get(i);
            x[1][i] = xValues[1].get(i);
            x[2][i] = xValues[2].get(i);
            y[i] = yValues.get(i);
            err[i] = errValues.get(i);
            //System.out.println(x[0][i]+", "+x[0][i]+", "+x[0][i]+", "+x[0][i]);
            fields[i] = fieldValues.get(i);
            idNums[i] = idValues.get(i);
        }
        calcR1Rho.setEquation(eqn);
        calcR1Rho.setXY(x, y);
        calcR1Rho.setIds(idNums);
        calcR1Rho.setErr(err);
        calcR1Rho.setFieldValues(fields);
        calcR1Rho.setFields(usedFields);
        calcR1Rho.setMap(stateCount, states);
    }

    @Override
    public List<ParValueInterface> guessPars(String eqn, boolean absMode) {
        setupFit(eqn, absMode);
        double[] guesses = calcR1Rho.guess();
        String[] parNames = calcR1Rho.getParNames();
        List<ParValueInterface> parValues = new ArrayList<>();
        int[][] map = calcR1Rho.getMap();
        for (int i = 0; i < parNames.length; i++) {
            ParValueInterface parValue = new ParValue(parNames[i], guesses[map[0][i]]);
            parValues.add(parValue);
        }
        return parValues;
    }

    @Override
    public double rms(double[] pars) {
        double rms = calcR1Rho.getRMS(pars);
        return rms;
    }

    @Override
    public CPMGFitResult doFit(String eqn, boolean absMode, boolean nonParBootStrap, double[] sliderguesses) {
        double[][] xvals = new double[xValues.length][xValues[0].size()];
        double[] yvals = new double[yValues.size()];
        for (int i = 0; i < xvals.length; i++) {
            for (int j = 0; j < xvals[0].length; j++) {
                xvals[i][j] = xValues[i].get(j);
            }
        }
        for (int i = 0; i < yvals.length; i++) {
            yvals[i] = yValues.get(i);
        }
        double[][] peaks = R1RhoEquations.r1rhoPeakGuess(xvals, yvals, fieldValues.get(0));
        if (peaks.length >= 1) {
            setupFit(eqn, absMode);
            int[][] map = calcR1Rho.getMap();
            double[] guesses;
            if (sliderguesses != null) {
                //fixme
                guesses = sliderguesses;
            } else {
                guesses = calcR1Rho.guess();
            }
            //        System.out.println("dofit guesses = " + guesses);
            //        double[] guesses = setupFit(eqn, absMode);
            double[][] boundaries = calcR1Rho.boundaries(guesses);
            double sigma = CoMDPreferences.getStartingRadius();
            for (int i = 0; i < guesses.length; i++) {
                System.out.println(i + " map " + map[0][i] + " bou0 " + boundaries[0][i] + " bou1 " + boundaries[1][i] + " gue " + guesses[i]);
            }
            PointValuePair result = calcR1Rho.refine(guesses, boundaries[0],
                    boundaries[1], sigma, CoMDPreferences.getOptimizer());
            double[] pars = result.getPoint();
            System.out.println(eqn);

            for (int[] map1 : map) {
                for (int j = 0; j < map1.length; j++) {
                    System.out.printf(" %3d", map1[j]);
                }
                System.out.println("");
            }

            System.out.print("Fit pars ");
            for (int i = 0; i < pars.length; i++) {
                System.out.printf(" %.3f", pars[i]);
            }
            System.out.println("");

            double aic = calcR1Rho.getAICc(pars);
            double rms = calcR1Rho.getRMS(pars);
            double rChiSq = calcR1Rho.getReducedChiSq(pars);
            System.out.println("rms " + rms);
            int nGroupPars = calcR1Rho.getNGroupPars();
            sigma /= 2.0;

            String[] parNames = calcR1Rho.getParNames();
            double[] errEstimates;
            double[][] simPars = null;
            if (FitModel.getCalcError()) {
                if (nonParBootStrap) {
                    errEstimates = calcR1Rho.simBoundsBootstrapStream(pars.clone(), boundaries[0], boundaries[1], sigma);
                } else {
                    errEstimates = calcR1Rho.simBoundsStream(pars.clone(), boundaries[0], boundaries[1], sigma);

                }
                simPars = calcR1Rho.getSimPars();
            } else {
                errEstimates = new double[pars.length];
            }
            // fixme
            double[] extras = new double[xValues.length];
            extras[0] = usedFields[0];
            for (int j = 1; j < extras.length; j++) {
                extras[j] = xValues[j].get(0);
            }
            return getResults(this, eqn, parNames, resNums, map, states, extras, nGroupPars, pars, errEstimates, aic, rms, rChiSq, simPars, true);
        } else {
            return null;
        }
    }

    @Override
    public double[] getSimX() {
        int nPoints = 100;
        double[] x = new double[nPoints];
        double firstValue = -5000.0;
        double value = firstValue;
        for (int i = 0; i < nPoints; i++) {
            x[i] = value;
            value += 2.0 * Math.abs(firstValue) / (nPoints - 1);

        }
        return x;
    }

}
