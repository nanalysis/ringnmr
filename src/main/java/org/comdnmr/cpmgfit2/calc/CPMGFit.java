package org.comdnmr.cpmgfit2.calc;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import org.apache.commons.math3.optim.PointValuePair;

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
    int nInGroup = 1;

    public void setData(Collection<ExperimentData> expDataList, String[] resNums) {
        nInGroup = resNums.length;
        int id = 0;
        for (String resNum : resNums) {
            for (ExperimentData expData : expDataList) {
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
            }
            id++;
        }
        usedFields = new double[expDataList.size()];
        int iExp = 0;
        for (ExperimentData expData : expDataList) {
            usedFields[iExp++] = expData.getField();
        }
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
        double[] guesses = calcR.guess();
        double[][] boundaries = calcR.boundaries();
        double[] sigma = new double[guesses.length];
        for (int i = 0; i < guesses.length; i++) {
            sigma[i] = (boundaries[1][i] - boundaries[0][i]) / 10.0;
        }
        PointValuePair result = calcR.refine(guesses, boundaries[0], boundaries[1], sigma);
        double[] pars = result.getPoint();
        double aic = calcR.getAICc(pars);
        double rms = calcR.getRMS(pars);
        int nGroupPars = calcR.getNGroupPars();

        String[] parNames = calcR.getParNames();
        double[] errEstimates;
        if (false) {
            errEstimates = calcR.simBoundsBootstrapStream(pars.clone(), boundaries[0], boundaries[1], sigma);
        } else {
            errEstimates = calcR.simBoundsStream(pars.clone(), boundaries[0], boundaries[1], sigma);

        }
        int nNonGroup = parNames.length - nGroupPars;
        List<List<ParValueInterface>> allParValues = new ArrayList<>();
        for (int iGroup = 0; iGroup < nInGroup; iGroup++) {
            List<ParValueInterface> parValues = new ArrayList<>();
            for (int i = 0; i < nGroupPars; i++) {
                ParValue parValue = new ParValue(parNames[i], pars[i], errEstimates[i]);
                parValues.add(parValue);
            }
            for (int j = 0; j < nNonGroup; j++) {
                int k = nGroupPars + iGroup * nNonGroup + j;
                ParValue parValue = new ParValue(parNames[nGroupPars + j], pars[k], errEstimates[k]);
                parValues.add(parValue);
            }
            allParValues.add(parValues);
        }
        CPMGFitResult fitResult = new CPMGFitResult(parNames, allParValues, eqn, nGroupPars, nInGroup, usedFields, aic, rms);
        return fitResult;
    }

}
