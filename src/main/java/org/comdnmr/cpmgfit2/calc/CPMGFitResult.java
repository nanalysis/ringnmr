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
     * @return the nInGroup
     */
    public int getnInGroup() {
        return nInGroup;
    }

    /**
     * @return the usedFields
     */
    public double[] getUsedFields() {
        return usedFields;
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
     * @return the parValues
     */
    public List<ParValueInterface> getParValues(int iGroup) {
        return allParValues.get(iGroup);
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

    public double[] getPars(int iGroup) {
        return pars[iGroup].clone();
    }

    private final String[] fitParNames;
    private final List<List<ParValueInterface>> allParValues;
    private final double aicc;
    private final double rms;
    private final double[] usedFields;
    private final String equationName;
    private final int nGroupPars;
    private final int nInGroup;
    private final double[][] pars;

    CPMGFitResult(String[] fitParNames, List<List<ParValueInterface>> allParValues, String equationName, int nGroupPars, int nInGroup, double[] usedFields, double aicc, double rms) {
        this.fitParNames = fitParNames.clone();
        this.allParValues = new ArrayList<>();
        this.allParValues.addAll(allParValues);
        this.usedFields = usedFields.clone();
        this.aicc = aicc;
        this.rms = rms;
        this.equationName = equationName;
        this.nGroupPars = nGroupPars;
        this.nInGroup = nInGroup;
        pars = new double[allParValues.size()][fitParNames.length];
        for (int i = 0; i < pars.length; i++) {
            for (int j = 0; j < fitParNames.length; j++) {
                pars[i][j] = allParValues.get(i).get(j).getValue();
            }
        }

    }
}
