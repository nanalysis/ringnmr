/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.data;

import java.util.HashMap;

/**
 *
 * @author brucejohnson
 */
public class OffsetExperiment extends DoubleArrayExperiment {

    double[] freqOffsets;
    double tau;
    double B1field;

    public OffsetExperiment(ExperimentSet experimentSet, String name, String nucleus, double field,
            double temperature, String expMode, double tau, double B1field) {
        super(experimentSet, name, nucleus, field, temperature, expMode);
        this.tau = tau;
        this.B1field = B1field;
    }

    public void setXVals(double[] xVals) {
        this.freqOffsets = xVals.clone();
    }

    public double[] getXVals() {
        return freqOffsets;
    }

    public double getTau() {
        return tau;
    }

    public double getB1Field() {
        return B1field;
    }

}
