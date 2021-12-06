/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.data;

/**
 *
 * @author brucejohnson
 */
public class DoubleArrayExperiment extends Experiment {

    public DoubleArrayExperiment(ExperimentSet experimentSet, String name, String nucleus, double field, double temperature, String expMode) {
        super(experimentSet, name, nucleus, field, temperature, expMode);
    }

    public void setXVals(double[] xVals) {
    }

    public double[] getXVals() {
        return null;
    }

}
