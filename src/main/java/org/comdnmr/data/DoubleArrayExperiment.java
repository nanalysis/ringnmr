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

    public DoubleArrayExperiment(String name, String nucleus, double field, double temperature, String expMode) {
        super(name, nucleus, field, temperature, expMode);
    }

    public void setXVals(double[] xVals) {
    }

    public double[] getXVals() {
        return null;
    }

}
