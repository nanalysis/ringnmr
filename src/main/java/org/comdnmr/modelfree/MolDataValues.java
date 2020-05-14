/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.modelfree;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author brucejohnson
 */
public class MolDataValues {

    final String specifier;
    final double[] vector = new double[3];
    final List<RelaxDataValue> dataValues = new ArrayList<>();
    int testModel = 1;

    public MolDataValues(String specifier, double[] vector) {
        this.specifier = specifier;
        System.arraycopy(vector, 0, this.vector, 0, 3);
    }

    public void addData(RelaxDataValue value) {
        dataValues.add(value);
    }

    public List<RelaxDataValue> getData() {
        return dataValues;
    }

    public void setTestModel(int model) {
        this.testModel = model;
    }

    public int getTestModel() {
        return testModel;
    }

}
