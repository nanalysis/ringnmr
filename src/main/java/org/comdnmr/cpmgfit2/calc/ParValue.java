/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.cpmgfit2.calc;

/**
 *
 * @author Bruce Johnson
 */
public class ParValue implements ParValueInterface {

    String name;
    double value;
    double err;

    public ParValue(String parName) {
        this.name = parName;
    }

    public ParValue(String parName, double value, double err) {
        this.name = parName;
        this.value = value;
        this.err = err;
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public double getValue() {
        return value;
    }

    @Override
    public double getError() {
        return err;
    }

}
