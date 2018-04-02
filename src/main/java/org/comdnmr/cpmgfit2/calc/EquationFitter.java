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
public interface EquationFitter {

    public CPMGFitResult doFit(ResidueProperties resProps, String eqn, boolean absMode, boolean nonParBootStrap);
    public void setData(ResidueProperties resProps, String[] resNums);

}
