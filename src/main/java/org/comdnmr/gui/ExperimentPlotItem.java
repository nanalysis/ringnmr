/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.gui;

import org.comdnmr.data.ExperimentResult;

/**
 *
 * @author brucejohnson
 */
public class ExperimentPlotItem implements BarPlotItem {
    ExperimentResult expResult;

    @Override
    public int getResidueNumber() {
        return expResult.getResNum();
    }


    @Override
    public double getError() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int getAtom() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double getValue() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
