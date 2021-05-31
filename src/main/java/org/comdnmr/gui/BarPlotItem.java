/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.gui;

/**
 *
 * @author brucejohnson
 */
public interface BarPlotItem {

    public int getResidueNumber();

    public int getAtom();

    public double getValue();

    public double getError();
}
