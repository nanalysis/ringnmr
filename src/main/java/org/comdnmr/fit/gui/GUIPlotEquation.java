/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.gui;

import javafx.scene.paint.Color;
import org.comdnmr.fit.calc.PlotEquation;

/**
 *
 * @author brucejohnson
 */
public class GUIPlotEquation extends PlotEquation {

    private Color color = Color.BLACK;
    private double scaleValue = 1.0;

    public GUIPlotEquation(String name, double[] pars, double[] errs, double[] extras) {
        super(name, pars, errs, extras);
    }

    public GUIPlotEquation(PlotEquation eqn) {
        super(eqn.getName(), eqn.getPars(), eqn.getErrs(), eqn.getExtras());
    }

    /**
     * @return the color
     */
    public Color getColor() {
        return color;
    }

    /**
     * @param color the color to set
     */
    public void setColor(Color color) {
        this.color = color;
    }

    /**
     * @return the scaleValue
     */
    public double getScaleValue() {
        return scaleValue;
    }

    /**
     * @param scaleValue the scaleValue to set
     */
    public void setScaleValue(double scaleValue) {
        this.scaleValue = scaleValue;
    }
    

}
