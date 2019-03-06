/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

/**
 *
 * @author Bruce Johnson
 */
public class PlotEquation {

    //        CalcRDisp.CPMGEquation equation = CalcRDisp.CPMGEquation.CPMGFAST;
    String name;
    double[] pars;
    double[] errs;
    double[] extras;

    public PlotEquation(String name, double[] pars, double[] errs, double[] extras) {
        this.name = name;
        this.pars = pars.clone();
        this.errs = errs.clone();
        //System.out.println("ploteq constructor " + extras.length);
        this.extras = extras.clone();
        //            equation = CalcRDisp.CPMGEquation.valueOf(name);
    }

    @Override
    public PlotEquation clone() {
        return new PlotEquation(name, pars, errs, extras);
    }

    public String getName() {
        return name;
    }

    public void setExtra(double[] extras) {
        this.extras = extras.clone();
    }

    public double getExtra(int i) {
        return extras[i];
    }

    public double[] getExtras() {
        return extras;
    }

    public double[] getPars() {
        return pars;
    }

    public double[] getErrs() {
        return pars;
    }

    public double calculate(double[] pars, double[] xValue, double field) {
        EquationType equationType = ResidueFitter.getEquationType(name);
        int[][] map = equationType.makeMap(1);
        return equationType.calculate(pars, map[0], xValue, 0, field);
    }

    public double calculate(double[] xValue, double field) {
        EquationType equationType = ResidueFitter.getEquationType(name);
        int[][] map = equationType.makeMap(1);

        double y = equationType.calculate(pars, map[0], xValue, 0, field);
        return y;
    }

    public double calculate(double xValue) {
        EquationType equationType = ResidueFitter.getEquationType(name);
        int[][] map = equationType.makeMap(1);
        double[] ax = new double[extras.length];
        ax[0] = xValue;
        for (int j = 1; j < extras.length; j++) {
            ax[j] = extras[j];
        }
        double y = calculate(ax, getExtra(0));
        return y;
    }

    public double getMinX() {
        EquationType equationType = ResidueFitter.getEquationType(name);
        return equationType.getMinX();
    }

    public double getMaxX() {
        EquationType equationType = ResidueFitter.getEquationType(name);
        return equationType.getMaxX();
    }

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append(name);
        for (double par : pars) {
            builder.append(" ").append(par);
        }
        for (double extra : extras) {
            builder.append(" ").append(extra);
        }
        return builder.toString();

    }

}
