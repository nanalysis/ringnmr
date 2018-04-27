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
public interface EquationType {

    double calculate(double[] par, int[] map, double[] x, int idNum, double field);

    default double[] calculate(double[] par, int[] map, double[][] x, int idNum, double field) {
        double[] yValues = new double[x[0].length];
        double[] ax = new double[x.length];
        for (int i = 0; i < x[0].length; i++) {
            for (int j = 0; j < ax.length; j++) {
                ax[j] = x[j][i];
            }
            yValues[i] = calculate(par, map, ax, idNum, field);
        }
        return yValues;
    }

    double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field);

    double[][] boundaries(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field);

    double getRex(double[] pars, int[] map);

    double getKex(double[] pars);

    double getKex(double[] pars, int id);

    int[][] makeMap(int n);

    int[][] makeMap(int n, int m);

    int[][] makeMap(int[] stateCount, int[][] states, int[] mask);

    public String[] getParNames();

    public int getNGroupPars();

    public void setFieldRef(double[] fields);

    public void setFieldRef(double field);

    public String getName();

}
