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
public interface EquationType {

    double calculate(double[] par, int[] map, double[] x, int idNum, double field);

    default double[] calculate(double[] par, int[] map, double[][] x, int idNum, double[] fields) {
        double[] yValues = new double[x[0].length];
        double[] ax = new double[x.length];
        for (int i = 0; i < x[0].length; i++) {
            for (int j = 0; j < ax.length; j++) {
                ax[j] = x[j][i];
            }
            yValues[i] = calculate(par, map, ax, idNum, fields[i]);
        }
        return yValues;
    }

    public default void constrain(String parName, double[] guesses, double[][] boundaries, int[][] map, int id, double lower, double upper) {
        String[] parNames = getParNames();
        int index = -1;
        for (int i = 0; i < parNames.length; i++) {
            if (parNames[i].equals(parName)) {
                index = i;
                break;
            }
        }
        if (index != -1) {
            int j = map[id][index];
            boundaries[0][j] = lower;
            boundaries[1][j] = upper;
            guesses[j] = (lower + upper) / 2;
        }
    }

    double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double[] fields);
    
    double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field);

    double getRex(double[] pars, int[] map, double field);

    double getKex(double[] pars);

    double getKex(double[] pars, int id);

    int[][] makeMap(int n);

    int[][] makeMap(int n, int m);

    int[][] makeMap(int[] stateCount, int[][] states, int[] mask);

    public String[] getParNames();

    public int getNGroupPars();

    public String getName();

    public default double getMinX() {
        return -10000.0;
    }

    public default double getMaxX() {
        return 10000.0;
    }

}
