/*
 * CoMD/NMR Software : A Program for Analyzing NMR Dynamics Data
 * Copyright (C) 2018-2019 Bruce A Johnson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.eqnfit;

/**
 *
 * @author Bruce Johnson
 */
public interface EquationType {

    double calculate(double[] par, int[] map, double[] x, int idNum);

    default double[] calculate(double[] par, int[] map, double[][] x, int idNum) {
        double[] yValues = new double[x[0].length];
        double[] ax = new double[x.length];
        for (int i = 0; i < x[0].length; i++) {
            for (int j = 0; j < ax.length; j++) {
                ax[j] = x[j][i];
            }
            yValues[i] = calculate(par, map, ax, idNum);
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

    double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID);

    double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID);

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
