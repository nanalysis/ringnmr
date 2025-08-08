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

import org.comdnmr.modelfree.RelaxEquations;
import org.comdnmr.util.CoMDPreferences;

/**
 *
 * @author Simon Hulse
 */
public enum SSR1RhoEquation implements EquationType {

    CSA("SSR1RhoCSA", 0, "tauc", "s2") {

        @Override
        // nu1 and nuR are both provided in kHz
        public double calculate(double[] par, int[] map, double[] x, int idNum) {
            // FIXME: During simulations, getting map = [1, 0]... not sure why
            // expected tauc = par[map[0]] and s2 = par[map[1]]
            double tauc = par[map[1]];
            double s2 = par[map[0]];
            double nu1kHz = x[0];
            double nuRkHz = x[1];
            double omega1 = 2.0e3 * Math.PI * nu1kHz;
            double omegaR = 2.0e3 * Math.PI * nuRkHz;

            // TODO: Check that this line is facilitated in the GUI, in order
            // to vary SIGMA.
            // N.B. Δσ = 3/2 δ
            RelaxEquations.setSigma("C", (3.0 / 2.0) * -36.77e-6);
            double b0 = 1.0e6 * CoMDPreferences.getRefField();
            RelaxEquations relaxEquations = new RelaxEquations(b0, "H", "C");
            double r1rho = relaxEquations.r1RhoCSA(omegaR, omega1, tauc, s2);
            return r1rho;
        }

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum) {
            // TODO: should this be X.length or X[0].length?
            int n = X[0].length;
            double[] yCalc = new double[n];
            for (int i = 0; i < n; i++) {
                double[] x = new double[] {X[0][i], X[1][i]};
                yCalc[i] = calculate(par, map, x, idNum);
            }
            return yCalc;
        }

        // TODO: integrate with `map`
        public double[][] boundaries(
            double[] guesses,
            double[][] xValues,
            double[] yValues,
            int[][] map,
            int[] idNums,
            int nID) {
            return new double[][]{{0.0, 0.0}, {1.0, 1.0}};
        }

        @Override
        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
            // TODO: Anything sensible possible here?
            return new double[2];
        }

        @Override
        public int getNGroupPars() {
            return nGroupPars;
        }

        @Override
        public double getRex(double[] pars, int[] map, double field) {
            return 0.0;
        }

        @Override
        public double getKex(double[] pars) {
            return pars[0];
        }

        @Override
        public double getKex(double[] pars, int id) {
            return pars[0];
        }

        // TODO: For Bruce
        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[][] map = makeMap(1);
            return map;
        }

        // TODO: For Bruce
        @Override
        public int[][] makeMap(int n) {
            return makeMap(n, 2);
        }

        // TODO: For Bruce
        @Override
        public int[][] makeMap(int n, int m) {
            int[][] map = new int[n][m];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    map[i][0] = m * i + j;
                }
            }
            return map;
        }
    };

    final String equationName;
    final int nGroupPars;
    String[] parNames;
    static final String[] equationNames = {"CSA"};

    SSR1RhoEquation(String equationName, int nGroupPars, String... parNames) {
        this.equationName = equationName;
        this.parNames = parNames;
        this.nGroupPars = nGroupPars;
    }

    public static String[] getAllEquationNames() {
        return equationNames;
    }

    @Override
    public String getName() {
        return equationName;
    }

    @Override
    public String[] getParNames() {
        return parNames;
    }

    @Override
    public int getNGroupPars() {
        return nGroupPars;
    }

    public double getMinX() {
        return 1.0;
    }

    public double getMaxX() {
        return 50.0;
    }
}
