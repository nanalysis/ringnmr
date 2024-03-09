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
package org.comdnmr.eqnfit;

import org.comdnmr.util.CoMDPreferences;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author Bruce Johnson
 */
public enum CESTEquation implements CESTEquationType {

    TROTT_PALMER("trott_Palmer", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            return CESTEquations.r1rhoApprox("trott", X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
        }
       
        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[] map0 = {2, 3, 4, 4, 5, 6};
            int n = states.length;
            int[][] map = new int[n][8];
            int offset = 0;
            int lastState = 0;
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
                if (states[i][0] != lastState) {
                    offset += 5;
                }
                for (int j = 0; j < map0.length; j++) {
                    map[i][j + 2] = map0[j] + offset;
                }
                lastState = states[i][0];
            }
            return map;
        }
    },
    NOEX("noEx", 0, "deltaA0", "R1A", "R2A") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum) {
            double deltaA0 = par[map[0]];
            double R1A = par[map[1]];
            double R2A = par[map[2]];
            return CESTEquations.noEx(X, deltaA0, R1A, R2A);
        }

        @Override
        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
            int nPars = CESTFitFunction.getNPars(map);
            double[] guesses = new double[nPars];
            for (int id = 0; id < map.length; id++) {
                int[] map1 = map[id];
                double[][] xy = CESTEquations.getXYValues(xValues, yValues, idNums, id);
                List<CESTPeak> peaks = CESTEquations.cestPeakGuess(xy, "cest");
                if (peaks.size() > 0) {
                    double tex = xValues[2][0];
                    double[] r1 = CESTEquations.cestR1Guess(yValues, tex, "cest");
                    double[][] r2 = CESTEquations.cestR2Guess(peaks, yValues, "cest");
                    guesses[map1[0]] = peaks.get(0).position; //400 * 2.0 * Math.PI; //deltaA
                    guesses[map1[1]] = r1[0]; //2.4; //R1A
                    guesses[map1[2]] = r2[0][0] / 2; //20.0; //R2A
                } else {
                    guesses = null;
                }
            }
            return guesses;

        }

        @Override
        public double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID) {
            double[][] boundaries = new double[2][guesses.length];
            int id = 0;
            int[] map1 = map[id];
            double[][] xy = CESTEquations.getXYValues(xValues, yValues, idNums, id);
            List<CESTPeak> peaks = CESTEquations.cestPeakGuess(xy, "cest");

            double field = xValues[xValues.length -1][0];
            double dAbound = (peaks.get(0).getWidths()[1] / field) / 2;
            double tex = xValues[2][0];
            double r1A = guesses[map[id][1]];
            double[] r1BouA = CESTEquations.r1Boundaries(r1A, tex, 0.1);

            boundaries[0][map[id][0]] = guesses[map[id][0]] - dAbound; //deltaA LB
            boundaries[1][map[id][0]] = guesses[map[id][0]] + dAbound; //deltaA UB
            boundaries[0][map[id][1]] = r1BouA[0]; //R1A LB
            boundaries[1][map[id][1]] = r1BouA[1]; //R1A UB
            boundaries[0][map[id][2]] = 0.1; //R2A LB
            boundaries[1][map[id][2]] = guesses[map[id][2]] * 4; //R2A UB
            if (boundaries[1][map[id][2]] < 200.0) {
                boundaries[1][map[id][2]] = 200.0;
            }
            return boundaries;
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            //int[][] map = {{0, 1, 2}};
            int[] map0 = {0, 1, 2};
            int n = states.length;
            int[][] map = new int[n][3];
            for (int i = 0; i < n; i++) {
                map[i] = map0;
            }
            return map;
        }

    },
    SD("sD", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            return CESTEquations.r1rhoApprox("sd", X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[] map0 = {2, 3, 4, 4, 5, 6};
            int n = states.length;
            int[][] map = new int[n][8];
            int offset = 0;
            int lastState = 0;
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
                if (states[i][0] != lastState) {
                    offset += 5;
                }
                for (int j = 0; j < map0.length; j++) {
                    map[i][j + 2] = map0[j] + offset;
                }
                lastState = states[i][0];
            }
            return map;
        }
    },
    BALDWINKAY("baldwinKay", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.r1rhoApprox("baldwinkay", X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
            return yCalc;
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[] map0 = {2, 3, 4, 4, 5, 6};
            int n = states.length;
            int[][] map = new int[n][8];
            int offset = 0;
            int lastState = 0;
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
                if (states[i][0] != lastState) {
                    offset += 5;
                }
                for (int j = 0; j < map0.length; j++) {
                    map[i][j + 2] = map0[j] + offset;
                }
                lastState = states[i][0];
            }
            return map;
        }
    },
    LAGUERRE("laguerre", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.r1rhoApprox("laguerre", X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
            return yCalc;
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[] map0 = {2, 3, 4, 4, 5, 5};
            int n = states.length;
            int[][] map = new int[n][8];
            int offset = 0;
            int lastState = 0;
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
                if (states[i][0] != lastState) {
                    offset += 4;
                }
                for (int j = 0; j < map0.length; j++) {
                    map[i][j + 2] = map0[j] + offset;
                }
                lastState = states[i][0];
            }
            return map;
        }
    },
    EIGENEXACT1("eigenExact1", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.eigenExact1(X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
            return yCalc;
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[] map0 = {2, 3, 4, 4, 5, 6};
            int n = states.length;
            int[][] map = new int[n][8];
            int offset = 0;
            int lastState = 0;
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
                if (states[i][0] != lastState) {
                    offset += 5;
                }
                for (int j = 0; j < map0.length; j++) {
                    map[i][j + 2] = map0[j] + offset;
                }
                lastState = states[i][0];
            }
            return map;
        }

    },
    EXACT0("exact0", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.exact0(X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
            return yCalc;
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[] map0 = {2, 3, 4, 5, 6, 7};
            int n = states.length;
            int[][] map = new int[n][8];
            int offset = 0;
            int lastState = 0;
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
                if (states[i][0] != lastState) {
                    offset += 6;
                }
                for (int j = 0; j < map0.length; j++) {
                    map[i][j + 2] = map0[j] + offset;
                }
                lastState = states[i][0];
            }
            return map;

        }

    },
    EXACT1("exact1", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.exact0(X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
            return yCalc;
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[] map0 = {2, 3, 4, 4, 5, 6};
            int n = states.length;
            int[][] map = new int[n][8];
            int offset = 0;
            int lastState = 0;
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
                if (states[i][0] != lastState) {
                    offset += 5;
                }
                for (int j = 0; j < map0.length; j++) {
                    map[i][j + 2] = map0[j] + offset;
                }
                lastState = states[i][0];
            }
            return map;
        }

    },
    EXACT2("exact2", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.exact0(X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
            return yCalc;
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[] map0 = {2, 3, 4, 5, 6, 6};
            int n = states.length;
            int[][] map = new int[n][8];
            int offset = 0;
            int lastState = 0;
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
                if (states[i][0] != lastState) {
                    offset += 5;
                }
                for (int j = 0; j < map0.length; j++) {
                    map[i][j + 2] = map0[j] + offset;
                }
                lastState = states[i][0];
            }
            return map;

        }

    };
    final String equationName;
    final int nGroupPars;
    String[] parNames;
    double fieldRef;

    @Override
    public String getName() {
        return equationName;
    }

    @Override
    public String[] getParNames() {
        return parNames;
    }

    public void setFieldRef(double[] fields) {
        fieldRef = Arrays.stream(fields).min().getAsDouble();
    }

    public void setFieldRef(double field) {
        fieldRef = field;
    }

    @Override
    public int getNGroupPars() {
        return nGroupPars;
    }

    CESTEquation(String equationName, int nGroupPars, String... parNames) {
        this.equationName = equationName;
        this.parNames = parNames;
        this.nGroupPars = nGroupPars;
    }

    public static String[] getAllEquationNames() {
        String[] equationNames = {"NOEX", "TROTT_PALMER", "SD", "BALDWINKAY", "LAGUERRE",
            "EIGENEXACT1", "EXACT0", "EXACT1", "EXACT2"};
        return equationNames;
    }

    public static String[] getEquationNames() {
        String[] equationNames = new String[CoMDPreferences.getActiveCESTEquations().size()];
        for (int i = 0; i < equationNames.length; i++) {
            equationNames[i] = CoMDPreferences.getActiveCESTEquations().get(i);
        }
        return equationNames;
    }

}
