/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

import java.util.Arrays;
import java.util.List;
import org.comdnmr.fit.calc.CESTEquations.Peak;

/**
 *
 * @author Bruce Johnson
 */
public enum R1RhoEquation implements R1RhoEquationType {
//     public double[] cestR1rhoPerturbation(double[][] X, double pb, double kex, double deltaA0, double deltaB0, double R10, double R2A, double R2B) {

    R1RHOPERTURBATION("R1rhoPerturbation", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double[] fields) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] omegarf = X[0];
            double[] b1Field = X[1];
            double[] deltaA = new double[omegarf.length];
            double[] deltaB = new double[omegarf.length];
            double[] omegaB1 = new double[omegarf.length];
            for (int i = 0; i < omegarf.length; i++) {
                deltaA[i] = (deltaA0 - omegarf[i]) * fields[i] * 2.0 * Math.PI;
                deltaB[i] = (deltaB0 - omegarf[i]) * fields[i] * 2.0 * Math.PI;
                omegaB1[i] = b1Field[i] * 2.0 * Math.PI;
            }
            double[] yCalc = R1RhoEquations.r1rhoPerturbation(omegaB1, pb, kex, deltaA, deltaB, R1A, R1B, R2A, R2B);
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
    R1RHOPERTURBATIONNOEX("R1rhoPerturbationNoEx", 0, "deltaA0", "R1A", "R2A") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double[] fields) {
            double deltaA0 = par[map[0]];
            double R1A = par[map[1]];
            double R2A = par[map[2]];
            double[] omegarf = X[0];
            double[] b1Field = X[1];
            double[] deltaA = new double[omegarf.length];
            double[] omegaB1 = new double[omegarf.length];

            for (int i = 0; i < omegarf.length; i++) {
                deltaA[i] = (deltaA0 - omegarf[i]) * fields[i] * 2.0 * Math.PI;
                omegaB1[i] = b1Field[i] * 2.0 * Math.PI;

            }
            double[] yCalc = R1RhoEquations.r1rhoPerturbationNoEx(omegaB1, deltaA, R1A, R2A);
            return yCalc;
        }

        @Override
        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            int nPars = CalcR1Rho.getNPars(map);
            double[] guesses = new double[nPars];
            for (int id = 0; id < map.length; id++) {
                int[] map1 = map[id];
                double[][] xy = CESTEquations.getXYValues(xValues, yValues, idNums, id);
                List<Peak> peaks = CESTEquations.cestPeakGuess(xy[0], xy[1], field, "r1rho");
                double tex = xValues[2][0];
                double[] r1 = CESTEquations.cestR1Guess(yValues, tex, "r1rho");
                double[][] r2 = CESTEquations.cestR2Guess(peaks, yValues, "r1rho");
                guesses[map1[0]] = peaks.get(peaks.size() - 1).position; //400 * 2.0 * Math.PI; //deltaA
                guesses[map1[1]] = r1[0]; //2.4; //R1A
                guesses[map1[2]] = 20.0; //20.0; //R2A
            }
            return guesses;

        }

        @Override
        public double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            double[][] boundaries = new double[2][guesses.length];
            int id = 0;
            int[] map1 = map[id];
            double[][] xy = CESTEquations.getXYValues(xValues, yValues, idNums, id);
            List<Peak> peaks = CESTEquations.cestPeakGuess(xy[0], xy[1], field, "r1rho");

            double dAbound = peaks.get(peaks.size() - 1).width / 2;
            double tex = xValues[2][0];
            double r1A = guesses[map[id][1]];
            double[] r1BouA = CESTEquations.r1Boundaries(r1A, tex, 0.1);

            boundaries[0][map[id][0]] = guesses[map[id][0]] - dAbound; //deltaA LB
            boundaries[1][map[id][0]] = guesses[map[id][0]] + dAbound; //deltaA UB
            boundaries[0][map[id][1]] = r1BouA[0]; //R1A LB
            boundaries[1][map[id][1]] = r1BouA[1]; //R1A UB
            boundaries[0][map[id][2]] = 0.1; //R2A LB
            boundaries[1][map[id][2]] = 250.0; //R2A UB
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
    R1RHOBALDWINKAY("R1rhoBaldwinKay", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double[] fields) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] omegarf = X[0];
            double[] b1Field = X[1];
            double[] deltaA = new double[omegarf.length];
            double[] deltaB = new double[omegarf.length];
            double[] omegaB1 = new double[omegarf.length];
            for (int i = 0; i < omegarf.length; i++) {
                deltaA[i] = (deltaA0 - omegarf[i]) * fields[i] * 2.0 * Math.PI;
                deltaB[i] = (deltaB0 - omegarf[i]) * fields[i] * 2.0 * Math.PI;
                omegaB1[i] = b1Field[i] * 2.0 * Math.PI;
            }
            double[] yCalc = R1RhoEquations.r1rhoBaldwinKay(omegaB1, pb, kex, deltaA, deltaB, R1A, R1B, R2A, R2B);
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
    R1RHOLAGUERRE("R1rhoLaguerre", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double[] fields) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] omegarf = X[0];
            double[] b1Field = X[1];
            double[] deltaA = new double[omegarf.length];
            double[] deltaB = new double[omegarf.length];
            double[] omegaB1 = new double[omegarf.length];
            for (int i = 0; i < omegarf.length; i++) {
                deltaA[i] = (deltaA0 - omegarf[i]) * fields[i] * 2.0 * Math.PI;
                deltaB[i] = (deltaB0 - omegarf[i]) * fields[i] * 2.0 * Math.PI;
                omegaB1[i] = b1Field[i] * 2.0 * Math.PI;
            }
            double[] yCalc = R1RhoEquations.r1rhoLaguerre(omegaB1, pb, kex, deltaA, deltaB, R1A, R1B, R2A, R2B);
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
    R1RHOEXACT("R1rhoExact", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double[] fields) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = new double[X[1].length];
            for (int i = 0; i < yCalc.length; i++) {
                double omegaB1 = X[1][i] * 2.0 * Math.PI;
                double deltaA = (deltaA0 - X[0][i]) * fields[i] * 2.0 * Math.PI;
                double deltaB = (deltaB0 - X[0][i]) * fields[i] * 2.0 * Math.PI;
                yCalc[i] = R1RhoEquations.r1rhoExact(omegaB1, pb, kex, deltaA, deltaB, R1A, R1B, R2A, R2B);
            }
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
    R1RHOEXACT0("R1rhoExact0", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double[] fields) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = new double[X[1].length];
            double delay = X[2][0];
            for (int i = 0; i < yCalc.length; i++) {
                double omegaB1 = X[1][i] * 2.0 * Math.PI;
                double deltaA = (deltaA0 - X[0][i]) * fields[i] * 2.0 * Math.PI;
                double deltaB = (deltaB0 - X[0][i]) * fields[i] * 2.0 * Math.PI;
                double r1rho = R1RhoEquations.r1rhoExact0(delay, omegaB1, pb, kex, deltaA, deltaB, R1A, R1B, R2A, R2B);
                yCalc[i] = r1rho;
            }
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

//        public double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
//            double[][] bounds = super.boundaries(guesses, xValues, yValues, map, idNums, nID, field);
//            for (int id = 0; id < map.length; id++) {
//                int[] map1 = map[id];
//                bounds[0][map1[4]] = 1.0; //R1A UB
//                bounds[1][map1[4]] = 2.5; //R1A UB
//            }
//            return bounds;
//        }
//
//        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
//            double[] guesses = super.guess(xValues, yValues, map, idNums, nID, field);
//            for (int id = 0; id < map.length; id++) {
//                int[] map1 = map[id];
//                guesses[map1[4]] = 2.0; //R1A UB
//            }
//            return guesses;
//        }
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

    R1RhoEquation(String equationName, int nGroupPars, String... parNames) {
        this.equationName = equationName;
        this.parNames = parNames;
        this.nGroupPars = nGroupPars;
    }

    public static String[] getAllEquationNames() {
        String[] equationNames = {"R1RHOPERTURBATIONNOEX", "R1RHOPERTURBATION", "R1RHOBALDWINKAY", "R1RHOLAGUERRE", "R1RHOEXACT", "R1RHOEXACT0"};
        return equationNames;
    }

    public static String[] getEquationNames() {
        String[] equationNames = new String[CoMDPreferences.getActiveR1RhoEquations().size()];
        for (int i = 0; i < equationNames.length; i++) {
            equationNames[i] = CoMDPreferences.getActiveR1RhoEquations().get(i);
        }
        return equationNames;
    }

}
