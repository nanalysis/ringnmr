/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.cpmgfit2.calc;

import java.util.Arrays;

/**
 *
 * @author Bruce Johnson
 */
public enum CESTEquation implements EquationType {
//     public double[] cestR1rhoPerturbation(double[][] X, double pb, double kex, double deltaA0, double deltaB0, double R10, double R2A, double R2B) {

    CESTR1RHO("cestR1rhoPerturbation", 0, "kEx", "pB", "deltaA0", "deltaB0", "R10", "R2A", "R2B") {
        @Override
        public double calculate(double[] par, int[] map, double[] X, int idNum, double field) {
            double[][] x = new double[2][1];
            x[0][0] = X[0];
            x[1][0] = X[1];
            double[] y = calculate(par, map, x, idNum, field);
            return y[0];
        }

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double field) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R10 = par[map[4]];
            double R2A = par[map[5]];
            double R2B = par[map[6]];
            double[] yCalc = CESTEquations.cestR1rhoPerturbation(X, pb, kex, deltaA0, deltaB0, R10, R2A, R2B);
            return yCalc;
        }

        @Override
        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            int nPars = CalcCEST.getNPars(map);
            double[] guesses = new double[nPars];
            for (int id = 0; id < map.length; id++) {
                double minY = DataUtil.getMinValue(yValues, idNums, id);
                double maxY = DataUtil.getMaxValue(yValues, idNums, id);
                double mean = DataUtil.getMeanValue(yValues, idNums, id);
                double vMid = DataUtil.getMidValue(yValues, xValues[0], idNums, id);
                System.out.println(minY + " " + maxY + " " + mean + " " + vMid);
                System.out.println(id + " " + map[id].length + " " + map[id][0] + " " + map[id][1]);
                guesses[map[id][0]] = 0.0;
                guesses[map[id][1]] = 500.0;
                guesses[map[id][2]] = 0.0;
                guesses[map[id][3]] = 0.0;
                guesses[map[id][4]] = 0.0;
                guesses[map[id][5]] = 0.0;
                guesses[map[id][6]] = 0.0;
            }
            return guesses;
        }

        @Override
        public double[][] boundaries(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            double[] guesses = guess(xValues, yValues, map, idNums, nID, field);
            double[][] boundaries = new double[2][guesses.length];
            for (int id = 0; id < map.length; id++) {
                int iPar = map[id][0];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = guesses[iPar] * 4;
                iPar = map[id][1];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = guesses[iPar] * 4;
            }
            return boundaries;
        }

        @Override
        public double getRex(double[] pars, int[] map) {
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

        @Override
        public int[][] makeMap(int n) {
            int nP = 7;
            int[][] map = new int[n][nP];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < nP; j++) {
                    map[i][0] = nP * i + j;
                }
            }
            return map;
        }

        @Override
        public int[][] makeMap(int n, int m) {
            int nP = m;
            int[][] map = new int[n][nP];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < nP; j++) {
                    map[i][0] = nP * i + j;
                }
            }
            return map;
        }

        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[][] map = makeMap(1);
            return map;
        }
    };
    final String equationName;
    final int nGroupPars;
    String[] parNames;
    double fieldRef;

    public String getName() {
        return equationName;
    }

    public String[] getParNames() {
        return parNames;
    }

    public void setFieldRef(double[] fields) {
        fieldRef = Arrays.stream(fields).min().getAsDouble();
    }

    public void setFieldRef(double field) {
        fieldRef = field;
    }

    public int getNGroupPars() {
        return nGroupPars;
    }

    CESTEquation(String equationName, int nGroupPars, String... parNames) {
        this.equationName = equationName;
        this.parNames = parNames;
        this.nGroupPars = nGroupPars;
    }

    public static String[] getEquationNames() {
        String[] equationNames = {"CESTR1RHO"};
        return equationNames;
    }

}
