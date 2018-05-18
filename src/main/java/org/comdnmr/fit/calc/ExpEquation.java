/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

import java.util.Arrays;

/**
 *
 * @author Bruce Johnson
 */
public enum ExpEquation implements EquationType {

    EXPAB("ExpAB", 0, "A", "R") {
        @Override
        public double calculate(double[] par, int[] map, double[] x, int idNum, double field) {
            double A = par[map[0]];
            double R = par[map[1]];
            double delay = x[0];
            double value = A * Math.exp(-R * delay);
//            System.out.println(delay + " " + A + " " + R + " " + value);
            return value;
        }

        @Override
        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            int nPars = CalcRDisp.getNPars(map);
            double[] guesses = new double[nPars];
            double kExSum = 0.0;
            for (int id = 0; id < map.length; id++) {
                double minY = DataUtil.getMinValue(yValues, idNums, id);
                double maxY = DataUtil.getMaxValue(yValues, idNums, id);
                double mean = DataUtil.getMeanValue(yValues, idNums, id);
                double vMid = DataUtil.getMidValue(yValues, xValues[0], idNums, id);
                System.out.println(minY + " " + maxY + " " + mean + " " + vMid);
                System.out.println(id + " " + map[id].length + " " + map[id][0] + " " + map[id][1]);
                guesses[map[id][0]] = maxY;
                guesses[map[id][1]] = 1.0 / vMid;

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
            return pars[map[2]];
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
            int[][] map = new int[n][2];
            for (int i = 0; i < n; i++) {
                map[i][0] = 2 * i + 0;
                map[i][1] = 2 * i + 1;
            }
            return map;
        }

        @Override
        public int[][] makeMap(int n, int m) {
            int[][] map = new int[n][3];
            for (int i = 0; i < n; i++) {
                map[i][0] = 2 * i + 0;
                map[i][1] = 2 * i + 1;
            }
            return map;
        }

        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int n = states.length;
            int[][] map = new int[n][3];
            int lastCount = 0;
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
            }
            lastCount++;
            int maxIndex = 0;
            for (int i = 0; i < n; i++) {
                map[i][1] = CPMGFit.getMapIndex(states[i], stateCount, r2Mask) + lastCount;
                maxIndex = Math.max(map[i][1], maxIndex);
            }
            lastCount = maxIndex + 1;
            for (int i = 0; i < n; i++) {
                map[i][2] = CPMGFit.getMapIndex(states[i], stateCount, 0, 3) + lastCount;
            }
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

    ExpEquation(String equationName, int nGroupPars, String... parNames) {
        this.equationName = equationName;
        this.parNames = parNames;
        this.nGroupPars = nGroupPars;
    }

    public static String[] getEquationNames() {
        String[] equationNames = {"EXPAB"};
        return equationNames;
    }

}
