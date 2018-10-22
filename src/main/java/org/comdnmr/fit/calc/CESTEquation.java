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
public enum CESTEquation implements CESTEquationType {
//     public double[] cestR1rhoPerturbation(double[][] X, double pb, double kex, double deltaA0, double deltaB0, double R10, double R2A, double R2B) {

    CESTR1RHOPERTURBATION("cestR1rhoPerturbation", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double field) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.cestR1rhoApprox("trott", X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
            return yCalc;
        }

        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            int nPars = CalcCEST.getNPars(map);
            double[] guesses = new double[nPars];
            for (int id = 0; id < map.length; id++) {
                double[][] peaks = CESTEquations.cestPeakGuess(xValues, yValues);
                double tex = xValues[2][0];
                double[] r1 = CESTEquations.cestR1Guess(yValues, tex);
                double[][] r2 = CESTEquations.cestR2Guess(peaks, yValues);
                guesses[map[id][0]] = CESTEquations.cestKexGuess(peaks); //112.0; //kex
                guesses[map[id][1]] = CESTEquations.cestPbGuess(peaks, yValues); //0.1; //pb
                guesses[map[id][2]] = peaks[0][0]; //-250 * 2.0 * Math.PI; //deltaB
                guesses[map[id][3]] = peaks[peaks.length - 1][0]; //400 * 2.0 * Math.PI; //deltaA
                guesses[map[id][4]] = r1[0]; //2.4; //R1A
                guesses[map[id][5]] = r1[1]; //2.4; //R1B
                guesses[map[id][6]] = r2[0][0]; //20.0; //R2A
                guesses[map[id][7]] = r2[1][0]; //100.0; //R2B
            }
            return guesses;

        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[][] map = {{0, 1, 2, 3, 4, 4, 5, 6}};

            return map;
        }

    },
    CESTR1RHOPERTURBATIONNOEX("cestR1rhoPerturbationNoEx", 0, "deltaA0", "R1A", "R2A") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double field) {
            double deltaA0 = par[map[0]];
            double R1A = par[map[1]];
            double R2A = par[map[2]];
            double[] yCalc = CESTEquations.cestR1rhoPerturbationNoEx(X, deltaA0, R1A, R2A);
            return yCalc;
        }

        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            int nPars = CalcCEST.getNPars(map);
            double[] guesses = new double[nPars];
            for (int id = 0; id < map.length; id++) {
                double[][] peaks = CESTEquations.cestPeakGuess(xValues, yValues);
                double tex = xValues[2][0];
                double[] r1 = CESTEquations.cestR1Guess(yValues, tex);
                double[][] r2 = CESTEquations.cestR2Guess(peaks, yValues);
                guesses[map[id][0]] = peaks[0][0]; //400 * 2.0 * Math.PI; //deltaA
                guesses[map[id][1]] = r1[0]; //2.4; //R1A
                guesses[map[id][2]] = r2[0][0] / 2; //20.0; //R2A
            }
            return guesses;

        }

        @Override
        public double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            double[][] boundaries = new double[2][guesses.length];
            int id = 0;
            double[][] peaks = CESTEquations.cestPeakGuess(xValues, yValues);
            double dAbound = 0;
            dAbound = peaks[0][2] / 2;
            double tex = xValues[2][0];
            double r1A = guesses[map[id][1]];
            double[] r1BouA = CESTEquations.r1Boundaries(r1A, tex, 0.1);

            boundaries[0][map[id][0]] = guesses[map[id][0]] - dAbound; //deltaA LB
            boundaries[1][map[id][0]] = guesses[map[id][0]] + dAbound; //deltaA UB
            boundaries[0][map[id][1]] = r1BouA[0]; //R1A LB
            boundaries[1][map[id][1]] = r1BouA[1]; //R1A UB
            boundaries[0][map[id][2]] = 0.1; //R2A LB
            boundaries[1][map[id][2]] = guesses[map[id][2]] * 4; //R2A UB
            return boundaries;
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[][] map = {{0, 1, 2}};

            return map;
        }

    },
    CESTR1RHOSD("cestR1rhoSD", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double field) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.cestR1rhoApprox("sd", X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
            return yCalc;
        }

        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            int nPars = CalcCEST.getNPars(map);
            double[] guesses = new double[nPars];
            for (int id = 0; id < map.length; id++) {
                double[][] peaks = CESTEquations.cestPeakGuess(xValues, yValues);
                double tex = xValues[2][0];
                double[] r1 = CESTEquations.cestR1Guess(yValues, tex);
                double[][] r2 = CESTEquations.cestR2Guess(peaks, yValues);
                guesses[map[id][0]] = CESTEquations.cestKexGuess(peaks); //112.0; //kex
                guesses[map[id][1]] = CESTEquations.cestPbGuess(peaks, yValues); //0.1; //pb
                guesses[map[id][2]] = peaks[0][0]; //-250 * 2.0 * Math.PI; //deltaB
                guesses[map[id][3]] = peaks[peaks.length - 1][0]; //400 * 2.0 * Math.PI; //deltaA
                guesses[map[id][4]] = r1[0]; //2.4; //R1A
                guesses[map[id][5]] = r1[1]; //2.4; //R1B
                guesses[map[id][6]] = r2[0][0]; //r2[1][r2[0].length-1]/4; //20.0; //R2A
                guesses[map[id][7]] = r2[1][0]; //100.0; //R2B
            }

            return guesses;
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[][] map = {{0, 1, 2, 3, 4, 4, 5, 6}};

            return map;
        }
    },
    CESTR1RHOBALDWINKAY("cestR1rhoBaldwinKay", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double field) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.cestR1rhoApprox("baldwinkay", X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
            return yCalc;
        }

        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            int nPars = CalcCEST.getNPars(map);
            double[] guesses = new double[nPars];
            for (int id = 0; id < map.length; id++) {
                double[][] peaks = CESTEquations.cestPeakGuess(xValues, yValues);
                double tex = xValues[2][0];
                double[] r1 = CESTEquations.cestR1Guess(yValues, tex);
                double[][] r2 = CESTEquations.cestR2Guess(peaks, yValues);
                guesses[map[id][0]] = CESTEquations.cestKexGuess(peaks); //112.0; //kex
                guesses[map[id][1]] = CESTEquations.cestPbGuess(peaks, yValues); //0.1; //pb
                guesses[map[id][2]] = peaks[0][0]; //-250 * 2.0 * Math.PI; //deltaB
                guesses[map[id][3]] = peaks[peaks.length - 1][0]; //400 * 2.0 * Math.PI; //deltaA
                guesses[map[id][4]] = r1[0]; //2.4; //R1A
                guesses[map[id][5]] = r1[1]; //2.4; //R1B
                guesses[map[id][6]] = r2[0][0]; //r2[1][r2[0].length-1]/4; //20.0; //R2A
                guesses[map[id][7]] = r2[1][0]; //100.0; //R2B
            }

            return guesses;
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[][] map = {{0, 1, 2, 3, 4, 4, 5, 6}};

            return map;
        }
    },
    CESTR1RHON("cestR1rhoN", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double field) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.cestR1rhoApprox("laguerre", X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
            return yCalc;
        }

        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            int nPars = CalcCEST.getNPars(map);
            double[] guesses = new double[nPars];
            for (int id = 0; id < map.length; id++) {
                double[][] peaks = CESTEquations.cestPeakGuess(xValues, yValues);
                double tex = xValues[2][0];
                double[] r1 = CESTEquations.cestR1Guess(yValues, tex);
                double[][] r2 = CESTEquations.cestR2Guess(peaks, yValues);
                guesses[map[id][0]] = CESTEquations.cestKexGuess(peaks); //112.0; //kex
                guesses[map[id][1]] = CESTEquations.cestPbGuess(peaks, yValues); //0.1; //pb
                guesses[map[id][2]] = peaks[0][0]; //-250 * 2.0 * Math.PI; //deltaB
                guesses[map[id][3]] = peaks[peaks.length - 1][0]; //400 * 2.0 * Math.PI; //deltaA
                guesses[map[id][4]] = r1[0]; //2.4; //R1A
                guesses[map[id][5]] = r1[1]; //2.4; //R1B
                guesses[map[id][6]] = r2[0][0]; //r2[1][r2[0].length-1]/4; //20.0; //R2A
                guesses[map[id][7]] = r2[0][0]; //100.0; //R2B
            }
            return guesses;
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[][] map = {{0, 1, 2, 3, 4, 4, 5, 5}};

            return map;
        }
    },
    CESTR1RHOEXACT1("cestR1rhoExact1", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double field) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.cestR1rhoExact1(X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
            return yCalc;
        }

        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            int nPars = CalcCEST.getNPars(map);
            double[] guesses = new double[nPars];
            for (int id = 0; id < map.length; id++) {
                double minY = DataUtil.getMinValue(yValues, idNums, id);
                double maxY = DataUtil.getMaxValue(yValues, idNums, id);
                double mean = DataUtil.getMeanValue(yValues, idNums, id);
                double vMid = DataUtil.getMidValue(yValues, xValues[0], idNums, id);
                double[][] peaks = CESTEquations.cestPeakGuess(xValues, yValues);
                double tex = xValues[2][0];
                double[] r1 = CESTEquations.cestR1Guess(yValues, tex);
                double[][] r2 = CESTEquations.cestR2Guess(peaks, yValues);
                guesses[map[id][0]] = CESTEquations.cestKexGuess(peaks); //112.0; //kex
                guesses[map[id][1]] = CESTEquations.cestPbGuess(peaks, yValues); //0.1; //pb
                guesses[map[id][2]] = peaks[0][0]; //-250 * 2.0 * Math.PI; //deltaB
                guesses[map[id][3]] = peaks[peaks.length - 1][0]; //400 * 2.0 * Math.PI; //deltaA
                guesses[map[id][4]] = r1[0]; //2.4; //R1A
                guesses[map[id][5]] = r1[1]; //2.4; //R1B
                guesses[map[id][6]] = r2[0][0]; //r2[1][r2[0].length-1]/4; //20.0; //R2A
                guesses[map[id][7]] = r2[1][0]; //100.0; //R2B
            }
            return guesses;
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[][] map = {{0, 1, 2, 3, 4, 4, 5, 6}};

            return map;
        }

    },
    CESTEXACT0("cestExact0", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double field) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.cestExact0(X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
            return yCalc;
        }

        @Override
        public double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            double[][] boundaries = new double[2][guesses.length];
            int id = 0;
            double[][] peaks = CESTEquations.cestPeakGuess(xValues, yValues);
            double dAbound = 0;
            double dBbound = 0;
            if (peaks.length > 1) {
                dAbound = peaks[0][2] / 2;
                dBbound = peaks[1][2] / 2;
            } else if (peaks.length == 1) {
                dAbound = peaks[0][2] / 2;
                dBbound = dAbound;
            }
            boundaries[0][map[id][0]] = 1.0; //kex LB
            boundaries[1][map[id][0]] = guesses[map[id][0]] * 4; //kex UB
            boundaries[0][map[id][1]] = 0.01; //pb LB
            boundaries[1][map[id][1]] = 0.25; //pb UB //guesses[1] * 4;
            boundaries[0][map[id][2]] = guesses[map[id][2]] - dAbound; //deltaA LB
            boundaries[1][map[id][2]] = guesses[map[id][2]] + dAbound; //deltaA UB
            boundaries[0][map[id][3]] = guesses[map[id][3]] - dBbound; //deltaB LB
            boundaries[1][map[id][3]] = guesses[map[id][3]] + dBbound; //deltaB UB
            boundaries[0][map[id][4]] = 0.1; //R1A LB
            boundaries[1][map[id][4]] = guesses[map[id][4]] * 4; //R1A UB
            boundaries[0][map[id][5]] = 0.1; //R1B LB
            boundaries[1][map[id][5]] = guesses[map[id][5]] * 4; //R1B UB
            boundaries[0][map[id][6]] = 0.1; //R2A LB
            boundaries[1][map[id][6]] = guesses[map[id][6]] * 4; //R2A UB
            boundaries[0][map[id][7]] = 0.1; //R2B LB
            boundaries[1][map[id][7]] = guesses[map[id][7]] * 4; //R2B UB
            return boundaries;
        }

    },
    CESTEXACT1("cestExact1", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double field) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.cestExact0(X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
            return yCalc;
        }

        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            int nPars = CalcCEST.getNPars(map);
            double[] guesses = new double[nPars];
            for (int id = 0; id < map.length; id++) {
                double[][] peaks = CESTEquations.cestPeakGuess(xValues, yValues);
                double tex = xValues[2][0];
                double[] r1 = CESTEquations.cestR1Guess(yValues, tex);
                double[][] r2 = CESTEquations.cestR2Guess(peaks, yValues);
                guesses[map[id][0]] = CESTEquations.cestKexGuess(peaks); //112.0; //kex
                guesses[map[id][1]] = CESTEquations.cestPbGuess(peaks, yValues); //0.1; //pb
                guesses[map[id][2]] = peaks[0][0]; //-250 * 2.0 * Math.PI; //deltaB
                guesses[map[id][3]] = peaks[peaks.length - 1][0]; //400 * 2.0 * Math.PI; //deltaA
                guesses[map[id][4]] = r1[0]; //2.4; //R1A
                guesses[map[id][5]] = r1[1]; //2.4; //R1B
                guesses[map[id][6]] = r2[0][0]; //r2[1][r2[0].length-1]/4; //20.0; //R2A
                guesses[map[id][7]] = r2[1][0]; //100.0; //R2B
            }
            return guesses;
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[][] map = {{0, 1, 2, 3, 4, 4, 5, 6}};

            return map;
        }

    },
    CESTEXACT2("cestExact2", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double field) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.cestExact0(X, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
            return yCalc;
        }

        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            int nPars = CalcCEST.getNPars(map);
            double[] guesses = new double[nPars];
            for (int id = 0; id < map.length; id++) {
                double[][] peaks = CESTEquations.cestPeakGuess(xValues, yValues);
                double tex = xValues[2][0];
                double[] r1 = CESTEquations.cestR1Guess(yValues, tex);
                double[][] r2 = CESTEquations.cestR2Guess(peaks, yValues);
                guesses[map[id][0]] = CESTEquations.cestKexGuess(peaks); //112.0; //kex
                guesses[map[id][1]] = CESTEquations.cestPbGuess(peaks, yValues); //0.1; //pb
                guesses[map[id][2]] = peaks[0][0]; //-250 * 2.0 * Math.PI; //deltaB
                guesses[map[id][3]] = peaks[peaks.length - 1][0]; //400 * 2.0 * Math.PI; //deltaA
                guesses[map[id][4]] = r1[0]; //2.4; //R1A
                guesses[map[id][5]] = r1[1]; //2.4; //R1B
                guesses[map[id][6]] = r2[0][0]; //r2[1][r2[0].length-1]/4; //20.0; //R2A
                guesses[map[id][7]] = r2[1][0]; //100.0; //R2B
            }
            return guesses;
        }

        @Override
        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[][] map = {{0, 1, 2, 3, 4, 5, 6, 6}};

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
        String[] equationNames = {"CESTR1RHOPERTURBATIONNOEX", "CESTR1RHOPERTURBATION", "CESTR1RHOSD", "CESTR1RHOBALDWINKAY", "CESTR1RHON", 
            "CESTR1RHOEXACT1","CESTEXACT0", "CESTEXACT1", "CESTEXACT2"};
        return equationNames;
    }

}
