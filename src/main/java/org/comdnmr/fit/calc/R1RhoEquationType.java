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
public interface R1RhoEquationType extends EquationType {

    @Override
    public default double calculate(double[] par, int[] map, double[] X, int idNum, double field) {
        double[][] x = new double[3][1];
        //System.out.println(x.length + " " + x[0].length + " X " + X.length);
        x[0][0] = X[0];
        x[1][0] = X[1];
        x[2][0] = X[2];
        double[] fields = {field};
        double[] y = calculate(par, map, x, idNum, fields);
        return y[0];
    }

    @Override
    public default double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
        int nPars = CalcR1Rho.getNPars(map);
        double[] guesses = new double[nPars];
        for (int id = 0; id < map.length; id++) {
            int[] map1 = map[id];
            double[][] xy = CESTEquations.getXYValues(xValues, yValues, idNums, id);
            double[][] peaks = R1RhoEquations.r1rhoPeakGuess(xy[0], xy[1], field);
            double tex = xValues[2][0];
            double[] r1 = R1RhoEquations.r1rhoR1Guess(xy[1], tex);
            double[][] r2 = R1RhoEquations.r1rhoR2Guess(peaks, xy[1]);
            guesses[map1[0]] = R1RhoEquations.r1rhoKexGuess(peaks); //112.0; //kex
            guesses[map1[1]] = R1RhoEquations.r1rhoPbGuess(peaks, xy[1]); //0.1; //pb
            guesses[map1[2]] = peaks[peaks.length - 1][0]; //-250 * 2.0 * Math.PI; //deltaA
            guesses[map1[3]] = peaks[0][0]; //400 * 2.0 * Math.PI; //deltaB
            guesses[map1[4]] = r1[0]; //2.4; //R1A
            guesses[map1[5]] = r1[1]; //2.4; //R1B
            guesses[map1[6]] = r2[0][0]; //20.0; //R2A
            guesses[map1[7]] = r2[1][0]; //100.0; //R2B
        }
//            for (int i=0; i<guesses.length; i++) {
//                System.out.println(guesses[i]);
//            }

        return guesses;
    }

    @Override
    public default double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
        double[][] boundaries = new double[2][guesses.length];
        for (int id = 0; id < map.length; id++) {
            int[] map1 = map[id];
            double[][] xy = CESTEquations.getXYValues(xValues, yValues, idNums, id);
            double[][] peaks = R1RhoEquations.r1rhoPeakGuess(xy[0], xy[1], field);
            double dAbound = 0;
            double dBbound = 0;
            if (peaks.length > 1) {
                dAbound = (peaks[0][2] / field) / 2.0;
                dBbound = (peaks[1][2] / field) / 2.0;
            } else if (peaks.length == 1) {
                dAbound = (peaks[0][2] / field) / 2.0;
                dBbound = dAbound;
            }
            double tex = xValues[2][0];
            double r1A = guesses[map1[4]];
            double[] r1BouA = R1RhoEquations.r1Boundaries(r1A, tex, 0.1);
            double r1B = guesses[map1[5]];
            double[] r1BouB = R1RhoEquations.r1Boundaries(r1B, tex, 0.1);

            boundaries[0][map1[0]] = 1.0; //kex LB
            // boundaries[1][map1[0]] = guesses[map1[0]] * 5; //kex UB
            boundaries[1][map1[0]] = 500.0; //kex UB
            boundaries[0][map1[1]] = 0.01; //pb LB
            boundaries[1][map1[1]] = 0.25; //pb UB //guesses[1] * 4;
            boundaries[0][map1[2]] = guesses[map1[2]] - dAbound; //deltaA LB
            boundaries[1][map1[2]] = guesses[map1[2]] + dAbound; //deltaA UB
            boundaries[0][map1[3]] = guesses[map1[3]] - dBbound; //deltaB LB
            boundaries[1][map1[3]] = guesses[map1[3]] + dBbound; //deltaB UB
            boundaries[0][map1[4]] = r1BouA[0]; //R1A LB
            boundaries[1][map1[4]] = r1BouA[1]; //R1A UB
            boundaries[0][map1[5]] = r1BouB[0]; //R1B LB
            boundaries[1][map1[5]] = r1BouB[1]; //R1B UB
            boundaries[0][map1[6]] = 1.0; //R2A LB
            boundaries[1][map1[6]] = 250.0; //R2A UB
            boundaries[0][map1[7]] = 1.0; //R2B LB
            boundaries[1][map1[7]] = 250.0; //R2B UB
        }
        return boundaries;
    }

    @Override
    public default double getRex(double[] pars, int[] map, double field) {
        return 0.0;
    }

    @Override
    public default double getKex(double[] pars) {
        return pars[0];
    }

    @Override
    public default double getKex(double[] pars, int id) {
        return pars[0];
    }

    @Override
    public default int[][] makeMap(int n) {
        int nP = 8;
        int[][] map = new int[n][nP];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < nP; j++) {
                map[i][j] = nP * i + j;
            }
        }
        return map;
    }

    @Override
    public default int[][] makeMap(int n, int m) {
        int nP = m;
        int[][] map = new int[n][nP];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < nP; j++) {
                map[i][0] = nP * i + j;
            }
        }
        return map;
    }

    @Override
    public default int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
        int[][] map = makeMap(1);
        return map;
    }
}
