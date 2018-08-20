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
public interface CESTEquationType extends EquationType {
    
        public default double calculate(double[] par, int[] map, double[] X, int idNum, double field) {
            double[][] x = new double[3][1];
            //System.out.println(x.length + " " + x[0].length + " X " + X.length);
            x[0][0] = X[0];
            x[1][0] = X[1];
            x[2][0] = X[2];
            double[] y = calculate(par, map, x, idNum, field);
            return y[0];
        }

        public default double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            int nPars = CalcCEST.getNPars(map);
            double[] guesses = new double[nPars];
            for (int id = 0; id < map.length; id++) {
                double minY = DataUtil.getMinValue(yValues, idNums, id);
                double maxY = DataUtil.getMaxValue(yValues, idNums, id);
                double mean = DataUtil.getMeanValue(yValues, idNums, id);
                double vMid = DataUtil.getMidValue(yValues, xValues[0], idNums, id);
                //System.out.println(minY + " " + maxY + " " + mean + " " + vMid);
                //System.out.println(id + " " + map[id].length + " " + map[id][0] + " " + map[id][1]);
                double[][] peaks = CESTEquations.cestPeakGuess(xValues[0], yValues);
//                for (int i=0; i<peaks.length; i++) {
//                    for (int j=0; j<peaks[i].length; j++) {
//                        System.out.println("peaks " + i + " " + j + " = " + peaks[i][j]);
//                    }
//                }
                double tex = xValues[2][0];
                double[] r1 = CESTEquations.cestR1Guess(yValues, tex);
                double[][] r2 = CESTEquations.cestR2Guess(peaks);
                guesses[map[id][0]] = CESTEquations.cestKexGuess(peaks); //112.0; //kex
                guesses[map[id][1]] = CESTEquations.cestPbGuess(peaks); //0.1; //pb
                guesses[map[id][2]] = peaks[peaks.length-1][0]; //400 * 2.0 * Math.PI; //deltaA
                guesses[map[id][3]] = peaks[0][0]; //-250 * 2.0 * Math.PI; //deltaB
                guesses[map[id][4]] = r1[0]; //2.4; //R1A
                guesses[map[id][5]] = r1[1]; //2.4; //R1B
                guesses[map[id][6]] = r2[0][0]; //20.0; //R2A
                guesses[map[id][7]] = r2[1][0]; //100.0; //R2B
            }
//            for (int i=0; i<guesses.length; i++) {
//                System.out.println(guesses[i]);
//            }
            
            return guesses;
        }

        
        public default double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            double[][] boundaries = new double[2][guesses.length];
            int id = 0;
                boundaries[0][map[id][0]] = 1.0; //kex LB
                boundaries[1][map[id][0]] = guesses[map[id][0]] * 4; //kex UB
                boundaries[0][map[id][1]] = 0.01; //pb LB
                boundaries[1][map[id][1]] = 0.25; //pb UB //guesses[1] * 4;
                boundaries[0][map[id][2]] = 0.0; //deltaA LB
                boundaries[1][map[id][2]] = guesses[map[id][2]] * 4; //deltaA UB
                boundaries[0][map[id][3]] = guesses[map[id][3]] * 4; //deltaB LB
                boundaries[1][map[id][3]] = 0.0; //deltaB UB
                boundaries[0][map[id][4]] = 0.0; //R1A LB
                boundaries[1][map[id][4]] = guesses[map[id][4]] * 4; //R1A UB
                boundaries[0][map[id][5]] = 0.0; //R1B LB
                boundaries[1][map[id][5]] = guesses[map[id][5]] * 4; //R1B UB
                boundaries[0][map[id][6]] = 0.0; //R2A LB
                boundaries[1][map[id][6]] = guesses[map[id][6]] * 4; //R2A UB
                boundaries[0][map[id][7]] = 0.0; //R2B LB
                boundaries[1][map[id][7]] = guesses[map[id][7]] * 4; //R2B UB
            return boundaries;
        }

        
        public default double getRex(double[] pars, int[] map) {
            return 0.0;
        }

        
        public default double getKex(double[] pars) {
            return pars[0];
        }

        
        public default double getKex(double[] pars, int id) {
            return pars[0];
        }

        
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

        public default int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int[][] map = makeMap(1);
            return map;
        }
}
