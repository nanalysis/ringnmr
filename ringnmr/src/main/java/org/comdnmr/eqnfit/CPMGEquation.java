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

import org.comdnmr.util.ANNLoader;
import org.comdnmr.util.CoMDPreferences;
import org.comdnmr.util.DataUtil;
import org.comdnmr.util.Utilities;

import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.util.FastMath;
import org.ojalgo.ann.ArtificialNeuralNetwork;
import org.ojalgo.matrix.store.MatrixStore;
import org.ojalgo.structure.Access1D;

/**
 *
 * @author Bruce Johnson
 */
public enum CPMGEquation implements EquationType {

    NOEX("noex", 0, "R2") {
        @Override
        public double calculate(double[] par, int[] map, double[] x, int idNum, double field) {
            double R2 = par[map[0]];
            double value = R2;
            return value;
        }

        @Override
        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double[] fields) {
            int nPars = CPMGFitFunction.getNPars(map);
            double[] guesses = new double[nPars];
            for (int id = 0; id < map.length; id++) {
                double mean = DataUtil.getMeanValue(yValues, idNums, id);
                guesses[map[id][0]] = mean;
            }
            return guesses;
        }

        @Override
        public double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            double[][] boundaries = new double[2][guesses.length];
            for (int id = 0; id < guesses.length; id++) {
                boundaries[0][id] = 0.0;
                boundaries[1][id] = guesses[id] * 4;
            }
            return boundaries;
        }

        @Override
        public double getKex(double[] pars) {
            return 0.0;
        }

        @Override
        public double getRex(double[] pars, int[] map, double field) {
            return 0.0;
        }

        @Override
        public double getKex(double[] pars, int id) {
            return 0.0;
        }

        @Override
        public int[][] makeMap(int n) {
            int[][] map = new int[n][1];
            for (int i = 0; i < n; i++) {
                map[i][0] = i;
            }
            return map;
        }

        @Override
        public int[][] makeMap(int n, int m) {
            int[][] map = new int[n][1];
            for (int i = 0; i < n; i++) {
                map[i][0] = i;
            }
            return map;
        }

        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int n = states.length;
            int[][] map = new int[n][1];
            for (int i = 0; i < n; i++) {
                map[i][0] = CPMGFitter.getMapIndex(states[i], stateCount, r2Mask);
            }
            return map;
        }
    }, CPMGFAST("cpmgfast", 1, "Kex", "R2", "dPPMmin") {
        @Override
        public double calculate(double[] par, int[] map, double[] x, int idNum, double field) {
            double kEx = par[map[0]];
            double R2 = par[map[1]];
            double dPPMmin = par[map[2]];
            double vu = x[0];
            double value;
            if (kEx <= 0.0) {
                value = R2;
            } else {
                double tauCP = 1.0 / (2.0 * vu);
                double dPPMMinRad = 2.0 * Math.PI * dPPMmin * field;
                double Rex = dPPMMinRad * dPPMMinRad / 4.0 / kEx;
                value = R2 + Rex * (1 - 2.0 * FastMath.tanh(0.5 * kEx * tauCP) / (kEx * tauCP));
            }
            return value;
        }

        @Override
        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double[] fields) {
            int nPars = CPMGFitFunction.getNPars(map);
            double[] guesses = new double[nPars];
            if (CoMDPreferences.getNeuralNetworkGuess()) {
                for (int id = 0; id < map.length; id++) {
                    double[] annGuess = new double[map[id].length];
                    double[][] xy = getXYValues(xValues, yValues, idNums, id);
                    double[] fields2 = {fields[id]};
                    try {
                        
                        annGuess = annCPMGGuesser("FAST", xy[0], xy[1], fields2);
                    } catch (Exception ex) {
                        Logger.getLogger(CPMGEquation.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    for (int k = 0; k < map[id].length; k++) {
                        guesses[map[id][k]] = annGuess[k];
                    }
                }
            } else {
                double kExSum = 0.0;
                for (int id = 0; id < map.length; id++) {
                    double minY = DataUtil.getMinValue(yValues, idNums, id);
                    double maxY = DataUtil.getMaxValue(yValues, idNums, id);
                    double mean = DataUtil.getMeanValue(yValues, idNums, id);
                    double vMid = DataUtil.getMidValue(yValues, xValues[0], idNums, id);
                    double r2 = minY * 0.95;
                    double rex = maxY - minY;
                    if (rex < 0.0) {
                        rex = 0.0;
                    }
                    double field = fields[id];
                    guesses[map[id][1]] = r2;
                    double tauMid = 1.0 / (2.0 * vMid);
                    double kEx = 1.915 / (0.5 * tauMid);
                    double dPPMMinRad = Math.sqrt(4.0 * rex / (field * field) * kEx);
                    double dPPMMin = dPPMMinRad / (2.0 * Math.PI);
                    guesses[map[id][2]] = dPPMMin;
                    if (rex >= 0) {
                        kExSum += kEx; // 1.915 comes from solving equation iteratively at tcp rex 0.5 half max
                    }
                }
                guesses[0] = kExSum / map.length;
                if (guesses[0] > CoMDPreferences.getCPMGMaxFreq()) {
                    guesses[0] = CoMDPreferences.getCPMGMaxFreq() * 0.9;
                }
            }

            return guesses;
        }

        @Override
        public double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            double[][] boundaries = new double[2][guesses.length];
            for (int id = 0; id < map.length; id++) {
                int iPar = map[id][0];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = Math.min(guesses[iPar] * 4, CoMDPreferences.getCPMGMaxFreq());
                iPar = map[id][1];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = guesses[iPar] * 4;
                iPar = map[id][2];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = guesses[iPar] * 4;
            }
            return boundaries;
        }

        @Override
        public double getRex(double[] pars, int[] map, double field) {
            double dPPMmin = pars[map[2]];

            double kEx = pars[0];
            double dPPMMinRad = 2.0 * Math.PI * dPPMmin * field;
            double Rex = dPPMMinRad * dPPMMinRad / 4.0 / kEx;
            return Rex;
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
            int[][] map = new int[n][3];
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 2 * i + 1;
                map[i][2] = 2 * i + 2;
            }
            return map;
        }

        @Override
        public int[][] makeMap(int n, int m) {
            int[][] map = new int[n][3];
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 2 * i + 1;
                map[i][2] = 2 * i + 2;
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
                map[i][1] = CPMGFitter.getMapIndex(states[i], stateCount, r2Mask) + lastCount;
                maxIndex = Math.max(map[i][1], maxIndex);
            }
            lastCount = maxIndex + 1;
            for (int i = 0; i < n; i++) {
                map[i][2] = CPMGFitter.getMapIndex(states[i], stateCount, 0, 3) + lastCount;
            }
            return map;
        }
    }, //        ISHIMA("isima", "R2", "Rex", "PaDw", "Tau") {
    //            double calculate(double[] par, double tcp, double field) {
    //                /*
    //                Ishima and Torchia approximation for skewed populations and all time scales
    //                R2(1/tcp)=R2+Rex/(1+Tau*sqrt(144/tcp**4+PaDw**4))
    //                 */
    //                double R2 = par[0];
    //                double Rex = par[1];
    //                double PaDw = par[2];
    //                double Tau = par[3];
    //                double value = R2 + Rex / (1 + Tau * FastMath.sqrt(144.0 / FastMath.pow(tcp, 4) + FastMath.pow(PaDw, 4)));
    //                return value;
    //            }
    //        },
    CPMGSLOW("cpmgslow", 2, "Kex", "pA", "R2", "dPPM") {
        @Override
        public double calculate(double[] par, int[] map, double[] x, int idNum, double field) {
            double kEx = par[map[0]];
            double pA = par[map[1]]; // p1-p2
            double r2 = par[map[2]];
            double dPPM = par[map[3]];
            double pB = 1.0 - pA;
            double pDelta = pA - pB;
            double dW = dPPM * field * 2.0 * Math.PI;
            double nu = x[0];
            double tauCP = 1.0 / (2.0 * nu);
            double psi = (pDelta * kEx) * (pDelta * kEx) - dW * dW + 4.0 * pA * pB * kEx * kEx;
            double zeta = -2.0 * dW * kEx * pDelta;
            double eta1 = Math.sqrt(psi * psi + zeta * zeta);
            double etaP = (1.0 / Math.sqrt(2.0)) * tauCP * Math.sqrt(eta1 + psi);
            double etaM = (1.0 / Math.sqrt(2.0)) * tauCP * Math.sqrt(eta1 - psi);
            double d1 = (psi + 2.0 * dW * dW) / Math.sqrt(psi * psi + zeta * zeta);
            double dP = 0.5 * (d1 + 1);
            double dM = 0.5 * (d1 - 1);
            double ch = dP * Math.cosh(etaP) - dM * Math.cos(etaM);
            if (Double.isInfinite(ch) || (ch > 1e40)) {
                ch = 1e40;
            }
            double rexContrib = 0.5 * (kEx - (1.0 / tauCP) * FastMath.acosh(ch));
            double dR = (1.0 / tauCP) * FastMath.acosh(ch);
            double value = r2 + rexContrib;
            if (Double.isInfinite(value)) {
                System.out.println("nu " + nu + "kex " + kEx + " pa " + pA + " r2 " + r2 + " dw " + dW + " psi " + psi + " rex " + rexContrib + " ch " + ch + " eta " + etaP + " " + etaM + " " + eta1 + " " + psi);
            }
            return value;
        }

        @Override
        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double[] fields) {
            int nPars = CPMGFitFunction.getNPars(map);
            double[] guesses = new double[nPars];
            if (CoMDPreferences.getNeuralNetworkGuess()) {
                for (int id = 0; id < map.length; id++) {
                    double[] annGuess = new double[map[id].length];
                    double[][] xy = getXYValues(xValues, yValues, idNums, id);
                    double[] fields2 = {fields[id]};
                    try {
                        
                        annGuess = annCPMGGuesser("SLOW", xy[0], xy[1], fields2);
                    } catch (Exception ex) {
                        Logger.getLogger(CPMGEquation.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    for (int k = 0; k < map[id].length; k++) {
                        guesses[map[id][k]] = annGuess[k];
                    }
                }
            } else {
                double kExSum = 0.0;
                double pa = 0.95;
                for (int id = 0; id < map.length; id++) {
                    double minY = DataUtil.getMinValue(yValues, idNums, id);
                    double maxY = DataUtil.getMaxValue(yValues, idNums, id);
                    double mean = DataUtil.getMeanValue(yValues, idNums, id);
                    double vMid = DataUtil.getMidValue(yValues, xValues[0], idNums, id);
                    double r2 = minY * 0.95;
                    double rex = maxY - r2;
                    double tauMid = 1.0 / (2.0 * vMid);
                    double kex = 1.915 / (0.5 * tauMid); // 1.915 comes from solving equation iteratively at tcp rex 0.5 half max
                    if (kex > CoMDPreferences.getCPMGMaxFreq()) {
                        kex = CoMDPreferences.getCPMGMaxFreq() * 0.9;
                    }
                    double field = fields[id];
                    double dw2 = rex / (pa * (1.0 - pa)) * kex;
                    double dPPM = Math.sqrt(dw2) / (2.0 * Math.PI) / field;
                    guesses[map[id][2]] = r2;
                    guesses[map[id][3]] = dPPM;
                    kExSum += kex;
                }
                guesses[0] = kExSum / map.length;
                guesses[1] = pa;
            }

            return guesses;
        }

        @Override
        public double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            double[][] boundaries = new double[2][guesses.length];
           // "Kex", "pA", "R2", "dPPM"
            for (int id = 0; id < map.length; id++) {
                int iPar = map[id][0];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = Math.min(guesses[iPar] * 4, CoMDPreferences.getCPMGMaxFreq());
                iPar = map[id][1];
                boundaries[0][iPar] = 0.5;
                boundaries[1][iPar] = 0.999;
                iPar = map[id][2];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = guesses[iPar] * 4;
                iPar = map[id][3];
                boundaries[0][iPar] = 0.0;
                boundaries[1][iPar] = guesses[iPar] * 4;
            }
            return boundaries;
        }
        //        CPMGSLOW("cpmgslow", 2, "Kex", "pA", "R2", "dW") {

        @Override
        public double getRex(double[] pars, int[] map, double field) {
            double[] x = new double[1];
            x[0] = 10.0;
            double y0 = calculate(pars, map, x, 0, field);
            x[0] = 1.0e4;
            double y1 = calculate(pars, map, x, 0, field);
            double rex = y0 - y1;
//            if (pars[map[3]] != 0.0) {
//                rex = pars[map[1]] * (1.0 - pars[map[1]]) * pars[map[0]] / (1.0 + Math.pow(pars[map[0]] / pars[map[3]], 2));
//            }
            return rex;
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
            int[][] map = new int[n][4];
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
                map[i][2] = 2 * i + 2;
                map[i][3] = 2 * i + 3;
            }
            return map;
        }

        @Override
        public int[][] makeMap(int n, int m) {
            int[][] map = new int[n][4];
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
                map[i][2] = 2 * i + 2;
                map[i][3] = 2 * i + 3;
            }
            return map;
        }

        public int[][] makeMap(int[] stateCount, int[][] states, int[] r2Mask) {
            int n = states.length;
            int[][] map = new int[n][4];
            int lastCount = 0;
            for (int i = 0; i < n; i++) {
                map[i][0] = 0;
                map[i][1] = 1;
            }
            int maxIndex = 0;
            lastCount = 2;
            for (int i = 0; i < n; i++) {
                map[i][2] = CPMGFitter.getMapIndex(states[i], stateCount, r2Mask) + lastCount;
                maxIndex = Math.max(map[i][2], maxIndex);
            }
            lastCount = maxIndex + 1;
            for (int i = 0; i < n; i++) {
                map[i][3] = CPMGFitter.getMapIndex(states[i], stateCount, 0, 3) + lastCount;
            }
            return map;
        }
    };
    final String equationName;
    final int nGroupPars;
    String[] parNames;

    public String getName() {
        return equationName;
    }

    public String[] getParNames() {
        return parNames;
    }

    public int getNGroupPars() {
        return nGroupPars;
    }

    @Override
    public double getMinX() {
        return 5.0;
    }

    CPMGEquation(String equationName, int nGroupPars, String... parNames) {
        this.equationName = equationName;
        this.parNames = parNames;
        this.nGroupPars = nGroupPars;
    }

    public static String[] getEquationNames() {
        String[] equationNames = {"NOEX", "CPMGFAST", "CPMGSLOW"};
        return equationNames;
    }

    public static double[] annCPMGGuesser(String exchangeType, double[] xValues, double[] yValues, double[] fields) throws Exception {
        //System.out.println("Using Neural Network!");
        double[] trainingX = {10.0, 20.0, 50.0, 100.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1100.0};
        double[] yANNInput = DataUtil.getCPMGInterpolation(trainingX, xValues, yValues);
        if (yANNInput != null) {
            int nFields = fields.length;
            exchangeType = exchangeType.toUpperCase();
            int nPars = (exchangeType.startsWith("FA")) ? 3 : (exchangeType.startsWith("SL")) ? 4 : 0;
            int annInputLength = nFields + yANNInput.length;
            double[] annInput = new double[annInputLength];
            double[] guesses = new double[nPars];
            String[] outputLabels = new String[nPars];
            String savedAnnFile = "";

            if (nPars == 3) {
                savedAnnFile = (nFields > 1) ? "data/ANNCPMGFast2.txt" : "data/ANNCPMGFast1.txt";
                outputLabels[0] = "kex";
                outputLabels[1] = "r2";
                outputLabels[2] = "dppmmin";
            } else if (nPars == 4) {
                savedAnnFile = (nFields > 1) ? "data/ANNCPMGSlow2.txt" : "data/ANNCPMGSlow1.txt";
                outputLabels[0] = "kex";
                outputLabels[1] = "pa";
                outputLabels[2] = "r2";
                outputLabels[3] = "dppm";
            } else {
                String msg = "Exchange type '{}' is not valid. Available options are 'FAST' or 'SLOW' exchange.".replace("{}", exchangeType);
                throw new Exception(msg);
            }

            ANNLoader ANN = ANNLoader.getInstance(savedAnnFile);
            ArtificialNeuralNetwork trainedNetwork = ANN.getTrainedNetwork();
            if (trainedNetwork != null) {
                HashMap scaleTracker = ANN.getScaleValues();
                if ((scaleTracker.containsKey("fields")) && (scaleTracker.containsKey("r2eff"))) {
                    for (int i = 0; i < nFields; i++) {
                        annInput[i] = Utilities.scale(fields[i], (double[]) scaleTracker.get("fields"), true);
                        //System.out.println("field " + fields[i] + " " + annInput[i]);
                    }
                    for (int j = 0; j < yANNInput.length; j++) {
                        annInput[j + nFields] = Utilities.scale(yANNInput[j], (double[]) scaleTracker.get("r2eff"), true);
                       // System.out.println("xy " + trainingX[j] + " " + yANNInput[j] + " " + annInput[j]);
                    }
                }
                Access1D wrappedInput = Access1D.wrap(annInput);
                MatrixStore<Double> result = trainedNetwork.invoke(wrappedInput);
                double[] resultArr = result.toRawCopy1D();

                int l = 0;
                for (String label : outputLabels) {
                    if (scaleTracker.containsKey(label)) {
                        guesses[l] = Utilities.scale(resultArr[l], (double[]) scaleTracker.get(label), false);
                        //System.out.println("guess " + guesses[l]);
                        l++;
                    }
                }
                return guesses;
            } else {
                throw new Exception("Could not load neural network.");
            }
        } else {
            return null;
        }
    }

    /**
     * Combines CPMG x and y values into a single matrix.
     *
     * @param xValues Matrix containing the offset values i.e. CEST/R1rho
     * irradiation frequency (X[0]), the B1 field values (X[1]), and the Tex
     * values (X[2]).
     * @param yValues Array of the CEST/R1rho intensities.
     * @param idNums Array of the ID numbers for the datasets.
     * @param id Integer ID number to retrieve the x and y values.
     * @return Matrix containing the offset (x) and intensity (y) values. X
     * values are in the first array, Y values are in the second.
     */
    public static double[][] getXYValues(double[][] xValues, double[] yValues, int[] idNums, int id) {
        int n = 0;
        for (int i = 0; i < idNums.length; i++) {
            if (idNums[i] == id) {
                n++;
            }
        }
        double[][] result = new double[2][n];
        int j = 0;
        for (int i = 0; i < idNums.length; i++) {
            if (idNums[i] == id) {
                result[0][j] = xValues[0][i];
                result[1][j] = yValues[i];
                j++;
            }
        }
        return result;

    }

}
