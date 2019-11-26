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

import org.comdnmr.util.CoMDPreferences;
import java.util.Arrays;
import java.util.List;
import static org.comdnmr.util.Utilities.TWO_PI;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;

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
                        deltaA[i] = (deltaA0 - omegarf[i]) * fields[i] * TWO_PI;
                        deltaB[i] = (deltaB0 - omegarf[i]) * fields[i] * TWO_PI;
                        omegaB1[i] = b1Field[i] * TWO_PI;
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
                        deltaA[i] = (deltaA0 - omegarf[i]) * fields[i] * TWO_PI;
                        omegaB1[i] = b1Field[i] * TWO_PI;

                    }
                    double[] yCalc = R1RhoEquations.r1rhoPerturbationNoEx(omegaB1, deltaA, R1A, R2A);
                    return yCalc;
                }

                @Override
                public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double[] fields) {
                    int nPars = R1RhoFitFunction.getNPars(map);
                    double[] guesses = new double[nPars];
                    for (int id = 0; id < map.length; id++) {
                        int[] map1 = map[id];
                        double[][] xy = CESTEquations.getXYValues(xValues, yValues, idNums, id);
                        double field = fields[id];
                        List<CESTPeak> peaks = CESTEquations.cestPeakGuess(xy[0], xy[1], field, "r1rho");
                        if (peaks.size() > 0) {
                            double tex = xValues[2][0];
                            double[] r1 = CESTEquations.cestR1Guess(yValues, tex, "r1rho");
                            double[][] r2 = CESTEquations.cestR2Guess(peaks, yValues, "r1rho");
                            guesses[map1[0]] = peaks.get(peaks.size() - 1).position; //400 * 2.0 * Math.PI; //deltaA
                            guesses[map1[1]] = r1[0]; //2.4; //R1A
                            guesses[map1[2]] = r2[0][0] / 2; //20.0; //R2A
                        } else {
                            guesses = null;
                        }
                    }
                    return guesses;

                }

                @Override
                public double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
                    double[][] boundaries = new double[2][guesses.length];
                    int id = 0;
                    int[] map1 = map[id];
                    double[][] xy = CESTEquations.getXYValues(xValues, yValues, idNums, id);
                    List<CESTPeak> peaks = CESTEquations.cestPeakGuess(xy[0], xy[1], field, "r1rho");

                    double dAbound = peaks.get(peaks.size() - 1).getWidths()[1] / 2;
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
                        deltaA[i] = (deltaA0 - omegarf[i]) * fields[i] * TWO_PI;
                        deltaB[i] = (deltaB0 - omegarf[i]) * fields[i] * TWO_PI;
                        omegaB1[i] = b1Field[i] * TWO_PI;
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
                        deltaA[i] = (deltaA0 - omegarf[i]) * fields[i] * TWO_PI;
                        deltaB[i] = (deltaB0 - omegarf[i]) * fields[i] * TWO_PI;
                        omegaB1[i] = b1Field[i] * TWO_PI;
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
                    double k1 = pb * kex;
                    double km1 = (1 - pb) * kex;
                    double[][] K = {{-k1, 0, 0, km1, 0, 0}, 
                        {0, -k1, 0, 0, km1, 0}, 
                        {0, 0, -k1, 0, 0, km1}, 
                        {k1, 0, 0, -km1, 0, 0}, 
                        {0, k1, 0, 0, -km1, 0}, 
                        {0, 0, k1, 0, 0, -km1}};
                    double[][] La = {{-R2A, 0, 0, 0, 0, 0}, //{-R2A, -deltaA, 0, 0, 0, 0}
                        {0, -R2A, 0, 0, 0, 0}, //{deltaA, -R2A, -omegaB1[i], 0, 0, 0}
                        {0, 0, -R1A, 0, 0, 0}, //{0, omegaB1[i], -R1A, 0, 0, 0}
                        {0, 0, 0, 0, 0, 0},
                        {0, 0, 0, 0, 0, 0},
                        {0, 0, 0, 0, 0, 0}};

                    double[][] Lb = {{0, 0, 0, 0, 0, 0},
                        {0, 0, 0, 0, 0, 0},
                        {0, 0, 0, 0, 0, 0},
                        {0, 0, 0, -R2B, 0, 0}, //{0, 0, 0, -R2B, -deltaB, 0}
                        {0, 0, 0, 0, -R2B, 0}, //{0, 0, 0, deltaB, -R2B, -omegaB1[i]}
                        {0, 0, 0, 0, 0, -R1B}}; //{0, 0, 0, 0, omegaB1[i], -R1B}

                    DMatrixRMaj La1 = new DMatrixRMaj(La);
                    DMatrixRMaj Lb1 = new DMatrixRMaj(Lb);
                    DMatrixRMaj K1 = new DMatrixRMaj(K);

                    DMatrixRMaj Z = new DMatrixRMaj(La.length, La[0].length);
                    for (int i = 0; i < yCalc.length; i++) {
                        double omegaB1 = X[1][i] * TWO_PI;
                        double deltaA = (deltaA0 - X[0][i]) * fields[i] * TWO_PI;
                        double deltaB = (deltaB0 - X[0][i]) * fields[i] * TWO_PI;
                        La1.set(0, 1, -deltaA); 
                        La1.set(1, 0, deltaA);
                        La1.set(1, 2, -omegaB1);
                        La1.set(2, 1, omegaB1);

                        Lb1.set(3, 4, -deltaB); 
                        Lb1.set(4, 3, deltaB);
                        Lb1.set(4, 5, -omegaB1);
                        Lb1.set(5, 4, omegaB1);         

                        CommonOps_DDRM.add(La1, Lb1, Z);
                        CommonOps_DDRM.addEquals(Z, K1);
//                        yCalc[i] = R1RhoEquations.r1rhoExact(omegaB1, pb, kex, deltaA, deltaB, R1A, R1B, R2A, R2B);
                        yCalc[i] = R1RhoEquations.r1rhoExact(Z);
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
                    double pA = 1.0 - pb;
                    double kAB = pb * kex;
                    double kBA = pA * kex;
            
                    double[][] K = {
                        {-kAB, 0, 0, kBA, 0, 0},
                        {0, -kAB, 0, 0, kBA, 0},
                        {0, 0, -kAB, 0, 0, kBA},
                        {kAB, 0, 0, -kBA, 0, 0},
                        {0, kAB, 0, 0, -kBA, 0},
                        {0, 0, kAB, 0, 0, -kBA}};

                    double[][] La = {{-R2A, 0, 0, 0, 0, 0}, //{-R2A, -deltaA, 0, 0, 0, 0}
                        {0, -R2A, 0, 0, 0, 0}, //{deltaA, -R2A, -omegaB1[i], 0, 0, 0}
                        {0, 0, -R1A, 0, 0, 0}, //{0, omegaB1[i], -R1A, 0, 0, 0}
                        {0, 0, 0, 0, 0, 0},
                        {0, 0, 0, 0, 0, 0},
                        {0, 0, 0, 0, 0, 0}};

                    double[][] Lb = {{0, 0, 0, 0, 0, 0},
                        {0, 0, 0, 0, 0, 0},
                        {0, 0, 0, 0, 0, 0},
                        {0, 0, 0, -R2B, 0, 0}, //{0, 0, 0, -R2B, -deltaB, 0}
                        {0, 0, 0, 0, -R2B, 0}, //{0, 0, 0, deltaB, -R2B, -omegaB1[i]}
                        {0, 0, 0, 0, 0, -R1B}}; //{0, 0, 0, 0, omegaB1[i], -R1B}
                    
                    DMatrixRMaj La1 = new DMatrixRMaj(La);
                    DMatrixRMaj Lb1 = new DMatrixRMaj(Lb);
                    DMatrixRMaj K1 = new DMatrixRMaj(K);

                    DMatrixRMaj Z = new DMatrixRMaj(La.length, La[0].length);
                    double[] m0 = new double[La.length];
                    double[] m1 = new double[La.length];
                    for (int i = 0; i < yCalc.length; i++) {
                        double omegaB1 = X[1][i] * TWO_PI;
                        double deltaA = (deltaA0 - X[0][i]) * fields[i] * TWO_PI;
                        double deltaB = (deltaB0 - X[0][i]) * fields[i] * TWO_PI;
                        double theta = Math.atan2(omegaB1, deltaA);
                        double cosA = Math.cos(theta);
                        double sinA = Math.sin(theta);

                        m0[0] = pA * sinA;
                        m0[2] = pA * cosA;
                        m1[0] = sinA;
                        m1[2] = cosA;
                        
                        La1.set(0, 1, -deltaA); 
                        La1.set(1, 0, deltaA);
                        La1.set(1, 2, -omegaB1);
                        La1.set(2, 1, omegaB1);

                        Lb1.set(3, 4, -deltaB); 
                        Lb1.set(4, 3, deltaB);
                        Lb1.set(4, 5, -omegaB1);
                        Lb1.set(5, 4, omegaB1);         

                        CommonOps_DDRM.add(La1, Lb1, Z);
                        CommonOps_DDRM.addEquals(Z, K1);
                        CommonOps_DDRM.scale(delay, Z);
//                        double r1rho = R1RhoEquations.r1rhoExact0(delay, omegaB1, pb, kex, deltaA, deltaB, R1A, R1B, R2A, R2B);
                        double r1rho = R1RhoEquations.r1rhoExact0(Z, m0, m1, delay);
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
