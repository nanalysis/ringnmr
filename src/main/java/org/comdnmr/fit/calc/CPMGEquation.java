/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

import org.apache.commons.math3.util.FastMath;

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
        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field, boolean neuralNetworkGuess) {
            int nPars = CalcRDisp.getNPars(map);
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
                map[i][0] = CPMGFit.getMapIndex(states[i], stateCount, r2Mask);
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
        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field, boolean neuralNetworkGuess) {
            int nPars = CalcRDisp.getNPars(map);
            double[] guesses = new double[nPars];
            double kExSum = 0.0;
            int nKex = 0;
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
                guesses[map[id][1]] = r2;
                double tauMid = 1.0 / (2.0 * vMid);
                double kEx = 1.915 / (0.5 * tauMid);
                double dPPMMinRad = Math.sqrt(4.0 * rex / (field * field) * kEx);
                double dPPMMin = dPPMMinRad / (2.0 * Math.PI);
                guesses[map[id][2]] = dPPMMin;
                if (rex >= 0) {
                    kExSum += kEx; // 1.915 comes from solving equation iteratively at tcp rex 0.5 half max
                    nKex++;
                }
            }
            guesses[0] = kExSum / nKex;
            if (guesses[0] > CoMDPreferences.getCPMGMaxFreq()) {
                guesses[0] = CoMDPreferences.getCPMGMaxFreq() * 0.9;
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
                map[i][1] = CPMGFit.getMapIndex(states[i], stateCount, r2Mask) + lastCount;
                maxIndex = Math.max(map[i][1], maxIndex);
            }
            lastCount = maxIndex + 1;
            for (int i = 0; i < n; i++) {
                map[i][2] = CPMGFit.getMapIndex(states[i], stateCount, 0, 3) + lastCount;
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
            double ch = dP * FastMath.cosh(etaP) - dM * FastMath.cos(etaM);
            double rexContrib = 0.5 * (kEx - (1.0 / tauCP) * FastMath.acosh(ch));
            double dR = (1.0 / tauCP) * FastMath.acosh(ch);
            double value = r2 + rexContrib;
            if (Double.isInfinite(value)) {
                System.out.println("infi");
                System.out.println(kEx + " pa " + pA + " r2 " + r2 + " dw " + dW + " psi " + psi + " rex " + rexContrib + " ch " + ch + " eta " + etaP + " " + etaM + " " + eta1 + " " + psi);
            }
            return value;
        }

        @Override
        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field, boolean neuralNetworkGuess) {
            int nPars = CalcRDisp.getNPars(map);
            double[] guesses = new double[nPars];
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
                double dw2 = rex / (pa * (1.0 - pa)) * kex;
                double dPPM = Math.sqrt(dw2) / (2.0 * Math.PI) / field;
                guesses[map[id][2]] = r2;
                guesses[map[id][3]] = dPPM;
                kExSum += kex;
            }
            guesses[0] = kExSum /= nID;
            guesses[1] = pa;
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
                boundaries[0][iPar] = 0.5;
                boundaries[1][iPar] = 0.99;
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
                map[i][2] = CPMGFit.getMapIndex(states[i], stateCount, r2Mask) + lastCount;
                maxIndex = Math.max(map[i][2], maxIndex);
            }
            lastCount = maxIndex + 1;
            for (int i = 0; i < n; i++) {
                map[i][3] = CPMGFit.getMapIndex(states[i], stateCount, 0, 3) + lastCount;
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

}
