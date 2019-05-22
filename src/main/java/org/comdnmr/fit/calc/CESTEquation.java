package org.comdnmr.fit.calc;

import java.util.Arrays;
import java.util.List;

/**
 *
 * @author Bruce Johnson
 */
public enum CESTEquation implements CESTEquationType {

    CESTR1RHOPERTURBATION("cestR1rhoPerturbation", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

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
            double[] yCalc = CESTEquations.cestR1rhoApprox("trott", X, fields, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
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
    CESTR1RHOPERTURBATIONNOEX("cestR1rhoPerturbationNoEx", 0, "deltaA0", "R1A", "R2A") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double[] fields) {
            double deltaA0 = par[map[0]];
            double R1A = par[map[1]];
            double R2A = par[map[2]];
            double[] yCalc = CESTEquations.cestR1rhoPerturbationNoEx(X, fields, deltaA0, R1A, R2A);
            return yCalc;
        }

        @Override
        public double[] guess(double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            int nPars = CalcCEST.getNPars(map);
            double[] guesses = new double[nPars];
            for (int id = 0; id < map.length; id++) {
                int[] map1 = map[id];
                double[][] xy = CESTEquations.getXYValues(xValues, yValues, idNums, id);
                List<CESTPeak> peaks = CESTEquations.cestPeakGuess(xy[0], xy[1], field, "cest");
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
        public double[][] boundaries(double[] guesses, double[][] xValues, double[] yValues, int[][] map, int[] idNums, int nID, double field) {
            double[][] boundaries = new double[2][guesses.length];
            int id = 0;
            int[] map1 = map[id];
            double[][] xy = CESTEquations.getXYValues(xValues, yValues, idNums, id);
            List<CESTPeak> peaks = CESTEquations.cestPeakGuess(xy[0], xy[1], field, "cest");

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
    CESTR1RHOSD("cestR1rhoSD", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

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
            double[] yCalc = CESTEquations.cestR1rhoApprox("sd", X, fields, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
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
    CESTR1RHOBALDWINKAY("cestR1rhoBaldwinKay", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double[] field) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.cestR1rhoApprox("baldwinkay", X, field, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
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
    CESTR1RHON("cestR1rhoN", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double[] field) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.cestR1rhoApprox("laguerre", X, field, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
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
    CESTR1RHOEXACT1("cestR1rhoExact1", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double[] field) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.cestR1rhoExact1(X, field, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
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
    CESTEXACT0("cestExact0", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double field[]) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.cestExact0(X, field, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
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
    CESTEXACT1("cestExact1", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double[] field) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.cestExact0(X, field, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
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
    CESTEXACT2("cestExact2", 0, "Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B") {

        @Override
        public double[] calculate(double[] par, int[] map, double[][] X, int idNum, double[] field) {
            double kex = par[map[0]];
            double pb = par[map[1]];
            double deltaA0 = par[map[2]];
            double deltaB0 = par[map[3]];
            double R1A = par[map[4]];
            double R1B = par[map[5]];
            double R2A = par[map[6]];
            double R2B = par[map[7]];
            double[] yCalc = CESTEquations.cestExact0(X, field, pb, kex, deltaA0, deltaB0, R1A, R1B, R2A, R2B);
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
        String[] equationNames = {"CESTR1RHOPERTURBATIONNOEX", "CESTR1RHOPERTURBATION", "CESTR1RHOSD", "CESTR1RHOBALDWINKAY", "CESTR1RHON",
            "CESTR1RHOEXACT1", "CESTEXACT0", "CESTEXACT1", "CESTEXACT2"};
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
