package org.comdnmr.modelfree;

public class BootstrapAggregator {
    int nFreq = 3;
    final int nExp;
    final int nJ;
    final int nSel;
    final int nBootStrap;
    final int[][] selections;
    static final int[][] SELECTIONS4 = {
            {0,1,2,3},

            {0,0,2,3},
            {0,1,0,3},
            {0,1,2,0},

            {1,1,2,3},
            {0,1,1,3},
            {0,1,2,1},

            {2,1,2,3},
            {0,2,2,3},
            {0,1,2,2},

            {3,1,2,3},
            {0,3,2,3},
            {0,1,3,3},

            {0,0,1,1},
            {0,0,2,2},
            {0,0,3,3},

            {1,1,2,2},
            {1,1,3,3},

            {2,2,3,3},
    };

    static final int[][] SELECTIONS3 = {
            {0, 1, 2},
            {1, 1, 2},
            {0, 0, 2},
            {2, 1, 2},
            {0, 1, 0},
            {0, 2, 2},
            {0, 1, 1}
    };

    static final int[][] SELECTIONS2 = {
            {0, 0},
            {0, 1},
            {1, 1}
    };
    public BootstrapAggregator(int nExp) throws IllegalArgumentException {
        if ((nExp != 2) && (nExp != 3) && (nExp != 4)) {
            throw new IllegalArgumentException("Invalid nFreq " + nExp);
        }
        this.nExp = nExp;
        this.nJ = nExp * nFreq;
        if (nExp == 2) {
            selections = SELECTIONS2;
        } else if (nExp == 3) {
                selections = SELECTIONS3;
        } else {
            selections = SELECTIONS4;
        }
        nSel = selections.length;
        nBootStrap = nSel * nSel * nSel;
    }

    public int getN() {
        return nBootStrap;
    }

    public  int[] getY( int iSelection) {
        int[] counts = new int[nJ];
        int[][] sel = getSelections(iSelection);
        for (int iFreq = 0; iFreq < nFreq; iFreq++) {
            for (int iExp = 0; iExp < nExp; iExp++) {
                int swapIndex = sel[iFreq][iExp] * nFreq + iFreq;
                counts[swapIndex]++;
            }
        }
        return counts;
    }

    public static void incrCounts(int[] totalCounts, int[] counts) {
        for (int i=0;i<counts.length;i++) {
            totalCounts[i] += counts[i];
        }
    }

    public  int[] getTotalUse() {
        int[] total = new int[nJ];
        for (int i=0;i<getN();i++) {
            int[] counts = getY( i);
            incrCounts(total, counts);
        }
        boolean ok = true;
        for (int i=0;i<nJ;i++) {
            if (total[i] != getN()) {
                ok = false;
                break;
            }
        }
        if (ok) {
            return new int[0];
        } else {
            return total;
        }
    }

    public int[][] getSelections(int iSelection) {
        int i0 = iSelection / (nSel * nSel);
        int r0 = iSelection % (nSel * nSel);
        int i1 = r0 / nSel;
        int i2 = r0 % nSel;
        int[] f0 = selections[i0];
        int[] f1 = selections[i1];
        int[] f2 = selections[i2];
        return new int[][]{f0, f1, f2};
    }

    public  double[][] getBootStrapJ(double[][] jValues, int iSelection) {
        int[][] sel = getSelections(iSelection);
        double[][] result = new double[3][nJ];
        for (int iFreq = 0; iFreq < nFreq; iFreq++) {
            for (int iExp = 0; iExp < nExp; iExp++) {
                int index = iExp * nFreq + iFreq;
                int swapIndex = sel[iFreq][iExp] * nFreq + iFreq;
                for (int jType = 0; jType < 3; jType++) {  // type is freq, value, error
                    result[jType][index] = jValues[jType][swapIndex];
                }
            }
        }
        return result;
    }
}
