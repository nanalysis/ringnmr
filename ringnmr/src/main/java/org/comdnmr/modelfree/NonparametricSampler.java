package org.comdnmr.modelfree;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class NonparametricSampler extends BootstrapWeightSampler {

    private final Iterator<Integer> iterator;

    private static final Map<Integer, int[][]> SELECTION_MAP;
    static {
        SELECTION_MAP = new HashMap<>();
        SELECTION_MAP.put(
            2,
            new int[][]{
                {0, 0},
                {0, 1},
                {1, 1}
            }
        );
        SELECTION_MAP.put(
            3,
            new int[][]{
                {0, 1, 2},
                {1, 1, 2},
                {0, 0, 2},
                {2, 1, 2},
                {0, 1, 0},
                {0, 2, 2},
                {0, 1, 1}
            }
        );
        SELECTION_MAP.put(
            4,
            new int[][]{
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
                {2,2,3,3}
            }
        );
    }

    public NonparametricSampler(int nExp) {
        super(nExp);
        iterator = generateIterator();
    }

    private Iterator<Integer> generateIterator() {
        // Create a list of integers from 0 to N-1, shuffle it, and create an iterator
        List<Integer> randomNumbers = IntStream.range(0, getNBootstraps())
            .boxed()
            .collect(Collectors.toList());
        Collections.shuffle(randomNumbers, new Random());
        return randomNumbers.iterator();
    }


    public int getNSelections() {
        return getSelections().length;
    }

    private int[][] getSelections() {
        return SELECTION_MAP.get(nExp);
    }

    private int[] getSelectionsRow(int index) {
        return getSelections()[index];
    }

    public int getNBootstraps() {
        int nSel = getNSelections();
        return nSel * nSel * nSel;
    }

    public double[] sample() {
        if (!iterator.hasNext()) {
            throw new NoSuchElementException("Maximum number of samples exceeded.");
        }
        return getWeights(iterator.next());
    }

    public double[] getWeights(int index) {
        int nSel = getNSelections();
        int nSelSq = nSel * nSel;

        // J(0) pointer
        int p0 = index / nSelSq;
        int r0 = index % nSelSq;
        // J(wN) pointer
        int pN = r0 / nSel;
        // J(0.87wH) pointer
        int pH = r0 % nSel;

        int[][] selections = {
            getSelectionsRow(p0),
            getSelectionsRow(pN),
            getSelectionsRow(pH)
        };

        double[] weights = new double[getNJ()];
        for (int i = 0; i < N_FREQS; i++) {
            for (int j = 0; j < getNExp(); j++) {
                weights[selections[i][j] * N_FREQS + i] += 1.0;
            }
        }

        return weights;
    }
}
