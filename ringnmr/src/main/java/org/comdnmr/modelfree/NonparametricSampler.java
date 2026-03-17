package org.comdnmr.modelfree;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Class representing a nonparametric bootstrap sampler that extends the {@link BootstrapWeightSampler}.
 * This sampler generates weights in accordance with the protocol outlined by
 * Crawley and Palmer in their
 * <a href="https://mr.copernicus.org/articles/2/251/2021/">2021 paper on bootstrap aggregation</a>.
 *
 * @author simonhulse
 */
public class NonparametricSampler extends BootstrapWeightSampler {

    // Iterator of randomly ordered ints which specifies the ordering of possible weight vectors.
    private final Iterator<Integer> iterator;

    // Static map containing permitted combinations of weight assignments for a
    // given frequency class. See Table 1 in Crawley and Palmer's paper.
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
                {0, 1, 2, 3},
                {0, 0, 2, 3},
                {0, 1, 0, 3},
                {0, 1, 2, 0},
                {1, 1, 2, 3},
                {0, 1, 1, 3},
                {0, 1, 2, 1},
                {2, 1, 2, 3},
                {0, 2, 2, 3},
                {0, 1, 2, 2},
                {3, 1, 2, 3},
                {0, 3, 2, 3},
                {0, 1, 3, 3},
                {0, 0, 1, 1},
                {0, 0, 2, 2},
                {0, 0, 3, 3},
                {1, 1, 2, 2},
                {1, 1, 3, 3},
                {2, 2, 3, 3}
            }
        );
    }

    /**
     * Constructs a NonparametricSampler with a specified number of static fields.
     *
     * @param nExp the number of static fields (must be between 2 and 4)
     */
    public NonparametricSampler(int nExp) {
        super(nExp);
        iterator = generateIterator();
    }

    /**
     * Generates a randomized iterator over bootstrap sample indices.
     *
     * @return an iterator over the randomized bootstrap indices
     */
    private Iterator<Integer> generateIterator() {
        // Create a list of integers from 0 to N-1, shuffle it, and create an iterator
        List<Integer> randomNumbers = IntStream.range(0, getNBootstraps())
            .boxed()
            .collect(Collectors.toList());
        Collections.shuffle(randomNumbers, new Random());
        return randomNumbers.iterator();
    }

    /**
     * Returns the number of selections based on the number of static fields.
     *
     * @return the number of selections
     */
    public int getNSelections() {
        return getSelections().length;
    }

    /**
     * Returns the selection patterns, given `nExp()` static fields.
     *
     * @return a 2D array of selection patterns
     */
    private int[][] getSelections() {
        return SELECTION_MAP.get(nExp);
    }

    /**
     * Retrieves a specific row of selections based on the provided index.
     *
     * @param index the index of the selection row
     * @return an array representing the selection row
     */
    private int[] getSelectionsRow(int index) {
        return getSelections()[index];
    }

    /**
     * Calculates the total number of bootstrap samples.
     *
     * @return the total number of bootstraps
     */
    public int getNBootstraps() {
        int nSel = getNSelections();
        return nSel * nSel * nSel;
    }

    /**
     * Samples a set of nonparametric bootstrap weights.
     * This method can be called a maximum of {@link #getNBootstraps()} times.
     *
     * @return an array of sampled weights
     * @throws NoSuchElementException if the maximum number of samples is exceeded
     */
    public double[] sample() {
        if (!iterator.hasNext()) {
            throw new NoSuchElementException(
                String.format(
                    "Maximum number of samples exceeded. " +
                    "This sampler cannot be sampled more than %d^3 = %d times.",
                    getNSelections(),
                    getNBootstraps()
                )
            );
        }
        return getWeights(iterator.next());
    }

    /**
     * Retrieves weights based on the index.
     *
     * @param index the index to retrieve weights for
     * @return an array of weights corresponding to the selection
     */
    private double[] getWeights(int index) {
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
