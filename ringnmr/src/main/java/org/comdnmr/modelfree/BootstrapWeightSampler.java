package org.comdnmr.modelfree;

import org.apache.commons.rng.sampling.ObjectSampler;


/**
 * Abstract class representing a bootstrap weight sampler.
 * This class implements the ObjectSampler interface to generate a double
 * arrays of weights for bootstrapping.
 *
 * @author simonhulse
 */
abstract class BootstrapWeightSampler implements ObjectSampler<double[]> {

    // Number of static fields R1, R2 and NOE data have been acquired with
    protected final int nExp;

    // Number of distinct frequency classes (0, N, 0.87H)
    protected static final int N_FREQS = 3;

    /**
     * Constructs a BootstrapWeightSampler with a specified number of experiments.
     *
     * @param nExp the number of static fields data was acquired at
     * @throws IllegalArgumentException if `nExp` is less than 2 or greater than 4
     */
    public BootstrapWeightSampler(int nExp) throws IllegalArgumentException {
        if (nExp < 2 || nExp > 4) {
            throw new IllegalArgumentException(String.format("Invalid nFreq: %d", nExp));
        }
        this.nExp = nExp;
    }

    /**
     * Abstract method for sampling bootstrap weights.
     * Concrete implementations must provide their own sampling logic.
     *
     * @return an array of double values representing sampled weights
     */
    abstract public double[] sample();

    /**
     * Gets the number of static fields used data has been acquired at.
     *
     * @return the number of experiments.
     */
    public int getNExp() {
        return nExp;
    }

    /**
     * Calculates and returns the total number of spectral density values.
     * Given by `3 * getNExp()`.
     *
     * @return the total number of parameters (NJ)
     */
    public int getNJ() {
        return N_FREQS * nExp;
    }
}
