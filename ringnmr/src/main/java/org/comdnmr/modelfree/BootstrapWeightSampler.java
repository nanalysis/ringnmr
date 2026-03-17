package org.comdnmr.modelfree;

import org.apache.commons.rng.sampling.ObjectSampler;

abstract class BootstrapWeightSampler implements ObjectSampler<double[]> {

    protected int nExp;
    protected static final int N_FREQS = 3;

    public BootstrapWeightSampler(int nExp) throws IllegalArgumentException {
        if (nExp < 2 || nExp > 4) {
            throw new IllegalArgumentException(String.format("Invalid nFreq: %d", nExp));
        }
        this.nExp = nExp;
    }

    abstract public double[] sample();

    public int getNExp() {
        return nExp;
    }

    public int getNJ() {
        return N_FREQS * nExp;
    }
}
