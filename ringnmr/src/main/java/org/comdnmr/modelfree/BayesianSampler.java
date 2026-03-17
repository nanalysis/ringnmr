package org.comdnmr.modelfree;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.DirichletSampler;

public class BayesianSampler extends BootstrapWeightSampler {

    private final DirichletSampler sampler;

    public BayesianSampler(int nExp, UniformRandomProvider rng) {
        super(nExp);
        sampler = DirichletSampler.symmetric(rng, getNJ(), 1.0);
    }

    public double[] sample() {
        int nJ = getNJ();
        double[] weights = new double[nJ];
        int index = 0;
        for (double w : sampler.sample()) {
            weights[index++] = nJ * w;
        }
        return weights;
    }
}
