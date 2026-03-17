package org.comdnmr.modelfree;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.DirichletSampler;

/**
 * Class representing a Bayesian sampler that extends the {@link BootstrapWeightSampler}.
 * This class utilizes a symmetric Dirichlet distribution to generate weights.
 *
 * @author simonhulse
 */
public class BayesianSampler extends BootstrapWeightSampler {

    // Concentration factor
    private static final double ALPHA = 1.0;

    // Dirichlet sampler instance
    private final DirichletSampler sampler;

    /**
     * Constructs a BayesianSampler for data with a specified number of static
     * fields and random number generator.
     *
     * @param nExp the number of static fields (must be between 2 and 4)
     * @param rng the random number generator used for sampling
     */
    public BayesianSampler(int nExp, UniformRandomProvider rng) {
        super(nExp);
        sampler = DirichletSampler.symmetric(rng, getNJ(), ALPHA);
    }

    /**
     * Generates a sample of bootstrapping weights.
     *
     * @return an array of double values representing the sampled weights.
     *         The sum of the weights is given by {@link #getNJ()}.
     */
    public double[] sample() {
        int nJ = getNJ();
        double[] weights = sampler.sample();
        // Weights for Dirichlet distribution sum to 1.
        // Scale these so that they sum to nJ.
        for (int i = 0; i < nJ; i++) weights[i] *= nJ;
        return weights;
    }
}
