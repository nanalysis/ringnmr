package org.comdnmr.modelfree;

import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;

public class SimonsScore extends Score {

    public CMAESOptimizer optimizer;
    private long runtime;

    public SimonsScore(
        double rss,
        int nValues,
        int nPars,
        boolean parsOK,
        double complexityS,
        double complexityTauF,
        double complexityTauS,
        double[] pars,
        CMAESOptimizer optimizer
    ) {
        super(rss, nValues, nPars, parsOK, complexityS, complexityTauF, complexityTauS, pars);
        this.optimizer = optimizer;
    }

    public CMAESOptimizer getOptimizer() { return optimizer; }

    public void setRuntime(long runtime) { this.runtime = runtime; }

    public long getRuntime() { return runtime; }
}
