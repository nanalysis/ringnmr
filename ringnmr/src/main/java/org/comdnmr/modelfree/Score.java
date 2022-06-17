package org.comdnmr.modelfree;

public class Score {

    final double rss;
    final int nValues;
    final int nPars;
    final boolean parsOK;
    final double complexity;
    final double[] pars;

    public Score(double rss, int nValues, int nPars, boolean parsOK, double complexity) {
        this(rss, nValues, nPars, parsOK, complexity, null);
    }

    public Score(double rss, int nValues, int nPars, boolean parsOK, double complexity, double[] pars) {
        this.rss = rss;
        this.nValues = nValues;
        this.nPars = nPars;
        this.parsOK = parsOK;
        this.complexity = complexity;
        this.pars = pars;
    }

    public double[] getPars() {
        return pars;
    }

    public double rms() {
        return Math.sqrt(rss / nValues);
    }

    public int getN() {
        return nValues;
    }

    public double value() {
        return value(0.0);
    }

    public double value(double lambda) {
        double score = rms();
        if (!parsOK) {
            score += nValues * 10.0;
        }
        score += complexity * lambda;
        return score;
    }

    public double complexity() {
        return complexity;
    }

    public boolean parsOK() {
        return parsOK;
    }

    public double aic() {
        return 2 * nPars + nValues * Math.log(rss);
    }

    public double aicc() {
        int k = nPars;
        return aic() + 2.0 * k * (k + 1) / (nValues - k - 1.0);
    }

    public double chiSq() {
        return rss;
    }

    public double reducedChiSq() {
        return rss / (nValues - nPars);
    }
}
