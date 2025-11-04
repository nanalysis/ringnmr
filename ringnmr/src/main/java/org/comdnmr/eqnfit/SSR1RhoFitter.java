package org.comdnmr.eqnfit;

import org.comdnmr.util.CoMDOptions;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;

public class SSR1RhoFitter implements EquationFitter {

    FitFunction fitFunc;
    CoMDOptions options;

    public SSR1RhoFitter(CoMDOptions options) {
        fitFunc = new SSR1RhoFitFunction(options);
        this.options = options;
    }

    @Override
    public List<String> getEquationNameList() {
        return Arrays.asList(SSR1RhoEquation.getAllEquationNames());
    }

    @Override
    public FitFunction getFitModel() {
        return fitFunc;
    }

    @Override
    public String getExpType() {
        return "";
    }

    @Override
    public Optional<FitResult> doFit(String eqn, double[] sliderGuesses, CoMDOptions options) {
        return Optional.empty();
    }

    @Override
    public int[] getStateCount() {
        return new int[0];
    }

    @Override
    public int[][] getStates() {
        return new int[0][];
    }

    // NOTE: Copied over from R1RhoFitter
    // TODO: Could make a utility method called linespace to do this in one line.
    // Probably useful elsewhere in the codebase too
    @Override
    public double[] getSimX(int nPts, double xLB, double xUB) {
        int nPoints = nPts;
        double[] x = new double[nPoints];
        double firstValue = xLB;
        double lastValue = xUB;
        double delta = (lastValue - firstValue) / (nPoints + 1);
        double value = firstValue;
        for (int i = 0; i < nPoints; i++) {
            x[i] = value;
            value += delta;

        }
        return x;
    }

    @Override
    public double rms(double[] pars) {
        return 0;
    }

    @Override
    public void setData(List<Double>[] allXValues, List<Double> yValues, List<Double> errValues) {

    }

    @Override
    public List<ParValueInterface> guessPars(String eqn) {
        return List.of();
    }

    @Override
    public double[] getSimXDefaults() {
        return new double[0];
    }
}
