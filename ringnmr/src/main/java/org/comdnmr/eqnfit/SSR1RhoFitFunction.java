package org.comdnmr.eqnfit;

import org.apache.commons.math3.util.FastMath;
import org.comdnmr.util.CoMDOptions;

import java.util.Optional;

public class SSR1RhoFitFunction extends FitFunction {

    public SSR1RhoFitFunction(CoMDOptions options) {
        super(options);
        this.equation = SSR1RhoEquation.CSA;
    }

    @Override
    public void setEquation(String eqName) {
        equation = SSR1RhoEquation.valueOf(eqName.toUpperCase());
    }

    // TODO
    @Override
    public int[] getMask() {
        return new int[0];
    }

    // TODO: simBounds has no usages. Delete from FitFunction
    @Override
    public Optional<double[]> simBounds(double[] start, double[] lowerBounds, double[] upperBounds, double inputSigma, CoMDOptions options) {
        return Optional.empty();
    }

    // TODO: Will need to be properly implemented for bootstrap results
    @Override
    public Optional<double[]> simBoundsStream(double[] start, double[] lowerBounds, double[] upperBounds, double inputSigma, CoMDOptions options) {
        return Optional.empty();
    }

    @Override
    public double[][] getSimPars() { return parValues; }

    // Currently only set up for simple case, without any mappings/IDs
    @Override
    public double value(double[] normPar) {
        double[] par = deNormalize(normPar);
        double[][] x = xValues;
        double[] ys = equation.calculate(par, map[0], x, 0);

        double result = 0.0;
        for (int i = 0; i < ys.length; i++) {
            double delta = ys[i] - yValues[i];
            if (weightFit) delta /= errValues[i];
            if (absMode) {
                result += FastMath.abs(delta);
            } else {
                result += FastMath.pow(delta, 2);
            }
        }
        return result / (yValues.length - par.length);
    }
}
