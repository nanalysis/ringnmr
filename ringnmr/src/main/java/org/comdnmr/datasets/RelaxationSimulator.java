package org.comdnmr.datasets;

import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

/**
 * Simulates "experimental" R1, R2 and NOE measurements from known ("true")
 * relaxation rates: R1 from a noisy inversion-recovery profile, R2 from a
 * noisy decaying exponential, and NOE from two noisy intensity measurements
 * (saturated and reference). R1 and R2 are recovered by nonlinear
 * least-squares fitting of the noisy profile, with uncertainties taken from
 * the fit's parameter covariance; NOE is recovered as an intensity ratio,
 * with its uncertainty from standard error propagation.
 */
class RelaxationSimulator {

    record FitResult(double value, double error) {}

    // Arbitrary reference signal intensity (units are irrelevant; only ratios/relative noise matter).
    private static final double REFERENCE_INTENSITY = 1.0;

    // Delay times, as fractions of 1/R1 or 1/R2, at which synthetic intensities are sampled.
    private static final double[] INV_DELAY_FRACTIONS = {0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.4, 2.0, 3.0};
    private static final double[] EXPDECAY_DELAY_FRACTIONS = {0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.4, 2.0};

    private static final double COV_SINGULARITY_THRESHOLD = 1.0e-14;
    private static final int MAX_EVALUATIONS = 1000;
    private static final int MAX_ITERATIONS = 1000;

    private RelaxationSimulator() {}

    static FitResult simulateInversionRecovery(double r1True, double noiseLevel, Random rng) {
        double[] delays = scale(INV_DELAY_FRACTIONS, 1.0 / r1True);
        double sigma = noiseLevel * REFERENCE_INTENSITY;
        double[] intensities = new double[delays.length];
        for (int i = 0; i < delays.length; i++) {
            double trueIntensity = REFERENCE_INTENSITY * (1.0 - 2.0 * Math.exp(-r1True * delays[i]));
            intensities[i] = trueIntensity + sigma * rng.nextGaussian();
        }
        double[] guess = {REFERENCE_INTENSITY, r1True * 1.2};
        return fitExponential(delays, intensities, guess, sigma, true);
    }

    static FitResult simulateExpDecay(double r2True, double noiseLevel, Random rng) {
        double[] delays = scale(EXPDECAY_DELAY_FRACTIONS, 1.0 / r2True);
        double sigma = noiseLevel * REFERENCE_INTENSITY;
        double[] intensities = new double[delays.length];
        for (int i = 0; i < delays.length; i++) {
            double trueIntensity = REFERENCE_INTENSITY * Math.exp(-r2True * delays[i]);
            intensities[i] = trueIntensity + sigma * rng.nextGaussian();
        }
        double[] guess = {REFERENCE_INTENSITY, r2True * 1.2};
        return fitExponential(delays, intensities, guess, sigma, false);
    }

    static FitResult simulateNOE(double noeTrue, double noiseLevel, Random rng) {
        double iSatTrue = noeTrue * REFERENCE_INTENSITY;
        double sigma = noiseLevel * REFERENCE_INTENSITY;

        double iUnsat = REFERENCE_INTENSITY + sigma * rng.nextGaussian();
        double iSat = iSatTrue + sigma * rng.nextGaussian();

        double noe = iSat / iUnsat;
        double error = Math.abs(noe) * Math.sqrt(
            (sigma / iSat) * (sigma / iSat) + (sigma / iUnsat) * (sigma / iUnsat)
        );
        return new FitResult(noe, error);
    }

    /**
     * Fits I(t) = I0 * (1 - 2*exp(-R*t)) (inversionRecovery=true) or
     * I(t) = I0*exp(-R*t) (inversionRecovery=false) to (delays, intensities)
     * via Levenberg-Marquardt, returning the fitted rate R and its standard
     * error.
     */
    // Must weight by 1/sigma^2: commons-math3's getSigma() reports the covariance implied by
    // the weights as-is (no residual/chi-square rescaling), so an unweighted fit would report
    // a fixed, noise-independent "error" even when the fit is exact.
    private static FitResult fitExponential(double[] delays, double[] intensities, double[] guess, double sigma, boolean inversionRecovery) {
        MultivariateJacobianFunction model = point -> {
            double i0 = point.getEntry(0);
            double rate = point.getEntry(1);
            RealVector value = new ArrayRealVector(delays.length);
            RealMatrix jacobian = new Array2DRowRealMatrix(delays.length, 2);
            for (int i = 0; i < delays.length; i++) {
                double t = delays[i];
                double expTerm = Math.exp(-rate * t);
                double predicted;
                double dI0;
                double dRate;
                if (inversionRecovery) {
                    predicted = i0 * (1.0 - 2.0 * expTerm);
                    dI0 = 1.0 - 2.0 * expTerm;
                    dRate = 2.0 * i0 * t * expTerm;
                } else {
                    predicted = i0 * expTerm;
                    dI0 = expTerm;
                    dRate = -i0 * t * expTerm;
                }
                value.setEntry(i, predicted);
                jacobian.setEntry(i, 0, dI0);
                jacobian.setEntry(i, 1, dRate);
            }
            return new Pair<>(value, jacobian);
        };

        LeastSquaresBuilder builder = new LeastSquaresBuilder()
            .start(guess)
            .model(model)
            .target(intensities)
            .maxEvaluations(MAX_EVALUATIONS)
            .maxIterations(MAX_ITERATIONS);
        if (sigma > 0.0) {
            double[] weights = new double[intensities.length];
            Arrays.fill(weights, 1.0 / (sigma * sigma));
            builder.weight(MatrixUtils.createRealDiagonalMatrix(weights));
        }

        LeastSquaresProblem.Evaluation evaluation = new LevenbergMarquardtOptimizer().optimize(builder.build());
        double rate = evaluation.getPoint().getEntry(1);
        double error = sigma > 0.0 ? evaluation.getSigma(COV_SINGULARITY_THRESHOLD).getEntry(1) : 0.0;
        return new FitResult(rate, error);
    }

    private static double[] scale(double[] fractions, double factor) {
        double[] result = new double[fractions.length];
        for (int i = 0; i < fractions.length; i++) {
            result[i] = fractions[i] * factor;
        }
        return result;
    }
}
