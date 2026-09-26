package org.comdnmr.data;

import java.util.Arrays;
import java.util.function.BiFunction;
import java.util.stream.IntStream;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NotPositiveException;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.DoubleAccumulator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.SynchronizedRandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;
import org.comdnmr.util.CoMDPreferences;

public class Fitter {

    RandomGenerator random = new Well19937c();

    boolean reportFitness = false;
    int reportAt = 10;
    double[][] parValues;
    double[][] xValues;
    double[] yValues;
    double[] errValues;

    double[] lowerBounds;
    double[] upperBounds;
    double[] start;
    double inputSigma;
    BiFunction<double[], double[][], Double> valuesFunction = null;

    private Fitter() {

    }

    public static Fitter getArrayFitter(BiFunction<double[], double[][], Double> function) {
        Fitter fitter = new Fitter();
        fitter.valuesFunction = function;
        return fitter;
    }

    public double value(double[] par) {
        // setup values array in case we've passed in a functin that uses it

        int nA = xValues.length + 2;
        double[][] values = new double[nA][];
        for (int i = 0; i < xValues.length; i++) {
            values[i] = xValues[i];
        }
        values[nA - 2] = yValues;
        values[nA - 1] = errValues;
        return valuesFunction.apply(par, values);
    }

    // ---- incumbent-guard diagnostics ------------------------------------
    // Atomic because FitModel fits residues on a parallelStream.
    private static final AtomicLong GUARD_CALLS = new AtomicLong();
    private static final AtomicLong GUARD_FIRES = new AtomicLong();
    private static final DoubleAccumulator GUARD_WORST_RATIO =
            new DoubleAccumulator(Math::max, 1.0);
    /** Set true to print a line every time the guard rejects a fit. */
    public static volatile boolean GUARD_VERBOSE = false;

    public static void guardReset() {
        GUARD_CALLS.set(0);
        GUARD_FIRES.set(0);
        GUARD_WORST_RATIO.reset();
    }

    /** Human-readable summary of how often a fit came back worse than its start. */
    public static String guardReport() {
        long c = GUARD_CALLS.get();
        long f = GUARD_FIRES.get();
        return String.format(
                "GUARD warm-started fits %d, rejected %d (%.1f%%), worst returned/start ratio %.2f",
                c, f, c > 0 ? 100.0 * f / c : 0.0, GUARD_WORST_RATIO.get());
    }

    public PointValuePair fit(double[] start, double[] lowerBounds, double[] upperBounds, double inputSigma, int nTry) throws Exception {
        random.setSeed(1);
        this.start = start;
        this.lowerBounds = lowerBounds.clone();
        this.upperBounds = upperBounds.clone();
        this.inputSigma = inputSigma;
        Optimizer opt = new Optimizer();
        if (xValues != null) {
            opt.setXYE(xValues, yValues, errValues);
        }
        PointValuePair bestPair = null;
        double[] bestStart = new double[start.length];
        double[] tryStart = new double[start.length];
        double range = 0.9;
        for (int iTry = 0; iTry < nTry; iTry++) {
            PointValuePair result;
            if (iTry > 0) {
                for (int j = 0; j < start.length; j++) {
                    double delta = (upperBounds[j] - lowerBounds[j]) * range;
                    tryStart[j] = lowerBounds[j] + delta * (1.0 - range) / 2.0 + random.nextDouble() * delta;
                }
            } else {
                System.arraycopy(start, 0, tryStart, 0, start.length);
            }
            if (CoMDPreferences.getOptimizer().equals("BOBYQA")) {
                result = opt.refineBOBYQA(tryStart, inputSigma);
            } else {
                result = opt.refineCMAES(tryStart, inputSigma);
            }
            if ((bestPair == null) || (result.getValue() < bestPair.getValue())) {
                bestPair = result;
                System.arraycopy(tryStart, 0, bestStart, 0, tryStart.length);
            }
        }
        // Incumbent guard.  CMA-ES evaluates its sampled population, never the
        // initial guess itself, so with a large start radius it can return a
        // point worse than where it began -- fatal for the warm-started
        // refinement passes, which then lock in the worse answer.  Warm starts
        // were observed going from rss 3.4 to 22.8 this way.  Never return
        // worse than the caller's starting point.
        if ((bestPair != null) && (start != null)) {
            GUARD_CALLS.incrementAndGet();
            double startValue = opt.value(opt.normalize(start));
            if (Double.isFinite(startValue) && (startValue < bestPair.getValue())) {
                GUARD_FIRES.incrementAndGet();
                if (startValue > 0.0) {
                    GUARD_WORST_RATIO.accumulate(bestPair.getValue() / startValue);
                }
                if (GUARD_VERBOSE) {
                    System.out.printf(
                            "GUARDFIRE nTry %d sigma %.1f returned %.6f start %.6f (%.2fx worse)%n",
                            nTry, inputSigma, bestPair.getValue(), startValue,
                            startValue > 0.0 ? bestPair.getValue() / startValue : Double.NaN);
                }
                bestPair = new PointValuePair(start.clone(), startValue);
            }
        }
        return bestPair;
    }

    public void setXYE(double[][] xValues, double[] yValues, double[] errValues) {
        this.xValues = xValues;
        this.yValues = yValues;
        this.errValues = errValues;
    }

    class Optimizer implements MultivariateFunction {

        RandomGenerator random = new SynchronizedRandomGenerator(new Well19937c());

        public class Checker extends SimpleValueChecker {

            public Checker(double relativeThreshold, double absoluteThreshold, int maxIter) {
                super(relativeThreshold, absoluteThreshold, maxIter);
            }

            @Override
            public boolean converged(final int iteration, final PointValuePair previous, final PointValuePair current) {
                boolean converged = super.converged(iteration, previous, current);
                if (reportFitness) {
                    if (converged) {
                        System.out.println(previous.getValue() + " " + current.getValue());
                    }
                    if (converged || (iteration == 1) || ((iteration % reportAt) == 0)) {
                        long time = System.currentTimeMillis();
                        long deltaTime = time - startTime;
                        System.out.println(deltaTime + " " + iteration + " " + current.getValue());
                    }
                }
                return converged;
            }
        }

        double[][] xValues;
        double[] yValues;
        double[] errValues;
        double[][] values;
        long startTime;
        long endTime;
        long fitTime;
        double tol = 1.0e-5;
        boolean absMode = false;
        boolean weightFit = false;

        @Override
        public double value(double[] normPar) {
            double[] par = deNormalize(normPar);
            return valuesFunction.apply(par, values);
        }

        void fixGuesses(double[] guesses) {
            for (int i = 0; i < guesses.length; i++) {
                if (guesses[i] > 98.0) {
                    guesses[i] = 98.0;
                } else if (guesses[i] < 2) {
                    guesses[i] = 2.0;
                }
            }
        }

        double[] normalize(double[] pars) {
            double[] normPars = new double[pars.length];
            for (int i = 0; i < pars.length; i++) {
                normPars[i] = 100.0 * (pars[i] - lowerBounds[i]) / (upperBounds[i] - lowerBounds[i]);
            }
            return normPars;
        }

        double[] deNormalize(double[] normPars) {
            double[] pars = new double[normPars.length];
            for (int i = 0; i < pars.length; i++) {
                pars[i] = normPars[i] / 100.0 * (upperBounds[i] - lowerBounds[i]) + lowerBounds[i];
            }
            return pars;
        }

        void setXYE(double[][] xValues, double[] yValues, double[] errValues) {
            this.xValues = xValues;
            this.yValues = yValues;
            this.errValues = errValues;
            // setup values array in case we've passed in a functin that uses it

            int nA = xValues.length + 2;
            this.values = new double[nA][];
            for (int i = 0; i < xValues.length; i++) {
                this.values[i] = xValues[i];
            }
            this.values[nA - 2] = yValues;
            this.values[nA - 1] = errValues;
        }

        public PointValuePair refineCMAES(double[] guess, double inputSigma) throws Exception {
            startTime = System.currentTimeMillis();
            random.setSeed(1);
            double lambdaMul = 10.0;
            int lambda = (int) (lambdaMul * FastMath.round(4 + 3 * FastMath.log(guess.length)));
            //int nSteps = guess.length*1000;
            int nSteps = 2000;
            double stopFitness = 0.0;
            int diagOnly = 0;
            double[] normLower = new double[guess.length];
            double[] normUpper = new double[guess.length];
            double[] sigma = new double[guess.length];
            Arrays.fill(normLower, 0.0);
            Arrays.fill(normUpper, 100.0);
            Arrays.fill(sigma, inputSigma);
            double[] normGuess = normalize(guess);
            fixGuesses(normGuess);

            //new Checker(100 * Precision.EPSILON, 100 * Precision.SAFE_MIN, nSteps));
            CMAESOptimizer cmaesOptimizer = new CMAESOptimizer(nSteps, stopFitness, true, diagOnly, 0,
                    random, true,
                    new Checker(tol, tol, nSteps));
            PointValuePair result = null;

            try {
                result = cmaesOptimizer.optimize(
                        new CMAESOptimizer.PopulationSize(lambda),
                        new CMAESOptimizer.Sigma(sigma),
                        new MaxEval(2000000),
                        new ObjectiveFunction(this), GoalType.MINIMIZE,
                        new SimpleBounds(normLower, normUpper),
                        new InitialGuess(normGuess));
            } catch (DimensionMismatchException | NotPositiveException | NotStrictlyPositiveException |
                     TooManyEvaluationsException e) {
                throw new Exception("failure to fit data " + e.getMessage());
            }
            endTime = System.currentTimeMillis();
            fitTime = endTime - startTime;
            PointValuePair deNormResult = new PointValuePair(deNormalize(result.getPoint()), result.getValue());

            return deNormResult;
        }

        public PointValuePair refineBOBYQA(double[] guess, double inputSigma) {
            startTime = System.currentTimeMillis();
            random.setSeed(1);
            double lambdaMul = 3.0;
            int lambda = (int) (lambdaMul * FastMath.round(4 + 3 * FastMath.log(guess.length)));
            //int nSteps = guess.length*1000;
            int nSteps = 2000;
            double stopFitness = 0.0;
            int diagOnly = 0;
            double[] normLower = new double[guess.length];
            double[] normUpper = new double[guess.length];
            Arrays.fill(normLower, 0.0);
            Arrays.fill(normUpper, 100.0);
            double[] normGuess = normalize(guess);
            fixGuesses(normGuess);
            //new Checker(100 * Precision.EPSILON, 100 * Precision.SAFE_MIN, nSteps));

            int n = guess.length;
            int nInterp = 2 * n + 1;
            double initialRadius = inputSigma;
            double stopRadius = CoMDPreferences.getFinalRadius();
            stopRadius = Math.pow(10.0, stopRadius);

            BOBYQAOptimizer optimizer = new BOBYQAOptimizer(nInterp, initialRadius, stopRadius);
            PointValuePair result = null;

            result = optimizer.optimize(
                    new MaxEval(2000000),
                    new ObjectiveFunction(this), GoalType.MINIMIZE,
                    new SimpleBounds(normLower, normUpper),
                    new InitialGuess(normGuess));
            endTime = System.currentTimeMillis();
            fitTime = endTime - startTime;
            PointValuePair deNormResult = new PointValuePair(deNormalize(result.getPoint()), result.getValue());
            return deNormResult;
        }
    }

    public double[] bootstrap(double[] guess, int nSim, boolean parametric, double[] yPred) {
        reportFitness = false;
        int nPar = start.length;
        parValues = new double[nPar + 1][nSim];

        IntStream.range(0, nSim).parallel().forEach(iSim -> {
            double[][] newX = new double[xValues.length][yValues.length];
            double[] newY = new double[yValues.length];
            double[] newErr = new double[yValues.length];
            Optimizer optimizer = new Optimizer();
            for (int iValue = 0; iValue < yValues.length; iValue++) {
                int rI = random.nextInt(yValues.length);
                for (int xIndex = 0; xIndex < newX.length; xIndex++) {
                    newX[xIndex][iValue] = xValues[xIndex][rI];
                }
                if (parametric) {
//                    System.out.println(iValue + ": yPred " + yPred[iValue] + " " + "errVal " + errValues[iValue]);
                    newY[iValue] = yPred[iValue] + errValues[iValue] * random.nextGaussian(); //parametric
                } else {
                    newY[iValue] = yValues[rI]; //non-parametric
                }
                newErr[iValue] = errValues[rI];
            }

            // fixme  idNum should be set in above loop
            optimizer.setXYE(newX, newY, newErr);

            PointValuePair result;
            try {
                if (CoMDPreferences.getOptimizer().equals("BOBYQA")) {
                    result = optimizer.refineBOBYQA(guess, inputSigma);
                } else {
                    result = optimizer.refineCMAES(guess, inputSigma);
                }
            } catch (Exception ex) {
                return;
            }
            double[] rPoint = result.getPoint();
            for (int j = 0; j < nPar; j++) {
                parValues[j][iSim] = rPoint[j];
            }
            parValues[nPar][iSim] = result.getValue();
        });

        double[] parSDev = new double[nPar];
        for (int i = 0; i < nPar; i++) {
            DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
            parSDev[i] = dStat.getStandardDeviation();
        }
        return parSDev;
    }
}
