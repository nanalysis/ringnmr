package org.comdnmr.fit.calc;

import java.util.Arrays;
import java.util.stream.IntStream;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;

/**
 *
 * @author Bruce Johnson
 */
public abstract class FitModel implements MultivariateFunction {

    // fixme is there a thread safe RandomGenerator
    public final RandomGenerator DEFAULT_RANDOMGENERATOR = new MersenneTwister(1);
    int reportAt = 10;

    EquationType equation;
    long startTime = 0;
    double[][] xValues;
    double[] fieldValues;
    double[] yValues;
    double[] errValues;
    double[] fields;
    int[] idNums;
    int[][] map;
    int nSim = 200;
    int nID = 1;
    boolean reportFitness = false;
    boolean absMode = false;
    double[][] parValues;

    public class Checker extends SimpleValueChecker {

        public Checker(double relativeThreshold, double absoluteThreshold, int maxIter) {
            super(relativeThreshold, absoluteThreshold, maxIter);
        }

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

    public abstract void setEquation(String eqName);

    public PointValuePair refine(double[] guess, double[] lowerBounds, double[] upperBounds, double[] inputSigma) {
        startTime = System.currentTimeMillis();
        DEFAULT_RANDOMGENERATOR.setSeed(1);
        double lambdaMul = 3.0;
        int lambda = (int) (lambdaMul * FastMath.round(4 + 3 * FastMath.log(guess.length)));
        //int nSteps = guess.length*1000;
        int nSteps = 2000;
        double stopFitness = 0.0;
        int diagOnly = 0;
        double tol = 1.0e-5;
        //new Checker(100 * Precision.EPSILON, 100 * Precision.SAFE_MIN, nSteps));
        CMAESOptimizer optimizer = new CMAESOptimizer(nSteps, stopFitness, true, diagOnly, 0,
                DEFAULT_RANDOMGENERATOR, true,
                new CalcRDisp.Checker(tol, tol, nSteps));
        PointValuePair result = null;

        try {
            result = optimizer.optimize(
                    new CMAESOptimizer.PopulationSize(lambda),
                    new CMAESOptimizer.Sigma(inputSigma),
                    new MaxEval(2000000),
                    new ObjectiveFunction(this), GoalType.MINIMIZE,
                    new SimpleBounds(lowerBounds, upperBounds),
                    new InitialGuess(guess));
        } catch (Exception e) {
            e.printStackTrace();
        }
        return result;

    }

    public double[] simY(double[] par) {
        equation.setFieldRef(fields[0]);
        double[] yCalc = equation.calculate(par, map[0], xValues, idNums[0], fields[0]);
        return yCalc;
    }

    public abstract int[] getMask();

    public abstract double[] simBounds(double[] start, double[] lowerBounds, double[] upperBounds, double[] inputSigma);

    public abstract double[] simBoundsStream(double[] start, double[] lowerBounds, double[] upperBounds, double[] inputSigma);

    public abstract double[] simBoundsBootstrapStream(double[] start, double[] lowerBounds, double[] upperBounds, double[] inputSigma);
    
    public abstract double[][] getSimPars();

    public void setXY(double[][] x, double[] y) {
        this.xValues = x;
        this.yValues = y;
        this.idNums = new int[x[0].length];
    }

    public void setErr(double[] err) {
        this.errValues = err;
    }

    public void setFieldValues(double[] fieldValues) {
        this.fieldValues = fieldValues;
    }

    boolean setNID() {
        nID = Arrays.stream(idNums).max().getAsInt() + 1;
        boolean[] checkIDs = new boolean[nID];
        Arrays.stream(idNums).forEach(id -> checkIDs[id] = true);
        return IntStream.range(0, checkIDs.length).anyMatch(id -> checkIDs[id] == false);
    }

    public void setIds(int[] idNums) throws IllegalArgumentException {
        this.idNums = idNums;
        if (setNID()) {
            for (int id : idNums) {
                System.out.print(id + " ");
            }
            System.out.println("");
            throw new IllegalArgumentException("Invalid idNums, some values not used in setIds");
        }
    }

    public void setFields(double[] fields) {
        this.fields = fields;
        this.equation.setFieldRef(fields);
    }

    public double[] guess() {
        double[] guess = equation.guess(xValues, yValues, map, idNums, nID, fieldValues[0]);
        return guess;
    }

    public double[][] boundaries() {
        //return equation.boundaries(xValues, yValues, fieldValues[0]);
        double[][] boundaries = equation.boundaries(xValues, yValues, map, idNums, nID, fieldValues[0]);
        return boundaries;
    }

    public double calculate(double[] par, int[] map, double[] x, int idNum, double field) {
        return equation.calculate(par, map, x, idNum, field);
    }

    public double getRSS(double[] par) {
        double rss = 0.0;
        double[] x = new double[xValues.length];
        for (int i = 0; i < yValues.length; i++) {
            for (int j = 0; j < x.length; j++) {
                x[j] = xValues[j][i];
            }
            final double value;
            value = calculate(par, map[idNums[i]], x, idNums[i], fieldValues[i]);
            double delta = value - yValues[i];
            rss += delta * delta;
        }
        return rss;
    }

    public double getRMS(double[] par) {
        double rss = 0.0;
        double[] x = new double[xValues.length];
        for (int i = 0; i < yValues.length; i++) {
            for (int j = 0; j < x.length; j++) {
                x[j] = xValues[j][i];
            }
            final double value;
            value = calculate(par, map[idNums[i]], x, idNums[i], fieldValues[i]);

            double delta = value - yValues[i];
            rss += delta * delta;
        }
        return Math.sqrt(rss / xValues.length);
    }

    public double getAICc(double[] par) {
        double rss = getRSS(par);
        int k = par.length;
        int n = xValues.length;
        double aic = 2 * k + n * Math.log(rss);
        double aicc = aic + 2 * k * (k + 1) / (n - k - 1);
        return aicc;
    }

    public void setAbsMode(boolean value) {
        this.absMode = value;
    }

    public void setNSim(int value) {
        this.nSim = value;
    }

    public double[] getPredicted(double[] par) {
        double[] yPred = simY(par);
        return yPred;
    }

    public void setMap(int[][] map) {
        this.map = map;
    }

    public void setMap(int[] stateCount, int[][] states) {
        this.map = new int[states.length][stateCount.length];
        this.map = equation.makeMap(stateCount, states, getMask());
        for (int i = 0; i < map.length; i++) {
//            for (int j = 0; j < map[i].length; j++) {
//                System.out.print(" " + map[i][j]);
//            }
//            System.out.println(" map");
        }
    }

    public int[][] getMap() {
        return map;
    }

    public String[] getParNames() {
        return equation.getParNames();
    }
        
    public int getNGroupPars() {
        return equation.getNGroupPars();
    }

    public boolean checkID(int[] idValues, int[] testValues, int checkNum) {
        int n = Arrays.stream(idValues).max().getAsInt() + 1;
        int[] count = new int[n];
        for (int i = 0; i < testValues.length; i++) {
            count[testValues[i]]++;
        }
        boolean ok = true;
        for (int i = 0; i < count.length; i++) {
            if (count[i] < checkNum) {
                ok = false;
                break;
            }
        }
        return ok;
    }

    public void dump(double[] par, double[][] xValues, double field) {
        equation.setFieldRef(field);
        double[] x = new double[xValues.length];
        for (int i = 0; i < yValues.length; i++) {
            for (int j = 0; j < x.length; j++) {
                x[j] = xValues[j][i];
            }
            double yCalc = equation.calculate(par, map[idNums[i]], x, idNums[i], field);
            System.out.printf("%8.5f %8.5f %8.5f\n", x[0], yCalc, field);
        }
    }

    public void dump(double[] par) {
        for (double field : fields) {
            System.out.println("# field " + field);
            for (int i = 0; i < xValues.length; i++) {
                if (FastMath.abs(fieldValues[i] - field) < 0.01) {
                    double yCalc = equation.calculate(par, map[idNums[i]], xValues[i], idNums[i], field);
                    System.out.printf("%8.5f %8.5f %8.5f\n", xValues[i], yValues[i], yCalc);
                }
            }
        }

        DescriptiveStatistics dstat = new DescriptiveStatistics(xValues[0]);
        double min = dstat.getMin();
        double max = dstat.getMax();
        int n = 20;
        double[] ax = new double[1];
        double delta = (max - min) / (n - 1);
        for (double field : fields) {
            System.out.println("# field " + field);
            for (int i = 0; i < n; i++) {
                double x = min + i * delta;
                ax[0] = x;
                double y = equation.calculate(par, map[idNums[i]], ax, idNums[i], field);
                System.out.printf("%8.5f %8.5f\n", x, y);
            }
        }
    }
}
