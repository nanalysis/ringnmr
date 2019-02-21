package org.comdnmr.fit.calc;

import java.util.ArrayList;
import java.util.stream.IntStream;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.SynchronizedRandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;

public class CalcCEST extends FitModel {

    static RandomGenerator random = new SynchronizedRandomGenerator(new Well19937c());
    int[] r2Mask = {0, 1, 3};
    double[] rexErrors = new double[nID];
    CESTEquations cestEq = new CESTEquations();
    double[][] parValues;

    public CalcCEST() {
        this.equation = CESTEquation.CESTR1RHOPERTURBATION;
    }

    public void setEquation(String eqName) {
        equation = CESTEquation.valueOf(eqName.toUpperCase());
    }

    public CalcCEST(double[][] x, double[] y, double[] err, double[] fieldValues, double[] fields) throws IllegalArgumentException {
        this(x, y, err, fieldValues, fields, new int[x.length]);
    }

    public CalcCEST(double[][] x, double[] y, double[] err, double[] fieldValues, double[] fields, int[] idNums) throws IllegalArgumentException {
        this.xValues = new double[x.length][];
        this.xValues[0] = x[0].clone();
        this.xValues[1] = x[1].clone();
        this.xValues[2] = x[2].clone();
        this.yValues = y.clone();
        this.errValues = err.clone();
        this.fieldValues = fieldValues.clone();
        this.fields = fields.clone();
        this.idNums = idNums.clone();
        this.idNums = new int[yValues.length];
        this.equation = CESTEquation.CESTR1RHOPERTURBATION;
        if (setNID()) {
            throw new IllegalArgumentException("Invalid idNums, some values not used");
        }
    }

    public static int getNPars(int[][] map) {
        int maxIndex = 0;
        for (int i = 0; i < map.length; i++) {
            for (int j = 0; j < map[i].length; j++) {
                maxIndex = Math.max(map[i][j], maxIndex);
            }
        }
        return maxIndex + 1;
    }

    public int[] getMask() {
        return r2Mask;
    }

    @Override
    public double value(double[] par) {
        double sumAbs = 0.0;
        double sumSq = 0.0;
        double[] yCalc = equation.calculate(par, map[0], xValues, 0, fields[0]);

        for (int i = 0; i < yValues.length; i++) {
            double delta = (yCalc[i] - yValues[i]);
            sumAbs += FastMath.abs(delta);
            sumSq += delta * delta;
        }
//        for (double p:par) {
//            System.out.print(p + " ");
//        }
//        System.out.println(Math.sqrt(sumSq/yValues.length));
        if (absMode) {
            return sumAbs;
        } else {
            return sumSq;
        }
    }

    public double[] getPredicted(double[] par) {
        double[] yPred = simY(par);
        return yPred;
    }

    @Override
    public double[][] getSimPars() {
        return parValues;
    }

    public double[] simBounds(double[] start, double[] lowerBounds, double[] upperBounds, double[] inputSigma) {
        reportFitness = false;
        int nPar = start.length;
        int nSim = CoMDPreferences.getSampleSize();
        parValues = new double[nPar][nSim];
        double[] yPred = getPredicted(start);
        double[] yValuesOrig = yValues.clone();
        double[][] rexValues = new double[nID][nSim];
        rexErrors = new double[nID];
        for (int i = 0; i < nSim; i++) {
            for (int k = 0; k < yValues.length; k++) {
                yValues[k] = yPred[k] + errValues[k] * random.nextGaussian();
            }
            PointValuePair result = refine(start, lowerBounds, upperBounds, inputSigma);
            double[] rPoint = result.getPoint();
            for (int j = 0; j < nPar; j++) {
                parValues[j][i] = rPoint[j];
            }
        }
        double[] parSDev = new double[nPar];
        for (int i = 0; i < nPar; i++) {
            DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
            double p5 = dStat.getPercentile(5.0);
            double p95 = dStat.getPercentile(95.0);
            parSDev[i] = dStat.getStandardDeviation();
        }
        yValues = yValuesOrig;
        return parSDev;
    }

    public double[] simBoundsStream(double[] start, double[] lowerBounds, double[] upperBounds, double[] inputSigma) {
        reportFitness = false;
        int nPar = start.length;
        int nSim = CoMDPreferences.getSampleSize();
        parValues = new double[nPar][nSim];
        double[][] rexValues = new double[nID][nSim];
        rexErrors = new double[nID];
        double[] yPred = getPredicted(start);
        IntStream.range(0, nSim).parallel().forEach(i -> {
//        IntStream.range(0, nSim).forEach(i -> {
            CalcCEST rDisp = new CalcCEST(xValues, yPred, errValues, fieldValues, fields, idNums);
            rDisp.setEquation(equation.getName());
            rDisp.setAbsMode(absMode);
            double[] newY = new double[yValues.length];
            for (int k = 0; k < yValues.length; k++) {
                newY[k] = yPred[k] + errValues[k] * random.nextGaussian();
            }
            rDisp.setXY(xValues, newY);
            rDisp.setIds(idNums);
            rDisp.setMap(map);

            PointValuePair result = rDisp.refine(start, lowerBounds, upperBounds, inputSigma);
            double[] rPoint = result.getPoint();
            for (int j = 0; j < nPar; j++) {
                parValues[j][i] = rPoint[j];
            }
        });

        double[] parSDev = new double[nPar];
        for (int i = 0; i < nPar; i++) {
            DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
            parSDev[i] = dStat.getStandardDeviation();
        }
        return parSDev;
    }

    public double[] simBoundsBootstrapStream(double[] start, double[] lowerBounds, double[] upperBounds, double[] inputSigma) {
        reportFitness = false;
        int nPar = start.length;
        int nSim = CoMDPreferences.getSampleSize();
        parValues = new double[nPar][nSim];
        double[][] rexValues = new double[nID][nSim];
        rexErrors = new double[nID];
        IntStream.range(0, nSim).parallel().forEach(i -> {
            CalcCEST rDisp = new CalcCEST(xValues, yValues, errValues, fieldValues, fields, idNums);
            rDisp.setEquation(equation.getName());
            rDisp.setAbsMode(absMode);
            double[][] newX = new double[3][yValues.length];
            double[] newY = new double[yValues.length];
            double[] newErr = new double[yValues.length];
            double[] newFieldValues = new double[yValues.length];
            int[] newID = new int[yValues.length];
            int iTry = 0;
            do {
                for (int k = 0; k < yValues.length; k++) {
                    int rI = random.nextInt(yValues.length);
                    newX[0][k] = xValues[0][rI];
                    newX[1][k] = xValues[1][rI];
                    newX[2][k] = xValues[2][rI];
                    newY[k] = yValues[rI];
                    newErr[k] = errValues[rI];
                    newFieldValues[k] = fieldValues[rI];
                    newID[k] = idNums[rI];
                }
                iTry++;
            } while (!checkID(idNums, newID, 2) && (iTry < 10));
            // fixme  idNum should be set in above loop
            rDisp.setXY(newX, newY);
            rDisp.setErr(newErr);
            rDisp.setFieldValues(newFieldValues);
            rDisp.setIds(newID);
            rDisp.setMap(map);

            PointValuePair result = rDisp.refine(start, lowerBounds, upperBounds, inputSigma);
            double[] rPoint = result.getPoint();
            for (int j = 0; j < nPar; j++) {
                parValues[j][i] = rPoint[j];
            }
        });

        double[] parSDev = new double[nPar];
        for (int i = 0; i < nPar; i++) {
            DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
            parSDev[i] = dStat.getStandardDeviation();
        }
        return parSDev;
    }

}
