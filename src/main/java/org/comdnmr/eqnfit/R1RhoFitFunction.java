/*
 * CoMD/NMR Software : A Program for Analyzing NMR Dynamics Data
 * Copyright (C) 2018-2019 Bruce A Johnson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
package org.comdnmr.eqnfit;

import org.comdnmr.util.CoMDPreferences;
import java.util.stream.IntStream;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.SynchronizedRandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;

public class R1RhoFitFunction extends FitFunction {

    static RandomGenerator random = new SynchronizedRandomGenerator(new Well19937c());
    int[] r2Mask = {0, 1, 3};
    double[] rexErrors = new double[nID];
    R1RhoEquations r1rhoEq = new R1RhoEquations();

    public R1RhoFitFunction() {
        this.equation = R1RhoEquation.R1RHOPERTURBATION;
    }

    @Override
    public void setEquation(String eqName) {
        equation = R1RhoEquation.valueOf(eqName.toUpperCase());
    }

    public R1RhoFitFunction(double[][] x, double[] y, double[] err, double[] fieldValues) throws IllegalArgumentException {
        this(x, y, err, fieldValues, new int[x.length]);
    }

    public R1RhoFitFunction(double[][] x, double[] y, double[] err, double[] fieldValues, int[] idNums) throws IllegalArgumentException {
        this.xValues = new double[x.length][];
        this.xValues[0] = x[0].clone();
        this.xValues[1] = x[1].clone();
        this.xValues[2] = x[2].clone();
        this.yValues = y.clone();
        this.errValues = err.clone();
        this.fieldValues = fieldValues.clone();
        this.idNums = idNums.clone();
        this.idNums = new int[yValues.length];
        this.equation = R1RhoEquation.R1RHOPERTURBATION;
        if (setNID()) {
            throw new IllegalArgumentException("Invalid idNums, some values not used");
        }
    }

    public static int getNPars(int[][] map) {
        int maxIndex = 0;
        for (int[] map1 : map) {
            for (int j = 0; j < map1.length; j++) {
                maxIndex = Math.max(map1[j], maxIndex);
            }
        }
        return maxIndex + 1;
    }

    @Override
    public int[] getMask() {
        return r2Mask;
    }

    @Override
    public double value(double[] normPar) {
        double[] par = deNormalize(normPar);

        double sumAbs = 0.0;
        double sumSq = 0.0;
        double[] yCalc = new double[yValues.length];
        for (int id = 0; id < map.length; id++) {
            double[][] x = CESTEquations.getXValues(xValues, idNums, id);
            double[] fields = CESTEquations.getValues(fieldValues, idNums, id);
            double[] yCalc1 = equation.calculate(par, map[id], x, id, fields);
            int[] indicies = CESTEquations.getIndicies(idNums, id);
            for (int i = 0; i < indicies.length; i++) {
                yCalc[indicies[i]] = yCalc1[i];
            }
        }

        for (int i = 0; i < yValues.length; i++) {
            double delta = (yCalc[i] - yValues[i]);
            if (weightFit) {
                delta /= errValues[i];
            }
            sumAbs += FastMath.abs(delta);
            sumSq += delta * delta;
        }
        if (absMode) {
            return sumAbs / (yValues.length - par.length);
        } else {
            return sumSq / (yValues.length - par.length);
        }
    }

    @Override
    public double[] getPredicted(double[] par) {
        double[] yPred = simY(par);
        return yPred;
    }

    @Override
    public double[][] getSimPars() {
        return parValues;
    }

    @Override
    public double[] simBounds(double[] start, double[] lowerBounds, double[] upperBounds, double inputSigma) {
        reportFitness = false;
        int nPar = start.length;
        int nSim = CoMDPreferences.getSampleSize();
        parValues = new double[nPar + 1][nSim];
        double[] yPred = getPredicted(start);
        double[] yValuesOrig = yValues.clone();
        double[][] rexValues = new double[nID][nSim];
        rexErrors = new double[nID];
        String optimizer = CoMDPreferences.getBootStrapOptimizer();

        for (int i = 0; i < nSim; i++) {
            for (int k = 0; k < yValues.length; k++) {
                yValues[k] = yPred[k] + errValues[k] * random.nextGaussian();
            }
            PointValuePair result = refine(start, lowerBounds, upperBounds,
                    inputSigma, optimizer);
            double[] rPoint = result.getPoint();
            for (int j = 0; j < nPar; j++) {
                parValues[j][i] = rPoint[j];
            }
            parValues[nPar][i] = result.getValue();
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

    @Override
    public double[] simBoundsStream(double[] start, double[] lowerBounds,
            double[] upperBounds, double inputSigma, boolean nonParametric) {
        if (nonParametric) {
            return simBoundsStreamNonParametric(start, lowerBounds, upperBounds, inputSigma);
        } else {
            return simBoundsStreamParametric(start, lowerBounds, upperBounds, inputSigma);
        }

    }

    public double[] simBoundsStreamParametric(double[] start, double[] lowerBounds, double[] upperBounds, double inputSigma) {
        reportFitness = false;
        int nPar = start.length;
        int nSim = CoMDPreferences.getSampleSize();
        parValues = new double[nPar + 1][nSim];
        double[][] rexValues = new double[nID][nSim];
        rexErrors = new double[nID];
        double[] yPred = getPredicted(start);
        String optimizer = CoMDPreferences.getBootStrapOptimizer();

        IntStream.range(0, nSim).parallel().forEach(i -> {
//        IntStream.range(0, nSim).forEach(i -> {
            R1RhoFitFunction rDisp = new R1RhoFitFunction(xValues, yPred, errValues, fieldValues, idNums);
            rDisp.setEquation(equation.getName());
            double[] newY = new double[yValues.length];
            for (int k = 0; k < yValues.length; k++) {
                newY[k] = yPred[k] + errValues[k] * random.nextGaussian();
            }
            rDisp.setXY(xValues, newY);
            rDisp.setIds(idNums);
            rDisp.setMap(map);

            PointValuePair result = rDisp.refine(start, lowerBounds, upperBounds,
                    inputSigma, optimizer);
            double[] rPoint = result.getPoint();
            for (int j = 0; j < nPar; j++) {
                parValues[j][i] = rPoint[j];
            }
            parValues[nPar][i] = result.getValue();
        });

        double[] parSDev = new double[nPar];
        for (int i = 0; i < nPar; i++) {
            DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
            parSDev[i] = dStat.getStandardDeviation();
        }
        return parSDev;
    }

    public double[] simBoundsStreamNonParametric(double[] start, double[] lowerBounds, double[] upperBounds, double inputSigma) {
        reportFitness = false;
        int nPar = start.length;
        int nSim = CoMDPreferences.getSampleSize();
        parValues = new double[nPar + 1][nSim];
        double[][] rexValues = new double[nID][nSim];
        rexErrors = new double[nID];
        String optimizer = CoMDPreferences.getBootStrapOptimizer();

        IntStream.range(0, nSim).parallel().forEach(i -> {
            R1RhoFitFunction rDisp = new R1RhoFitFunction(xValues, yValues, errValues, fieldValues, idNums);
            rDisp.setEquation(equation.getName());
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

            PointValuePair result = rDisp.refine(start, lowerBounds, upperBounds,
                    inputSigma, optimizer);
            double[] rPoint = result.getPoint();
            for (int j = 0; j < nPar; j++) {
                parValues[j][i] = rPoint[j];
            }
            parValues[nPar][i] = result.getValue();

        });

        double[] parSDev = new double[nPar];
        for (int i = 0; i < nPar; i++) {
            DescriptiveStatistics dStat = new DescriptiveStatistics(parValues[i]);
            parSDev[i] = dStat.getStandardDeviation();
        }
        return parSDev;
    }

}
