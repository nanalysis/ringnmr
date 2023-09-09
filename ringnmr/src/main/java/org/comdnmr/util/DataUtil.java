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
package org.comdnmr.util;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.linear.*;

import java.util.Arrays;
import java.util.DoubleSummaryStatistics;
import java.util.stream.IntStream;

public class DataUtil {

    public static double[][] clone2DArray(double[][] values) {
        double[][] result = new double[values.length][];
        for (int i = 0; i < result.length; i++) {
            result[i] = values[i].clone();
        }
        return result;
    }

    public static double getMeanValue(double[] values, int[] idValues, int id) {
        DoubleSummaryStatistics dStat = IntStream.range(0, values.length).filter(i -> idValues[i] == id).
                mapToDouble(i -> values[i]).summaryStatistics();
        return dStat.getAverage();
    }

    public static double getMinValue(double[] values) {
        double minVal = Arrays.stream(values).min().getAsDouble();
        return minVal;
    }

    public static double getMinValue(double[] values, int[] idValues, int id) {
        double minVal = IntStream.range(0, values.length).filter(i -> idValues[i] == id).mapToDouble(i -> values[i]).min().getAsDouble();
        return minVal;
    }

    public static double getMaxValue(double[] values) {
        double maxVal = Arrays.stream(values).max().getAsDouble();
        return maxVal;
    }

    public static double getMaxValue(double[] values, int[] idValues, int id) {
        double maxVal = IntStream.range(0, values.length).filter(i -> idValues[i] == id).mapToDouble(i -> values[i]).max().getAsDouble();
        return maxVal;
    }

    public static double getMidValue(double[] yValues, double[] xValues) {
        double hh = (DataUtil.getMaxValue(yValues) + DataUtil.getMinValue(yValues)) / 2.0;
        double deltaUp = Double.MAX_VALUE;
        double deltaDown = Double.MAX_VALUE;
        double dUp, iUp;
        double dDown, iDown;
        dUp = iUp = dDown = iDown = 0;

        for (int i = 0; i < xValues.length; i++) {
            double dvar = yValues[i];
            double ivar = xValues[i];
            double ddvar = dvar - hh;

            if (ddvar >= 0 && ddvar < deltaUp) {
                deltaUp = ddvar;
                dUp = dvar;
                iUp = ivar;
            } else if (ddvar < 0 && -ddvar < deltaDown) {
                deltaDown = -ddvar;
                dDown = dvar;
                iDown = ivar;
            }
        }

        double mid;

        if (dUp == dDown) {
            mid = (iUp + iDown) / 2.0;
        } else {
            mid = ((hh - dDown) / (dUp - dDown) * (iUp - iDown)) + iDown;
        }

        return mid;
    }

    public static double getMidValue(double[] yValues, double[] xValues, int[] idValues, int id) {
        double hh = (DataUtil.getMaxValue(yValues) + DataUtil.getMinValue(yValues)) / 2.0;
        double deltaUp = Double.MAX_VALUE;
        double deltaDown = Double.MAX_VALUE;
        double dUp, iUp;
        double dDown, iDown;
        dUp = iUp = dDown = iDown = 0;

        for (int i = 0; i < xValues.length; i++) {
            if (idValues[i] != id) {
                continue;
            }
            double dvar = yValues[i];
            double ivar = xValues[i];
            double ddvar = dvar - hh;

            if (ddvar >= 0 && ddvar < deltaUp) {
                deltaUp = ddvar;
                dUp = dvar;
                iUp = ivar;
            } else if (ddvar < 0 && -ddvar < deltaDown) {
                deltaDown = -ddvar;
                dDown = dvar;
                iDown = ivar;
            }
        }

        double mid;

        if (dUp == dDown) {
            mid = (iUp + iDown) / 2.0;
        } else {
            mid = ((hh - dDown) / (dUp - dDown) * (iUp - iDown)) + iDown;
        }

        return mid;
    }

    public static double getMidValueZero(double[] yValues, double[] xValues, int[] idValues, int id) {
        double hh = DataUtil.getMaxValue(yValues) / 2.0;
        double deltaUp = Double.MAX_VALUE;
        double deltaDown = Double.MAX_VALUE;
        double dUp, iUp;
        double dDown, iDown;
        dUp = iUp = dDown = iDown = 0;

        for (int i = 0; i < xValues.length; i++) {
            if (idValues[i] != id) {
                continue;
            }
            double dvar = yValues[i];
            double ivar = xValues[i];
            double ddvar = dvar - hh;

            if (ddvar >= 0 && ddvar < deltaUp) {
                deltaUp = ddvar;
                dUp = dvar;
                iUp = ivar;
            } else if (ddvar < 0 && -ddvar < deltaDown) {
                deltaDown = -ddvar;
                dDown = dvar;
                iDown = ivar;
            }
        }

        double mid;

        if (dUp == dDown) {
            mid = (iUp + iDown) / 2.0;
        } else {
            mid = ((hh - dDown) / (dUp - dDown) * (iUp - iDown)) + iDown;
        }

        return mid;
    }

    public static double getMidValue1Side(double[] yValues, double[] xValues) {
        double hh = DataUtil.getMaxValue(yValues) / 2.0;
        double deltaMin = Double.MAX_VALUE;
        double iMid = 0;

        for (int i = 0; i < xValues.length; i++) {
            double dvar = yValues[i];
            double ivar = xValues[i];
            double ddvar = Math.abs(dvar - hh);

            if (ddvar < deltaMin) {
                deltaMin = ddvar;
                iMid = ivar;
            }
        }

        return iMid;
    }

    public static double getXAtMinValue(double[] yValues, double[] xValues) {
        double minVal = yValues[0];
        double curVal;
        double minX = 0;
        for (int i = 0; i < xValues.length; i++) {
            curVal = yValues[i];
            if (curVal < minVal) {
                minVal = curVal;
                minX = xValues[i];
            }
        }

        return minX;
    }

    public static double getYAtMaxX(double[] yValues, double[] xValues) {
        double maxVal = Double.NEGATIVE_INFINITY;
        double curVal;
        double maxY = 0;
        for (int i = 0; i < xValues.length; i++) {
            curVal = xValues[i];
            if (curVal > maxVal) {
                maxVal = curVal;
                maxY = yValues[i];
            }
        }
        return maxY;
    }

    public static double getYAtMinX(double[] yValues, double[] xValues) {
        double minVal = Double.MAX_VALUE;
        double curVal;
        double minY = 0;
        for (int i = 0; i < xValues.length; i++) {
            curVal = xValues[i];
            if (curVal < minVal) {
                minVal = curVal;
                minY = yValues[i];
            }
        }
        return minY;
    }

    public static double[] getCPMGInterpolation(double[] fixedXValues, double[] xDataPoints, double[] yDataPoints) {
        int lenOfFixedX = fixedXValues.length;
        int lenOfXDP = xDataPoints.length;
        double[] x = new double[lenOfXDP +1 ];
        double[] y = new double[lenOfXDP +1 ];
        if ((lenOfFixedX > 3) && (lenOfXDP > 3)) {
            double min = Double.MAX_VALUE;
            for (int i=0;i<xDataPoints.length;i++) {
                x[i] = xDataPoints[i];
                y[i] = yDataPoints[i];
                min = Math.min(min,yDataPoints[i]);
            }
            x[lenOfXDP] = fixedXValues[lenOfFixedX - 1];
            y[lenOfXDP] = min;
            var fitPars = fitPoly(x, y, 5);
            var poly = new PolynomialFunction(fitPars.toArray());
            double[] newArray = new double[lenOfFixedX];

            for (int i = 0; i < lenOfFixedX; i++) {
                newArray[i] = poly.value(fixedXValues[i]);
            }
            return newArray;
        } else {
            return null;
        }
    }

    public static RealVector fitPoly(double[] xValues, double[] yValues, int order) {
        int nRows = xValues.length;
        RealMatrix A = new Array2DRowRealMatrix(nRows, order);
        RealVector B = new ArrayRealVector(nRows);

        for (int i = 0; i < nRows; i++) {
            A.setEntry(i, 0, 1.0);

            for (int j = 1; j < order; j++) {
                A.setEntry(i, j, A.getEntry(i, j - 1) * (xValues[i]));
            }

            B.setEntry(i, yValues[i]);
        }

        SingularValueDecomposition svd = new SingularValueDecomposition(A);
        RealMatrix U = svd.getU();
        RealMatrix V = svd.getV();
        double[] svs = svd.getSingularValues();
        double coef;
        DecompositionSolver solver = svd.getSolver();
        RealVector X = solver.solve(B);
        return X;
    }

}
