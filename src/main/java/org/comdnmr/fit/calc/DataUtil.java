package org.comdnmr.fit.calc;

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
}
