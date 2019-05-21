package org.comdnmr.fit.calc;

import java.util.Arrays;

/**
 *
 * @author brucejohnson
 */
public class Utilities {

    public static double[][] copy2DArray(double[][] inData) {
        double[][] outData = new double[inData.length][];
        for (int i = 0; i < inData.length; i++) {
            outData[i] = inData[i].clone();
        }
        return outData;
    }

    /**
     * Purpose - This method scales (or unscales) a "value" between [0,1] using
     * the maximum and minimum values specified in "bounds"
     *
     * @param value - the value that is to be scaled between 0 and 1
     * @param bounds - an array of double values whose first index is the
     * maximum bound and the second index is the minimum bound.
     * @param scaling - parameter to indicate whether we are scaling or
     * unscaling the value provided. (true if scaling; false otherwise)
     *
     * @return A scaled value between 0 and 1, or null.
     */
    public static Double scale(double value, double[] bounds, boolean scaling) {
        if (bounds.length == 2) {
            double maxValue = bounds[0], minValue = bounds[1];
            if (scaling) {
                // scaling
                return (value - minValue) / (maxValue - minValue);
            } else { 
                // unscaling
                return ((value * maxValue) - (value * minValue)) + minValue;
            }
        } else {
            return null;
        }
    }
    
    
}
