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

/**
 *
 * @author brucejohnson
 */
public class Utilities {
    public static final double TWO_PI = 2.0 * Math.PI;

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
