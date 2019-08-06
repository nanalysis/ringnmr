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
 * @author Bruce Johnson
 */
public class SavitzkyGolay { 
    
    //quadratic/cubic Coefficient:
    final double M2L02R02[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 12, 17, 12, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    final double M2L03R03[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 3, 6, 7, 6, 3, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    final double M2L04R04[] = {0, 0, 0, 0, 0, 0, 0, 0, -21, 14, 39, 54, 59, 54, 39, 14, -21, 0, 0, 0, 0, 0, 0, 0, 0};
    final double M2L05R05[] = {0, 0, 0, 0, 0, 0, 0, -36, 9, 44, 69, 84, 89, 84, 69, 44, 9, -36, 0, 0, 0, 0, 0, 0, 0};
    final double M2L06R06[] = {0, 0, 0, 0, 0, 0, -11, 0, 9, 16, 21, 24, 25, 24, 21, 16, 9, 0, -11, 0, 0, 0, 0, 0, 0};
    final double M2L07R07[] = {0, 0, 0, 0, 0, -78, -13, 42, 87, 122, 147, 162, 167, 162, 147, 122, 87, 42, -13, -78, 0, 0, 0, 0, 0};
    final double M2L08R08[] = {0, 0, 0, 0, -21, -6, 7, 18, 27, 34, 39, 42, 43, 42, 39, 34, 27, 18, 7, -6, -21, 0, 0, 0, 0};
    final double M2L09R09[] = {0, 0, 0, -136, -51, 24, 89, 144, 189, 224, 249, 264, 269, 264, 249, 224, 189, 144, 89, 24, -51, -136, 0, 0, 0};
    final double M2L10R10[] = {0, 0, -171, -76, 9, 84, 149, 204, 249, 284, 309, 324, 329, 324, 309, 284, 249, 204, 149, 84, 9, -76, -171, 0, 0};
    final double M2L11R11[] = {0, -42, -21, -2, 15, 30, 43, 54, 63, 70, 75, 78, 79, 78, 75, 70, 63, 54, 43, 30, 15, -2, -21, -42, 0};
    final double M2L12R12[] = {-253, -138, -33, 62, 147, 222, 287, 322, 387, 422, 447, 462, 467, 462, 447, 422, 387, 322, 287, 222, 147, 62, -33, -138, -253};

    // quartic/quintic Coefficient:
    final double M4L03R03[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 5, -30, 75, 131, 75, -30, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    final double M4L04R04[] = {0, 0, 0, 0, 0, 0, 0, 0, 15, -55, 30, 135, 179, 135, 30, -55, 15, 0, 0, 0, 0, 0, 0, 0, 0};
    final double M4L05R05[] = {0, 0, 0, 0, 0, 0, 0, 18, -45, -10, 60, 120, 143, 120, 60, -10, -45, 18, 0, 0, 0, 0, 0, 0, 0};
    final double M4L06R06[] = {0, 0, 0, 0, 0, 0, 110, -198, -160, 110, 390, 600, 677, 600, 390, 110, -160, -198, 110, 0, 0, 0, 0, 0, 0};
    final double M4L07R07[] = {0, 0, 0, 0, 0, 2145, -2860, -2937, -165, 3755, 7500, 10125, 11053, 10125, 7500, 3755, -165, -2937, -2860, 2145, 0, 0, 0, 0, 0};
    final double M4L08R08[] = {0, 0, 0, 0, 195, -195, -260, -117, 135, 415, 660, 825, 883, 825, 660, 415, 135, -117, -260, -195, 195, 0, 0, 0, 0};
    final double M4L09R09[] = {0, 0, 0, 340, -255, -420, -290, 18, 405, 790, 1110, 1320, 1393, 1320, 1110, 790, 405, 18, -290, -420, -255, 340, 0, 0, 0};
    final double M4L10R10[] = {0, 0, 11628, -6460, -13005, -11220, -3940, 6378, 17655, 28190, 36660, 42120, 44003, 42120, 36660, 28190, 17655, 6378, -3940, -11220, -13005, -6460, 11628, 0, 0};
    final double M4L11R11[] = {0, 285, -114, -285, -285, -165, 30, 261, 495, 705, 870, 975, 1011, 975, 870, 705, 495, 261, 30, -165, -285, -285, -114, 285, 0};
    final double M4L12R12[] = {1265, -345, -1122, -1255, -915, -255, 590, 1503, 2385, 3155, 3750, 4125, 4253, 4125, 3750, 3155, 2385, 1503, 590, -255, -915, -1255, -1122, -345, 1265};

    double norm = 999999.9; //Normalization factor for smoothing coefficients.
        
    /*
     Implementation of the Savitzky-Golay smoothing filter.
     *
     The function input involves the number of points to use (eg.5,7,9,...,25)
     and order of polynomial (range from 2 to 5). Please note that for either a
     quadratic (2nd order) or cubic(3rd order) function, the set of coefficients
     is the same; so is that of a quartic function (4th order) and a quintic
     function (5th order).
     */
    public double[] runningSavitzkyGolaySmoothing(double[] vec, int j1, int j2, int order, int smoothSize, double[] X) throws Exception {
        /* Get - Number of points used      smoothSize
         - Order of polynomial        order
         - region of smoothing        j1(lower bound) amd j2(upper bound)
         - 1-D spectrum data          vec */

        if ((smoothSize > 25) || (smoothSize < 5) || (smoothSize / 2 * 2 == smoothSize)) {
            throw new Exception("Smooth-size should be one of 5,7,9,...,25");
        }

        double[] smoothingCoefficient;

        //smoothingCoefficient=M4L12R12;
        if (order == 2 || order == 3) //quadratic/cubic function:
        {
            switch (smoothSize) {
                case 5: {
                    smoothingCoefficient = M2L02R02;
                    norm = 35;
                    break;
                }
                case 7: {
                    smoothingCoefficient = M2L03R03;
                    norm = 21;
                    break;
                }
                case 9: {
                    smoothingCoefficient = M2L04R04;
                    norm = 231;
                    break;
                }
                case 11: {
                    smoothingCoefficient = M2L05R05;
                    norm = 429;
                    break;
                }
                case 13: {
                    smoothingCoefficient = M2L06R06;
                    norm = 143;
                    break;
                }
                case 15: {
                    smoothingCoefficient = M2L07R07;
                    norm = 1105;
                    break;
                }
                case 17: {
                    smoothingCoefficient = M2L08R08;
                    norm = 323;
                    break;
                }
                case 19: {
                    smoothingCoefficient = M2L09R09;
                    norm = 2261;
                    break;
                }
                case 21: {
                    smoothingCoefficient = M2L10R10;
                    norm = 3059;
                    break;
                }
                case 23: {
                    smoothingCoefficient = M2L11R11;
                    norm = 805;
                    break;
                }
                case 25: {
                    smoothingCoefficient = M2L12R12;
                    norm = 5175;
                    break;
                }
                default: {
                    smoothingCoefficient = M2L02R02;
                    norm = 35;
                    break;
                }
            }
        } else if (order == 4 || order == 5) // quartic/quintic function
        {
            switch (smoothSize) {
                case 7: {
                    smoothingCoefficient = M4L03R03;
                    norm = 231;
                    break;
                }
                case 9: {
                    smoothingCoefficient = M4L04R04;
                    norm = 429;
                    break;
                }
                case 11: {
                    smoothingCoefficient = M4L05R05;
                    norm = 429;
                    break;
                }
                case 13: {
                    smoothingCoefficient = M4L06R06;
                    norm = 2431;
                    break;
                }
                case 15: {
                    smoothingCoefficient = M4L07R07;
                    norm = 46189;
                    break;
                }
                case 17: {
                    smoothingCoefficient = M4L08R08;
                    norm = 4199;
                    break;
                }
                case 19: {
                    smoothingCoefficient = M4L09R09;
                    norm = 7429;
                    break;
                }
                case 21: {
                    smoothingCoefficient = M4L10R10;
                    norm = 260015;
                    break;
                }
                case 23: {
                    smoothingCoefficient = M4L11R11;
                    norm = 6555;
                    break;
                }
                case 25: {
                    smoothingCoefficient = M4L12R12;
                    norm = 30015;
                    break;
                }
                default: {
                    smoothingCoefficient = M4L03R03;
                    norm = 231;
                    break;
                }
            }
        } else {
            throw new Exception("Polynomial functions of degree 5 or higher haven't been implemented yet.");
        }
        for (int i = 0; i < smoothingCoefficient.length; i++) {
            smoothingCoefficient[i] = smoothingCoefficient[i] / (double) norm;
        }

        for (int j = j1; j < j2; j++) {
            int k1 = j - (smoothSize / 2);
            int k2 = j + (smoothSize / 2);
            double sum = 0.0;
            for (int k = k1; k <= k2; k++) {
                int position = k;
                if (k < 0) {
                    position = 0;
                }
                if (k >= vec.length) {
                    position = vec.length - 1;
                }
                sum += vec[position] * smoothingCoefficient[k - j + 12];
            }
//            X.set(j, sum);
            X[j] = sum;
        }
        return X;
    }
}
