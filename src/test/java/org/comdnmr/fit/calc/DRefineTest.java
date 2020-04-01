/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.optim.PointValuePair;
import org.comdnmr.data.Fitter;
import org.comdnmr.data.RelaxFit;
import org.comdnmr.data.RelaxFit.DiffusionType;
import static org.comdnmr.data.RelaxFit.DiffusionType.ANISOTROPIC;
import org.junit.Test;
import org.junit.Assert;

/**
 *
 * @author Martha
 */
public class DRefineTest {

    //Values for 1P7F.pdb residue 55
    private final double[][] xVals = {{0, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 2, 0}, //iField, iRes, eType, iRes1
                                        {1, 0, 0, 0}, {1, 0, 1, 0}, {1, 0, 2, 0},
                                        {2, 0, 0, 0}, {2, 0, 1, 0}, {2, 0, 2, 0},
                                        {0, 1, 0, 1}, {0, 1, 1, 1}, {0, 1, 2, 1},
                                        {1, 1, 0, 1}, {1, 1, 1, 1}, {1, 1, 2, 1},
                                        {2, 1, 0, 1}, {2, 1, 1, 1}, {2, 1, 2, 1},
                                        {0, 2, 0, 2}, {0, 2, 1, 2}, {0, 2, 2, 2},
                                        {1, 2, 0, 2}, {1, 2, 1, 2}, {1, 2, 2, 2},
                                        {2, 2, 0, 2}, {2, 2, 1, 2}, {2, 2, 2, 2}};
    private final double[] yVals = {2.984570, 4.713244, 0.515370, 2.276911, 4.778002, 0.695735, 1.853511, 5.417860, 0.763318, //res 3 
                                    2.932273, 4.964714, 0.592637, 2.187923, 5.240804, 0.738716, 1.785262, 6.078894, 0.793904, //res 24
                                    3.014416, 4.793736, 0.526066, 2.323677, 5.001800, 0.715022, 1.852569, 5.4888, 0.781323}; //res 55

    private final double[] errVals = {0.039921, 0.102811, 0.006740, 0.024362, 0.042812, 0.009900, 0.033970, 0.056956, 0.009744,
                                        0.032936, 0.064766, 0.004377, 0.048618, 0.052214, 0.006000, 0.020917, 0.158875, 0.006845, 
                                        0.028600, 0.061819, 0.006449, 0.058305, 0.034659, 0.009300, 0.019226, 0.024974, 0.009085};
    private final double[][] coords = {{-0.840520439988, 0.251236972341, -0.480005597563}, //res 3
                                        {0.423475119964, -0.834705176219, -0.352045581657}, //res 24
                                        {0.688828228668, -0.0591882033522, -0.722504275402}}; //normalized unit vector for residue 55 of 1P7F.pdb

    private final int[] residueModels = {1, 1, 1};
    private final DiffusionType diffusionType = ANISOTROPIC;

    private final double tauC = 3.3e-9;
    private final double Di = 1.0 / (6.0 * tauC);
    private final double[] guesses = {0.75 * Di, Di, 1.25 * Di, 0, 0, 0};
    private final double[] fields = {400.0e6, 600e6, 800.0e6};

    @Test
    public void testValueDMat1() {
        System.out.println("testValueDMat1");
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.makeRelaxObjs(fields, "H", "N");

        relaxFit.setXYE(xVals, yVals, errVals);
        relaxFit.setCoords(coords);
        relaxFit.setDiffusionType(diffusionType);
        relaxFit.setResidueModels(residueModels, true);
        double[] pars = {4.4170 * 1e7, 4.5832 * 1e7, 6.0129 * 1e7, Math.toRadians(98.06),
            Math.toRadians(68.64), Math.toRadians(77.42)};

        double value = relaxFit.calcD(fields[0], pars);
        System.out.println("value " + value);
        Assert.assertEquals(0.0, value, 4.0e-2);
    }
    
    @Test
    public void testValueDMat() {
        double[] guesses = {4.4170 * 1e7, 4.5832 * 1e7, 6.0129 * 1e7, Math.toRadians(98.06),
            Math.toRadians(68.64), Math.toRadians(77.42)};
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.makeRelaxObjs(fields, "H", "N");

        relaxFit.setXYE(xVals, yVals, errVals);
        relaxFit.setCoords(coords);
        relaxFit.setDiffusionType(diffusionType);
        relaxFit.setResidueModels(residueModels, true);

        PointValuePair fitResult = relaxFit.fitD(fields[0], guesses);
        Fitter fitter = relaxFit.getFitter();
        double[] fitPars = relaxFit.getPars();
  //      relaxFit.calcParErrs(fitResult, fitter, fitPars, true, true, true);
//        double[] fitErrs = relaxFit.getParErrs();
        double fitRMS = fitResult.getValue();

        double[][][] rotResults = relaxFit.rotateD(fitPars);
        double[][] D = rotResults[0];
        double[][] VT = rotResults[1];

        System.out.println("Fit RMS = " + fitRMS);
        System.out.println("Fit Pars: ");
        for (double par : fitPars) {
            System.out.print(par + " ");
        }
        System.out.println("\nFitPars Scaled: ");
        for (int i = 0; i < fitPars.length; i++) {
            double par = fitPars[i];
            if (i < 3) {
                par /= 1.0e7;
            } else if (i >= 3 && i <= 6) {
                par *= 180.0 / Math.PI;
            }
            System.out.print(par + " ");
        }
        double[] rotDifPars = {4.4170, 4.5832, 6.0129, 98.06, 68.64, 77.42};
        System.out.println("\nRotDif Pars: ");
        for (double par : rotDifPars) {
            System.out.print(par + " ");
        }
//        System.out.println("\nFitErrs Scaled: ");
//        for (int i = 0; i < fitErrs.length; i++) {
//            double par = fitErrs[i];
//            if (i < 3) {
//                par /= 1.0e7;
//            } else if (i >= 3 && i <= 6) {
//                par *= 180.0 / Math.PI;
//            }
//            System.out.print(par + " ");
//        }
        System.out.println();
        System.out.println("Fit D = " + new Array2DRowRealMatrix(D).toString());
        System.out.println("Fit VT = " + new Array2DRowRealMatrix(VT).toString());
    }
}
