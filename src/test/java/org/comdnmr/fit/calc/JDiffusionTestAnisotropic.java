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
package org.comdnmr.fit.calc;

import java.util.ArrayList;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealMatrixFormat;
import org.comdnmr.data.RelaxFit;
import org.comdnmr.data.RelaxFit.DiffusionType;
import static org.comdnmr.data.RelaxFit.DiffusionType.ANISOTROPIC;
import org.comdnmr.modelfree.RelaxEquations;
import org.junit.Assert;
import org.junit.Test;

public class JDiffusionTestAnisotropic {

    // private final double[] resPars = {4.417e7, 4.5832e7, 6.0129e7, 98.06 / 180. * Math.PI, 68.64 / 180. * Math.PI, 77.42 / 180. * Math.PI};//Dxx, Dyy, Dzz, alpha, beta, gamma
    // private final double[] v = {0.688828228668, -0.0591882033522, -0.722504275402}; //normalized unit vector for residue 55 of 1P7F.pdb
    private final DiffusionType dType = ANISOTROPIC;
    // line 745674 of output
    private final double[] resPars = {4.4183e7, 4.5887e7, 6.0036e7, 98.17 / 180. * Math.PI, 68.88 / 180. * Math.PI, 76.68 / 180. * Math.PI};//Dxx, Dyy, Dzz, alpha, beta, gamma

    double[][] D = {{resPars[0], 0.0, 0.0}, {0.0, resPars[1], 0.0}, {0.0, 0.0, resPars[2]}};
    //double[][] VT = {{-0.9775, -0.0583, -0.2029}, {-0.1659, -0.3825, 0.9089}, {-0.1306, 0.9221, 0.3642}};
    double[][] VT = {{-0.9750, -0.0561, -0.2149}, {-0.1783, -0.3797, 0.9078}, {-0.1325, 0.9234, 0.3602}};
//3       TYR     -0.840520439988 0.251236972341  -0.480005597563

    //private final double[] v = {-0.840520439988, 0.251236972341, -0.480005597563};
   // LYS 4
    private final double[] v = {0.981, -0.192, 0.001};

    public double[] genJ(RelaxFit relaxFit, RelaxEquations relaxObj, double[] vec, int modelNum, double[] modelPars) {
        double[] pars = new double[resPars.length + 1];
        double s2 = modelPars[0];
        switch (modelNum) {
            case 1:
                pars[6] = s2;
                break;
            case 2:
                pars = new double[resPars.length + 2];
                System.arraycopy(resPars, 0, pars, 0, resPars.length);
                double tau = modelPars[1];
                pars[6] = s2;
                pars[7] = tau;
                break;
            case 5:
                pars = new double[resPars.length + 3];
                System.arraycopy(resPars, 0, pars, 0, resPars.length);
                tau = modelPars[1];
                double sf2 = modelPars[2];
                pars[6] = s2;
                pars[7] = tau;
                pars[8] = sf2;
                break;
            case 6:
                pars = new double[resPars.length + 4];
                System.arraycopy(resPars, 0, pars, 0, resPars.length);
                tau = modelPars[1];
                sf2 = modelPars[2];
                double tauS = modelPars[3];
                pars[6] = s2;
                pars[7] = tau;
                pars[8] = sf2;
                pars[9] = tauS;
                break;
            default:
                break;
        }
        double[] J = relaxFit.getJDiffusion(pars, relaxObj, modelNum, vec, dType, D, VT, false);
        return J;
    }

    @Test
    public void testGetRotVT() {
        RelaxFit relaxFit = new RelaxFit();
        Rotation rot = relaxFit.getDRotation(resPars, dType);
        double[][] VTres = new Array2DRowRealMatrix(rot.getMatrix()).transpose().getData();
        RealMatrixFormat formatter = new RealMatrixFormat();
        RealMatrix VTR = MatrixUtils.createRealMatrix(VTres);
        RealMatrix validVTR = MatrixUtils.createRealMatrix(VT);
        System.out.println("result VT = " + formatter.format(VTR));
        System.out.println("valid VT = " + formatter.format(validVTR));
        double norm = VTR.subtract(validVTR).getNorm();
        Assert.assertEquals(0, norm, 1.0e-3);

    }

    @Test
    public void testGetJDiffusionModel0() {

        double RotDifJ0 = 1.1378239504681047E-9;
        double RotDifJp = 6.681298484151205E-10;

        RelaxFit relaxFit = new RelaxFit();
        double[] fields = {400.0e6};
        relaxFit.makeRelaxObjs(fields, "H", "N");
        ArrayList<RelaxEquations> relaxObjs = relaxFit.getRelaxObjs();
        RelaxEquations relaxObj = relaxObjs.get(0);

        double[] modelPars = {1.0};
        double[] J = genJ(relaxFit, relaxObj, v, 1, modelPars); //relaxFit.getJDiffusion(resPars, relaxObj, 1, v, dType, D, VT, false);
        double[] validJ = {1.3107950303472567E-9, 7.71213958846025E-10, 1.5842901000533268E-11, 1.915852662686483E-11, 2.3627355958876023E-11};

        System.out.println("Model 0:");
        System.out.println("J = ");
        for (double val : J) {
            System.out.print(val + " ");
        }
        System.out.println("\nvalid J = ");
        for (double val : validJ) {
            System.out.print(val + " ");
        }
        System.out.println("\nRotDif J0 = " + RotDifJ0 + " RotDif Jp = " + RotDifJp);

        double[] diff = new double[J.length];
        Assert.assertArrayEquals(validJ, J, 1.0e-11);
        for (int i = 0; i < diff.length; i++) {
            diff[i] = validJ[i] - J[i];
            Assert.assertEquals(0, diff[i], 1.0e-3);
        }

//        Assert.assertEquals(1.73e-10, J[0] - RotDifJ0, 1.0e-3);
//        Assert.assertEquals(1.03e-10, J[1] - RotDifJp, 1.0e-3);
        Assert.assertEquals(J[0] * 1.0e9, RotDifJ0 * 1.0e9, 1.0e-0);
        Assert.assertEquals(J[1] * 1.0e9, RotDifJp * 1.0e9, 1.0e-0);
    }

    @Test
    public void testRhoCalcs() {
        double r1 = 3.014416;
        double r2 = 4.793736;
        double noe = 0.526066;

        RelaxFit relaxFit = new RelaxFit();
        double[] fields = {400.0e6};
        relaxFit.makeRelaxObjs(fields, "H", "N");
        ArrayList<RelaxEquations> relaxObjs = relaxFit.getRelaxObjs();
        RelaxEquations relaxObj = relaxObjs.get(0);

        double[] modelPars = {1.0};
        double[] J = genJ(relaxFit, relaxObj, v, 1, modelPars); //relaxFit.getJDiffusion(resPars, relaxObj, modelNum, v, dType, D, VT, false);

        double rhoExp = relaxObj.calcRhoExp(r1, r2, noe, J);
        double rhoPred = relaxObj.calcRhoPred(J);

        double validRhoExp = 2.273;
        double validRhoPred = 2.266;

        double RotDifRhoExp = 2.273;
        double RotDifRhoPred = 2.245;

        System.out.println("Model 0: ");
        System.out.println("RhoExp, valid RhoExp, RotDif RhoExp = " + rhoExp + ", " + validRhoExp + ", " + RotDifRhoExp);
        System.out.println("RhoPred, valid RhoPred, RotDif RhoPred = " + rhoPred + ", " + validRhoPred + ", " + RotDifRhoPred);

        Assert.assertEquals(0, validRhoExp - rhoExp, 1.0e-3);
        Assert.assertEquals(0, validRhoPred - rhoPred, 1.0e-3);
        // Assert.assertEquals(-1.42e-4, RotDifRhoExp-rhoExp, 1.0e-3);
        Assert.assertEquals(RotDifRhoExp, rhoExp, 1.0e-3);
        // Assert.assertEquals(-0.021, RotDifRhoPred - rhoPred, 1.0e-3);
        Assert.assertEquals(RotDifRhoPred, rhoPred, 3.0e-2);

    }

    @Test
    public void testGetJDiffusionModel1() {

        double r1 = 3.031857;
        double r2 = 4.987638;
        double noe = 0.566376;

        RelaxFit relaxFit = new RelaxFit();
        double[] fields = {400.0e6};
        relaxFit.makeRelaxObjs(fields, "H", "N");
        ArrayList<RelaxEquations> relaxObjs = relaxFit.getRelaxObjs();
        RelaxEquations relaxObj = relaxObjs.get(0);

        double[] vec = {0.693167203159, -0.720730552917, 0.00816691844665}; //res 26 of 1P7F.pdb

        double[] modelPars = {0.858};
        double[] J = genJ(relaxFit, relaxObj, vec, 1, modelPars); //relaxFit.getJDiffusion(resPars, relaxObj, 1, v, dType, D, VT, false);

        double rhoExp = relaxObj.calcRhoExp(r1, r2, noe, J);
        double rhoPred = relaxObj.calcRhoPred(J);

        double validRhoExp = 2.381;
        double validRhoPred = 2.363;

        double RotDifRhoExp = 2.381;
        double RotDifRhoPred = 2.368;

        System.out.println("Model 1:");
        System.out.println("RhoExp, valid RhoExp, RotDif RhoExp = " + rhoExp + ", " + validRhoExp + ", " + RotDifRhoExp);
        System.out.println("RhoPred, valid RhoPred, RotDif RhoPred = " + rhoPred + ", " + validRhoPred + ", " + RotDifRhoPred);

        Assert.assertEquals(0, validRhoExp - rhoExp, 1.0e-3);
        Assert.assertEquals(0, validRhoPred - rhoPred, 1.0e-3);
        Assert.assertEquals(-1.42e-4, RotDifRhoExp - rhoExp, 1.0e-3);
        Assert.assertEquals(5.36e-3, RotDifRhoPred - rhoPred, 1.0e-3);
    }

    @Test
    public void testGetJDiffusionModel2() {

        double r1 = 3.165870;
        double r2 = 4.990445;
        double noe = 0.569069;

        RelaxFit relaxFit = new RelaxFit();
        double[] fields = {400.0e6};
        relaxFit.makeRelaxObjs(fields, "H", "N");
        ArrayList<RelaxEquations> relaxObjs = relaxFit.getRelaxObjs();
        RelaxEquations relaxObj = relaxObjs.get(0);

        double[] vec = {-0.978994554933, -0.186815436447, 0.0816679503594}; //res 52 of 1P7F.pdb

        double[] modelPars = {0.806, 1.25e-11};
        double[] J = genJ(relaxFit, relaxObj, vec, 2, modelPars); //relaxFit.getJDiffusion(resPars, relaxObj, 1, v, dType, D, VT, false);

        double rhoExp = relaxObj.calcRhoExp(r1, r2, noe, J);
        double rhoPred = relaxObj.calcRhoPred(J);

        double validRhoExp = 2.235;
        double validRhoPred = 2.208;

        double RotDifRhoExp = 2.235;
        double RotDifRhoPred = 2.212;

        System.out.println("Model 2:");
        System.out.println("RhoExp, valid RhoExp, RotDif RhoExp = " + rhoExp + ", " + validRhoExp + ", " + RotDifRhoExp);
        System.out.println("RhoPred, valid RhoPred, RotDif RhoPred = " + rhoPred + ", " + validRhoPred + ", " + RotDifRhoPred);

        Assert.assertEquals(0, validRhoExp - rhoExp, 1.0e-3);
        Assert.assertEquals(0, validRhoPred - rhoPred, 1.0e-3);
        Assert.assertEquals(-1.42e-4, RotDifRhoExp - rhoExp, 1.0e-3);
        Assert.assertEquals(4.16e-3, RotDifRhoPred - rhoPred, 1.0e-3);
    }

    @Test
    public void testGetJDiffusionModel5() {

        double r1 = 3.014416;
        double r2 = 4.793736;
        double noe = 0.526066;

        RelaxFit relaxFit = new RelaxFit();
        double[] fields = {400.0e6};
        relaxFit.makeRelaxObjs(fields, "H", "N");
        ArrayList<RelaxEquations> relaxObjs = relaxFit.getRelaxObjs();
        RelaxEquations relaxObj = relaxObjs.get(0);

        double[] modelPars = {0.822, 2.70e-10, 0.845};
        double[] J = genJ(relaxFit, relaxObj, v, 5, modelPars); //relaxFit.getJDiffusion(resPars, relaxObj, 1, v, dType, D, VT, false);

        double rhoExp = relaxObj.calcRhoExp(r1, r2, noe, J);
        double rhoPred = relaxObj.calcRhoPred(J);

        double validRhoExp = 2.273;
        double validRhoPred = 2.263;

        double RotDifRhoExp = 2.273;
        double RotDifRhoPred = 2.245;

        System.out.println("Model 5:");
        System.out.println("RhoExp, valid RhoExp, RotDif RhoExp = " + rhoExp + ", " + validRhoExp + ", " + RotDifRhoExp);
        System.out.println("RhoPred, valid RhoPred, RotDif RhoPred = " + rhoPred + ", " + validRhoPred + ", " + RotDifRhoPred);

        Assert.assertEquals(0, validRhoExp - rhoExp, 1.0e-3);
        Assert.assertEquals(0, validRhoPred - rhoPred, 1.0e-3);
        Assert.assertEquals(-1.42e-4, RotDifRhoExp - rhoExp, 1.0e-3);
        Assert.assertEquals(-0.018, RotDifRhoPred - rhoPred, 1.0e-3);
    }

    @Test
    public void testGetJDiffusionModel6() {

        double r1 = 2.996007;
        double r2 = 4.966411;
        double noe = 0.490225;

        RelaxFit relaxFit = new RelaxFit();
        double[] fields = {400.0e6};
        relaxFit.makeRelaxObjs(fields, "H", "N");
        ArrayList<RelaxEquations> relaxObjs = relaxFit.getRelaxObjs();
        RelaxEquations relaxObj = relaxObjs.get(0);

        double[] vec = {-0.640451991371, -0.275353563169, 0.716939092252}; //residue 56 of 1P7F.pdb

        double[] modelPars = {0.831, 8.95e-9, 0.858, 1.16e-11};
        double[] J = genJ(relaxFit, relaxObj, vec, 6, modelPars); //relaxFit.getJDiffusion(resPars, relaxObj, 1, v, dType, D, VT, false);

        double rhoExp = relaxObj.calcRhoExp(r1, r2, noe, J);
        double rhoPred = relaxObj.calcRhoPred(J);

        double validRhoExp = 2.425;
        double validRhoPred = 2.177;

        double RotDifRhoExp = 2.425;
        double RotDifRhoPred = 2.179;

        System.out.println("Model 6:");
        System.out.println("RhoExp, valid RhoExp, RotDif RhoExp = " + rhoExp + ", " + validRhoExp + ", " + RotDifRhoExp);
        System.out.println("RhoPred, valid RhoPred, RotDif RhoPred = " + rhoPred + ", " + validRhoPred + ", " + RotDifRhoPred);

        Assert.assertEquals(0, validRhoExp - rhoExp, 1.0e-3);
        Assert.assertEquals(0, validRhoPred - rhoPred, 1.0e-3);
        Assert.assertEquals(-1.42e-4, RotDifRhoExp - rhoExp, 1.0e-3);
        Assert.assertEquals(2.16e-3, RotDifRhoPred - rhoPred, 1.0e-3);

    }
}
