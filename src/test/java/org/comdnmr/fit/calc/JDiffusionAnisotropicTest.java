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

public class JDiffusionAnisotropicTest {
    
    private final double Dxx = 4.4170*1e7; //44169507.013
    private final double Dyy = 4.5832*1e7; //45831876.669
    private final double Dzz = 6.0129*1e7; //60129330.314
    private final double[] resPars = {Dxx, Dyy, Dzz, 98.06/180.*Math.PI, 68.64/180.*Math.PI, 77.42/180.*Math.PI};//Dxx, Dyy, Dzz, alpha, beta, gamma
    private final double[] v = {0.688828228668, -0.0591882033522, -0.722504275402}; //normalized unit vector for residue 55 of 1P7F.pdb
    private final DiffusionType dType = ANISOTROPIC;
    
    double[][] D = {{resPars[0], 0.0, 0.0}, {0.0, resPars[1], 0.0}, {0.0, 0.0, resPars[2]}}; 
    double[][] VT = {{-0.9775, -0.0583, -0.2029}, {-0.1659, -0.3825, 0.9089}, {-0.1306, 0.9221, 0.3642}};

    
    public double[] genJ(RelaxFit relaxFit, RelaxEquations relaxObj, double[] vec, int modelNum, double[] modelPars) {
        double[] pars = new double[resPars.length + 1];
        double s2 = modelPars[0];
        switch (modelNum) {
            case 1:
                for (int p=0; p<resPars.length; p++) {
                    pars[p] = resPars[p];
                }
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
//        for (double par : pars) {
//            System.out.print(par + " ");
//        }
//        System.out.println();
        double[] J = relaxFit.getJDiffusion(pars, relaxObj, modelNum, vec, dType, D, VT);
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
        System.out.println("testGetJDiffModel0:");
        double[] fields = {400.0e6, 800e6};
        double[] modelPars = {1.0};
        
//        double[] RotDifJ0 = {1.1378239504681047E-9, 1.1378239504681047E-9}; //Model 1: S2 = 0.868
//        double[] RotDifJp = {6.681298484151205E-10, 3.003341664262763E-10}; //Model 1: S2 = 0.868
        
        double[] RotDifJ0 = {1.3107654670742526E-9, 1.3107654670742526E-9};
        double[] RotDifJp = {7.696810536343281E-10, 3.4598292374172923E-10};
        
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.makeRelaxObjs(fields, "H", "N");
        ArrayList<RelaxEquations> relaxObjs = relaxFit.getRelaxObjs();
        for (int r=0; r<relaxObjs.size(); r++) {
            System.out.println("Model 0: " + fields[r]/1.e6 + " MHz");
            RelaxEquations relaxObj = relaxObjs.get(r);
        
            double[] J = genJ(relaxFit, relaxObj, v, 1, modelPars); //relaxFit.getJDiffusion(resPars, relaxObj, 1, v, dType, D, VT, false);

            System.out.println("J = ");
            for (double val : J) {
                System.out.print(val + " ");
            }
            
            System.out.println("\nRotDif J0 = " + RotDifJ0[r] + " RotDif Jp = " + RotDifJp[r]);

            Assert.assertEquals(RotDifJ0[r] * 1e9, J[0] * 1e9, 1.0e-3);
            Assert.assertEquals(RotDifJp[r] * 1e10, J[1] * 1e10, 1.0e-3);
        }
        System.out.println();
    }
    
    @Test
    public void testRhoCalcs() {
        System.out.println("testRhoCalcs:");
        double[] fields = {400.0e6, 800e6};
        double[] modelPars = {1.0};
         
        double[] r1 = {3.014416, 1.852569};
        double[] r2 = {4.793736, 5.488800};
        double[] noe = {0.526066, 0.781323};
        
        double[] RotDifRhoExps = {2.273, 5.045};
//        double[] RotDifRhoPreds = {2.245, 4.837}; //Model 1: S2 = 0.868
        double[] RotDifRhoPreds = {2.271, 5.051};
        
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.makeRelaxObjs(fields, "H", "N");
        ArrayList<RelaxEquations> relaxObjs = relaxFit.getRelaxObjs();
        for (int r=0; r<relaxObjs.size(); r++) {
            System.out.println("Model 0: " + fields[r]/1.e6 + " MHz");
            RelaxEquations relaxObj = relaxObjs.get(r);
        
            double[] J = genJ(relaxFit, relaxObj, v, 1, modelPars); //relaxFit.getJDiffusion(resPars, relaxObj, modelNum, v, dType, D, VT, false);

            double rhoExp = relaxObj.calcRhoExp(r1[r], r2[r], noe[r], J);
            double rhoPred = relaxObj.calcRhoPred(J);

            double RotDifRhoExp = RotDifRhoExps[r];
            double RotDifRhoPred = RotDifRhoPreds[r];

            System.out.println("RhoExp, RotDif RhoExp = " + rhoExp + ", " + RotDifRhoExp);
            System.out.println("RhoPred, RotDif RhoPred = " + rhoPred + ", " + RotDifRhoPred);

            Assert.assertEquals(RotDifRhoExp, rhoExp, 1.0e-3);
            Assert.assertEquals(RotDifRhoPred, rhoPred, 1.0e-3);
        }
        System.out.println();
        
    }
    
    }
