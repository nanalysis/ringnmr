package org.comdnmr.data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import org.apache.commons.math3.geometry.euclidean.threed.NotARotationMatrixException;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.RotationConvention;
import org.apache.commons.math3.geometry.euclidean.threed.RotationOrder;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.optim.PointValuePair;
import static org.comdnmr.data.RelaxFit.ExptType.NOE;
import static org.comdnmr.data.RelaxFit.ExptType.R1;
import static org.comdnmr.data.RelaxFit.ExptType.R2;
import static org.comdnmr.data.RelaxFit.DiffusionType.ANISOTROPIC;
import static org.comdnmr.data.RelaxFit.DiffusionType.OBLATE;
import static org.comdnmr.data.RelaxFit.DiffusionType.PROLATE;
import org.comdnmr.modelfree.RelaxEquations;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import static org.ejml.dense.row.CommonOps_DDRM.transpose;

/**
 *
 * @author brucejohnson
 */
public class RelaxFit {

    boolean reportFitness = true;
    int reportAt = 10;
    long startTime = 0;
    double B0;
    double[][] xValues;
    double[] yValues;
    double[] errValues;
    double[] bestPars;
    double[] parErrs;
    double bestAIC;
    double bestChiSq;
    ArrayList<RelaxEquations> relaxObjs = new ArrayList<>();
    Fitter bestFitter;
    int[] residueModels;
    int[] resParStart;
    int[] nParsPerModel = {0, 1, 2, 0, 0, 3, 4};
    static ExptType[] expTypeArray = {R1, R2, NOE};
    double[] tauMBounds = {3.8e-9, 4.5e-9};
    double[][] coords;
    DiffusionType diffusionType;
    
    public void setXYE(double[][] xValues, double[] yValues, double[] errValues) {
        this.xValues = xValues;
        this.yValues = yValues;
        this.errValues = errValues;
    }
    
    public void setCoords(double[][] coords) {
        this.coords = coords;
    }
    
    public void setDiffusionType(DiffusionType type) {
        this.diffusionType = type;
    }

    public enum ExptType {
        R1, R2, NOE;
    }
    
    public enum DiffusionType {
        ISOTROPIC, PROLATE, OBLATE, ANISOTROPIC;
    }

    public void makeRelaxObjs(double[] fields, String elem1, String elem2) {
        for (int i = 0; i < fields.length; i++) {
            relaxObjs.add(new RelaxEquations(fields[i], elem1, elem2));
        }
    }

    public double[] getJ(double[] pars, RelaxEquations relaxObj, int modelNum) {
        double tauM = pars[0];//4.5e-9;
        double s2 = pars[1];
        double[] J = new double[5];
        switch (modelNum) {
            case 1:
                J = relaxObj.getJModelFree(tauM, s2);
                break;
            case 2:
                double tau = pars[2];
                J = relaxObj.getJModelFree(tau, tauM, s2);
                break;
            case 5:
                tau = pars[2];
                double sf2 = pars[3];
                J = relaxObj.getJModelFree(tau, tauM, s2, sf2);
                break;
            case 6:
                tau = pars[2];
                sf2 = pars[3];
                double tauS = pars[4];
//                System.out.println("tau, sf2, tauS = " + tau + " " + sf2 + " " + tauS);
                J = relaxObj.getJModelFree(tau, tauM, tauS, s2, sf2);
                break;
            default:
                break;
        }
//        for (double Jval : J) {
//            System.out.println("J: " + Jval);
//        }
        return J;
    }
    
    public double[] getJDiffusion(double[] pars, RelaxEquations relaxObj, int modelNum, double[] v, DiffusionType dType, double[][] D) {
        int nEqlDiffPars = 0;
        int nNullAngles = 0;
        if (dType == PROLATE || dType == OBLATE) { //Dxx = Dyy or Dyy = Dzz
            nEqlDiffPars = 1;
            nNullAngles = 1;
        }
        
        //fixme only do EigenDecomposition of D instead of rotating D then EigenDecomposition??
        if (D == null) {
            D = rotateD(pars, dType, null);
        }
        
        EigenDecomposition eig = new EigenDecomposition(new Array2DRowRealMatrix(D));
        D = eig.getD().getData();
        double[][] VT = eig.getVT().getData();
    
        double s2 = pars[6-nEqlDiffPars-nNullAngles];
        double[] J = new double[5];
        switch (modelNum) {
            case 1:
                J = relaxObj.getJDiffusion(dType, D, VT, v, s2, null, null, null);
                break;
            case 2:
                double tau = pars[7-nEqlDiffPars-nNullAngles];
                J = relaxObj.getJDiffusion(dType, D, VT, v, s2, tau, null, null);
                break;
            case 5:
                tau = pars[7-nEqlDiffPars-nNullAngles];
                double sf2 = pars[8-nEqlDiffPars-nNullAngles];
                J = relaxObj.getJDiffusion(dType, D, VT, v, s2, tau, sf2, null);
                break;
            case 6:
                tau = pars[7-nEqlDiffPars-nNullAngles];
                sf2 = pars[8-nEqlDiffPars-nNullAngles];
                double tauS = pars[9-nEqlDiffPars-nNullAngles];
//                System.out.println("tau, sf2, tauS = " + tau + " " + sf2 + " " + tauS);
                J = relaxObj.getJDiffusion(dType, D, VT, v, s2, tau, sf2, tauS);
                break;
            default:
                break;
        }
//        for (double Jval : J) {
//            System.out.println("J: " + Jval);
//        }
        return J;
    }
    
    public double[][] rotateD(double[] pars, DiffusionType dType, double[][] initD) {
        int nEqlDiffPars = 0;
        double Dxx = pars[0];
        double Dyy = pars[1];
        double Dzz = 0.0;
        double alpha = 0.0;
        double beta = 0.0;
        double gamma = 0.0;
        switch (dType) {
            case PROLATE:
                //Dxx = Dyy
                nEqlDiffPars = 1;
                Dyy = pars[1-nEqlDiffPars];
                Dzz = pars[2-nEqlDiffPars];
                alpha = pars[3-nEqlDiffPars];
                beta = pars[4-nEqlDiffPars];
                gamma = 0;
                break;
            case OBLATE:
                //Dyy = Dzz
                nEqlDiffPars = 1;
                Dzz = pars[2-nEqlDiffPars];
                alpha = pars[3-nEqlDiffPars];
                beta = pars[4-nEqlDiffPars];
                gamma = 0;
                break;
            case ANISOTROPIC:
                nEqlDiffPars = 0;
                Dxx = pars[0];
                Dyy = pars[1];
                Dzz = pars[2];
                alpha = pars[3];
                beta = pars[4];
                gamma = pars[5];
                break;
        }
        double[][] D = {{Dxx, 0.0, 0.0}, 
                        {0.0, Dyy, 0.0}, 
                        {0.0, 0.0, Dzz}};
        if (initD != null) {
            D = initD;
        }
        Rotation rot = null;
        try {
            rot = new Rotation(RotationOrder.ZYZ, RotationConvention.VECTOR_OPERATOR, alpha, beta, gamma);
        } catch (NotARotationMatrixException nE) {
            System.out.println("Can't create rot mat:" + nE.getMessage());
            double[][] rotMatCatch = rot.getMatrix();
            for (int i = 0; i < 3; i++) {
                rotMatCatch[1][i] = -rotMatCatch[1][i];
            }
            try {
                rot = new Rotation(rotMatCatch, 1e-6);
            } catch (NotARotationMatrixException nE2) {
                System.out.println("Can't create rot mat 2nd try:" + nE.getMessage());
                rot = null;
            }
        }
        if (rot != null) {
            double[][] rotMat = rot.getMatrix();
            DMatrixRMaj D1 = new DMatrixRMaj(D);
            DMatrixRMaj rotMat1 = new DMatrixRMaj(rotMat);
            DMatrixRMaj rotMat1T = new DMatrixRMaj(rotMat1.numRows, rotMat1.numCols);
            DMatrixRMaj res1 = new DMatrixRMaj(D1.numRows, D1.numCols);
            DMatrixRMaj res2 = new DMatrixRMaj(res1.numRows, res1.numCols);
            transpose(rotMat1, rotMat1T);
            CommonOps_DDRM.mult(rotMat1, D1, res1);
            CommonOps_DDRM.mult(res1, rotMat1T, res2);
//            DMatrixRMaj DScaled = new DMatrixRMaj(D1.numRows, D1.numCols);
            for (int i=0; i<res2.numRows; i++) {
                for (int j=0; j<res2.numCols; j++) {
                    D[i][j] = res2.get(i, j);
//                    DScaled.set(i, j, res2.get(i, j)/1e7);
                } 
            }
//            System.out.println(res2.toString());
//            System.out.println("scaled D = " + DScaled.toString());
        }
        return D;
    }

    public double getYVal(double[] pars, RelaxEquations relaxObj, double[] J, ExptType type) {
        double y = 0.0;
        switch (type) {
            case R1:
                y = relaxObj.R1(J);
//                System.out.print("R1 = " + y + " ");
                break;
            case R2:
//                double Rex = pars[0];
                y = relaxObj.R2(J, 0.0);
//                System.out.print("R2 = " + y + " ");
                break;
            case NOE:
                y = relaxObj.NOE(J);
//                System.out.print("NOE = " + y + " ");
                break;
            default:
                break;
        }
        return y;
    }

    public double[][] getSimValues(double first, double last, int n, boolean adjust) {
        double[][] result = new double[2][n];
        double delta = (last - first) / (n - 1);
//        int i = 0;
        for (int i=0; i<n; i++) {
            int iField = (int) xValues[i][0];
            RelaxEquations relaxObj = relaxObjs.get(iField);
            int iRes = (int) xValues[i][1];
            int iExpType = (int) xValues[i][2];
            int modelNum = residueModels[iRes];
            double[] J = getJ(bestPars, relaxObj, modelNum);
            ExptType type = expTypeArray[iExpType];
            double x = first + delta * i;
            double y = getYVal(bestPars, relaxObj, J, type);
            result[0][i] = x;
            result[1][i] = y;
        }
        return result;
    }

    public double value(double[] pars, double[][] values) {
        int n = values[0].length;
        double sum = 0.0;
        for (int i=0; i<n; i++) {
            int iField = (int) xValues[i][0];
            RelaxEquations relaxObj = relaxObjs.get(iField);
            int iRes = (int) xValues[i][1];
            int iExpType = (int) xValues[i][2];
            int modelNum = residueModels[iRes];
            double[] J = getJ(pars, relaxObj, modelNum);
            ExptType type = expTypeArray[iExpType];
            double y = getYVal(pars, relaxObj, J, type);
            double delta = y - values[2][i];
            sum += (delta * delta) / (errValues[i] * errValues[i]);
        }
        double rms = Math.sqrt(sum / n);

        return rms;

    }

    public double valueMultiResidue(double[] pars, double[][] values) {
        int n = values[0].length;
        double sum = 0.0;
        for (int i=0; i<values[2].length; i++) {
            int iField = (int) xValues[i][0];
            RelaxEquations relaxObj = relaxObjs.get(iField);
            int iRes = (int) xValues[i][1];
            int iExpType = (int) xValues[i][2];
            int modelNum = residueModels[iRes];
            double[] resPars = new double[nParsPerModel[modelNum] + 1];
            resPars[0] = pars[0];
            int parStart = resParStart[iRes] + 1;
            System.arraycopy(pars, parStart, resPars, 1, nParsPerModel[modelNum]);
//            parStart += resParStart[i];
            double[] J = getJ(resPars, relaxObj, modelNum);
            ExptType type = expTypeArray[iExpType];
            double y = getYVal(resPars, relaxObj, J, type);
            double delta = y - values[2][i];
            sum += (delta * delta) / (errValues[i] * errValues[i]);
        }

        double rms = Math.sqrt(sum / n);

        return rms;

    }
    
    public double diffSqDiffusion(int iRes, int iExpType, int i, double[][] values, 
            RelaxEquations relaxObj, double[] J, Map<Integer, double[]> expValMap, Map<Integer, double[]> expErrMap) {
        
        double val = 0.0;
        if (!expValMap.containsKey(iRes)) {
            expValMap.put(iRes, new double[3]);
            expErrMap.put(iRes, new double[3]);
        }
        if (expValMap.get(iRes)[iExpType] == 0.0) {
            expValMap.get(iRes)[iExpType] = values[2][i];
            expErrMap.get(iRes)[iExpType] = errValues[i];
        }
        double r1 = expValMap.get(iRes)[0];
        double r2 = expValMap.get(iRes)[1];
        double noe = expValMap.get(iRes)[2];
        if (r1 != 0.0 && r2 != 0.0 && noe != 0.0) {
            double r1Err = expErrMap.get(iRes)[0];
            double r2Err = expErrMap.get(iRes)[1];
            double noeErr = expErrMap.get(iRes)[2];
            double rhoExp = relaxObj.calcRhoExp(r1, r2, noe, J);
            double rhoPred = relaxObj.calcRhoPred(J);
            double delta = rhoPred - rhoExp; 
            double error = relaxObj.calcRhoExpError(r1, r2, noe, J, r1Err, r2Err, noeErr, rhoExp);//(r1Err + r2Err + noeErr) / 3;//
            val= (delta * delta) / (error * error); 
        }
        return val;
    }
    
    public double valueDiffusion(double[] pars, double[][] values) {
        int n = values[0].length;
        double sum = 0.0;
        Map<Integer, double[]> expValMap = new HashMap<>();
        Map<Integer, double[]> expErrMap = new HashMap<>();
        for (int i=0; i<n; i++) {
            int iField = (int) xValues[i][0];
            RelaxEquations relaxObj = relaxObjs.get(iField);
            int iRes = (int) xValues[i][1];
            int iExpType = (int) xValues[i][2];
            int iCoord = (int) xValues[i][3];
            int modelNum = residueModels[iRes];
            double[] v = coords[iCoord];
            double[] J = getJDiffusion(pars, relaxObj, modelNum, v, diffusionType, null);
            sum += diffSqDiffusion(iRes, iExpType, i, values, relaxObj, J, expValMap, expErrMap);
        }
        double rms = Math.sqrt(sum / n / 3);

        return rms;

    }

    public double valueDiffusionMultiResidue(double[] pars, double[][] values) {
        int n = values[0].length;
        double sum = 0.0;
        Map<Integer, double[]> expValMap = new HashMap<>();
        Map<Integer, double[]> expErrMap = new HashMap<>();
        for (int i=0; i<values[2].length; i++) {
            int iField = (int) xValues[i][0];
            RelaxEquations relaxObj = relaxObjs.get(iField);
            int iRes = (int) xValues[i][1];
            int iExpType = (int) xValues[i][2];
            int iCoord = (int) xValues[i][3];
            int modelNum = residueModels[iRes];
            int nDiffPars = 6;
            if (diffusionType == OBLATE || diffusionType == PROLATE) {
                nDiffPars = 4;
            }
            double[] resPars = new double[nParsPerModel[modelNum] + nDiffPars];
            for (int j=0; j<nDiffPars; j++) {
                resPars[j] = pars[j];
            }
            int parStart = resParStart[iRes] + nDiffPars;
            System.arraycopy(pars, parStart, resPars, nDiffPars, nParsPerModel[modelNum]);
//            parStart += resParStart[i];
            double[] v = coords[iCoord];
            double[] J = getJDiffusion(resPars, relaxObj, modelNum, v, diffusionType, null);
            sum += diffSqDiffusion(iRes, iExpType, i, values, relaxObj, J, expValMap, expErrMap);
        }

        double rms = Math.sqrt(sum / n);

        return rms;

    }
    
    public double getRMSDIso(double tauC) {
        int n = yValues.length;
        double sumIso = 0.0;
        Map<Integer, double[]> isoValMap = new HashMap<>();
        Map<Integer, double[]> expErrMap = new HashMap<>();
        for (int i=0; i<n; i++) {
            int iField = (int) xValues[i][0];
            RelaxEquations relaxObj = relaxObjs.get(iField);
            int iRes = (int) xValues[i][1];
            int iExpType = (int) xValues[i][2];
            double[] JIso = relaxObj.getJ(tauC);
            if (!isoValMap.containsKey(iRes)) {
                isoValMap.put(iRes, new double[3]);
                expErrMap.put(iRes, new double[3]);
            }
            if (isoValMap.get(iRes)[iExpType] == 0.0) {
                isoValMap.get(iRes)[iExpType] = yValues[i];
                expErrMap.get(iRes)[iExpType] = errValues[i];
            }
            double r1Iso = isoValMap.get(iRes)[0];
            double r2Iso = isoValMap.get(iRes)[1];
            double noeIso = isoValMap.get(iRes)[2];
            if (r1Iso != 0.0 && r2Iso != 0.0 && noeIso != 0.0) {
                double r1Err = expErrMap.get(iRes)[0];
                double r2Err = expErrMap.get(iRes)[1];
                double noeErr = expErrMap.get(iRes)[2];
                double rhoExpIso = relaxObj.calcRhoExp(r1Iso, r2Iso, noeIso, JIso);
                double rhoPredIso = relaxObj.calcRhoPred(JIso);
                double deltaIso = rhoPredIso - rhoExpIso; 
                double errorIso = relaxObj.calcRhoExpError(r1Iso, r2Iso, noeIso, JIso, r1Err, r2Err, noeErr, rhoExpIso);//(r1Err + r2Err + noeErr) / 3;//
                sumIso += (deltaIso * deltaIso) / (errorIso * errorIso); 
            }
        }
        double rms = Math.sqrt(sumIso / n / 3);

        return rms;

    }
    
    public double valueDMat(double[] pars, double[][] values) {
        double tauC = 4.2e-9;
        double Di = 1/(6*tauC);
        double[][] D = {{0.75*Di, 0, 0}, {0, Di, 0}, {0, 0, 1.25*Di}};
        if (diffusionType == PROLATE) {
            D[1][1] = 0.75*Di;
        } else if (diffusionType == OBLATE) {
            D[1][1] = 1.25*Di;
        } 
        double[][] x0 = rotateD(pars, diffusionType, D);
        int n = values[0].length;
        double sum = 0.0;
        Map<Integer, double[]> expValMap = new HashMap<>();
        Map<Integer, double[]> expErrMap = new HashMap<>();
        for (int i=0; i<n; i++) {
            int iField = (int) xValues[i][0];
            RelaxEquations relaxObj = relaxObjs.get(iField);
            int iRes = (int) xValues[i][1];
            int iExpType = (int) xValues[i][2];
            int iCoord = (int) xValues[i][3];
            int modelNum = residueModels[iRes];
            double[] v = coords[iCoord];
            double[] J = getJDiffusion(pars, relaxObj, modelNum, v, diffusionType, x0);
            sum += diffSqDiffusion(iRes, iExpType, i, values, relaxObj, J, expValMap, expErrMap);
        }
        double rms = Math.sqrt(sum / n / 3);

        return rms;

    }

    public double[] calcYVals(double[] pars, boolean multiRes, boolean diffusion) {
        int n = yValues.length;
        double[] calcY = new double[n];
        for (int i=0; i<n; i++) {
            int iField = (int) xValues[i][0];
            RelaxEquations relaxObj = relaxObjs.get(iField);
            int iRes = (int) xValues[i][1];
            int iExpType = (int) xValues[i][2];
            int iCoord = (int) xValues[i][3];
            int modelNum = residueModels[iRes];
            double[] J;
            if (diffusion) {
                double[] v = coords[iCoord];
                J = getJDiffusion(pars, relaxObj, modelNum, v, diffusionType, null);
            } else {
                J = getJ(pars, relaxObj, modelNum);
            }
            ExptType type = expTypeArray[iExpType];
            double y = getYVal(pars, relaxObj, J, type);
            if (multiRes) {
                double[] resPars;
                if (diffusion) {
                    int nDiffPars = 6;
                    if (diffusionType == OBLATE || diffusionType == PROLATE) {
                        nDiffPars = 4;
                    }
                    resPars = new double[nParsPerModel[modelNum] + nDiffPars];
                    for (int j=0; j<nDiffPars; j++) {
                        resPars[j] = pars[j];
                    }
                    int parStart = resParStart[iRes] + nDiffPars;
                    System.arraycopy(pars, parStart, resPars, nDiffPars, nParsPerModel[modelNum]);
                    double[] v = coords[iCoord];
                    J = getJDiffusion(resPars, relaxObj, modelNum, v, diffusionType, null);
                } else {
                    resPars = new double[nParsPerModel[modelNum] + 1];
                    resPars[0] = pars[0];
                    int parStart = resParStart[iRes] + 1;
                    System.arraycopy(pars, parStart, resPars, 1, nParsPerModel[modelNum]);
        //            parStart += resParStart[i];
                    J = getJ(resPars, relaxObj, modelNum);
                }
                y = getYVal(resPars, relaxObj, J, type);
            }
            calcY[i] = y;
        }

        return calcY;
    }

    public double[] getPars() {
        return bestPars;
    }

    public double getAIC() {
        return bestAIC;
    }

    public double getChiSq() {
        return bestChiSq;
    }

    public double[] getParErrs() {
        return parErrs;
    }
    
    public double[][] getXVals() {
        return xValues;
    }

    public double[] getYVals() {
        return yValues;
    }

    public double[] getErrVals() {
        return errValues;
    }
    
    public double[][] getCoords() {
        return coords;
    }

    public Fitter getFitter() {
        return bestFitter;
    }
    
    public int[] getResParStart() {
        return resParStart;
    }
    
    public int[] getResidueModels() {
        return residueModels;
    }
    
    public ArrayList<RelaxEquations> getRelaxObjs() {
        return relaxObjs;
    }
    
    public void setResidueModels(int[] bestModels, boolean diffusion) {
        residueModels = bestModels;
        resParStart = new int[residueModels.length];
        int nParsPrev = 0;
        resParStart[0] = 0;
        for (int i=1; i<residueModels.length; i++) {
            resParStart[i] = nParsPerModel[residueModels[i-1]] + nParsPrev;
            nParsPrev = resParStart[i];
//            System.out.println(i + " resModel " + residueModels[i] + " resParStart " + resParStart[i]);
        }
    }
    
    public void setTauMBounds(double[] bounds) {
        tauMBounds[0] = bounds[0];
        tauMBounds[1] = bounds[1];
    }

    public PointValuePair fit(double sf, double[] guesses, boolean globalFit, boolean diffusion) {
        Fitter fitter;
        this.B0 = sf * 2.0 * Math.PI / RelaxEquations.GAMMA_H;
        if (globalFit && diffusion) {
            fitter = Fitter.getArrayFitter(this::valueDiffusionMultiResidue);
        } else if (globalFit && !diffusion) {
            fitter = Fitter.getArrayFitter(this::valueMultiResidue);
        } else if (!globalFit && diffusion) {
            fitter = Fitter.getArrayFitter(this::valueDiffusion);
        } else {
            fitter = Fitter.getArrayFitter(this::value);
        }
        double[] xVals0 = new double[yValues.length];
        double[][] xValues2 = {xVals0, xVals0};//{xValues[0], xValues[1]};
        fitter.setXYE(xValues2, yValues, errValues);
        double[] start = guesses; //{max0, r0, max1, 10.0};
        double[] lower = new double[guesses.length]; //{max0 / 2.0, r0 / 2.0, max1 / 2.0, 1.0};
        double[] upper = new double[guesses.length]; //{max0 * 2.0, r0 * 2.0, max1 * 2.0, 100.0};
        for (int i = 0; i < lower.length; i++) {
            lower[i] = guesses[i] / 20.0;
            upper[i] = guesses[i] * 20.0;
        }
        int iS2diff = 6;
        //bounds for the global fit parameters (e.g. tauM, Dxx, Dyy, Dzz, alpha, 
        //beta, gamma), where the indices in the parameter array are always the
        //same regardless of the model used.
        if (diffusion) {
            if (diffusionType == OBLATE || diffusionType == PROLATE) {
                int nEqlDiffPars = 1;
                int nNullAngles = 1;
                int iAlpha = 3 - nEqlDiffPars;
                int iBeta = 4 - nEqlDiffPars;
                lower[iAlpha] = 0.0;
                upper[iAlpha] = Math.PI/2;
                lower[iBeta] = 0.0;
                upper[iBeta] = Math.PI/2;
                iS2diff -= (nEqlDiffPars + nNullAngles);
            } else { //anisotropic
                int iAlpha = 3;
                int iBeta = 4;
                int iGamma = 5;
                lower[iAlpha] = 0.0;
                upper[iAlpha] = Math.PI/2;
                lower[iBeta] = 0.0;
                upper[iBeta] = Math.PI/2;
                lower[iGamma] = 0.0;
                upper[iGamma] = Math.PI/2;
            }
        } else {
            int iTauM = 0;
            lower[iTauM] = tauMBounds[0];
            upper[iTauM] = tauMBounds[1];
        }
        //bounds for the other fit parameters (e.g. S2, tau, Sf2, tauS), which 
        //can have different indices in the parameter array depending on the 
        //model used and individual vs. global fitting.
        int iS2;
        if (globalFit) {
            for (int i=0; i<resParStart.length; i++) {
                int modelNum = residueModels[i];
                if (diffusion) {
                    iS2 = resParStart[i] + iS2diff;
                } else {
                    iS2 = resParStart[i] + 1; 
                }
                lower[iS2] = 0.0;
                upper[iS2] = 1.0;
                if (modelNum >= 5) {
                    int iS2f = iS2 + 2;
                    lower[iS2f] = 0.0;
                    upper[iS2f] = 1.0;
                }
            }
        } else {
            if (diffusion) {
                iS2 = iS2diff;
            } else {
                iS2 = 1;
            }
            lower[iS2] = 0.0;
            upper[iS2] = 1.0;
            int iS2f = iS2 + 2;
            if (lower.length > iS2f) {
                lower[iS2f] = 0.0;
                upper[iS2f] = 1.0;
            }
        }
//        for (int i=0; i<guesses.length; i++) {
//            System.out.println("guess lower upper " + i + " " + guesses[i] + " " + lower[i] + " " + upper[i]);
//        }
        try {
            PointValuePair result = fitter.fit(start, lower, upper, 10.0);
            bestPars = result.getPoint();
            bestChiSq = result.getValue();
            int k = bestPars.length;
            int n = yValues.length;
            bestAIC = n * Math.log(bestChiSq) + 2 * k;
//            bestAIC += (2*k*k + 2*k)/(n - k - 1); //corrected AIC for small sample sizes
            bestFitter = fitter;
            return result;
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }

    }

    public void calcParErrs(PointValuePair result, Fitter fitter, double[] fitPars, boolean globalFit, boolean diffusion) {
        double[] yCalc;
        yCalc = calcYVals(fitPars, globalFit, diffusion);
        parErrs = fitter.bootstrap(result.getPoint(), 300, true, yCalc);
    }
    
    public PointValuePair fitD(double sf, double[] guesses) {
        Fitter fitter = Fitter.getArrayFitter(this::valueDMat);
        this.B0 = sf * 2.0 * Math.PI / RelaxEquations.GAMMA_H;
        double[] xVals0 = new double[yValues.length];
        double[][] xValues2 = {xVals0, xVals0};//{xValues[0], xValues[1]};
        fitter.setXYE(xValues2, yValues, errValues);
        double[] start = guesses; //{max0, r0, max1, 10.0};
        double[] lower = new double[guesses.length]; //{max0 / 2.0, r0 / 2.0, max1 / 2.0, 1.0};
        double[] upper = new double[guesses.length]; //{max0 * 2.0, r0 * 2.0, max1 * 2.0, 100.0};
        for (int i = 0; i < lower.length; i++) {
            lower[i] = guesses[i] / 20.0;
            upper[i] = guesses[i] * 20.0;
        }
        //bounds for the global fit parameters (e.g. tauM, Dxx, Dyy, Dzz, alpha, 
        //beta, gamma), where the indices in the parameter array are always the
        //same regardless of the model used.
        if (diffusionType == OBLATE || diffusionType == PROLATE) {
            int iAlpha = 2;
            int iBeta = 3;
            lower[iAlpha] = 0.0;
            upper[iAlpha] = Math.PI/2;
            lower[iBeta] = 0.0;
            upper[iBeta] = Math.PI/2;
        } else { //anisotropic
            int iAlpha = 3;
            int iBeta = 4;
            int iGamma = 5;
            lower[iAlpha] = 0.0;
            upper[iAlpha] = Math.PI/2;
            lower[iBeta] = 0.0;
            upper[iBeta] = Math.PI/2;
            lower[iGamma] = 0.0;
            upper[iGamma] = Math.PI/2;
        }
//        for (int i=0; i<guesses.length; i++) {
//            System.out.println("guess lower upper " + i + " " + guesses[i] + " " + lower[i] + " " + upper[i]);
//        }
        try {
            PointValuePair result = fitter.fit(start, lower, upper, 10.0);
            bestPars = result.getPoint();
            bestChiSq = result.getValue();
            int k = bestPars.length;
            int n = yValues.length;
            bestAIC = n * Math.log(bestChiSq) + 2 * k;
//            bestAIC += (2*k*k + 2*k)/(n - k - 1); //corrected AIC for small sample sizes
            bestFitter = fitter;
            return result;
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }

    }
}
