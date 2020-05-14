package org.comdnmr.modelfree;

import java.util.Arrays;
import java.util.Map;
import org.apache.commons.math3.geometry.euclidean.threed.NotARotationMatrixException;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.RotationConvention;
import org.apache.commons.math3.geometry.euclidean.threed.RotationOrder;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.optim.PointValuePair;
import org.comdnmr.data.Fitter;
import static org.comdnmr.modelfree.RelaxFit.DiffusionType.ANISOTROPIC;
import static org.comdnmr.modelfree.RelaxFit.DiffusionType.OBLATE;
import static org.comdnmr.modelfree.RelaxFit.DiffusionType.PROLATE;

/**
 *
 * @author brucejohnson
 */
public class RelaxFit {

    boolean reportFitness = true;
    int reportAt = 10;
    long startTime = 0;
    double B0;
    Map<String, MolDataValues> molDataValues;
    double[] bestPars;
    double[] parErrs;
    double bestAIC;
    double bestChiSq;
    Fitter bestFitter;
    int[] nParsPerModel = {0, 1, 2, 0, 0, 3, 4};
    DiffusionType diffusionType;
    static int ANISO_ANGLE_STARTS[][] = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {1, 1, 0},
        {0, 0, 1},
        {1, 0, 1},
        {0, 1, 1},
        {1, 1, 1}
    };
    static int ANGLE_STARTS[][] = {
        {0, 0},
        {1, 0},
        {0, 1},
        {1, 1}
    };

    double globalTau = 4.0e-9;
    boolean useGlobalTau = false;

    public double getGlobalTau() {
        return globalTau;
    }

    public void setGlobalTau(double globalTau) {
        this.globalTau = globalTau;
    }

    public boolean isUseGlobalTau() {
        return useGlobalTau;
    }

    public void setUseGlobalTau(boolean useGlobalTau) {
        this.useGlobalTau = useGlobalTau;
    }

    public static double[] getDValues(double isoD) {
        double[] anisoD = {0.75 * isoD, isoD, 1.25 * isoD};
        return anisoD;
    }

    public void setRelaxData(Map<String, MolDataValues> molDataValues) {
        this.molDataValues = molDataValues;
    }

    public void setDiffusionType(DiffusionType type) {
        this.diffusionType = type;
    }

    public enum DiffusionType {
        ISOTROPIC(1, 0) {
            @Override
            public double[] getGuess(double isoD) {
                double[] guess = {isoD};
                return guess;
            }

            public double[] getAngles(int iStart) {
                return new double[0];
            }
        },
        PROLATE(2, 2) {
            @Override
            public double[] getGuess(double isoD) {
                double[] guess = {0.75 * isoD, 1.25 * isoD};
                return guess;
            }

            public double[] getAngles(int iStart) {
                int[] jAng = ANGLE_STARTS[iStart];
                double[] angles = new double[2];
                for (int i = 0; i < 2; i++) {
                    angles[i] = (Math.PI * (jAng[i] * 2 + 1)) / 4;
                }
                return angles;
            }
        }, OBLATE(2, 2) {
            @Override
            public double[] getGuess(double isoD) {
                double[] guess = {0.75 * isoD, 1.25 * isoD};
                return guess;
            }

            public double[] getAngles(int iStart) {
                int[] jAng = ANGLE_STARTS[iStart];
                double[] angles = new double[2];
                for (int i = 0; i < 2; i++) {
                    angles[i] = (Math.PI * (jAng[i] * 2 + 1)) / 4;
                }
                return angles;
            }
        }, ANISOTROPIC(3, 3) {
            @Override
            public double[] getGuess(double isoD) {
                double[] anisoD = {0.75 * isoD, isoD, 1.25 * isoD};
                return anisoD;
            }

            public double[] getAngles(int iStart) {
                int[] jAng = ANISO_ANGLE_STARTS[iStart];
                double[] angles = new double[3];
                for (int i = 0; i < 3; i++) {
                    angles[i] = (Math.PI * (jAng[i] * 2 + 1)) / 4;
                }
                return angles;
            }
        };

        int nAnglePars;
        int nDiffPars;
        int nAngleGuesses;

        DiffusionType(int nDiffPars, int nAnglePars) {
            this.nAnglePars = nAnglePars;
            this.nDiffPars = nDiffPars;
            this.nAngleGuesses = (int) Math.pow(2, nAnglePars);
        }

        public abstract double[] getGuess(double isoD);

        public abstract double[] getAngles(int iStart);

        public int getNAnglePars() {
            return nAnglePars;
        }

        public int getNDiffusionPars() {
            return nDiffPars;
        }

        public int getNAngleGuesses() {
            return nAngleGuesses;
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

    public double[] getJDiffusion(double[] pars, RelaxEquations relaxObj, int modelNum, double[] v, DiffusionType dType, double[][] D, double[][] VT) {
        int nEqlDiffPars = 0;
        int nNullAngles = 0;
        if (dType == PROLATE || dType == OBLATE) { //Dxx = Dyy or Dyy = Dzz
            nEqlDiffPars = 1;
            nNullAngles = 1;
        }

        if (D == null && VT == null) {
            double Dxx = pars[0];
            double Dyy = pars[1];
            double Dzz = 0.0;
            switch (dType) {
                case PROLATE:
                    //Dxx = Dyy
                    nEqlDiffPars = 1;
                    Dyy = pars[1 - nEqlDiffPars];
                    Dzz = pars[2 - nEqlDiffPars];
                    break;
                case OBLATE:
                    //Dyy = Dzz
                    nEqlDiffPars = 1;
                    Dzz = pars[2 - nEqlDiffPars];
                    break;
                case ANISOTROPIC:
                    nEqlDiffPars = 0;
                    Dxx = pars[0];
                    Dyy = pars[1];
                    Dzz = pars[2];
                    break;
            }
            double[][] D1 = {{Dxx, 0.0, 0.0},
            {0.0, Dyy, 0.0},
            {0.0, 0.0, Dzz}};

            Rotation rot = getDRotation(pars, dType);
            D = D1;
            VT = getRotationMatrix(rot);
//            System.out.println("Rotated Deig = " + new Array2DRowRealMatrix(D).toString());
//            System.out.println("Rotated VTeig = " + new Array2DRowRealMatrix(VT).toString());
        }
        int nDiffPars = 6;
        double s2 = pars[nDiffPars - nEqlDiffPars - nNullAngles];
        double[] J = new double[5];
        switch (modelNum) {
            case 1:
                J = relaxObj.getJDiffusion(dType, D, VT, v, s2, null, null, null);
                break;
            case 2:
                double tau = pars[nDiffPars + 1 - nEqlDiffPars - nNullAngles];
                J = relaxObj.getJDiffusion(dType, D, VT, v, s2, tau, null, null);
                break;
            case 5:
                tau = pars[nDiffPars + 1 - nEqlDiffPars - nNullAngles];
                double sf2 = pars[nDiffPars + 2 - nEqlDiffPars - nNullAngles];
                J = relaxObj.getJDiffusion(dType, D, VT, v, s2, tau, sf2, null);
                break;
            case 6:
                tau = pars[nDiffPars + 1 - nEqlDiffPars - nNullAngles];
                sf2 = pars[nDiffPars + 2 - nEqlDiffPars - nNullAngles];
                double tauS = pars[nDiffPars + 3 - nEqlDiffPars - nNullAngles];
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

    public Rotation getDRotation(double[] pars, DiffusionType dType) {
        int nEqlDiffPars = 0;
        double alpha = 0.0;
        double beta = 0.0;
        double gamma = 0.0;
        switch (dType) {
            case PROLATE:
                //Dxx = Dyy
                nEqlDiffPars = 1;
                alpha = pars[3 - nEqlDiffPars];
                beta = pars[4 - nEqlDiffPars];
                gamma = 0;
                break;
            case OBLATE:
                //Dyy = Dzz
                nEqlDiffPars = 1;
                alpha = pars[3 - nEqlDiffPars];
                beta = pars[4 - nEqlDiffPars];
                gamma = 0;
                break;
            case ANISOTROPIC:
                nEqlDiffPars = 0;
                alpha = pars[3];
                beta = pars[4];
                gamma = pars[5];
                break;
        }
//        System.out.println("alpha = " + alpha*180./Math.PI + " beta = " + beta*180./Math.PI + " gamma = " + gamma*180./Math.PI);
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
        return rot;
    }

    public double[][] getRotationMatrix(Rotation rot) {
        return new Array2DRowRealMatrix(rot.getMatrix()).transpose().getData();
    }

    public double[][][] rotateD(double[] resPars) {
        double[][][] resList = new double[2][3][3];
        double dx = resPars[0];
        double d1 = resPars[1];
        double d2 = 0.0;
        double dy = 0.0;
        double dz = 0.0;
        if (diffusionType == PROLATE) {
            dy = dx;
            dz = d1;
        } else if (diffusionType == OBLATE) {
            dy = d1;
            dz = d1;
        } else if (diffusionType == ANISOTROPIC) {
            dy = d1;
            dz = resPars[2];
        }
        double[][] D = {{dx, 0.0, 0.0}, {0.0, dy, 0.0}, {0.0, 0.0, dz}};
        Rotation rot = getDRotation(resPars, diffusionType);
        double[][] VT = getRotationMatrix(rot);//getVT();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                resList[0][i][j] = D[i][j];
                resList[1][i][j] = VT[i][j];
            }
        }

        return resList;
    }

    public double value(double[] pars, double[][] values) {
        double sumSq = 0.0;
        int n = 0;
        boolean parsOK = true;
        for (MolDataValues molData : molDataValues.values()) {
            int testModel = molData.getTestModel();
            double[] resPars;
            if (useGlobalTau) {
                int nResPars = nParsPerModel[testModel];
                resPars = new double[nResPars + 1];
                resPars[0] = globalTau;
                System.arraycopy(pars, 0, resPars, 1, nResPars);
            } else {
                resPars = pars;
            }
            if (!checkParConstraint(resPars, testModel)) {
                parsOK = false;
            }

            for (RelaxDataValue dValue : molData.getData()) {
                RelaxEquations relaxObj = dValue.relaxObj;
                double[] J = getJ(resPars, relaxObj, testModel);
                double r1 = relaxObj.R1(J);
                double r2 = relaxObj.R2(J, 0.0);
                double noe = relaxObj.NOE(J);
                double delta2 = dValue.score2(r1, r2, noe);
                sumSq += delta2;
                n += 3;
            }
        }
        double rms = Math.sqrt(sumSq / n);
        if (!parsOK) {
            rms += n * 10.0;
        }
        return rms;
    }

    public double valueMultiResidue(double[] pars, double[][] values) {
        double sumSq = 0.0;
        int n = 0;
        boolean parsOK = true;
        int parStart = 1;
        for (MolDataValues molData : molDataValues.values()) {
            int testModel = molData.getTestModel();
            int nResPars = nParsPerModel[testModel];
            double[] resPars = new double[nResPars + 1];
            resPars[0] = pars[0];
            System.arraycopy(pars, parStart, resPars, 1, nResPars);
            parStart += nResPars;
            if (!checkParConstraint(resPars, testModel)) {
                parsOK = false;
            }

            for (RelaxDataValue dValue : molData.getData()) {
                RelaxEquations relaxObj = dValue.relaxObj;
                double[] J = getJ(resPars, relaxObj, testModel);
                double r1 = relaxObj.R1(J);
                double r2 = relaxObj.R2(J, 0.0);
                double noe = relaxObj.NOE(J);
                double delta2 = dValue.score2(r1, r2, noe);
                sumSq += delta2;
                n += 3;
            }
        }
        double rms = Math.sqrt(sumSq / n);
        if (!parsOK) {
            rms += n * 10.0;
        }
        return rms;
    }

    public double valueDMat(double[] pars, double[][] values) {
        if (diffusionType == OBLATE || diffusionType == PROLATE) {
            Arrays.sort(pars, 0, 2);
        } else if (diffusionType == ANISOTROPIC) {
            Arrays.sort(pars, 0, 3);
        }
        double[][][] rotResults = rotateD(pars);
        double[][] D = rotResults[0];
        double[][] VT = rotResults[1];
        double sumSq = 0.0;
        int modelNum = 1;
        int nDiffPars = diffusionType.getNDiffusionPars() + diffusionType.getNAnglePars();
        int n = 0;
        for (MolDataValues molData : molDataValues.values()) {
            for (RelaxDataValue dValue : molData.getData()) {
                RelaxEquations relaxObj = dValue.relaxObj;
                double[] v = molData.vector;
                double[] resPars = new double[nParsPerModel[modelNum] + nDiffPars];
                System.arraycopy(pars, 0, resPars, 0, nDiffPars);
                resPars[resPars.length - 1] = 1.0; //Model 0: S2 = 1.0, all others null.
                double[] J = getJDiffusion(resPars, relaxObj, modelNum, v, diffusionType, D, VT);
                double rhoExp = dValue.calcExpRho(J);
                double rhoPred = dValue.calcPredRho(J);
                double delta = rhoPred - rhoExp;
                sumSq += delta * delta;
                n++;
            }
        }
        double rms = Math.sqrt(sumSq / n);

        return rms;

    }

    public void dumpValues(double[] pars) {
        if (diffusionType == OBLATE || diffusionType == PROLATE) {
            Arrays.sort(pars, 0, 2);
        } else if (diffusionType == ANISOTROPIC) {
            Arrays.sort(pars, 0, 3);
        }
        double[][][] rotResults = rotateD(pars);
        double[][] D = rotResults[0];
        double[][] VT = rotResults[1];
        int modelNum = 1;
        int nDiffPars = diffusionType.getNDiffusionPars() + diffusionType.getNAnglePars();
        for (MolDataValues molData : molDataValues.values()) {
            for (RelaxDataValue dValue : molData.getData()) {
                RelaxEquations relaxObj = dValue.relaxObj;
                double[] v = molData.vector;
                double[] resPars = new double[nParsPerModel[modelNum] + nDiffPars];
                System.arraycopy(pars, 0, resPars, 0, nDiffPars);
                resPars[resPars.length - 1] = 1.0; //Model 0: S2 = 1.0, all others null.
                double[] J = getJDiffusion(resPars, relaxObj, modelNum, v, diffusionType, D, VT);
                double rhoExp = dValue.calcExpRho(J);
                double rhoPred = dValue.calcPredRho(J);
                System.out.println(rhoExp + " " + rhoPred + " " + (rhoExp - rhoPred) + " " + molData.specifier);
            }
        }
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

    public Fitter getFitter() {
        return bestFitter;
    }

    private boolean checkParConstraint(double[] pars, int modelNum) {
        boolean result;
        double scale = 1.0;
        switch (modelNum) {
            case 2:
                result = pars[2] < scale * pars[0];
                break;
            case 5:
                result = pars[2] < scale * pars[0];
                break;
            case 6:
                result = pars[2] < scale * pars[0] && pars[4] < scale * pars[0];
                break;
            default:
                result = true;

        }
        return result;

    }

    public static double[][] guesses(int model, double tau) {
        double[] guess1 = {tau, 0.9};
        double[] lower1 = {tau / 10.0, 0.0};
        double[] upper1 = {tau * 10.0, 1.0};
        double[] guess2 = {tau, 0.9, tau / 10.0};
        double[] lower2 = {tau / 10.0, 0.0};
        double[] upper2 = {tau * 10.0, 1.0};
        double[] guess5 = {tau, 0.9, tau / 10.0, 0.9};
        double[] lower5 = {tau / 10.0, 0.0, tau / 100.0, 0.0};
        double[] upper5 = {tau * 10.0, 1.0, tau * 2.0, 1.0};
        double[] guess6 = {tau, 0.9, tau / 10.0, 0.9, tau / 4.0, 0.9};
        double[] lower6 = {tau / 10.0, 0.0, tau / 100.0, tau / 100.0, 0.0};
        double[] upper6 = {tau * 10.0, 1.0, tau * 2.0, tau * 2.0, 1.0};
        double[][] guesses = {null, guess1, guess2, null, null, guess5, guess6};
        double[][] lower = {null, lower1, lower2, null, null, lower5, lower6};
        double[][] upper = {null, upper1, upper2, null, null, upper5, upper6};
        double[][] result = new double[3][];
        result[0] = guesses[model];
        result[1] = lower[model];
        result[2] = upper[model];
        return result;

    }

    public PointValuePair fitResidueToModel(double[] start, double[] lower, double[] upper) {
        Fitter fitter = Fitter.getArrayFitter(this::value);
        try {
            PointValuePair result = fitter.fit(start, lower, upper, 10.0);
            return result;
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }
    }

    public PointValuePair fitMultiResidueToModel(double[] start, double[] lower, double[] upper) {
        Fitter fitter = Fitter.getArrayFitter(this::valueMultiResidue);
        try {
            PointValuePair result = fitter.fit(start, lower, upper, 10.0);
            return result;
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }
    }

    public PointValuePair fitDiffusion(double[] guesses) {
        Fitter fitter = Fitter.getArrayFitter(this::valueDMat);
        double[] start = guesses;
        double[] lower = new double[guesses.length];
        double[] upper = new double[guesses.length];
        int nDiffPars = diffusionType.getNDiffusionPars();
        int nAnglePars = diffusionType.getNAnglePars();
//        System.out.println(diffusionType + " " + nDiffPars + " " + nAnglePars);
        for (int i = 0; i < nDiffPars; i++) {
            lower[i] = guesses[i] / 2;
            upper[i] = guesses[i] * 2;
        }
        for (int i = nDiffPars; i < nDiffPars + nAnglePars; i++) {
            lower[i] = guesses[i] - Math.PI / 4.0;
            upper[i] = guesses[i] + Math.PI / 4.0;
        }
        try {
            PointValuePair result = fitter.fit(start, lower, upper, 10.0);
//            System.out.println("Scaled guess, bounds:");
            for (int i = 0; i < lower.length; i++) {
                double lb = lower[i];
                double ub = upper[i];
                double guess = guesses[i];
                if (i < nDiffPars) {
                    lb /= 1e7;
                    ub /= 1e7;
                    guess /= 1e7;
                } else {
                    lb = Math.toDegrees(lb);
                    ub = Math.toDegrees(ub);
                    guess = Math.toDegrees(guess);
                }
//                System.out.printf("guess %7.3f LB %7.3f UB %7.3f\n", guess, lb, ub);
            }
            bestPars = result.getPoint();
            bestChiSq = result.getValue();
            return result;
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }
    }
}
