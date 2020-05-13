/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.optim.PointValuePair;
import org.comdnmr.data.BondVectorData;
import org.comdnmr.modelfree.RelaxDataValue;
import org.comdnmr.modelfree.RelaxEquations;
import org.comdnmr.modelfree.RelaxFit;
import org.comdnmr.modelfree.RelaxFit.DiffusionType;
import org.junit.Test;
import org.junit.Assert;

/**
 *
 * @author Martha
 */
public class DRefineTest {

//    private final double[] guesses = {0.75 * Di, Di, 1.25 * Di, Math.PI / 2, Math.PI / 2, Math.PI / 2};
//    private final double[] guesses = {4.4170 * 1e7, 4.5832 * 1e7, 6.0129 * 1e7, Math.toRadians(98.06),
//            Math.toRadians(68.64), Math.toRadians(77.42)};
    private final double[] fields = {400.0e6, 500e6, 600e6, 700e6, 800.0e6};

    private BondVectorData loadVectors(String vectorFileName) throws IOException {
        File file = new File(vectorFileName);
        Path path = file.toPath();
        Map<String, double[]> coordMap = new HashMap<>();

        Files.lines(path).forEach(line -> {
            String[] fields = line.split("\t");
            if (fields.length == 5) {
                int iRes = Integer.valueOf(fields[0]);
                String chain = "A";
                String resSpecifier = chain + iRes;
                double[] v = new double[3];
                v[0] = Double.valueOf(fields[2]);
                v[1] = Double.valueOf(fields[3]);
                v[2] = Double.valueOf(fields[4]);
                coordMap.put(resSpecifier, v);
            }
        });
        Map<String, Integer> coordIndices = new HashMap<>();
        double[][] vectors = new double[coordMap.size()][3];
        List<String> names = new ArrayList<>();
        coordMap.keySet().stream().sorted().forEach(k -> {

            int i = coordIndices.size();
            double[] v = coordMap.get(k);
            vectors[i][0] = v[0];
            vectors[i][1] = v[1];
            vectors[i][2] = v[2];
            coordIndices.put(k, i);
            names.add(k);
        });
        BondVectorData bData = new BondVectorData(vectors, names, coordIndices);
        return bData;
    }

    private List<RelaxDataValue> loadData(String relaxDataFileName, BondVectorData bData) throws IOException {

        File relaxFile = new File(relaxDataFileName);
        Path relaxPath = relaxFile.toPath();
        Map<Integer, Integer> fieldMap = new HashMap<>();
        List<RelaxDataValue> result = new ArrayList<>();
        Files.lines(relaxPath).forEach(line -> {
            String[] fields = line.split("\\s+");
            if (fields.length == 11) {
                int iRes = Integer.valueOf(fields[0]);
                String chain = fields[1];
                String resSpecifier = chain + iRes;
                Integer resIndex = bData.getIndex(resSpecifier);
                if (resIndex == null) {
                    System.out.println("no res " + resSpecifier);
                } else {
                    double field = Double.valueOf(fields[4]);
                    String e1 = fields[3];
                    String e2 = fields[2];
                    RelaxEquations relaxObj = RelaxEquations.getRelaxEquations(field * 1e6, e1, e2);

                    double r1 = Double.valueOf(fields[5]);
                    double r1Error = Double.valueOf(fields[6]);
                    double r2 = Double.valueOf(fields[7]);
                    double r2Error = Double.valueOf(fields[8]);
                    double noe = Double.valueOf(fields[9]);
                    double noeError = Double.valueOf(fields[10]);
                    RelaxDataValue dValue = new RelaxDataValue(resIndex, resSpecifier, r1, r1Error, r2, r2Error, noe, noeError, relaxObj);
                    result.add(dValue);
                }
            }
        });
        return result;
    }

    public double[] scalePars(double[] pars, int nDiffPars) {
        double[] scaledPars = new double[pars.length];
        for (int p = 0; p < pars.length; p++) {
            double newPar = pars[p];
            if (p < nDiffPars) {
                newPar /= 1.0e7;
            } else if (p >= nDiffPars) {
                newPar *= 180.0 / Math.PI;
            }
            scaledPars[p] = newPar;
        }
        return scaledPars;
    }

    private void loadTestData(RelaxFit relaxFit) {
        try {
            BondVectorData bData = loadVectors("src/test/data/1P7F_A_vectors.txt");
            System.out.println(bData);
            List<RelaxDataValue> rData = loadData("src/test/data/1P7F_data.txt", bData);
            relaxFit.setBondVectorData(bData);
            relaxFit.setRelaxData(rData);
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
            return;
        }

    }

    @Test
    public void testValueDMat1() {
        RelaxFit relaxFit = new RelaxFit();
        loadTestData(relaxFit);
        relaxFit.setDiffusionType(DiffusionType.ANISOTROPIC);
        double[] pars = {4.4170 * 1e7, 4.5832 * 1e7, 6.0129 * 1e7, Math.toRadians(98.06),
            Math.toRadians(68.64), Math.toRadians(77.42)};

        double value = relaxFit.valueDMat(pars, null);
        System.out.println("value " + value);
        Assert.assertEquals(0.0, value, 1.2e-1);
    }

    @Test
    public void testValueDMatFile() {
        RelaxFit relaxFit = new RelaxFit();
        loadTestData(relaxFit);
        double tauCGuess = 3.3e-9;

        double[] rotDifPars = {4.4170, 4.5832, 6.0129, 98.06, 68.64, 77.42};
        double bestFitRMS = Double.MAX_VALUE;
        double[] bestGuesses = null;
        double[] bestFitPars = null;
        DiffusionType bestType = null;
        double isoD = 1.0 / (6.0 * tauCGuess);
        for (DiffusionType diffType : DiffusionType.values()) {
            if (diffType == DiffusionType.ISOTROPIC) {
                continue;
            }
            int nPars = diffType.getNAnglePars() + diffType.getNDiffusionPars();
            double[] guess = new double[nPars];
            relaxFit.setDiffusionType(diffType);
            System.arraycopy(diffType.getGuess(isoD), 0, guess, 0, diffType.getNDiffusionPars());
            int nGuesses = diffType.getNAngleGuesses();
            for (int iGuess = 0; iGuess < nGuesses; iGuess++) {
                double[] angleGuess = diffType.getAngles(iGuess);
                System.arraycopy(angleGuess, 0, guess, diffType.getNDiffusionPars(), diffType.getNAnglePars());
                PointValuePair fitResult = relaxFit.fitDiffusion(guess);
                double fitRMS = fitResult.getValue();
                double[] scaledGuesses = scalePars(guess, diffType.getNDiffusionPars());
                double[] fitPars = relaxFit.getPars();
                Arrays.sort(fitPars, 0, diffType.getNDiffusionPars());
                double[] scaledPars = scalePars(fitPars, diffType.getNDiffusionPars());

//                System.out.println("RMS: " + fitRMS);
//                System.out.println("Scaled guesses, Fit Pars, RotDif Pars: ");
//                for (int i = 0; i < scaledGuesses.length; i++) {
//                    System.out.printf("guess %7.3f best Fit %7.3f RotDif %7.3f\n",
//                            scaledGuesses[i], scaledPars[i], rotDifPars[i]);
//                }
//                System.out.println();
                if (fitRMS < bestFitRMS) {
                    bestGuesses = guess.clone();
                    bestFitRMS = fitRMS;
                    bestFitPars = fitPars.clone();
                    bestType = diffType;
                }
            }
        }
        relaxFit.setDiffusionType(bestType);
        double[][][] rotResults = relaxFit.rotateD(bestFitPars);
        double[][] D = rotResults[0];
        double[][] VT = rotResults[1];

        double[] scaledBestGuesses = scalePars(bestGuesses, bestType.getNDiffusionPars());
        double[] scaledBestFitPars = scalePars(bestFitPars, bestType.getNDiffusionPars());

        System.out.println("\nbest Fit RMS: " + bestFitRMS + " " + bestType.toString());
        System.out.println("Scaled best guesses, best Fit Pars, RotDif Pars: ");
        for (int i = 0; i < scaledBestGuesses.length; i++) {
            System.out.printf("guess %7.3f best Fit %7.3f RotDif %7.3f\n",
                    scaledBestGuesses[i], scaledBestFitPars[i], rotDifPars[i]);
        }
        System.out.println();
        System.out.println("best Fit D = " + new Array2DRowRealMatrix(D).toString());
        System.out.println("best Fit VT = " + new Array2DRowRealMatrix(VT).toString());
        //  relaxFit.dumpValues(bestFitPars);

    }

}
