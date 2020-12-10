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
package org.comdnmr.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.comdnmr.eqnfit.CESTEquation;
import org.comdnmr.eqnfit.CESTFitter;
import org.comdnmr.eqnfit.CPMGEquation;
import org.comdnmr.eqnfit.CPMGFitter;
import org.comdnmr.eqnfit.CurveFit;
import org.comdnmr.eqnfit.EquationType;
import org.comdnmr.eqnfit.ExpEquation;
import org.comdnmr.eqnfit.ExpFitter;
import org.comdnmr.eqnfit.PlotEquation;
import org.comdnmr.eqnfit.R1RhoEquation;
import org.comdnmr.eqnfit.R1RhoFitter;
import org.nmrfx.chemistry.Compound;
import org.nmrfx.chemistry.Entity;
import org.nmrfx.chemistry.InvalidMoleculeException;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.Polymer;
import org.nmrfx.chemistry.Residue;
import org.nmrfx.chemistry.io.MMcifReader;
import org.nmrfx.chemistry.io.MoleculeIOException;
import org.nmrfx.chemistry.io.NMRNEFReader;
import org.nmrfx.chemistry.io.NMRStarReader;
import org.nmrfx.chemistry.io.NMRStarWriter;
import org.nmrfx.chemistry.io.PDBFile;
import org.nmrfx.datasets.DatasetBase;
import org.nmrfx.peaks.InvalidPeakException;
import org.nmrfx.peaks.PeakList;
import org.nmrfx.star.ParseException;
import org.yaml.snakeyaml.Yaml;

/**
 *
 * @author Bruce Johnson
 */
public class DataIO {

    enum XCONV {
        IDENTITY() {
        },
        TAU2() {
            @Override
            double convert(double value, double[] pars, ExperimentData expData) {
                return 1000.0 / (2.0 * value);
            }
        },
        PPMTOHZ() {
            @Override
            double convert(double value, double[] pars, ExperimentData expData) {
                return value * expData.getNucleusField();
            }
        },
        HZTOPPM() {
            @Override
            double convert(double value, double[] pars, ExperimentData expData) {
                return value / expData.getNucleusField();
            }
        },
        CALC() {
            @Override
            double convert(double value, double[] pars, ExperimentData expData) {
                return pars[0] + pars[2] * (value + pars[1]);
            }
        };

        double convert(double value, double[] pars, ExperimentData expData) {
            return value;
        }
    }

    enum YCONV {
        IDENTITY() {
        },
        RATE() {
            @Override
            Double convert(double value, double refValue, double tau) {
                Double result = null;
                if ((value < refValue) && (value > 0.0)) {
                    result = -Math.log(value / refValue) / tau;

                }
                return result;
            }
        },
        NORMALIZE() {
            @Override
            Double convert(double value, double refValue, double tau) {
                return value / refValue;
            }
        };

        Double convert(double value, double refValue, double tau) {
            return value;
        }
    }

    static Pattern resPatter = Pattern.compile("[^0-9]*([0-9]+)[^0-9]*");
    //Nested Map for T1/T2 fit results: <expType, List<frameName, fields, fitResultsMap <resNum, <field, <parName, parValue>>>>>
    static Map<String, List<Object>> expFitResults = new TreeMap<>();

    public static void loadFromPeakList(PeakList peakList, ExperimentData expData,
            ResidueProperties resProp, String xConvStr, String yConvStr) {
        loadFromPeakList(peakList, expData, resProp, XCONV.valueOf(xConvStr.toUpperCase()),
                YCONV.valueOf(yConvStr.toUpperCase()));

    }

    public static void loadFromPeakList(PeakList peakList, ExperimentData expData,
            ResidueProperties resProp, XCONV xConv, YCONV yConv) {
        String expMode = expData.getExpMode();
        DatasetBase dataset = DatasetBase.getDataset(peakList.fileName);
        final double[] xValues;
        if (dataset != null) {
            xValues = dataset.getValues(2);
        } else {
            xValues = null;
        }
        resProp.addExperimentData(expData.getName(), expData);
        boolean eSet = false;
        peakList.peaks().stream().filter(peak -> peak.getStatus() >= 0).forEach(peak -> {
            String resName = "";
            int peakField = -1;
            String residueNum;
            int peakNum = peak.getIdNum();
            String label = peak.getPeakDim(0).getLabel();
            if (label.contains(".")) {
                residueNum = label.substring(0, label.indexOf("."));
            } else {
                residueNum = String.valueOf(peak.getIndex());
            }
            int residueNumInt = 0;
            if (residueNum.length() > 0) {
                try {
                    residueNumInt = Integer.parseInt(residueNum);
                } catch (NumberFormatException nfE) {
                    Matcher matcher = resPatter.matcher(residueNum);
                    if (matcher.matches()) {
                        residueNum = matcher.group(1);
                        residueNumInt = Integer.parseInt(residueNum);
                    }
                }
            }
            double refIntensity = 1.0;
            double refError = 0.0;
            List<Double> xValueList = new ArrayList<>();
            List<Double> yValueList = new ArrayList<>();
            List<Double> errValueList = new ArrayList<>();
            List<String> peakRefList = new ArrayList<>();

            if (peak.getMeasures().isPresent()) {
                double[][] v = peak.getMeasures().get();
                for (int i = 0; i < v[0].length; i++) {
                    xValueList.add(xValues[i]);
                    yValueList.add(v[0][i]);
                    errValueList.add(v[1][i]);
                    peakRefList.add(String.valueOf(i));

                }
            }
            if (expMode.equals("cest")) {
                processCESTData(expData, residueNum, xValueList, yValueList, errValueList, peakNum);
            } else {
                DynamicsSource dynSource = DynamicsSource.createFromPeak(peak);
                System.out.println(dynSource.toString());
                ResidueData residueData = new ResidueData(expData, dynSource, xValueList, yValueList, errValueList);
                expData.addResidueData(residueNum, residueData);
            }
            boolean ok = true;

            ResidueInfo residueInfo = resProp.getResidueInfo(residueNum);

            if (residueInfo == null) {
                residueInfo = new ResidueInfo(resProp, residueNumInt, 0, 0, 0);
                residueInfo.setResName(resName);
                resProp.addResidueInfo(residueNum, residueInfo);
            } else {
                residueInfo.setResName(resName);
            }
            if (expMode.equals("noe")) {
                residueInfo.value = yValueList.get(0);
                residueInfo.err = errValueList.get(0);
            }
        });
        if (!eSet) {
            double errValue = estimateErrors(expData);
            setErrors(expData, errValue);
        }

    }

    public static void loadPeakFile(String fileName, ExperimentData expData,
            ResidueProperties resProp, XCONV xConv, YCONV yConv)
            throws IOException, IllegalArgumentException { //(String fileName, ResidueProperties resProp, String nucleus,
//            double temperature, double field, double tau, double[] vcpmgs, String expMode,
//            HashMap<String, Object> errorPars, double[] delayCalc) throws IOException, IllegalArgumentException {
        System.out.println("load peak file");
//        Path path = Paths.get(fileName);
//        String fileTail = path.getFileName().toString();
//        fileTail = fileTail.substring(0, fileTail.indexOf('.'));

//        ExperimentData expData = new ExperimentData(fileTail, nucleus, field, temperature, tau, xvals, expMode, errorPars, delayCalc, B1field);
//        String fileName = expData.getName();
        Path path = Paths.get(fileName);
        if (Files.notExists(path, LinkOption.NOFOLLOW_LINKS)) {
            throw new FileNotFoundException(fileName);
        }
        HashMap<String, Object> errorPars = expData.getErrorPars();
        String expMode = expData.getExpMode();
        double[] xVals = expData.getXVals();
        double[] delayCalc = expData.getDelayCalc();
        double tau = expData.getTau();
        resProp.addExperimentData(expData.getName(), expData);
        boolean eSet = false;
        double errF = 0.05;
        double noise = 1.0;
        String errorMode = "";
        if (errorPars != null) {
            if (errorPars.containsKey("mode")) {
                errorMode = errorPars.get("mode").toString();
                if (errorMode.equals("percent")) {
                    if (errorPars.containsKey("value")) {
                        String percentValue = errorPars.get("value").toString();
                        errF = Double.parseDouble(percentValue) * 0.01;
                        eSet = true;
                    }
                } else if (errorMode.equals("noise")) {
                    if (errorPars.containsKey("value")) {
                        String noiseValueStr = errorPars.get("value").toString();
                        noise = Double.parseDouble(noiseValueStr);
                        eSet = true;
                    }
                }
            }
        }

        boolean gotHeader = false;
        boolean hasErrColumns = false;
        String[] peakRefs = null;
        double[] xValues = null;
        int offset = 0;
        int residueField = -1;
        int peakField = -1;
        String header = "";
        List<String> peakRefList = new ArrayList<>();
        //  Peak       Residue N       T1      T2      T11     T3      T4      T9      T5      T10     T12     T6      T7      T8
        int fakeRes = 1;
        try (BufferedReader fileReader = Files.newBufferedReader(path)) {
            while (true) {
                String line = fileReader.readLine();
                if (line == null) {
                    break;
                }
//            String line = fileReader.readLine();
                String sline = line.trim();
                if (sline.length() == 0) {
                    continue;
                }
                if (sline.charAt(0) == '#') {
                    continue;
                }
                String[] sfields = line.split("\t", -1);
                if (!gotHeader) {
                    int nfields = sfields.length;
                    if (fileName.endsWith(".mpk2")) {
                        // find last field that starts with "lab"
                        //   .mpk2 files have peak labels in columns like "lab1", "lab2"
                        //   intensities start in next column

                        for (int i = nfields - 1; i >= 0; i--) {
                            if (sfields[i].startsWith("lab")) {
                                offset = i + 1;
                                break;
                            }
                        }
                        // find first lab1,lab2... field to get residue number from
                        for (int i = 0; i < nfields; i++) {
                            if (sfields[i].startsWith("lab")) {
                                residueField = i;
                                break;
                            }
                        }
                        // find id label to get peak number
                        for (int i = 0; i < nfields; i++) {
                            if (sfields[i].startsWith("id")) {
                                peakField = i;
                                break;
                            }
                        }
                        // check to see if file has err values
                        for (int i = 0; i < nfields; i++) {
                            if (sfields[i].startsWith("err")) {
                                hasErrColumns = true;
                                break;
                            }
                        }
                    } else {
                        for (int i = 0; i < nfields; i++) {
                            if (sfields[i].startsWith("Res")) {
                                residueField = i;
                                offset = i + 1;
                                break;
                            }
                        }
                    }
                    if (expMode.equals("cest") || expMode.equals("cpmg") || expMode.equals("noe")) {
                        offset++;
                        if (hasErrColumns) {
                            offset++;
                        }
                    }
                    int nValues = nfields - offset;
                    if (hasErrColumns) {
                        nValues /= 2;
                    }
                    xValues = new double[nValues];
                    peakRefs = new String[nValues];
                    System.out.println("off " + offset + " " + nfields + " " + hasErrColumns);
                    int iX = 0;
                    for (int i = offset; i < nfields; i++) {
                        // fixme assumes first vcpmg is the 0 ref   
                        // fixme. need to explicitly account for alternating x-value, errorHeader fields
                        if (xVals == null) {
                            try {
                                double x = Double.parseDouble(sfields[i].trim());
                                if (hasErrColumns) {
                                    i++;
                                }
                                xValues[iX] = xConv.convert(x, delayCalc, expData);

                            } catch (NumberFormatException nFE) {
                            }
                        } else {
                            xValues[iX] = xVals[iX];
                        }
                        peakRefs[iX] = String.valueOf(iX);
                        peakRefList.add(peakRefs[iX]);
                        iX++;
                    }
                    gotHeader = true;
                    header = sline;
                } else {
                    String residueNum = "";
                    int peakNum = -1;
                    if (residueField != -1) {
                        residueNum = sfields[residueField].trim();
                    }
                    if (peakField != -1) {
                        peakNum = (int) Double.parseDouble(sfields[peakField].trim());
                    }
                    if (residueNum.equals("")) {
                        residueNum = String.valueOf(fakeRes);
                    } else if (residueNum.indexOf('.') != -1) {
                        int dotIndex = residueNum.indexOf('.');
                        residueNum = residueNum.substring(0, dotIndex);
                    }
                    int residueNumInt = 0;
                    if (residueNum.length() > 0) {
                        try {
                            residueNumInt = Integer.parseInt(residueNum);
                        } catch (NumberFormatException nfE) {
                            Matcher matcher = resPatter.matcher(residueNum);
                            if (matcher.matches()) {
                                residueNum = matcher.group(1);
                                residueNumInt = Integer.parseInt(residueNum);
                            }
                        }
                    }
                    double refIntensity = 1.0;
                    double refError = 0.0;
                    if (expMode.equals("cest") || expMode.equals("cpmg") || expMode.equals("noe")) {
                        if (hasErrColumns) {
                            refIntensity = Double.parseDouble(sfields[offset - 2].trim());
                            refError = Double.parseDouble(sfields[offset - 1].trim());
                        } else {
                            refIntensity = Double.parseDouble(sfields[offset - 1].trim());
                        }

                    }
                    List<Double> xValueList = new ArrayList<>();
                    List<Double> yValueList = new ArrayList<>();
                    List<Double> errValueList = new ArrayList<>();
                    boolean ok = true;
                    int iX = 0;
                    for (int i = offset; i < sfields.length; i++) {
                        String valueStr = sfields[i].trim();
                        if ((valueStr.length() == 0) || valueStr.equals("NA")) {
                            continue;
                        }
                        Double yValue;
                        double intensity;
                        try {
                            intensity = Double.parseDouble(sfields[i].trim());
                            yValue = yConv.convert(intensity, refIntensity, tau);
                        } catch (NumberFormatException nFE) {
                            ok = false;
                            continue;
                        }
                        if ((yValue == null) || Double.isNaN(yValue) || Double.isInfinite(yValue)) {
                            ok = false;
                            continue;
                        }
                        xValueList.add(xValues[iX]);
                        yValueList.add(yValue);

                        double eValue = 0.0;
                        if (hasErrColumns) {
                            i++;
                        }
                        if (hasErrColumns && errorMode.equals("measured")) {
                            eValue = Double.parseDouble(sfields[i].trim());
                            eSet = true;
                            if (expMode.equals("cpmg")) {
                                eValue = eValue / (intensity * tau);
                            } else if (expMode.equals("noe") && (yConv == YCONV.NORMALIZE)) {
                                double r1 = refError / refIntensity;
                                double r2 = eValue / intensity;
                                eValue = Math.abs(yValue) * Math.sqrt(r1 * r1 + r2 * r2);
                            }
                        } else if (expMode.equals("cpmg")) {
                            if (errorMode.equals("noise")) {
                                eValue = noise / (intensity * tau);
                            } else if (errorMode.equals("percent")) {
                                eValue = (errF * refIntensity) / (intensity * tau);
                            }
                        } else if (expMode.equals("noe") && (yConv == YCONV.NORMALIZE)) {
                            double r1 = noise / refIntensity;
                            double r2 = noise / intensity;
                            eValue = Math.abs(yValue) * Math.sqrt(r1 * r1 + r2 * r2);
                        } else {
                            if (errorMode.equals("percent")) {
                                eValue = Math.abs(yValue) * errF;
                            } else if (errorMode.equals("noise")) {
                                eValue = noise;
                            }
                        }
                        errValueList.add(eValue);
                        iX++;
                    }
                    if (!ok) {
                        continue;
                    }
                    String resName = "";
                    if (expMode.equals("cest")) {
                        processCESTData(expData, residueNum, xValueList, yValueList, errValueList, peakNum);
                    } else {
                        DynamicsSource dynSource = DynamicsSource.createFromSpecifiers(expMode + "." + peakNum, residueNum, "H", "N");
                        ResidueData residueData = new ResidueData(expData, dynSource, xValueList, yValueList, errValueList);
                        expData.addResidueData(residueNum, residueData);
                    }
                    ResidueInfo residueInfo = resProp.getResidueInfo(residueNum);

                    if (residueInfo == null) {
                        residueInfo = new ResidueInfo(resProp, residueNumInt, 0, 0, 0);
                        residueInfo.setResName(resName);
                        resProp.addResidueInfo(residueNum, residueInfo);
                    } else {
                        residueInfo.setResName(resName);
                    }
                    if (expMode.equals("noe")) {
                        residueInfo.value = yValueList.get(0);
                        residueInfo.err = errValueList.get(0);
                    }

                    fakeRes++;
                }
            }
        }
        if (!eSet) {
            double errValue = estimateErrors(expData);
            setErrors(expData, errValue);
        }
    }

    public static void processCESTData(ExperimentData expData, String residueNum,
            List<Double> xValueList, List<Double> yValueList, List<Double> errValueList, int peakNum) {
        Double B1field = expData.getB1Field();
        List<Double> B1fieldList = new ArrayList<>();
        for (int i = 0; i < xValueList.size(); i++) {
            B1fieldList.add(B1field);
        }
        double tau = expData.getTau();
        List<Double> tauList = new ArrayList<>();
        for (int i = 0; i < xValueList.size(); i++) {
            tauList.add(tau);
        }
        List<Double> bFieldUniqueValue = new ArrayList<>();
        bFieldUniqueValue.add(B1fieldList.get(0));
        List<Double> tauList1 = new ArrayList<>();
        tauList1.add(tauList.get(0));
        List<Double>[] xValueLists = new ArrayList[3];
        xValueLists[0] = xValueList;
        xValueLists[1] = B1fieldList;
        xValueLists[2] = tauList;
        DynamicsSource dynSource = DynamicsSource.createFromSpecifiers("cest." + peakNum, residueNum, "H", "N");

        ResidueData residueData = new ResidueData(expData, dynSource, xValueLists, yValueList, errValueList);
        expData.addResidueData(residueNum, residueData);
        expData.getExtras().clear();
        expData.setExtras(bFieldUniqueValue);
        expData.setExtras(tauList1);

    }

    /**
     * Get the residue name from a molecule, if present.
     *
     * @param residueNum int. The residue number.
     * @return name String. The name of the corresponding residue in the
     * molecule.
     */
    public static String getResidueName(int residueNum) {
        MoleculeBase mol = MoleculeFactory.getActive();
        String name = "";
        if (mol != null) {
            Iterator entityIterator = mol.entityLabels.values().iterator();
            while (entityIterator.hasNext()) {
                Entity entity = (Entity) entityIterator.next();
                if (entity instanceof Polymer) {
                    List<Residue> resList = ((Polymer) entity).getResidues();
                    for (Residue res : resList) {
                        int num = res.getIDNum();
                        if (num == residueNum) {
                            name = res.getName();
                            break;
                        }
                    }
                }
            }
        }
        return name;
    }

    public static void loadResidueDataFile(String fileName, ExperimentData expData,
            String residueNum, ResidueProperties resProp, String nucleus,
            double temperature, double field,
            HashMap<String, Object> errorPars, XCONV xConv, YCONV yConv, double refIntensity)
            throws IOException, IllegalArgumentException {
        boolean gotHeader = false;
        int nValues = 0;
        List<Double> xValueList = new ArrayList<>();
        List<Double> yValueList = new ArrayList<>();
        List<Double> errValueList = new ArrayList<>();
        System.out.println("Load XY file " + fileName);

        resProp.addExperimentData(expData.getName(), expData);
        String splitPattern = "\t";
        String[] sfields;

        try (BufferedReader fileReader = new BufferedReader(new FileReader(fileName))) {
            while (true) {
                String line = fileReader.readLine();
                if (line == null) {
                    break;
                }

//            String line = fileReader.readLine();
                String sline = line.trim();
                if (sline.length() == 0) {
                    continue;
                }
                if (sline.charAt(0) == '#') {
                    continue;
                }
                if (!gotHeader) {
                    if (!line.contains("\t")) {
                        splitPattern = ",";
                        if (!line.contains(",")) {
                            splitPattern = " +";
                        }
                    }
                    sfields = line.split(splitPattern, -1);
                    nValues = sfields.length;
                    gotHeader = true;
                } else {
                    sfields = line.split(splitPattern, -1);
                    try {

                        double offsetFreq = Double.parseDouble(sfields[0].trim());
                        double intensity = Double.parseDouble(sfields[1].trim());
                        double error = 0.01;
                        if (nValues > 2) {
                            error = Double.parseDouble(sfields[2].trim());
                        }
                        offsetFreq = xConv.convert(offsetFreq, null, expData);
                        intensity = yConv.convert(intensity, refIntensity, 0.0);
                        xValueList.add(offsetFreq);
                        yValueList.add(intensity);
                        errValueList.add(error);
                    } catch (NumberFormatException nFE) {
                        System.out.println(nFE.getMessage());
                        continue;
                    }
                }
            }
        } catch (IOException ex) {
            System.out.println(ex.getMessage());

            Logger.getLogger(DataIO.class.getName()).log(Level.SEVERE, null, ex);
        }
        processCESTData(expData, residueNum, xValueList, yValueList, errValueList, 0);
        ResidueInfo residueInfo = resProp.getResidueInfo(residueNum);
        if (residueInfo == null) {
            residueInfo = new ResidueInfo(resProp, Integer.parseInt(residueNum), 0, 0, 0);
            resProp.addResidueInfo(residueNum, residueInfo);
        }

    }

    public static void loadTextFile(String fileName, ResidueProperties resProp,
            String nucleus, double temperature, double field, XCONV xConv, String expMode)
            throws IOException, IllegalArgumentException {
        Path path = Paths.get(fileName);
        String fileTail = path.getFileName().toString();
        fileTail = fileTail.substring(0, fileTail.indexOf('.'));

//        ExperimentData expData = new ExperimentData(fileTail, nucleus, field, temperature);
        ExperimentData expData = new ExperimentData(fileTail, nucleus, field, temperature, null, null, expMode, null, null, null);
        resProp.addExperimentData(fileTail, expData);
        boolean gotHeader = false;
        String[] peakRefs = null;
        double[][] xValues = null;
//        List<Double> xValues = new ArrayList<>();
        int peakNum = 0;
        try (BufferedReader fileReader = Files.newBufferedReader(path)) {
            while (true) {
                String line = fileReader.readLine();
                if (line == null) {
                    break;
                }
//            String line = fileReader.readLine();
                String sline = line.trim();
                if (sline.length() == 0) {
                    continue;
                }
                String[] sfields = line.split("\t", -1);
                if (!gotHeader) {
                    int j = 0;
                    int nValues = (sfields.length - 1) / 2;
                    xValues = new double[1][nValues];
                    peakRefs = new String[nValues];
                    for (int i = 1; i < sfields.length - 1; i += 2) {
                        double xValue = Double.parseDouble(sfields[i].trim());
                        xValues[0][j] = xConv.convert(xValue, null, expData);
                        peakRefs[j] = String.valueOf(j);
                        j++;
                    }
                    gotHeader = true;
                } else {
                    int nValues = (sfields.length - 1) / 2;
                    String residueNum = sfields[0].trim();
                    double[] yValues = new double[nValues];
                    double[] errValues = new double[nValues];
                    int j = 0;
                    for (int i = 1; i < sfields.length - 1; i += 2) {
                        double r2Eff = Double.parseDouble(sfields[i].trim());
                        double r2EffErr = Double.parseDouble(sfields[i + 1].trim());
                        yValues[j] = r2Eff;
                        errValues[j] = r2EffErr;
                        j++;
                    }

                    DynamicsSource dynSource = DynamicsSource.createFromSpecifiers(expMode + "." + peakNum, residueNum, "H", "N");
                    ResidueData residueData = new ResidueData(expData, dynSource, xValues, yValues, errValues);
                    expData.addResidueData(residueNum, residueData);

                    ResidueInfo residueInfo = resProp.getResidueInfo(residueNum);
                    if (residueInfo == null) {
                        residueInfo = new ResidueInfo(resProp, Integer.parseInt(residueNum), 0, 0, 0);
                        resProp.addResidueInfo(residueNum, residueInfo);
                    }
                    peakNum++;
                }
            }
        }
    }

    public static void setPercentileErrors(ExperimentData expData, double fraction) {
        for (ResidueData residueData : expData.residueData.values()) {
            double[] yValues = residueData.getYValues();
            double[] errValues = residueData.getErrValues();
            for (int i = 0; i < yValues.length; i++) {
                errValues[i] = yValues[i] * fraction;
            }
        }
    }

    public static void setErrors(ExperimentData expData, double error) {
        for (ResidueData residueData : expData.residueData.values()) {
            double[] errValues = residueData.getErrValues();
            Arrays.fill(errValues, error);
        }

    }

    public static double estimateErrors(ExperimentData expData) {
        int nDups = 0;
        double sumDelta2 = 0.0;
        double sumAbs = 0.0;
        for (ResidueData residueData : expData.residueData.values()) {
            double[][] xValues = residueData.getXValues();
            double[] yValues = residueData.getYValues();
            for (int i = 0; i < yValues.length - 1; i++) {
                for (int j = (i + 1); j < yValues.length; j++) {
                    if (xValues[0][i] == xValues[0][j]) {
                        double delta = yValues[i] - yValues[j];
                        sumDelta2 += delta * delta;
                        sumAbs += Math.abs(delta);
                        nDups++;
                    }

                }
            }
        }
        double errorA = 0.0;
        double error2 = 0.0;
        if (nDups > 0) {
            error2 = Math.sqrt(sumDelta2 / (2.0 * nDups));
            errorA = Math.sqrt(Math.PI / 2.0) * sumAbs / (2.0 * nDups);
        }
        System.out.println("data " + expData.name + " errors " + error2 + "errorA " + errorA + " ndup " + nDups);
        return error2;
    }

    public static ResidueProperties loadResultsFile(String fitMode, String fileName) throws IOException {
        Path path = Paths.get(fileName);
        String fileTail = path.getFileName().toString();
        fileTail = fileTail.substring(0, fileTail.indexOf('.'));
        ResidueProperties resProp = new ResidueProperties(fileTail, fileName);
        File file = new File(fileName);
        if (!file.exists()) {
            return resProp;
        }

        /*
Residue	 Peak	GrpSz	Group	Equation	   RMS	   AIC	Best	     R2	  R2.sd	    Rex	 Rex.sd	    Kex	 Kex.sd	     pA	  pA.sd	     dPPM	  dPPM.sd
36	    26	    1	    0	    NOEX	  1.28	 49.95		   9.22	   0.09								
36	    26	    1	    0	CPMGFAST	  0.25	  7.44	best	   8.88	   0.09	   2.94	   0.08	 259.33	  17.97				
36	    26	    1	    0	CPMGSLOW	  0.28	 14.05		   8.89	   0.09			 164.36	  53.17	   0.96	   0.14	  24.51	   6.76
        
         */
        boolean gotHeader = false;
        String[] header = null;
        HashMap<String, Integer> headerMap = new HashMap<>();
        try (BufferedReader fileReader = Files.newBufferedReader(path)) {
            while (true) {
                String line = fileReader.readLine();
                if (line == null) {
                    break;
                }
//            String line = fileReader.readLine();
                String sline = line.trim();
                if (sline.length() == 0) {
                    continue;
                }
                if (sline.startsWith("#")) {
                    continue;
                }
                if (!gotHeader) {
                    header = line.split("\t");
                    for (int i = 0; i < header.length; i++) {
                        headerMap.put(header[i].trim(), i);
                    }
                    gotHeader = true;
                    continue;
                }

                String[] sfields = line.split("\t", -1);
                String[] stdHeaders = {"Best", "Residue", "Equation", "RMS", "GrpSz", "Group", "Peak", "State"};
                for (String stdHeader : stdHeaders) {
                    int index = headerMap.get(stdHeader);
                    if (index == -1) {
                        throw new IOException("Missing header " + stdHeader);
                    }
                }
                String bestValue = sfields[headerMap.get("Best")];
                String residueNumber = sfields[headerMap.get("Residue")].trim();
                String equationName = sfields[headerMap.get("Equation")].trim().toUpperCase();
                String rmsString = sfields[headerMap.get("RMS")].trim();
                double rms = Double.parseDouble(rmsString);
                String groupSzStr = sfields[headerMap.get("GrpSz")].trim();
                int groupSize = Integer.parseInt(groupSzStr);
                String groupIdStr = sfields[headerMap.get("Group")].trim();
                int groupId = Integer.parseInt(groupIdStr);
                String peakIdStr = sfields[headerMap.get("Peak")].trim();
                int peakNum = Integer.parseInt(peakIdStr);
                String stateString = sfields[headerMap.get("State")].trim();
                ResidueInfo residueInfo = resProp.getResidueInfo(residueNumber);

                if (residueInfo == null) {
                    residueInfo = new ResidueInfo(resProp, Integer.parseInt(residueNumber), groupId, groupSize, peakNum);
                    resProp.addResidueInfo(residueNumber, residueInfo);
                }
                double[] fields = new double[1];
                fields[0] = 1.0;
                int parStart = headerMap.get("Best") + 1;
                HashMap<String, Double> parMap = new HashMap<>();
                EquationType equationType = null;
                List<String> equationNames;
                switch (fitMode) {
                    case "exp":
                        equationType = ExpEquation.valueOf(equationName);
                        equationNames = ExpFitter.getEquationNames();
                        break;
                    case "cest":
                        equationName = equationName.replace("CEST", "");
                        if (equationName.startsWith("R1RHO")) {
                            equationName = equationName.replace("R1RHO", "");
                            switch (equationName) {
                                case "EXACT1":
                                    equationName = "EIGENEXACT1";
                                    break;
                                case "N":
                                    equationName = "LAGUERRE";
                                    break;
                                case "PERTURBATION":
                                    equationName = "TROTT_PALMER";
                                    break;
                                case "PERTURBATIONNOEX":
                                    equationName = "NOEX";
                                    break;
                                default:
                                    break;
                            }
                        }
                        equationType = CESTEquation.valueOf(equationName);
                        equationNames = CESTFitter.getEquationNames();
                        break;
                    case "r1rho":
                        equationName = equationName.replace("R1RHO", "");
                        if (equationName.equals("PERTURBATION")) {
                            equationName = "TROTT_PALMER";
                        } else if (equationName.equals("PERTURBATIONNOEX")) {
                            equationName = "NOEX";
                        }
                        equationType = R1RhoEquation.valueOf(equationName);
                        equationNames = R1RhoFitter.getEquationNames();
                        break;
                    default:
                        equationType = CPMGEquation.valueOf(equationName);
                        equationNames = CPMGFitter.getEquationNames();
                        break;
                }
                String[] eqnParNames = equationType.getParNames();
                int nPar = eqnParNames.length;

                double[] pars = new double[nPar];
                double[] errs = new double[nPar];
                int iPar = 0;
                for (String parName : eqnParNames) {
                    int index = headerMap.get(parName);
                    pars[iPar] = Double.parseDouble(sfields[index].trim());
                    index = headerMap.get(parName + ".sd");
                    errs[iPar++] = Double.parseDouble(sfields[index].trim());
                }
                //parMap.put("Kex", 0.0); // Set Kex to 0.0.  Overwritten below if exists in output file

                for (int i = parStart; i < header.length; i++) {
                    if (sfields[i].trim().length() > 0) {
                        String parName = header[i].trim();
                        double parValue = Double.parseDouble(sfields[i].trim());
                        parMap.put(parName, parValue);
                    }
                }
                String[] extraPars = {"AIC", "RMS"};
                for (String extraPar : extraPars) {
                    int parIndex = headerMap.get(extraPar);
                    if (parIndex != -1) {
                        parMap.put(extraPar, Double.parseDouble(sfields[parIndex].trim()));
                    }
                }
                parMap.put("Equation", 1.0 + equationNames.indexOf(equationName));
                PlotEquation plotEquation = new PlotEquation(fitMode, equationName, pars, errs, fields);
                CurveFit curveSet = new CurveFit(stateString, residueNumber, parMap, plotEquation);
                residueInfo.addCurveSet(curveSet, bestValue.equals("best"));
            }

        }
        return resProp;
    }

    static void getFitParameters(ResidueProperties resProp, Map<String, Object> dataMap2) {
        HashMap<String, Object> fitParMap = (HashMap<String, Object>) dataMap2.get("parameters");
        if (fitParMap != null) {
            Boolean absValueMode = (Boolean) fitParMap.get("absValue");
            if (absValueMode != null) {
                resProp.setAbsValueMode(absValueMode);
            }
            String bootStrapMode = (String) fitParMap.get("bootStrap");
            if (bootStrapMode != null) {
                resProp.setBootStrapMode(bootStrapMode);
            }

            System.out.println("absmode " + absValueMode + " bootstrap " + bootStrapMode);
        }
    }

    static void processYAMLDataSections(ResidueProperties resProp, Map<String, Object> dataMap2, Path dirPath, String expMode) throws IOException {

        ArrayList<HashMap<String, Object>> dataList = (ArrayList<HashMap<String, Object>>) dataMap2.get("data");
        DataIO.processYAMLDataSections(resProp, dirPath, expMode, dataList);
    }

    static Map<String, List<Double>> getConstraints(Map... maps) {

        Map<String, List<Double>> constraintMap = null;
        for (Map map : maps) {
            if (map.containsKey("constraints")) {
                constraintMap = (Map<String, List<Double>>) map.get("constraints");
                if (constraintMap != null) {
                    System.out.println(constraintMap.toString());
                }
            }
        }
        return constraintMap;
    }

    static Double getDoubleValue(Map<String, Object> dataMap, String key, Double defaultValue) {
        Double value = defaultValue;
        if (dataMap.containsKey(key)) {
            Object oValue = dataMap.get(key);
            if (oValue instanceof Number) {
                value = ((Number) oValue).doubleValue();
            }
        }
        return value;
    }

    public static ResidueProperties processPeakList(PeakList peakList) {
        ResidueProperties resProp = null;
        if (peakList != null) {
            String peakListName = peakList.getName();
            List<Number> vcpmgList = null;
            DatasetBase dataset = DatasetBase.getDataset(peakList.fileName);
            String expMode = "exp";
            resProp = new ResidueProperties(peakListName, peakListName);
            expMode = expMode.toLowerCase();
            resProp.setExpMode(expMode);
            if (vcpmgList == null) {
                String nucleus = dataset.getNucleus(1).getName();
                double B0field = dataset.getSf(1);
                double temperature = dataset.getTempK();
                Double tau = null;
                Double B1field = null;
                double[] delayCalc = {0.0, 0.0, 1.0};
                HashMap<String, Object> errorPars = new HashMap<>();
                ExperimentData expData = new ExperimentData(peakListName,
                        nucleus, B0field, temperature, tau, null, expMode,
                        errorPars, delayCalc, B1field);

                loadFromPeakList(peakList, expData, resProp, XCONV.IDENTITY, YCONV.IDENTITY);
            }
        }
        return resProp;

    }

    public static void processYAMLDataSections(ResidueProperties resProp, Path dirPath, String expMode, ArrayList<HashMap<String, Object>> dataList) throws IOException {
        for (HashMap<String, Object> dataMap3 : dataList) {
            Double temperature = getDoubleValue(dataMap3, "temperature", null);
            if (temperature != null) {
                temperature += 273.15;
            } else {
                temperature = getDoubleValue(dataMap3, "temperatureK", null);
            }

            XCONV xConv;
            YCONV yConv;
            switch (expMode) {
                case "cpmg":
                    xConv = XCONV.TAU2;
                    yConv = YCONV.RATE;
                    break;
                case "cest":
                    xConv = XCONV.IDENTITY;
                    yConv = YCONV.NORMALIZE;
                    break;
                default:
                    xConv = XCONV.IDENTITY;
                    yConv = YCONV.IDENTITY;
            }
            if (dataMap3.containsKey("xconv")) {
                String xConvStr = dataMap3.get("xconv").toString();
                xConv = XCONV.valueOf(xConvStr.toUpperCase());
                if (xConv == null) {
                    throw new IOException("Bad xconversion type");
                }
            }
            if (dataMap3.containsKey("yconv")) {
                String xConvStr = dataMap3.get("yconv").toString();
                yConv = YCONV.valueOf(xConvStr.toUpperCase());
                if (yConv == null) {
                    throw new IOException("Bad yconversion type");
                }
            }

            Double B0field = ((Number) dataMap3.get("B0")).doubleValue();
            String nucleus = (String) dataMap3.get("nucleus");
            List<Number> vcpmgList = (List<Number>) dataMap3.get("vcpmg");
            Double tau = getDoubleValue(dataMap3, "tau", 1.0);
            Double B1field = getDoubleValue(dataMap3, "B1", null);

            String fileFormat = (String) dataMap3.get("format");

//            String dataFileName = (String) dataMap3.get("file");
//            File file = new File(dataFileName).getAbsoluteFile();
//            dataFileName = file.getName();
//            String fileTail = dataFileName.substring(0, dataFileName.indexOf('.'));
//            String textFileName = FileSystems.getDefault().getPath(dirPath.toString(), dataFileName).toString();
            HashMap<String, Object> errorPars = (HashMap<String, Object>) dataMap3.get("error");
            Object delayField = dataMap3.get("delays");
            System.out.println("delays " + delayField);
            double[] delayCalc = {0.0, 0.0, 1.0};
            if (delayField instanceof Map) {
                Map<String, Object> delayMap = (Map<String, Object>) delayField;
                delayCalc[0] = getDoubleValue(delayMap, "delta0", 0.0);
                delayCalc[1] = getDoubleValue(delayMap, "c0", 0.0);
                delayCalc[2] = getDoubleValue(delayMap, "delta", 1.0);
            }
            System.out.println("err " + errorPars);

            if ((fileFormat != null) && fileFormat.equals("mpk2")) {

                String dataFileName = (String) dataMap3.get("file");
                File dataFile = new File(dataFileName);
                String fileTail = dataFile.getName();
                int dot = fileTail.lastIndexOf(".");
                fileTail = fileTail.substring(0, dot);
                if (!dataFile.isAbsolute()) {
                    dataFileName = FileSystems.getDefault().getPath(dirPath.toString(), dataFileName).toString();
                }
                if (vcpmgList == null) {
                    ExperimentData expData = new ExperimentData(fileTail,
                            nucleus, B0field, temperature, tau, null, expMode,
                            errorPars, delayCalc, B1field);
//                    loadPeakFile(textFileName, resProp, nucleus, temperature, field, tau, null, expMode, errorPars, delayCalc);
                    loadPeakFile(dataFileName, expData, resProp, xConv, yConv);
                } else {
                    double[] vcpmgs = new double[vcpmgList.size()];
                    for (int i = 0; i < vcpmgs.length; i++) {
                        vcpmgs[i] = vcpmgList.get(i).doubleValue();
                    }
                    ExperimentData expData;
                    try {
                        expData = new ExperimentData(fileTail,
                                nucleus, B0field, temperature, tau, vcpmgs, expMode,
                                errorPars, delayCalc, B1field);
                    } catch (Exception ex) {
                        ex.printStackTrace();
                        return;
                    }
//                    loadPeakFile(textFileName, resProp, nucleus, temperature, field, tau, vcpmgs, expMode, errorPars, delayCalc);
                    loadPeakFile(dataFileName, expData, resProp, xConv, yConv);
                }

            } else if ((fileFormat != null) && fileFormat.equals("ires")) {
                List<Map<String, Object>> filesMaps = (List<Map<String, Object>>) dataMap3.get("files");
                String expName = (String) dataMap3.get("name").toString();
                ExperimentData expData = new ExperimentData(expName,
                        nucleus, B0field, temperature, tau, null, expMode,
                        errorPars, delayCalc, B1field);
                for (Map<String, Object> filesMap : filesMaps) {
                    Map<String, List<Double>> constraintMap = getConstraints(dataMap3, filesMap);
                    expData.setConstraints(constraintMap);
                    String residueNum = filesMap.get("residue").toString();
                    String dataFileName = (String) filesMap.get("file");
                    double refIntensity = 1.0;
                    if (filesMap.containsKey("refIntensity")) {
                        refIntensity = ((Number) filesMap.get("refIntensity")).doubleValue();
                    }
                    File file = new File(dataFileName).getAbsoluteFile();
                    dataFileName = file.getName();
                    String fileTail = dataFileName.substring(0, dataFileName.indexOf('.'));
                    String textFileName = FileSystems.getDefault().getPath(dirPath.toString(), dataFileName).toString();
                    loadResidueDataFile(textFileName, expData, residueNum, resProp, nucleus,
                            temperature, B0field, errorPars, xConv, yConv, refIntensity);
                }
            } else if (vcpmgList == null) {
                String dataFileName = (String) dataMap3.get("file");
                File file = new File(dataFileName).getAbsoluteFile();
                dataFileName = file.getName();
                String textFileName = FileSystems.getDefault().getPath(dirPath.toString(), dataFileName).toString();
                loadTextFile(textFileName, resProp, nucleus, temperature, B0field, xConv, expMode);
            } else {
                String dataFileName = (String) dataMap3.get("file");
                String textFileName = FileSystems.getDefault().getPath(dirPath.toString(), dataFileName).toString();
                double[] vcpmgs = new double[vcpmgList.size()];
                for (int i = 0; i < vcpmgs.length; i++) {
                    vcpmgs[i] = vcpmgList.get(i).doubleValue();
                }
                String fileTail = dataFileName;
                int dot = dataFileName.lastIndexOf(".");
                fileTail = dataFileName.substring(0, dot);

                ExperimentData expData = new ExperimentData(fileTail,
                        nucleus, B0field, temperature, tau, vcpmgs, expMode,
                        errorPars, delayCalc, B1field);
                loadPeakFile(textFileName, expData, resProp, xConv, yConv);
            }
        }

    }

    public static ResidueProperties loadYAMLFile(String fileName) throws FileNotFoundException, IOException {
        File yamlFile = new File(fileName).getAbsoluteFile();
        ResidueProperties resProp = null;
        try (InputStream input = new FileInputStream(yamlFile)) {
            Path path = yamlFile.toPath();
            Path dirPath = path.getParent();

            Yaml yaml = new Yaml();
            for (Object data : yaml.loadAll(input)) {
                Map dataMap = (HashMap<String, Object>) data;
                Map dataMap2 = (HashMap<String, Object>) dataMap.get("fit");
                if (dataMap2 != null) {
                    String expMode = (String) dataMap2.get("mode");
                    expMode = expMode == null ? "cpmg" : expMode;
                    String parName = (String) dataMap2.get("file");
                    if (parName == null) {
                        String yamlName = yamlFile.getName();
                        parName = yamlName.substring(0, yamlName.length() - 5) + "_out.txt";
                    }
                    String parFileName = FileSystems.getDefault().getPath(dirPath.toString(), parName).toString();
                    resProp = DataIO.loadResultsFile(expMode, parFileName);

                    resProp.setExpMode(expMode);
                    getFitParameters(resProp, dataMap2);
                    processYAMLDataSections(resProp, dataMap2, dirPath, expMode);
                }
            }
        }
        //        residueProperties.put(resProp.getName(), resProp);
        //        residueProperties.put(resProp.getName(), resProp);
        //        residueProperties.put(resProp.getName(), resProp);
        //        residueProperties.put(resProp.getName(), resProp);

        return resProp;

    }

    public static void saveResultsFile(String fileName, ResidueProperties resProp, boolean saveStats) {
        String[] headerFields = {"Residue", "Peak", "GrpSz", "Group", "State", "Equation", "RMS", "AIC", "Best"};
        String[] headerFields2 = {"Residue", "Peak", "GrpSz", "Group", "State", "RefineOpt", "RefineTime",
            "BootstrapOpt", "BootstrapTime", "Samples", "AbsMode", "NonParametricMode", "StartRadius", "FinalRadius",
            "Tolerance", "Weight", "Equation", "RMS", "AIC", "Best"};
        if (saveStats) {
            headerFields = headerFields2;
        }
        StringBuilder headerBuilder = new StringBuilder();
        String[] cpmgFields = {"R2", "dPPMmin", "Kex", "pA", "dPPM"};
        String[] expFields = {"A", "R"};
        String[] cestFields = {"Kex", "Pb", "deltaA0", "deltaB0", "R1A", "R1B", "R2A", "R2B"};
        String[] parFields;
        if (resProp.getExpMode().equals("cpmg")) {
            parFields = cpmgFields;
        } else if (resProp.getExpMode().equals("cest") || resProp.getExpMode().equals("r1rho")) {
            parFields = cestFields;
        } else {
            parFields = expFields;
        }

        for (String field : headerFields) {
            if (headerBuilder.length() > 0) {
                headerBuilder.append('\t');
            }
            headerBuilder.append(field);
        }
        for (String field : parFields) {
            headerBuilder.append('\t').append(field).append('\t').append(field).append(".sd");
        }

        try (FileWriter writer = new FileWriter(fileName)) {
            writer.write(headerBuilder.toString());
            resProp.getResidueValues().stream().
                    sorted((a, b) -> Integer.compare(a.getResNum(), b.getResNum())).
                    forEach(resInfo -> {
                        String outputLine = resInfo.toOutputString(parFields, saveStats);
                        try {
                            writer.write(outputLine);
                        } catch (IOException ex) {
                            Logger.getLogger(DataIO.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    });
        } catch (IOException ex) {
            Logger.getLogger(DataIO.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public static List<Double>[] loadSimDataFile(File file) throws IOException, IllegalArgumentException {
        Path path = file.toPath();
        List<Double> xValues = new ArrayList<>();
        List<Double> yValues = new ArrayList<>();
        List<Double> errValues = new ArrayList<>();
        try (BufferedReader fileReader = Files.newBufferedReader(path)) {
            while (true) {
                String line = fileReader.readLine();
                if (line == null) {
                    break;
                }

                String sline = line.trim();
                if (sline.length() == 0) {
                    continue;
                }
                if (sline.contains("#")) {
                    continue;
                }
                String[] sfields = line.split("\t", -1);

                xValues.add(Double.parseDouble(sfields[0]));
                yValues.add(Double.parseDouble(sfields[1]));
                if (sfields.length > 2) {
                    errValues.add(Double.parseDouble(sfields[2]));
                }
            }
        }
        List[] result = {xValues, yValues, errValues};
        return result;
    }

    /**
     * Load molecular information from a file.
     *
     * @param fileName String. The file to read.
     * @param type String. The filetype: pdb, star, nef, or cif.
     * @throws ParseException
     * @throws MoleculeIOException
     * @throws IOException
     */
    public static void readMoleculeFile(String fileName, String type) throws ParseException, MoleculeIOException, IOException {
        switch (type) {
            case "pdb":
                PDBFile pdb = new PDBFile();
                MoleculeBase mol = pdb.read(fileName);
                mol.updateAtomArray();
                break;
            case "star":
                File starFile = new File(fileName).getAbsoluteFile();
                NMRStarReader.read(starFile);
                break;
            case "nef":
                NMRNEFReader.read(fileName);
                break;
            case "cif":
                File cifFile = new File(fileName).getAbsoluteFile();
                MMcifReader.read(cifFile);
                break;
            default:
                break;
        }
        System.out.println("loaded molecule " + MoleculeFactory.getActive().getName() + " from " + fileName);

    }

    /**
     * Add the T1/T2 fit results to a map, and to a molecule, if present.
     *
     * @param resProp ResidueProperties. The residue properties of the fit.
     * @param expType String. The experiment type, T1 or T2.
     */
    public static void addT1T2FitResults(ResidueProperties resProp, String expType) {
        String frameName = resProp.getName();
        ExperimentData expData = resProp.getExperimentData(resProp.getName());
        String nucName = expData.getNucleusName();
        String expName = "EXPAB";
        double[] fields = resProp.getFields();
        List<Object> expResults = new ArrayList<>();
        expResults.add(frameName);
        expResults.add(nucName);
        expResults.add(fields);
//        Collections.sort(fields, (a, b) -> Double.compare(a, b));
        List<Object> states = new ArrayList<>();
        //Nested Map for fit results: <resNum, <field, <parName, parValue>>>
        Map<Integer, Map<Integer, Map<String, Double>>> allFitResults = new TreeMap<>();
        for (int f = 0; f < fields.length; f++) {
            int field = (int) fields[f];
            List<ResidueInfo> resInfoList = resProp.getResidueValues();
            for (ResidueInfo resInfo : resInfoList) {
                Map<Integer, Map<String, Double>> fitResFieldMap;
                Map<String, CurveFit> parMap = resInfo.curveSets.get(expName);
                if (resInfoList.indexOf(resInfo) == 0) {
                    states = Arrays.asList(parMap.keySet().toArray());
                    Collections.sort(states, (a, b) -> a.toString().compareTo(b.toString()));
                }
                CurveFit curveFit = parMap.get(states.get(f).toString());
                Map<String, Double> fitPars = curveFit.getParMap();
//                    System.out.println(resInfo.resNum + " " + polymer.getResidue(resInfo.resNum - 1));
                if (!allFitResults.containsKey(resInfo.resNum)) {
                    fitResFieldMap = new HashMap<>();
                } else {
                    fitResFieldMap = allFitResults.get(resInfo.resNum);
                }
                fitResFieldMap.put(field, fitPars);
                allFitResults.put(resInfo.resNum, fitResFieldMap);
            }
        }
        expResults.add(allFitResults);
        expFitResults.put(expType, expResults);
        System.out.println(expType + " fit results added to map");

        MoleculeBase mol = MoleculeFactory.getActive();
        if (mol != null) {
            mol.setProperty(expType + "frameName", frameName);
            mol.setProperty(expType + "nucName", nucName);
            Iterator entityIterator = mol.entityLabels.values().iterator();
            List<Integer> addedRes = new ArrayList<>();
            int iPolyRes = 0;
            while (entityIterator.hasNext()) {
                Entity entity = (Entity) entityIterator.next();
                TreeSet<String> expTypes = (TreeSet<String>) entity.getPropertyObject("expTypes");
                if (expTypes == null) {
                    expTypes = new TreeSet<>();
                }
                expTypes.add(expType);
                entity.setPropertyObject("expTypes", expTypes);
                entity.setPropertyObject(expType + "fields", fields);
                List<ResidueInfo> resInfoList = resProp.getResidueValues();
                Collections.sort(resInfoList, (a, b) -> Integer.compare(a.resNum, b.resNum));
                if (entity instanceof Polymer) {
                    iPolyRes = 0;
                }
                for (ResidueInfo resInfo : resInfoList) {
                    if (addedRes.contains(resInfo.resNum)) {
                        continue;
                    }
                    if (entity instanceof Polymer) {
                        Polymer polymer = (Polymer) entity;
                        if (resInfo.resNum <= Integer.valueOf(polymer.getLastResidue().getNumber())) {
                            polymer.getResidue(iPolyRes).setPropertyObject(expType + "fitResults", allFitResults.get(resInfo.resNum));
                            addedRes.add(resInfo.resNum);
                            iPolyRes++;
                        }
                    } else if (entity instanceof Compound) {
                        Compound compound = (Compound) entity;
                        compound.setNumber(String.valueOf(resInfo.resNum));
                        compound.setPropertyObject(expType + String.valueOf(resInfo.resNum) + "fitResults", allFitResults.get(resInfo.resNum));
                        addedRes.add(resInfo.resNum);
                        break;
                    }
                }
            }
            System.out.println(expType + " fit results added to molecule " + mol.getName());
        }
    }

    public static void writeSTAR3File(String fileName, ResidueProperties resProp) throws IOException, ParseException, InvalidMoleculeException, InvalidPeakException {
        try (FileWriter writer = new FileWriter(fileName)) {
            writeSTAR3File(writer, resProp);
        }
    }

    public static void writeSTAR3File(File file, ResidueProperties resProp) throws IOException, ParseException, InvalidMoleculeException, InvalidPeakException {
        try (FileWriter writer = new FileWriter(file)) {
            writeSTAR3File(writer, resProp);
        }
    }

    public static void writeSTAR3File(FileWriter chan, ResidueProperties resProp) throws IOException, ParseException, InvalidMoleculeException, InvalidPeakException {
        if (MoleculeFactory.getActive() != null) {
            NMRStarWriter.writeAll(chan);
        } else {
            NMRStarWriter.writeAll(chan, expFitResults);
        }
    }
}
