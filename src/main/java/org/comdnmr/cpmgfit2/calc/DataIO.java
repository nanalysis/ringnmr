package org.comdnmr.cpmgfit2.calc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.yaml.snakeyaml.Yaml;

/**
 *
 * @author Bruce Johnson
 */
public class DataIO {

    public static void loadPeakFile(String fileName, ResidueProperties resProp, String nucleus, double temperature, double field, double tauCPMG, double[] vcpmgs) throws IOException, IllegalArgumentException {
        Path path = Paths.get(fileName);
        String fileTail = path.getFileName().toString();
        fileTail = fileTail.substring(0, fileTail.indexOf('.'));

        ExperimentData expData = new ExperimentData(fileTail, nucleus, field, temperature);
        resProp.addExperimentData(fileTail, expData);
        boolean gotHeader = false;
        String[] peakRefs = null;
        double[] xValues = null;
        List<String> peakRefList = new ArrayList<>();
        //  Peak       Residue N       T1      T2      T11     T3      T4      T9      T5      T10     T12     T6      T7      T8

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
                int offset = 4;
                if (!gotHeader) {
                    int nValues = sfields.length - offset;
                    xValues = new double[nValues];
                    peakRefs = new String[nValues];
                    for (int i = offset; i < sfields.length; i++) {
                        int j = i - offset;
                        // fixme assumes first vcpmg is the 0 ref 
                        xValues[j] = vcpmgs[j + 1];
                        peakRefs[j] = sfields[i];
                        peakRefList.add(peakRefs[j]);
                    }
                    gotHeader = true;
                } else {
                    String residueNum = sfields[1].trim();
                    double refIntensity = Double.parseDouble(sfields[offset - 1].trim());
                    List<Double> xValueList = new ArrayList<>();
                    List<Double> yValueList = new ArrayList<>();
                    List<Double> errValueList = new ArrayList<>();
                    boolean ok = true;
                    for (int i = offset; i < sfields.length; i++) {
                        int j = i - offset;
                        String valueStr = sfields[i].trim();
                        if ((valueStr.length() == 0) || valueStr.equals("NA")) {
                            continue;
                        }
                        double r2Eff;
                        try {
                            double intensity = Double.parseDouble(sfields[i].trim());
                            if (intensity > refIntensity) {
                                ok = false;
                                continue;
                            }
                            r2Eff = -Math.log(intensity / refIntensity) / tauCPMG;
                        } catch (NumberFormatException nFE) {
                            continue;
                        }
                        xValueList.add(xValues[j]);
                        yValueList.add(r2Eff);
                        errValueList.add(0.0);
                    }
                    if (!ok) {
                        continue;
                    }
                    ResidueData residueData = new ResidueData(expData, residueNum, xValueList, yValueList, errValueList, peakRefList);
                    expData.addResidueData(residueNum, residueData);
                    ResidueInfo residueInfo = resProp.getResidueInfo(residueNum);
                    if (residueInfo == null) {
                        residueInfo = new ResidueInfo(Integer.parseInt(residueNum), 0, 0, 0);
                        resProp.addResidueInfo(residueNum, residueInfo);
                    }

                }
            }
        }
        double errValue = estimateErrors(expData);
        setErrors(expData, errValue);
    }

    public static void loadTextFile(String fileName, ResidueProperties resProp, String nucleus, double temperature, double field) throws IOException, IllegalArgumentException {
        Path path = Paths.get(fileName);
        String fileTail = path.getFileName().toString();
        fileTail = fileTail.substring(0, fileTail.indexOf('.'));

        ExperimentData expData = new ExperimentData(fileTail, nucleus, field, temperature);
        resProp.addExperimentData(fileTail, expData);
        boolean gotHeader = false;
        String[] peakRefs = null;
        double[] xValues = null;
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
                    xValues = new double[nValues];
                    peakRefs = new String[nValues];
                    for (int i = 1; i < sfields.length - 1; i += 2) {
                        double tau = Double.parseDouble(sfields[i].trim());
                        double vcpmg = 1000.0 / (2.0 * tau);
                        xValues[j] = vcpmg;
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
                    ResidueData residueData = new ResidueData(expData, residueNum, xValues, yValues, errValues);
                    expData.addResidueData(residueNum, residueData);

                    ResidueInfo residueInfo = resProp.getResidueInfo(residueNum);
                    if (residueInfo == null) {
                        residueInfo = new ResidueInfo(Integer.parseInt(residueNum), 0, 0, 0);
                        resProp.addResidueInfo(residueNum, residueInfo);
                    }
                }
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
            double[] xValues = residueData.getXValues();
            double[] yValues = residueData.getYValues();
            for (int i = 0; i < xValues.length - 1; i++) {
                for (int j = (i + 1); j < xValues.length; j++) {
                    if (xValues[i] == xValues[j]) {
                        double delta = yValues[i] - yValues[j];
                        sumDelta2 += delta * delta;
                        sumAbs += Math.abs(delta);
                        nDups++;
                    }

                }
            }
        }
        double error2 = Math.sqrt(sumDelta2 / (2.0 * nDups));
        double errorA = Math.sqrt(Math.PI / 2.0) * sumAbs / (2.0 * nDups);
        System.out.println("data " + expData.name + " errors " + error2 + "errorA " + errorA + " ndup " + nDups);
        return errorA;
    }

    public static ResidueProperties loadParametersFromFile(String fileName) throws IOException {
        Path path = Paths.get(fileName);
        String fileTail = path.getFileName().toString();
        fileTail = fileTail.substring(0, fileTail.indexOf('.'));
        ResidueProperties resProp = new ResidueProperties(fileTail, fileName);
        File file = new File(fileName);
        if (!file.exists()) {
            return resProp;
        }

        /*
Residue	 Peak	GrpSz	Group	Equation	   RMS	   AIC	Best	     R2	  R2.sd	    Rex	 Rex.sd	    Kex	 Kex.sd	     pA	  pA.sd	     dW	  dW.sd
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
                    residueInfo = new ResidueInfo(Integer.parseInt(residueNumber), groupId, groupSize, peakNum);
                    resProp.addResidueInfo(residueNumber, residueInfo);
                }
                double[] fields = new double[1];
                fields[0] = 1.0;
                int parStart = headerMap.get("Best") + 1;
                HashMap<String, Double> parMap = new HashMap<>();
                CalcRDisp.CPMGEquation cpmgEquation = CalcRDisp.CPMGEquation.valueOf(equationName);
                String[] eqnParNames = cpmgEquation.getParNames();
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
                List<String> equationNames = Arrays.asList("NOEX", "CPMGFAST", "CPMGSLOW");
                parMap.put("Equation", 1.0 + equationNames.indexOf(equationName));
                PlotEquation plotEquation = new PlotEquation(equationName, pars, errs, fields);
                CurveFit curveSet = new CurveFit(stateString, residueNumber, parMap, plotEquation);
                residueInfo.addCurveSet(curveSet, bestValue.equals("best"));
            }

        }
        return resProp;
    }

    public static ResidueProperties loadParameters(String fileName) {
        File yamlFile = new File(fileName).getAbsoluteFile();
        ResidueProperties resProp = null;
        try (InputStream input = new FileInputStream(yamlFile)) {
            Path path = yamlFile.toPath();
            Path dirPath = path.getParent();

            Yaml yaml = new Yaml();
            int counter = 0;
            for (Object data : yaml.loadAll(input)) {
                Map dataMap = (HashMap<String, Object>) data;
                Map dataMap2 = (HashMap<String, Object>) dataMap.get("fit");
                String parName = (String) dataMap2.get("file");
                String parFileName = FileSystems.getDefault().getPath(dirPath.toString(), parName).toString();

                resProp = DataIO.loadParametersFromFile(parFileName);
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

                ArrayList<HashMap<String, Object>> dataList = (ArrayList<HashMap<String, Object>>) dataMap2.get("data");
                for (HashMap<String, Object> dataMap3 : dataList) {
                    String dataFileName = (String) dataMap3.get("file");
                    Double temperature = (Double) dataMap3.get("temperature");
                    Double field = (Double) dataMap3.get("field");
                    String nucleus = (String) dataMap3.get("nucleus");
                    List<Number> vcpmgList = (List<Number>) dataMap3.get("vcpmg");
                    Double tauCPMG = (Double) dataMap3.get("tau");
                    String textFileName = FileSystems.getDefault().getPath(dirPath.toString(), dataFileName).toString();
                    if (vcpmgList == null) {
                        loadTextFile(textFileName, resProp, nucleus, temperature, field);
                    } else {
                        double[] vcpmgs = new double[vcpmgList.size()];
                        for (int i = 0; i < vcpmgs.length; i++) {
                            vcpmgs[i] = vcpmgList.get(i).doubleValue();
                        }
                        loadPeakFile(textFileName, resProp, nucleus, temperature, field, tauCPMG, vcpmgs);
                    }
                }
                counter++;
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(DataIO.class.getName()).log(Level.SEVERE, null, ex);
            System.out.println(ex.getMessage());
            return null;
        } catch (IOException ex) {
            Logger.getLogger(DataIO.class.getName()).log(Level.SEVERE, null, ex);
            System.out.println(ex.getMessage());
            return null;
        }
        //        residueProperties.put(resProp.getName(), resProp);

        return resProp;

    }

    public static void saveParametersToFile(String fileName, ResidueProperties resProp) {
        String header = "Residue	Peak	GrpSz	Group	State	Equation	RMS	   AIC	Best	     R2	  R2.sd	    Rex	 Rex.sd	    Kex	 Kex.sd	     pA	  pA.sd	     dW	  dW.sd";
        try (FileWriter writer = new FileWriter(fileName)) {
            writer.write(header);
            resProp.getResidueMap().values().stream().
                    sorted((a, b) -> Integer.compare(a.getResNum(), b.getResNum())).
                    forEach(resInfo -> {
                        String outputLine = resInfo.toOutputString();
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

}
