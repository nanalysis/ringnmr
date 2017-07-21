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

    public static void loadTextFile(String fileName, ResidueProperties resProp, double temperature, double field) throws IOException, IllegalArgumentException {
        Path path = Paths.get(fileName);
        String fileTail = path.getFileName().toString();
        fileTail = fileTail.substring(0, fileTail.indexOf('.'));

        ExperimentData expData = new ExperimentData(fileTail, field, temperature);
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
                    ResidueData residueData = new ResidueData(xValues, yValues, errValues);
                    expData.addResidueData(residueNum, residueData);
                }

            }
        }
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
                int iPar = 0;
                for (String parName : eqnParNames) {
                    int index = headerMap.get(parName);
                    pars[iPar++] = Double.parseDouble(sfields[index].trim());
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
                PlotEquation plotEquation = new PlotEquation(equationName, pars, fields);
                residueInfo.addEquation(plotEquation, bestValue.equals("best"));
                residueInfo.addParMap(parMap, equationName);
            }

        }
        return resProp;
    }

    public static ResidueProperties loadParameters(String fileName) {
        File yamlFile = new File(fileName);
        ResidueProperties resProp = null;
        try (InputStream input = new FileInputStream(yamlFile)) {
            Path path = Paths.get(fileName);
            String fileTail = path.getFileName().toString();
            Path dirPath = path.getParent();

            Yaml yaml = new Yaml();
            int counter = 0;
            for (Object data : yaml.loadAll(input)) {
                Map dataMap = (HashMap<String, Object>) data;
                Map dataMap2 = (HashMap<String, Object>) dataMap.get("fit");
                String parName = (String) dataMap2.get("file");
                String parFileName = FileSystems.getDefault().getPath(dirPath.toString(), parName).toString();

                resProp = DataIO.loadParametersFromFile(parFileName);

                ArrayList<HashMap<String, Object>> dataList = (ArrayList<HashMap<String, Object>>) dataMap2.get("data");
                for (HashMap<String, Object> dataMap3 : dataList) {
                    String dataFileName = (String) dataMap3.get("file");
                    Double temperature = (Double) dataMap3.get("temperature");
                    Double field = (Double) dataMap3.get("field");
                    String textFileName = FileSystems.getDefault().getPath(dirPath.toString(), dataFileName).toString();
                    DataIO.loadTextFile(textFileName, resProp, temperature, field);
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
        String header = "Residue	Peak	GrpSz	Group	Equation	   RMS	   AIC	Best	     R2	  R2.sd	    Rex	 Rex.sd	    Kex	 Kex.sd	     pA	  pA.sd	     dW	  dW.sd";
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
