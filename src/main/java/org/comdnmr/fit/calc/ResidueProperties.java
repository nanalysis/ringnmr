package org.comdnmr.fit.calc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 *
 * @author Bruce Johnson
 */
public class ResidueProperties {

    Map<String, ExperimentData> expMaps = new HashMap<>();
    private HashMap<String, ResidueInfo> residueMap = new HashMap<>();
    final String name;
    String fileName = null;
    Map<Double, Integer> fieldMap = new LinkedHashMap();
    Map<Double, Integer> tempMap = new LinkedHashMap();
    Map<String, Integer> nucMap = new LinkedHashMap();
    List<Double> fieldList = new ArrayList<>();
    List<Double> tempList = new ArrayList<>();
    List<String> nucList = new ArrayList<>();
    private boolean absValueMode = false;
    private String bootStrapMode = "parametric";
    private String expMode = "cpmg";

    public ResidueProperties(String name, String fileName) {
        this.name = name;
        this.fileName = fileName;
    }

    @Override
    public String toString() {
        StringBuilder sBuilder = new StringBuilder();
        sBuilder.append(name).append(" ").append(fileName).append(" resmap ").append(residueMap.size());
        return sBuilder.toString();
    }

    public String getName() {
        return name;
    }

    /**
     * @return the absValueMode
     */
    public boolean isAbsValueMode() {
        return absValueMode;
    }

    /**
     * @param absValueMode the absValueMode to set
     */
    public void setAbsValueMode(boolean absValueMode) {
        this.absValueMode = absValueMode;
    }

    /**
     * @return the bootStrapMode
     */
    public String getBootStrapMode() {
        return bootStrapMode;
    }

    /**
     * @param bootStrapMode the bootStrapMode to set
     */
    public void setBootStrapMode(String bootStrapMode) {
        this.bootStrapMode = bootStrapMode;
    }

    /**
     * @return the expMode
     */
    public String getExpMode() {
        return expMode;
    }

    /**
     * @param expMode the expMode to set
     */
    public void setExpMode(String expMode) {
        this.expMode = expMode;
    }

    public void addExperimentData(String name, ExperimentData data) {
        expMaps.put(name, data);
    }

    public ExperimentData getExperimentData(String name) {
        return expMaps.get(name);
    }

    public Collection<ExperimentData> getExperimentData() {
        return expMaps.values();
    }

    public List<String> getEquationNames() {
        Set<String> equationSet = new HashSet<>();
        residueMap.values().stream().forEach(rInfo -> {
            equationSet.addAll(rInfo.curveSets.keySet());
        });
        List<String> equationNames = equationSet.stream().sorted().collect(Collectors.toList());
        return equationNames;
    }

    public static boolean matchStateString(String target, String state) {
        boolean match = true;
        String[] targetElems = target.split(":");
        String[] stateElems = state.split(":");
        for (int i = 0; i < targetElems.length; i++) {
            if (!targetElems[i].equals("*") && !targetElems[i].equals(stateElems[i])) {
                match = false;
                break;
            }
        }
        return match;
    }

    public List<String> getStateStrings() {
        if (fieldMap.isEmpty()) {
            setupMaps();
        }
        // List<String> stateStrings = 
        Set<String> states = expMaps.values().stream().map(e -> e.getState()).map(s
                -> "0:" + s.substring(2)).collect(Collectors.toSet());
        List<String> stateStrings = states.stream().sorted().collect(Collectors.toList());
        return stateStrings;
    }

    public static String getStateString(int[] state) {
        StringBuilder builder = new StringBuilder();
        builder.append(state[1]);
        for (int i = 2; i < state.length; i++) {
            builder.append(':');
            builder.append(state[i]);
        }
        return builder.toString();
    }

    public void setupMaps() {
        fieldMap.clear();
        tempMap.clear();
        nucMap.clear();
        for (ExperimentData expData : expMaps.values()) {
            if (!fieldMap.containsKey(Math.floor(expData.field))) {
                fieldMap.put(Math.floor(expData.field), fieldMap.size());
                fieldList.add(expData.field);
            }
            if (!tempMap.containsKey(Math.floor(expData.temperature))) {
                tempMap.put(Math.floor(expData.temperature), tempMap.size());
                tempList.add(expData.temperature);
            }
            if (!nucMap.containsKey(expData.nucleus)) {
                nucMap.put(expData.nucleus, nucMap.size());
                nucList.add(expData.nucleus);
            }
        }
        for (ExperimentData expData : expMaps.values()) {
            int[] state = getStateIndices(0, expData);
            expData.setState(getStateString(state));
        }
    }

    public int[] getStateIndices(int resIndex, ExperimentData expData) {
        if (fieldMap.size() == 0) {
            setupMaps();
        }
        int[] state = new int[4];
        state[0] = resIndex;
        state[1] = fieldMap.get(Math.floor(expData.field));
        state[2] = tempMap.get(Math.floor(expData.temperature));
        state[3] = nucMap.get(expData.nucleus);
//        System.out.println(resIndex + " " + expData.field + " " + expData.temperature + " " + expData.nucleus);
//        System.out.println("state index residue:" + state[0] + " field:" + state[1] + " temp:" + state[2] + " nuc:" + state[3]);

        return state;
    }

    int[] getStateCount(int nResidues) {
        if (fieldMap.size() == 0) {
            setupMaps();
        }
        int[] state = new int[4];
        state[0] = nResidues;
        state[1] = fieldMap.size();
        state[2] = tempMap.size();
        state[3] = nucMap.size();
//        System.out.println("state count residues:" + state[0] + " fields:" + state[1] + " temps:" + state[2] + " nucs:" + state[3]);
        return state;
    }

    public void addResidueInfo(String resNum, ResidueInfo value) {
        residueMap.put(resNum, value);
    }

    public void clearResidueMap() {
        residueMap.clear();
    }

    public List<ResidueInfo> getResidueValues() {
        List<ResidueInfo> values = new ArrayList<>();
        values.addAll(residueMap.values());
        return values;
    }

    public ResidueInfo getResidueInfo(String resNum) {
        return residueMap.get(resNum);
    }

    public double[] getFields() {
        double[] fields = new double[fieldList.size()];
        int i = 0;
        for (Double field : fieldList) {
            fields[i++] = field;
        }
        return fields;
    }

    public double[] getTemperatures() {
        double[] temperatures = new double[tempList.size()];
        int i = 0;
        for (Double temperature : tempList) {
            temperatures[i++] = temperature;
        }
        return temperatures;
    }

    public String[] getNuclei() {
        String[] nuclei = new String[nucList.size()];
        int i = 0;
        for (String nucleus : nucList) {
            nuclei[i++] = nucleus;
        }
        return nuclei;
    }

}
