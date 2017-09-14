package org.comdnmr.cpmgfit2.calc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author Bruce Johnson
 */
public class ResidueProperties {

    Map<String, ExperimentData> expMaps = new HashMap<>();
    HashMap<String, ResidueInfo> residueMap = new HashMap<>();
    final String name;
    String fileName = null;

    public ResidueProperties(String name, String fileName) {
        this.name = name;
        this.fileName = fileName;
    }

    public String getName() {
        return name;
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
        List<String> equationNames = new ArrayList<>();
       // fixme  equationNames.addAll(residueMap.values().stream().findFirst().get().curveSets.parMaps.keySet());
        return equationNames;
    }

    public void addResidueInfo(String resNum, ResidueInfo value) {
        residueMap.put(resNum, value);
    }

    public HashMap<String, ResidueInfo> getResidueMap() {
        return residueMap;
    }

    public ResidueInfo getResidueInfo(String resNum) {
        return residueMap.get(resNum);
    }

    public double[] getFields() {
        Set<Double> fieldSet = new HashSet<>();
        for (ExperimentData expData : expMaps.values()) {
            fieldSet.add(expData.getField());
        }
        double[] fields = new double[fieldSet.size()];
        int i = 0;
        for (Double field : fieldSet) {
            fields[i++] = field;
        }
        Arrays.sort(fields);
        return fields;
    }

}
