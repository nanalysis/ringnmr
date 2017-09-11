package org.comdnmr.cpmgfit2.calc;

import java.util.HashMap;
import java.util.Set;

/**
 *
 * @author Bruce Johnson
 */
public class ExperimentData {

    String name;
    HashMap<String, ResidueData> residueData;
    final double field;
    final double temperature;
    final String nucleus;

    public ExperimentData(String name, String nucleus, double field, double temperature) {
        this.name = name;
        this.field = field;
        this.temperature = temperature;
        if (nucleus == null) {
            nucleus = "H";
        }
        this.nucleus = nucleus;
        residueData = new HashMap<>();
    }

    public String getName() {
        return name;
    }

    public void addResidueData(String resNum, ResidueData data) {
        residueData.put(resNum, data);
    }

    public ResidueData getResidueData(String resNum) {
        return residueData.get(resNum);
    }

    public Set<String> getResidues() {
        return residueData.keySet();
    }

    public double getField() {
        return field;
    }

}
