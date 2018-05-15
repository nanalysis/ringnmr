package org.comdnmr.cpmgfit2.calc;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Bruce Johnson
 */
public class ExperimentData {

    String name;
    HashMap<String, ResidueData> residueData;
    final double field;
    final double temperature;
    double errFraction = 0.05;
    final String nucleus;
    List<Double> extras = new ArrayList<>();
    private String state = "";

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

    public void setErrFraction(double errFraction) {
        this.errFraction = errFraction;
    }

    public double getErrFraction() {
        return (errFraction);
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

    /**
     * @return the state
     */
    public String getState() {
        return state;
    }

    /**
     * @param state the state to set
     */
    public void setState(String state) {
        this.state = state;
    }

    public void setExtras(List extravals) {
        this.extras.addAll(extravals);
    }

    public List<Double> getExtras() {
        return extras;
    }

}
