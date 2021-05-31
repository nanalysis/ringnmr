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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.nmrfx.datasets.Nuclei;

/**
 *
 * @author brucejohnson
 */
public class Experiment {

    ExperimentSet experimentSet;
    String name;
    HashMap<String, ExperimentData> experimentalDataSets;
    protected final double B0field;
    protected final String nucleusName;
    protected final double nucleusField;
    protected final double temperature;
    final String expMode;
    HashMap<String, Object> errorPars;
    double errFraction = 0.05;
    List<Double> extras = new ArrayList<>();
    protected String state = "";
    Map<String, List<Double>> constraints = null;

    public Experiment(ExperimentSet experimentSet, String name, String nucleus, double B0field, double temperature,
            String expMode) {
        this.experimentSet = experimentSet;
        this.name = name;
        this.B0field = B0field;
        this.temperature = temperature;
        this.expMode = expMode;
        if (nucleus == null) {
            nucleus = "H1";
        }
        this.nucleusName = nucleus;
        double ratio = 1.0;
        Nuclei nuc = Nuclei.findNuclei(nucleusName);
        if (nuc != null) {
            ratio = nuc.getFreqRatio();
        }
        nucleusField = B0field * ratio;
        experimentalDataSets = new HashMap<>();
    }

    @Override
    public String toString() {
        StringBuilder sBuilder = new StringBuilder();
        sBuilder.append(name).append(" ").append(B0field).append(" ").append(" nres ").append(experimentalDataSets.size());
        return sBuilder.toString();
    }

    public String getName() {
        return name;
    }

    public double getErrFraction() {
        return errFraction;
    }

    public void setErrFraction(double errFraction) {
        this.errFraction = errFraction;
    }

    public double getTemperature() {
        return temperature;
    }

    public void addResidueData(String resNum, ExperimentData data) {
        experimentalDataSets.put(resNum, data);
    }

    public ExperimentData getResidueData(String resNum) {
        return experimentalDataSets.get(resNum);
    }

    public Set<String> getResidues() {
        return experimentalDataSets.keySet();
    }

    public double getB0Field() {
        return B0field;
    }

    public double getNucleusField() {
        return nucleusField;
    }

    public String getNucleusName() {
        return nucleusName;
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

    public List<Double> getExtras() {
        return extras;
    }

    public String getExpMode() {
        return expMode;
    }

    public HashMap<String, Object> getErrorPars() {
        return errorPars;
    }

    public double[] getConstraint(String key) {
        double[] result = null;
        if (constraints != null) {
            List<Double> values = constraints.get(key);
            if (values != null) {
                result = new double[2];
                result[0] = values.get(0);
                result[1] = values.get(1);
            }
        }
        return result;
    }

    public Map<String, List<Double>> getConstraints() {
        return constraints;
    }

    public void setConstraints(Map<String, List<Double>> constraints) {
        this.constraints = constraints;
    }

    public void setExtras(List extravals) {
        this.extras.addAll(extravals);
    }

}
