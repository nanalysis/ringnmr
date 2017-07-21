/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.cpmgfit2.calc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author Bruce Johnson
 */
public class ResidueInfo {

    int resNum;
    HashMap<String, PlotEquation> plotEquations = new HashMap<>();
    PlotEquation bestEquation = null;
    ResidueData resData = null;
    double[] extraValues;
    String[] peakRefs;
    int groupSize = 1;
    int groupId = 0;
    int peakNum = 0;
    HashMap<String, HashMap<String, Double>> parMaps = new HashMap<>();
    static final String[] parNames = {"R2", "Rex", "Kex", "pA", "dW"};

    public class DataValue {

        int index;
        ResidueInfo resInfo;

        public DataValue(ResidueInfo resInfo, int index) {
            this.resInfo = resInfo;
            this.index = index;
        }

        public double getX() {
            return resInfo.resData.getXValues()[index];
        }

        public double getY() {
            return resInfo.resData.getYValues()[index];
        }

        public String getPeak() {
            String peak = "";
            if (resInfo.peakRefs != null) {
                peak = resInfo.peakRefs[index];
            }
            return peak;
        }
    }

    public class ParValue implements ParValueInterface {

        String parName;
        ResidueInfo resInfo;
        String equationName;

        public ParValue(ResidueInfo resInfo, String parName, String equationName) {
            this.parName = parName;
            this.resInfo = resInfo;
            this.parName = parName;
            this.equationName = equationName;
        }

        public String getName() {
            return parName;
        }

        public double getValue(String equationName) {
            Double value = resInfo.parMaps.get(equationName).get(parName);
            if (value == null) {
                return 0.0;
            } else {
                return value;
            }
        }

        public double getValue() {
            return getValue(equationName);
        }

        public double getError(String equationName) {
            Double value = resInfo.parMaps.get(equationName).get(parName + ".sd");
            if (value == null) {
                return 0.0;
            } else {
                return value;
            }
        }

        public double getError() {
            return getError(equationName);
        }
    }

    public ResidueInfo(int resNum, int groupId, int groupSize) {
        this(resNum, groupId, groupSize, 0);
    }

    public ResidueInfo(int resNum, int groupId, int groupSize, int peakNum) {
        this.resNum = resNum;
        this.groupId = groupId;
        this.groupSize = groupSize;
        this.peakNum = peakNum;
    }

    public void addEquation(PlotEquation equation, boolean best) {
        plotEquations.put(equation.name, equation);
        if (best) {
            bestEquation = equation;
        }
    }

    public String getBestEquationName() {
        return bestEquation == null ? "" : bestEquation.name;
    }

    public int getResNum() {
        return resNum;
    }

    public PlotEquation getEquation(String name) {
        return plotEquations.get(name);
    }

    public void addParMap(HashMap<String, Double> parMap, String equationName) {
        parMaps.put(equationName, parMap);
    }

    public Double getParValue(String equationName, String parName) {
        return parMaps.get(equationName).get(parName);
    }

    public List<ParValueInterface> getParValues(String equationName) {
        final String useEquationName;
        if (equationName.equals("best")) {
            useEquationName = bestEquation.name;
        } else {
            useEquationName = equationName;
        }
        List<ParValueInterface> dataValues = new ArrayList<>();
        HashMap<String, Double> parMap = parMaps.get(useEquationName);
        parMap.keySet().stream().sorted().filter((parName) -> (parMap.containsKey(parName + ".sd"))).forEachOrdered((parName) -> {
            dataValues.add(new ParValue(this, parName, useEquationName));
        });
        return dataValues;
    }

    ArrayList<DataValue> getDataValues() {
        ArrayList<DataValue> dataValues = new ArrayList<>();
        double[] xValues = resData.getXValues();

        for (int i = 0; i < xValues.length; i++) {
            dataValues.add(new DataValue(this, i));
        }
        return dataValues;
    }

    public String toString() {
        StringBuilder sBuilder = new StringBuilder();
        PlotEquation plotEquation = bestEquation;
        sBuilder.append(plotEquation.name);
        sBuilder.append('\n');
        double[] xValues = null;
        double[] yValues = null;
        if (resData != null) {
            xValues = resData.getXValues();
            yValues = resData.getYValues();
        }
        if (xValues != null) {
            for (double value : xValues) {
                sBuilder.append(value);
                sBuilder.append(" ");
            }

            sBuilder.append('\n');
        }
        if (yValues != null) {
            for (double value : yValues) {
                sBuilder.append(value);
                sBuilder.append(" ");
            }
        }
        if (parMaps != null) {
            sBuilder.append('\n');
            sBuilder.append(parMaps.toString());
        }
        return sBuilder.toString();
    }

    //         String header = "Residue	Peak	GrpSz	Group	Equation	   RMS	   AIC	Best	     R2	  R2.sd	    Rex	 Rex.sd	    Kex	 Kex.sd	     pA	  pA.sd	     dW	  dW.sd";
    public String toOutputString() {
        char sep = '\t';
        StringBuilder sBuilder = new StringBuilder();
        sBuilder.append(resNum).append(sep);
        sBuilder.append(peakNum).append(sep);
        sBuilder.append(groupSize).append(sep);
        sBuilder.append(groupId).append(sep);
        String commonString = sBuilder.toString();
        StringBuilder allBuilder = new StringBuilder();
        Double parValue;
        for (PlotEquation plotEquation : plotEquations.values()) {
            sBuilder.setLength(0);
            sBuilder.append('\n');
            sBuilder.append(commonString);
            sBuilder.append(plotEquation.name).append(sep);
            Map<String, Double> parMap = parMaps.get(plotEquation.name);
            parValue = parMap.get("RMS");
            sBuilder.append(String.format("%.2f", parValue)).append(sep);
            parValue = parMap.get("AIC");
            sBuilder.append(String.format("%.1f", parValue)).append(sep);
            if (bestEquation == plotEquation) {
                sBuilder.append("best");
            }

            for (String parName : parNames) {
                sBuilder.append(sep);
                parValue = parMap.get(parName);
                if (parValue != null) {
                    sBuilder.append(String.format("%.3f", parValue));
                }
                sBuilder.append(sep);
                parValue = parMap.get(parName + ".sd");
                if (parValue != null) {
                    sBuilder.append(String.format("%.3f", parValue));
                }
            }
            allBuilder.append(sBuilder.toString());

        }

        return allBuilder.toString();

    }
}
