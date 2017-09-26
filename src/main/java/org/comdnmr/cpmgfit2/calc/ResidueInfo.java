/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.cpmgfit2.calc;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author Bruce Johnson
 */
public class ResidueInfo {

    int resNum;
    HashMap<String, CurveFit> curveSets = new LinkedHashMap<>();
    PlotEquation bestEquation = null;
    double[] extraValues;
    String[] peakRefs;
    int groupSize = 1;
    int groupId = 0;
    int peakNum = 0;
    static final String[] parNames = {"R2", "Rex", "Kex", "pA", "dW"};

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

        public double getValue(String key) {
            Double value = resInfo.curveSets.get(equationName).parMap.get(parName);
            if (value == null) {
                return 0.0;
            } else {
                return value;
            }
        }

        public double getValue() {
            return getValue(equationName);
        }

        public double getError(String key) {
            Double value = resInfo.curveSets.get(key).parMap.get(parName + ".sd");
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

    public Collection<CurveFit> getCurveSets() {
        return curveSets.values();
    }

    public void addCurveSet(CurveFit curveSet, boolean best) {
        String key = curveSet.plotEquation.getName() + "." + curveSet.state;
        curveSets.put(key, curveSet);
        if (best) {
            bestEquation = curveSet.plotEquation;
        }
    }

    public String checkKey(String key) {
        if (key.startsWith("best")) {
            if (key.indexOf('.') == -1) {
                key = bestEquation.name + ".0:0:0:0";
            } else {
                key = bestEquation.name + key.substring(4);
            }
        } else {
            if (key.indexOf('.') == -1) {
                key = key + ".0:0:0:0";
            }
        }
        return key;
    }

    public CurveFit getCurveSet(String key) {
        return curveSets.get(key);
    }

    public String getBestEquationName() {
        return bestEquation == null ? "" : bestEquation.name;
    }

    public int getResNum() {
        return resNum;
    }

    public Double getParValue(String key, String parName) {
        System.out.println("key 1: " + key + " " + bestEquation.name + " " + parName);
        key = checkKey(key);
        System.out.println("key 2: " + key);
        CurveFit curveFit = curveSets.get(key);
        if (curveFit == null) {
            for (String keyName : curveSets.keySet()) {
                System.out.println("keyname " + keyName);
            }
            return null;
        } else {
            return curveFit.parMap.get(parName);
        }
    }

    public List<ParValueInterface> getParValues(String key) {
        System.out.println("key 1: " + key + " " + bestEquation.name);
        key = checkKey(key);
        System.out.println("key 2: " + key);
        // fixme
        final String useEquationName;
        if (key.startsWith("best")) {
            useEquationName = bestEquation.name;
        } else {
            useEquationName = key;
        }
        List<ParValueInterface> dataValues = new ArrayList<>();
        Map<String, Double> parMap = curveSets.get(key).getParMap();
        parMap.keySet().stream().sorted().filter((parName) -> (parMap.containsKey(parName + ".sd"))).forEachOrdered((parName) -> {
            dataValues.add(new ParValue(this, parName, useEquationName));
        });
        return dataValues;
    }

    public String toString() {
        StringBuilder sBuilder = new StringBuilder();
        PlotEquation plotEquation = bestEquation;
        sBuilder.append(plotEquation.name);
        sBuilder.append(" ");
        sBuilder.append(resNum);
        sBuilder.append('\n');
        double[] xValues = null;
        double[] yValues = null;
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
        for (CurveFit curveSet : curveSets.values()) {
            if (curveSet.parMap != null) {
                sBuilder.append('\n');
                sBuilder.append(curveSet.parMap.toString());
            }
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
        for (CurveFit curveSet : curveSets.values()) {
            sBuilder.setLength(0);
            sBuilder.append('\n');
            sBuilder.append(commonString);
            PlotEquation plotEquation = curveSet.plotEquation;
            sBuilder.append(curveSet.state).append(sep);// fixme
            sBuilder.append(plotEquation.name).append(sep);
            Map<String, Double> parMap = curveSet.getParMap();
            parValue = parMap.get("RMS");
            sBuilder.append(String.format("%.2f", parValue)).append(sep);
            parValue = parMap.get("AIC");
            sBuilder.append(String.format("%.1f", parValue)).append(sep);
            if (bestEquation.getName().equals(plotEquation.getName())) {
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
