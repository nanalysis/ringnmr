/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

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

    ResidueProperties resProps;

    int resNum;
    Map<String, Map<String, CurveFit>> curveSets = new LinkedHashMap<>();
    Map<String, CPMGFitResult> fitResults = new LinkedHashMap<>();
    PlotEquation bestEquation = null;
    double[] extraValues;
    String[] peakRefs;
    int groupSize = 1;
    int groupId = 0;
    int peakNum = 0;

    public class ParValue implements ParValueInterface {

        String parName;
        ResidueInfo resInfo;
        String equationName;
        String state;

        public ParValue(ResidueInfo resInfo, String parName, String equationName, String state) {
            this.parName = parName;
            this.resInfo = resInfo;
            this.parName = parName;
            this.equationName = equationName;
            this.state = state;
        }

        @Override
        public String getName() {
            return parName;
        }

        public double getValue(String equationName, String state) {
            Double value = resInfo.curveSets.get(equationName).get(state).parMap.get(parName);
            if (value == null) {
                return 0.0;
            } else {
                return value;
            }
        }

        @Override
        public double getValue() {
            Double value = resInfo.curveSets.get(equationName).get(state).parMap.get(parName);
            return value;
        }

        public double getError(String key) {
            Double value = resInfo.curveSets.get(equationName).get(state).parMap.get(parName + ".sd");
            if (value == null) {
                return 0.0;
            } else {
                return value;
            }
        }

        @Override
        public double getError() {
            return getError(equationName);
        }

        @Override
        public String getResidue() {
            return String.valueOf(resInfo.resNum);
        }

        @Override
        public String getState() {
            return state;
        }
    }

    public ResidueInfo(ResidueProperties resProps, int resNum, int groupId, int groupSize) {
        this(resProps, resNum, groupId, groupSize, 0);
    }

    public ResidueInfo(ResidueProperties resProps, int resNum, int groupId, int groupSize, int peakNum) {
        this.resProps = resProps;
        this.resNum = resNum;
        this.groupId = groupId;
        this.groupSize = groupSize;
        this.peakNum = peakNum;
    }

    public void addFitResult(CPMGFitResult fitResult) {
        fitResults.put(fitResult.getEquationName(), fitResult);
    }

    public CPMGFitResult getFitResult(String equationName) {
        String useEquationName;
        if (equationName.startsWith("best")) {
            useEquationName = bestEquation.name;
        } else {
            useEquationName = equationName;
        }
        return fitResults.get(useEquationName);
    }

    public Collection<CurveFit> getCurveSets(String equationName) {
        return curveSets.get(equationName).values();
    }

    public void addCurveSet(CurveFit curveSet, boolean best) {
        Map<String, CurveFit> fitMap = curveSets.get(curveSet.plotEquation.getName());
        if (fitMap == null) {
            fitMap = new HashMap<>();
            curveSets.put(curveSet.plotEquation.getName(), fitMap);
        }
        fitMap.put(curveSet.state, curveSet);
        if (best) {
            bestEquation = curveSet.plotEquation;
        }
    }

    public CurveFit getCurveSet(String equationName, String state) {
        CurveFit curveFit = null;
        if (curveSets != null) {
            Map<String, CurveFit> curveFits = curveSets.get(equationName);
            if (curveFits != null) {
                curveFit = curveSets.get(equationName).get(state);
            }
        }
        return curveFit;
    }

    public String getBestEquationName() {
        return bestEquation == null ? "" : bestEquation.name;
    }

    public void setBestEquationName(String equationName) {
        Map<String, CurveFit> curveFits = curveSets.get(equationName);
        for (CurveFit curveFit : curveFits.values()) {
            if (curveFit.plotEquation.getName().equals(equationName)) {
                bestEquation = curveFit.plotEquation;
            }
        }
    }

    public int getResNum() {
        return resNum;
    }

    public Double getParValue(String equationName, String state, String parName) {
        Map<String, CurveFit> curveFits = curveSets.get(equationName);
        if (curveFits == null) {
            return null;
        }
        CurveFit curveFit = curveFits.get(state);
        if (curveFit == null) {
            return null;
        } else {
            return curveFit.parMap.get(parName);
        }
    }

    public List<ParValueInterface> getParValues(String equationName, String state) {
        final String useEquationName;
        if (equationName.startsWith("best")) {
            useEquationName = bestEquation.name;
        } else {
            useEquationName = equationName;
        }
        List<ParValueInterface> dataValues = new ArrayList<>();
        Map<String, CurveFit> curveFits = curveSets.get(useEquationName);
        if (curveFits != null) {
            curveFits.values().stream().forEach(cf -> {
                if (ResidueProperties.matchStateString(state, cf.getState())) {
                    Map<String, Double> parMap = cf.getParMap();
                    parMap.keySet().stream().sorted().filter((parName) -> (parMap.containsKey(parName + ".sd"))).forEachOrdered((parName) -> {
                        dataValues.add(new ParValue(this, parName, useEquationName, cf.getState()));
                    });
                }

            });
        }
        return dataValues;
    }

    @Override
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
        for (Map<String, CurveFit> fitMap : curveSets.values()) {
            for (CurveFit curveFit : fitMap.values()) {
                if (curveFit.parMap != null) {
                    sBuilder.append('\n');
                    sBuilder.append(curveFit.parMap.toString());
                }
            }
        }
        return sBuilder.toString();
    }

    //         String header = "Residue	Peak	GrpSz	Group	Equation	   RMS	   AIC	Best	     R2	  R2.sd	    Rex	 Rex.sd	    Kex	 Kex.sd	     pA	  pA.sd	     dW	  dW.sd";
    public String toOutputString(String[] parNames, boolean saveStats) {
        char sep = '\t';
        StringBuilder sBuilder = new StringBuilder();
        sBuilder.append(resNum).append(sep);
        sBuilder.append(peakNum).append(sep);
        sBuilder.append(groupSize).append(sep);
        sBuilder.append(groupId).append(sep);
        String commonString = sBuilder.toString();
        StringBuilder allBuilder = new StringBuilder();
        Double parValue;
        for (Map<String, CurveFit> fitMap : curveSets.values()) {
            for (CurveFit curveFit : fitMap.values()) {
                sBuilder.setLength(0);
                sBuilder.append('\n');
                sBuilder.append(commonString);
                PlotEquation plotEquation = curveFit.plotEquation;
                if (saveStats) {
                    CurveFit.CurveFitStats curveStats = getFitResult(plotEquation.name).getCurveFitStats();
                    String statString = curveStats.toString();
                    sBuilder.append(statString);
                }
                sBuilder.append(curveFit.state).append(sep);// fixme
                sBuilder.append(plotEquation.name).append(sep);
                Map<String, Double> parMap = curveFit.getParMap();
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

        }

        return allBuilder.toString();

    }
}
