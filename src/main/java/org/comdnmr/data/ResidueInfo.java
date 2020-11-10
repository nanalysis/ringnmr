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
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.data;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.comdnmr.eqnfit.FitResult;
import org.comdnmr.eqnfit.CurveFit;
import org.comdnmr.eqnfit.ParValueInterface;
import org.comdnmr.eqnfit.PlotEquation;

/**
 *
 * @author Bruce Johnson
 */
public class ResidueInfo {

    ResidueProperties resProps;

    int resNum;
    String resName;
    Map<String, Map<String, CurveFit>> curveSets = new LinkedHashMap<>();
    Map<String, FitResult> fitResults = new LinkedHashMap<>();
    PlotEquation bestEquation = null;
    double[] extraValues;
    String[] peakRefs;
    int groupSize = 1;
    int groupId = 0;
    int peakNum = 0;
    Double value = null;
    Double err = null;

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
            Double value = resInfo.curveSets.get(equationName).get(state).getParMap().get(parName);
            if (value == null) {
                return 0.0;
            } else {
                return value;
            }
        }

        @Override
        public double getValue() {
            Double value = resInfo.curveSets.get(equationName).get(state).getParMap().get(parName);
            return value;
        }

        public double getError(String key) {
            Double value = resInfo.curveSets.get(equationName).get(state).getParMap().get(parName + ".sd");
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
        public String getResidueName() {
            return String.valueOf(resInfo.resName);
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

    public void addFitResult(FitResult fitResult) {
        if (fitResult != null) {
            fitResults.put(fitResult.getEquationName(), fitResult);
        }
    }

    public FitResult getFitResult(String equationName) {
        String useEquationName;
        if (equationName.startsWith("best")) {
            useEquationName = bestEquation.getName();
        } else {
            useEquationName = equationName;
        }
        return fitResults.get(useEquationName);
    }

    public Collection<CurveFit> getCurveSets(String equationName) {
        return curveSets.get(equationName).values();
    }

    public void addCurveSet(CurveFit curveSet, boolean best) {
        Map<String, CurveFit> fitMap = curveSets.get(curveSet.getEquation().getName());
        if (fitMap == null) {
            fitMap = new HashMap<>();
            curveSets.put(curveSet.getEquation().getName(), fitMap);
        }
        fitMap.put(curveSet.getState(), curveSet);
        if (best) {
            bestEquation = curveSet.getEquation();
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
        return bestEquation == null ? "" : bestEquation.getName();
    }

    public void setBestEquationName(String equationName) {
        Map<String, CurveFit> curveFits = curveSets.get(equationName);
        for (CurveFit curveFit : curveFits.values()) {
            if (curveFit.getEquation().getName().equals(equationName)) {
                bestEquation = curveFit.getEquation();
            }
        }
    }

    public int getResNum() {
        return resNum;
    }
    
    public String getResName() {
        return resName;
    }
    
    public void setResName(String name) {
        this.resName = name;
    }

    public Double getParValue(String equationName, String state, String parName) {
        Map<String, CurveFit> curveFits = curveSets.get(equationName);
        if (curveFits == null) {
            if (parName.endsWith(".sd")) {
                return err;
            } else {
                return value;
            }
        }
        CurveFit curveFit = curveFits.get(state);
        if (curveFit == null) {
            return null;
        } else {
            return curveFit.getParMap().get(parName);
        }
    }

    public List<ParValueInterface> getParValues(String equationName, String state) {
        final String useEquationName;
        if (equationName.startsWith("best")) {
            useEquationName = bestEquation.getName();
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
        if (plotEquation == null) {
            return sBuilder.toString();
        }
        sBuilder.append(plotEquation.getName());
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
                if (curveFit.getParMap() != null) {
                    sBuilder.append('\n');
                    sBuilder.append(" eqn " + curveFit.getEquation().getName() + " " + curveFit.getState());
                    sBuilder.append(curveFit.getParMap().toString());
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
                PlotEquation plotEquation = curveFit.getEquation();
                if (saveStats) {
                    CurveFit.CurveFitStats curveStats = getFitResult(plotEquation.getName()).getCurveFitStats();
                    String statString = curveStats.toString();
                    sBuilder.append(statString);
                }
                sBuilder.append(curveFit.getState()).append(sep);// fixme
                sBuilder.append(plotEquation.getName()).append(sep);
                Map<String, Double> parMap = curveFit.getParMap();
                parValue = parMap.get("RMS");
                sBuilder.append(String.format("%.3f", parValue)).append(sep);
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
