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
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.relax.ResonanceSource;

/**
 *
 * @author Bruce Johnson
 */
public class ExperimentResult {

    ExperimentSet experimentSet;
    ResonanceSource resonanceSource;
    Map<String, Map<String, CurveFit>> curveFits = new LinkedHashMap<>();
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
        ExperimentResult resInfo;
        String equationName;
        String state;

        public ParValue(ExperimentResult resInfo, String parName, String equationName, String state) {
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
            Double value = resInfo.curveFits.get(equationName).get(state).getParMap().get(parName);
            if (value == null) {
                return 0.0;
            } else {
                return value;
            }
        }

        @Override
        public double getValue() {
            Double value = resInfo.curveFits.get(equationName).get(state).getParMap().get(parName);
            return value;
        }

        public double getError(String key) {
            Double value = resInfo.curveFits.get(equationName).get(state).getParMap().get(parName + ".sd");
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
        public ResonanceSource getDynamicsSource() {
            return resInfo.resonanceSource;
        }

        @Override
        public String getAtomName() {
            return resInfo.resonanceSource.getAtoms()[0].getName();
        }

        @Override
        public String getResidue() {
            return String.valueOf(resInfo.resonanceSource.getAtoms()[0].getResidueNumber());
        }

        @Override
        public String getResName() {
            return resInfo.resonanceSource.getAtoms()[0].getResidueName();
        }

        @Override
        public String getState() {
            return state;
        }
    }

    public ExperimentResult(ExperimentSet experimentSet, ResonanceSource resonanceSource, int groupId, int groupSize) {
        this(experimentSet, resonanceSource, groupId, groupSize, 0);
    }

    public ExperimentResult(ExperimentSet experimentSet, ResonanceSource resonanceSource, int groupId, int groupSize, int peakNum) {
        this.experimentSet = experimentSet;
        this.resonanceSource = resonanceSource;
        this.groupId = groupId;
        this.groupSize = groupSize;
        this.peakNum = peakNum;
    }

    public ResonanceSource getResonanceSource() {
        return resonanceSource;
    }

    public Atom getAtom() {
        return resonanceSource.getAtoms()[0];
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
        return curveFits.get(equationName).values();
    }

    public void addCurveFit(CurveFit curveFit, boolean best) {
        Map<String, CurveFit> fitMap = curveFits.get(curveFit.getEquation().getName());
        if (fitMap == null) {
            fitMap = new HashMap<>();
            curveFits.put(curveFit.getEquation().getName(), fitMap);
        }
        fitMap.put(curveFit.getState(), curveFit);
        if (best) {
            bestEquation = curveFit.getEquation();
        }
    }

    public CurveFit getCurveFit(String equationName, String state) {
        CurveFit curveFit = null;
        if (curveFits != null) {
            Map<String, CurveFit> curveFits = this.curveFits.get(equationName);
            if (curveFits != null) {
                curveFit = this.curveFits.get(equationName).get(state);
            }
        }
        return curveFit;
    }

    public String getBestEquationName() {
        return bestEquation == null ? "" : bestEquation.getName();
    }

    public void setBestEquationName(String equationName) {
        Map<String, CurveFit> curveFits = this.curveFits.get(equationName);
        for (CurveFit curveFit : curveFits.values()) {
            if (curveFit.getEquation().getName().equals(equationName)) {
                bestEquation = curveFit.getEquation();
            }
        }
    }

    public Double getParValue(String equationName, String state, String parName) {
        Map<String, CurveFit> curveFits = this.curveFits.get(equationName);
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
        Map<String, CurveFit> curveFits = this.curveFits.get(useEquationName);
        if (curveFits != null) {
            curveFits.values().stream().forEach(cf -> {
                if (ExperimentSet.matchStateString(state, cf.getState())) {
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
        sBuilder.append(getAtom().getShortName());
        sBuilder.append('\n');
        double[] xValues = null;
        double[] yValues = null;
        if (xValues != null) {
            for (double xValue : xValues) {
                sBuilder.append(xValue);
                sBuilder.append(" ");
            }

            sBuilder.append('\n');
        }
        if (yValues != null) {
            for (double yValue : yValues) {
                sBuilder.append(yValue);
                sBuilder.append(" ");
            }
        }
        for (Map<String, CurveFit> fitMap : curveFits.values()) {
            for (CurveFit curveFit : fitMap.values()) {
                if (curveFit.getParMap() != null) {
                    sBuilder.append('\n');
                    sBuilder.append(" eqn ").append(curveFit.getEquation().
                            getName()).append(" ").append(curveFit.getState());
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
        sBuilder.append(getAtom().getShortName()).append(sep);
        sBuilder.append(peakNum).append(sep);
        sBuilder.append(groupSize).append(sep);
        sBuilder.append(groupId).append(sep);
        String commonString = sBuilder.toString();
        StringBuilder allBuilder = new StringBuilder();
        Double parValue;
        for (Map<String, CurveFit> fitMap : curveFits.values()) {
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
