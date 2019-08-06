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
package org.comdnmr.eqnfit;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 *
 * @author Bruce Johnson
 */
public class CurveFit {

    final String state;
    final String resNum;
    final Map<String, Double> parMap;
    final PlotEquation plotEquation;

    public CurveFit(String state, String resNum, Map<String, Double> parMap, PlotEquation plotEquation) {
        this.state = state;
        this.resNum = resNum;
        this.parMap = parMap;
        this.plotEquation = plotEquation;
    }

    public Map<String, Double> getParMap() {
        return parMap;
    }

    public PlotEquation getEquation() {
        return plotEquation;
    }

    public String getResNum() {
        return resNum;
    }

    public String getState() {
        return state;
    }

    public List<ParValueInterface> getParValues() {
        List<ParValueInterface> dataValues = new ArrayList<>();
        parMap.keySet().stream().sorted().filter((parName) -> (parMap.containsKey(parName + ".sd"))).forEachOrdered((parName) -> {
            double value = parMap.get(parName);
            double err = parMap.get(parName + ".sd");
            dataValues.add(new ParValue(resNum, state, parName, value, err));
        });
        return dataValues;
    }

    public static class CurveFitStats {

        final String refineOptimizer;
        final String bootstrapOptimizer;
        final long refineTime;
        final long bootstrapTime;
        final int nBootstrapSamples;
        final boolean absMode;
        final boolean nonParametricMode;
        final double startRadius;
        final double finalRadius;
        final double tolerance;
        final boolean weight;

        public CurveFitStats(String refineOpt, String bootstrapOpt, long fitTime, long bootTime, int nSamples, boolean useAbs, 
                boolean useNonParametric, double sRadius, double fRadius, double tol, boolean useWeight) {
            this.refineOptimizer = refineOpt;
            this.bootstrapOptimizer = bootstrapOpt;
            this.refineTime = fitTime;
            this.bootstrapTime = bootTime;
            this.nBootstrapSamples = nSamples;
            this.absMode = useAbs;
            this.nonParametricMode = useNonParametric;
            this.startRadius = sRadius;
            this.finalRadius = fRadius;
            this.tolerance = tol;
            this.weight = useWeight;
        }

        public String getRefineOptimizer() {
            return refineOptimizer;
        }

        public String getBootstrapOptimizer() {
            return bootstrapOptimizer;
        }

        public long getRefineTime() {
            return refineTime;
        }

        public long getBootstrapTime() {
            return bootstrapTime;
        }

        public int nSamples() {
            return nBootstrapSamples;
        }

        @Override
        public String toString() {
            char sep = '\t';
            StringBuilder sBuilder = new StringBuilder();
            sBuilder.append(refineOptimizer).append(sep);
            sBuilder.append(refineTime).append(sep);
            sBuilder.append(bootstrapOptimizer).append(sep);
            sBuilder.append(bootstrapTime).append(sep);
            sBuilder.append(nBootstrapSamples).append(sep);
            sBuilder.append(absMode).append(sep);
            sBuilder.append(nonParametricMode).append(sep);
            sBuilder.append(startRadius).append(sep);
            sBuilder.append(finalRadius).append(sep);
            sBuilder.append(tolerance).append(sep);
            sBuilder.append(weight).append(sep);
            String statString = sBuilder.toString();
            return statString;
        }
    }

}
