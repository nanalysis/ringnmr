package org.comdnmr.fit.calc;

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

        public CurveFitStats(String refineOpt, String bootstrapOpt, long fitTime, long bootTime, int nSamples) {
            this.refineOptimizer = refineOpt;
            this.bootstrapOptimizer = bootstrapOpt;
            this.refineTime = fitTime;
            this.bootstrapTime = bootTime;
            this.nBootstrapSamples = nSamples;
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
            String statString = sBuilder.toString();
            return statString;
        }
    }

}
