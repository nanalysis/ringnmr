/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.calc;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

/**
 *
 * @author brucejohnson
 */
public class CorrelationTime {

    public static double estimateTau(Map<String, ResidueProperties> residueProps) {
        ResidueProperties resPropsR1 = residueProps.get("t1_output");
        ResidueProperties resPropsR2 = residueProps.get("t2_output");
        double tau = 0.0;
        if ((resPropsR1 != null) && (resPropsR2 != null)) {
            Map<Integer, Double> r1Map = resPropsR1.getParMapData(
                    "best", "0:0:0", "R");
            Map<Integer, Double> r2Map = resPropsR2.getParMapData(
                    "best", "0:0:0", "R");

            List<Double> t1List = new ArrayList<>();
            List<Double> t2List = new ArrayList<>();
            SummaryStatistics stats = new SummaryStatistics();
            r2Map.values().stream().forEach(v -> stats.addValue(v));
            double mean = stats.getMean();
            double sdev = stats.getStandardDeviation();
            stats.clear();
            for (Integer res : r1Map.keySet()) {
                Double r1 = r1Map.get(res);
                Double r2 = r2Map.get(res);
                if ((r1 != null) && (r2 != null)) {
                    if (r2 > (mean - sdev)) {
                        stats.addValue(r2 / r1);
                    }
                }
            }
            double ratio = stats.getMean();
            double sf = 600e6;
            tau = fit(sf, ratio);
            System.out.println("tau " + ratio + " " + tau);
        }
        return tau;
    }

    static class MatchFunction implements UnivariateFunction {

        double ratio;
        RelaxEquations r;

        MatchFunction(double sf, double ratio) {
            this.ratio = ratio;
            r = new RelaxEquations(sf, "H", "N");
        }

        @Override
        public double value(double tauC) {

            double calcRatio = r.r2r1Ratio(tauC);
            double value = Math.abs(ratio - calcRatio);
            return value;
        }
    }

    public static double fit(double sf, double ratio) {
        double tolAbs = 1E-12;
        double min = 1.0e-9;
        double max = 100.0e-9;
        MatchFunction f = new MatchFunction(sf, ratio);
        UnivariateObjectiveFunction fOpt = new UnivariateObjectiveFunction(f);
        SearchInterval searchInterval = new SearchInterval(min, max);
        MaxEval maxEval = new MaxEval(100);
        double best;

        BrentOptimizer brentOptimizer = new BrentOptimizer(tolAbs * 10.0, tolAbs);
        UnivariatePointValuePair optValue = brentOptimizer.optimize(fOpt, GoalType.MINIMIZE, searchInterval, maxEval);

        best = optValue.getPoint();
        return best;
    }

}
