/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.modelfree;

import org.comdnmr.data.ResidueProperties;
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

    public static double estimateTau(Map<String, ResidueProperties> residueProps,
            String r1SetName, String r2SetName, String noeSetName) {
        ResidueProperties resPropsR1 = residueProps.get(r1SetName);
        ResidueProperties resPropsR2 = residueProps.get(r2SetName);
        ResidueProperties resPropsNOE = residueProps.get(noeSetName);
        double tau = 0.0;
        if ((resPropsR1 != null) && (resPropsR2 != null)) {
            Map<Integer, Double> r1Map = resPropsR1.getParMapData(
                    "best", "0:0:0", "R");
            Map<Integer, Double> r2Map = resPropsR2.getParMapData(
                    "best", "0:0:0", "R");

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

    static class R1R2NOEDelta implements UnivariateFunction {

        double tau;
        double[] targets;
        double[] sigmas;
        RelaxEquations r;
        double[] calcVals;

        R1R2NOEDelta(double sf, double tau, double[] targets, double[] sigmas) {
            this.tau = tau;
            r = new RelaxEquations(sf, "H", "N");
            calcVals = new double[targets.length];
            this.targets = targets;
            this.sigmas = sigmas;
        }

        @Override
        public double value(double S2) {
            double[] jVals = r.getJ(tau, S2);

            double r1 = r.R1(jVals);
            double r2 = r.R2(jVals, 0.0);
            double noe = r.NOE(jVals);
            calcVals[0] = r1;
            calcVals[1] = r2;
            calcVals[2] = noe;

            double sumsq = 0.0;
            for (int i = 0; i < calcVals.length; i++) {
                double delta = (calcVals[i] - targets[i]) / sigmas[i];
                sumsq += delta * delta;
            }

            return sumsq;
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

    public static double fitS(double sf, double tau, double[] targets, double[] sigmas) {
        double tolAbs = 1E-12;
        double min = 0.1;
        double max = 1.0;
        R1R2NOEDelta f = new R1R2NOEDelta(sf, tau, targets, sigmas);
        UnivariateObjectiveFunction fOpt = new UnivariateObjectiveFunction(f);
        SearchInterval searchInterval = new SearchInterval(min, max);
        MaxEval maxEval = new MaxEval(100);
        double best;

        BrentOptimizer brentOptimizer = new BrentOptimizer(tolAbs * 10.0, tolAbs);
        UnivariatePointValuePair optValue = brentOptimizer.optimize(fOpt, GoalType.MINIMIZE, searchInterval, maxEval);

        best = optValue.getPoint();
        return best;
    }

    public static double fitS(Map<String, ResidueProperties> residueProps,
            String r1SetName, String r2SetName, String noeSetName) {
        ResidueProperties resPropsR1 = residueProps.get(r1SetName);
        ResidueProperties resPropsR2 = residueProps.get(r2SetName);
        ResidueProperties resPropsNOE = residueProps.get(noeSetName);
        System.out.println("fit S " + resPropsR1 + " " +resPropsR2 + " " + 
                resPropsNOE);
        double tau = 0.0;
        if ((resPropsR1 != null) && (resPropsR2 != null)) {
            Map<Integer, Double> r1Map = resPropsR1.getParMapData(
                    "best", "0:0:0", "R");
            Map<Integer, Double> r2Map = resPropsR2.getParMapData(
                    "best", "0:0:0", "R");
            Map<Integer, Double> noeMap = resPropsNOE.getParMapData(
                    "best", "0:0:0", "NOE");
            Map<Integer, Double> r1ErrMap = resPropsR1.getParMapData(
                    "best", "0:0:0", "R.sd");
            Map<Integer, Double> r2ErrMap = resPropsR2.getParMapData(
                    "best", "0:0:0", "R.sd");
            Map<Integer, Double> noeErrMap = resPropsNOE.getParMapData(
                    "best", "0:0:0", "NOE.sd");

            double sf = 700.0e6;
            double tauc = 4.4e-9;
            for (Integer res : r1Map.keySet()) {
                Double r1 = r1Map.get(res);
                Double r2 = r2Map.get(res);
                Double noe = noeMap.get(res);
                    System.out.print("res " + res + " " + r1  + " " + r2 + " " + noe);
                Double r1Err = r1ErrMap.get(res);
                Double r2Err = r2ErrMap.get(res);
                Double noeErr = noeErrMap.get(res);
                    System.out.println(" " + r1Err  + " " + r2Err + " " + noeErr);
                double[] targets = {r1, r2, noe};
                double[] sigmas = {r1Err, r2Err, noeErr};
                if ((r1 != null) && (r2 != null) && (noe != null)) {
                    double S2 = fitS(sf, tauc, targets, sigmas);
                    System.out.println("res " + res + " " + S2);
                }
            }
        }
        return tau;
    }

}
