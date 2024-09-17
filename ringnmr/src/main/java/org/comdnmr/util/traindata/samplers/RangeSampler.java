package org.comdnmr.util.traindata.samplers;

import java.util.*;

public class RangeSampler extends MultiSampler<Double> {
    double minStart;
    double maxStart;
    double minStop;
    double maxStop;
    ArrayList<Double> stepOptions;
    boolean wholeNumbers;

    public RangeSampler(String name, double minStart, double maxStart, double minStop, double maxStop, ArrayList<Double> stepOptions) {
        this(minStart, maxStart, minStop, maxStop, stepOptions);
        this.name = name;
    }

    public RangeSampler(double minStart, double maxStart, double minStop, double maxStop, ArrayList<Double> stepOptions) {
        this.minStart = minStart;
        this.maxStart = maxStart;
        this.minStop = minStop;
        this.maxStop = maxStop;
        this.stepOptions = stepOptions;
        this.wholeNumbers = true;
    }

    @Override
    public ArrayList<Double> sample() {
        double start = new RandomDoubleSampler(minStart, maxStart).sample();
        double stop = new RandomDoubleSampler(minStop, maxStop).sample();
        double step = new ChooseSampler<Double>(stepOptions, 1, 1).sample().get(0);

        if (wholeNumbers) {
            start = Math.rint(start);
            stop = Math.rint(stop);
        }

        double curr = start;
        ArrayList<Double> result = new ArrayList<Double>();
        while (curr < stop) {
            result.add(curr);
            curr += step;
        }

        return result;
    }
}
