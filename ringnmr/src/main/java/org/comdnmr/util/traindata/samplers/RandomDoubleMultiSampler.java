package org.comdnmr.util.traindata.samplers;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

public class RandomDoubleMultiSampler extends MultiSampler<Double> {
    public String name;
    public double min;
    public double max;
    public int minSamples;
    public int maxSamples;

    public RandomDoubleMultiSampler(double min, double max, int minSamples, int maxSamples) {
        this.min = min;
        this.max = max;
        this.minSamples = minSamples;
        this.maxSamples = maxSamples;
    }

    public RandomDoubleMultiSampler(String name, double min, double max, int minSamples, int maxSamples) {
        this.name = name;
        this.min = min;
        this.max = max;
        this.minSamples = minSamples;
        this.maxSamples = maxSamples;
    }

    @Override
    public ArrayList<Double> sample() {
        ArrayList<Double> value = new ArrayList<Double>();
        int nSamples = ThreadLocalRandom.current().nextInt(minSamples, maxSamples + 1);
        for (int i = 0; i < nSamples; i++) {
            value.add(ThreadLocalRandom.current().nextDouble(min, max));
        }

        Collections.sort(value);
        return value;
    }
}
