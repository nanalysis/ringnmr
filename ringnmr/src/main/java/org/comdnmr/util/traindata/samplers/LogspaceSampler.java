package org.comdnmr.util.traindata.samplers;

import java.util.*;

public class LogspaceSampler extends MultiSampler<Double> {
    public double minFirst;
    public double maxFirst;
    public double scale;
    public int minSamples;
    public int maxSamples;

    public LogspaceSampler(String name, double minFirst, double maxFirst, double scale, int minSamples, int maxSamples) {
        this.name = name;
        this.minFirst = minFirst;
        this.maxFirst = maxFirst;
        this.scale = scale;
        this.minSamples = minSamples;
        this.maxSamples = maxSamples;
    }

    public LogspaceSampler(double minFirst, double maxFirst, double scale, int minSamples, int maxSamples) {
        this.minFirst = minFirst;
        this.maxFirst = maxFirst;
        this.scale = scale;
        this.minSamples = minSamples;
        this.maxSamples = maxSamples;
    }

    @Override
    public ArrayList<Double> sample() {
        RandomDoubleSampler firstSampler = new RandomDoubleSampler(minFirst, maxFirst);
        double first = firstSampler.sample();

        RandomIntegerSampler nSamplesSampler = new RandomIntegerSampler(minSamples, maxSamples);
        int nSamples = nSamplesSampler.sample();

        double curr = first;

        ArrayList<Double> result = new ArrayList<Double>();
        for (int i = 0; i < nSamples; i++) {
            result.add(curr);
            curr = curr * scale;
        }

        return result;
    }
}
