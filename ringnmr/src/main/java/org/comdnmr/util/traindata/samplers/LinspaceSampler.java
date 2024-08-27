package org.comdnmr.util.traindata.samplers;

import java.util.*;

public class LinspaceSampler extends MultiSampler<Double> {
    double minFirst;
    double maxFirst;
    double minLast;
    double maxLast;
    int minSamples;
    int maxSamples;

    public LinspaceSampler(String name, double minFirst, double maxFirst, double minLast, double maxLast, int minSamples, int maxSamples) {
        this.name = name;
        this.minFirst = minFirst;
        this.maxFirst = maxFirst;
        this.minLast = minLast;
        this.maxLast = maxLast;
        this.minSamples = minSamples;
        this.maxSamples = maxSamples;
    }

    public LinspaceSampler(double minFirst, double maxFirst, double minLast, double maxLast, int minSamples, int maxSamples) {
        this.minFirst = minFirst;
        this.maxFirst = maxFirst;
        this.minLast = minLast;
        this.maxLast = maxLast;
        this.minSamples = minSamples;
        this.maxSamples = maxSamples;
    }

    @Override
    public ArrayList<Double> sample() {
        RandomDoubleSampler firstSampler = new RandomDoubleSampler(minFirst, maxFirst);
        double first = firstSampler.sample();

        RandomDoubleSampler lastSampler = new RandomDoubleSampler(minLast, maxLast);
        double last = lastSampler.sample();

        RandomIntegerSampler nSamplesSampler = new RandomIntegerSampler(minSamples, maxSamples);
        int nSamples = nSamplesSampler.sample();

        double incr = (last - first) / (nSamples - 1);
        double curr = first;

        ArrayList<Double> result = new ArrayList<Double>();
        for (int i = 0; i < nSamples; i++) {
            result.add(curr);
            curr = curr + incr;
        }

        return result;
    }
}
