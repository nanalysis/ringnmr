package org.comdnmr.util.traindata.samplers;

import java.util.*;

public class ChooseSampler<T extends Comparable<T>> extends MultiSampler<T> {
    public ArrayList<T> options;
    public int minSamples;
    public int maxSamples;

    public ChooseSampler(String name, ArrayList<T> options, int minSamples, int maxSamples) {
        this.name = name;
        this.options = options;
        this.minSamples = minSamples;
        this.maxSamples = maxSamples;
    }

    public ChooseSampler(ArrayList<T> options, int minSamples, int maxSamples) {
        this.options = options;
        this.minSamples = minSamples;
        this.maxSamples = maxSamples;
    }

    @Override
    public ArrayList<T> sample() {
        Collections.shuffle(options);

        RandomIntegerSampler nSamplesSampler = new RandomIntegerSampler(minSamples, maxSamples);
        int nSamples = nSamplesSampler.sample();

        ArrayList<T> result = new ArrayList<T>();
        for (int i = 0; i < nSamples; i++) {
            result.add(options.get(i));
        }

        Collections.sort(result);
        return result;
    }
}
