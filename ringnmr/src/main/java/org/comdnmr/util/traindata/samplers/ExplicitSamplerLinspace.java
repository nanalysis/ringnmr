package org.comdnmr.util.traindata.samplers;

import java.util.*;

public class ExplicitSamplerLinspace extends MultiSampler<Double> {
    String name;
    double first;
    double last;
    int nSamples;

    public ExplicitSamplerLinspace(String name, double first, double last, int nSamples) {
        this.name = name;
        this.first = first;
        this.last = last;
        this.nSamples = nSamples;
    }

    public ExplicitSamplerLinspace(double first, double last, int nSamples) {
        this.first = first;
        this.last = last;
        this.nSamples = nSamples;
    }

    @Override
    public ArrayList<Double> sample() {
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
