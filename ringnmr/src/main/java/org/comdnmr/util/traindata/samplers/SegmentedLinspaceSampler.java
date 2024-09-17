package org.comdnmr.util.traindata.samplers;

import java.util.*;

public class SegmentedLinspaceSampler extends MultiSampler<Double> {
    ArrayList<Double> minFirsts;
    ArrayList<Double> maxFirsts;
    ArrayList<Double> minLasts;
    ArrayList<Double> maxLasts;
    ArrayList<Integer> minSampless;
    ArrayList<Integer> maxSampless;

    public SegmentedLinspaceSampler(
        String name,
        ArrayList<Double> minFirsts,
        ArrayList<Double> maxFirsts,
        ArrayList<Double> minLasts,
        ArrayList<Double> maxLasts,
        ArrayList<Integer> minSampless,
        ArrayList<Integer> maxSampless
    ) {
        this.name = name;
        this.minFirsts = minFirsts;
        this.maxFirsts = maxFirsts;
        this.minLasts = minLasts;
        this.maxLasts = maxLasts;
        this.minSampless = minSampless;
        this.maxSampless = maxSampless;
    }

    public SegmentedLinspaceSampler(
        ArrayList<Double> minFirsts,
        ArrayList<Double> maxFirsts,
        ArrayList<Double> minLasts,
        ArrayList<Double> maxLasts,
        ArrayList<Integer> minSampless,
        ArrayList<Integer> maxSampless
    ) {
        this.minFirsts = minFirsts;
        this.maxFirsts = maxFirsts;
        this.minLasts = minLasts;
        this.maxLasts = maxLasts;
        this.minSampless = minSampless;
        this.maxSampless = maxSampless;
    }

    @Override
    public ArrayList<Double> sample() {
        ArrayList<Double> segmentFirsts = new ArrayList<>();
        for (int i = 0; i < minFirsts.size(); i++) {
            segmentFirsts.add(new RandomDoubleSampler(
                minFirsts.get(i),
                maxFirsts.get(i)
            ).sample());
        }

        ArrayList<Double> segmentLasts = new ArrayList<>();
        for (int i = 0; i < minLasts.size(); i++) {
            segmentLasts.add(new RandomDoubleSampler(
                minLasts.get(i),
                maxLasts.get(i)
            ).sample());
        }

        ArrayList<Integer> segmentSamples = new ArrayList<>();
        for (int i = 0; i < minSampless.size(); i++) {
            segmentSamples.add(new RandomIntegerSampler(
                minSampless.get(i),
                maxSampless.get(i)
            ).sample());
        }

        ArrayList<Double> values = new ArrayList<>();

        for (int i = 0; i < segmentFirsts.size(); i++) {
            double first = segmentFirsts.get(i);
            double last = segmentLasts.get(i);
            int nSamples = segmentSamples.get(i);

            double incr = (last - first) / (nSamples - 1);
            double curr = first;

            ArrayList<Double> result = new ArrayList<Double>();
            for (int j = 0; j < nSamples; j++) {
                values.add(curr);
                curr = curr + incr;
            }
        }

        return values;
    }
}
