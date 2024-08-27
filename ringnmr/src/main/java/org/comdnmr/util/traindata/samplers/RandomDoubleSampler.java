package org.comdnmr.util.traindata.samplers;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

public class RandomDoubleSampler extends SingleSampler<Double> {
    public double min;
    public double max;

    public RandomDoubleSampler(String name, double min, double max) {
        this.name = name;
        this.min = min;
        this.max = max;
    }

    public RandomDoubleSampler(double min, double max) {
        this.min = min;
        this.max = max;
    }

    @Override
    public Double sample() {
        double value = ThreadLocalRandom.current().nextDouble(min, max);
        return value;
    }
}
