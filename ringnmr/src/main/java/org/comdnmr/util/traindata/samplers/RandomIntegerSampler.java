package org.comdnmr.util.traindata.samplers;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

public class RandomIntegerSampler extends SingleSampler<Integer> {
    public int minInt;
    public int maxInt;
    public int minSamples;
    public int maxSamples;

    public RandomIntegerSampler(String name, int minInt, int maxInt) {
        this.name = name;
        this.minInt = minInt;
        this.maxInt = maxInt;
    }

    public RandomIntegerSampler(int minInt, int maxInt) {
        this.minInt = minInt;
        this.maxInt = maxInt;
    }

    @Override
    public Integer sample() {
        int randomNum = ThreadLocalRandom.current().nextInt(minInt, maxInt + 1);
        return randomNum;
    }
}
