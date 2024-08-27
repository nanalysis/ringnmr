package org.comdnmr.util.traindata.samplers;

import java.util.ArrayList;

public class MultiSampler <T extends Comparable<T>> {
    public String name;
    public ArrayList<T> sample() {
        return new ArrayList<T>();
    };
}
