package org.comdnmr.util.traindata.samplers;

import java.util.*;

public class ExplicitSampler<T extends Comparable<T>> extends MultiSampler<T> {
    public T value;

    public ExplicitSampler(String name, T value) {
        this.name = name;
        this.value = value;
    }

    public ExplicitSampler(T value) {
        this.value = value;
    }

    @Override
    public ArrayList<T> sample() {
        ArrayList<T> result = new ArrayList<T>();
        result.add(this.value);
        return result;
    }
}
