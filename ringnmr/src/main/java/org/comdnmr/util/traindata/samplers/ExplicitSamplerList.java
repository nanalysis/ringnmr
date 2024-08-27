package org.comdnmr.util.traindata.samplers;

import java.util.*;

public class ExplicitSamplerList<T extends Comparable<T>> extends MultiSampler<T> {
    public String name;
    public ArrayList<T> value;

    public ExplicitSamplerList(String name, ArrayList<T> value) {
        this.name = name;
        this.value = value;
    }

    public ExplicitSamplerList(ArrayList<T> value) {
        this.value = value;
    }

    @Override
    public ArrayList<T> sample() {
        return this.value;
    }
}
