package org.comdnmr.util.traindata.dependants;

import org.comdnmr.util.traindata.samplers.Sampler;

public class Multiplier extends DependantGenerator {
    Sampler sampler;

    public Multiplier(Sampler sampler) {
        this.sampler = sampler;
    }

    @Override
    public double fetch(double value) {
        double mul = sampler.sample().get(0);
        return value * mul;
    }
}
