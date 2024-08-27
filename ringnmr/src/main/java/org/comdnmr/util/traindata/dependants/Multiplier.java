package org.comdnmr.util.traindata.dependants;

public class Multiplier extends DependantGenerator {
    public double mul;

    public Multiplier(double mul) {
        this.mul = mul;
    }
    @Override
    public double fetch(double value) {
        return value * mul;
    }
}
