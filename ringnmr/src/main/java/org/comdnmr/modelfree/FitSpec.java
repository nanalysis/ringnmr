package org.comdnmr.modelfree;

import java.util.*;


// TODO: if tauM is not specified, it must be estimated using the data.
// Including the Map<String, MolDataValues> into the build method will allow
// this to be computed without the risk if an IllegalStateException
// (see getTauM below).
// estimateTau() is currently defined in FitR1R2NOEModel.
// Makes more sense to me to make this a method of this class
abstract class FitSpec{

    protected double tauM;
    protected boolean tauMNeedsComputing;
    protected final boolean fitTauM;
    protected final boolean fitJ;
    protected final BootstrapMode bootstrapMode;
    protected final boolean fitExchange; // TODO: Currently not supported
    protected final double tauFraction;
    protected final double t2Limit;
    protected final int nReplicates;

    protected FitSpec(Builder<?> builder) {
        if (builder.tauM.isEmpty()) {
            this.tauM = 0.0;
            tauMNeedsComputing = true;
        } else {
            this.tauM = builder.tauM.get();
            tauMNeedsComputing = false;
        }
        this.fitTauM = builder.fitTauM;
        this.fitJ = builder.fitJ;
        this.bootstrapMode = builder.bootstrapMode;
        this.fitExchange = builder.fitExchange;
        this.tauFraction = builder.tauFraction;
        this.t2Limit = builder.t2Limit;
        this.nReplicates = builder.nReplicates;
    }

    public void setTauM(double tauM) {
        this.tauM = tauM;
        tauMNeedsComputing = false;
    }

    public double getTauM() {
        if (tauMNeedsComputing) {
            throw new IllegalStateException("tauM has not be set!");
        }
        return tauM;
    }

    public boolean localfitTauM(Map<String, MolDataValues> data) {
        return (
            t2Limit < 1.0e-6 ||
            data
                .values()
                .stream()
                .anyMatch(
                    value ->
                    value
                        .getData()
                        .stream()
                        .anyMatch(d -> d.R2 > t2Limit)
                )
        );
    }

    public boolean tauMNeedsComputing() {
        return tauMNeedsComputing;
    }

    public enum BootstrapMode {
        PARAMETRIC,
        AGGREGATE,
        BAYESIAN
    }

    private Map<String, Double> estimateTau(Map<String, MolDataValues> molData) {
        Map<Long, List<RelaxDataValue>> map = new HashMap<>();
        for (var entry : molData.entrySet()) {
            MolDataValues value = entry.getValue();
            for (RelaxDataValue rlxValue : value.getData()) {
                long sfMHz = Math.round(rlxValue.relaxObj.getSF() / 1.0e6);
                List<RelaxDataValue> values = map.computeIfAbsent(sfMHz, k -> new ArrayList<>());
                values.add(rlxValue);
            }
        }
        int max = 0;
        long maxMHz = 0;
        for (long sfMHz : map.keySet()) {
            int size = map.get(sfMHz).size();
            if (size > max) {
                maxMHz = sfMHz;
                max = size;
            }
        }
        if (max > 0) {
            return CorrelationTime.estimateTau(maxMHz, "N", map.get(maxMHz));
        } else {
            return Collections.emptyMap();
        }
    }

    public static abstract class Builder<T extends Builder<T>> {
        private Optional<Double> tauM = Optional.empty();
        private boolean fitTauM = false;
        private boolean fitJ = false;
        private BootstrapMode bootstrapMode = BootstrapMode.PARAMETRIC;
        private boolean fitExchange = false; // Currently not supported
        private double tauFraction = 0.25;
        private double t2Limit = 0.0;
        private int nReplicates = 0;

        @SuppressWarnings("unchecked")
        private T self() {
            return (T) this;
        }

        public T tauM(double tauM) {
            this.tauM = Optional.ofNullable(tauM);
            return self();
        }

        public T fitTauM(boolean fitTauM) {
            this.fitTauM = fitTauM;
            return self();
        }

        public T fitJ(boolean fitJ) {
            this.fitJ = fitJ;
            return self();
        }

        public T bootstrapMode(BootstrapMode bootstrapMode) {
            this.bootstrapMode = bootstrapMode;
            return self();
        }

        public T fitExchange(boolean fitExchange) {
            this.fitExchange = fitExchange;
            return self();
        }

        public T tauFraction(double tauFraction) {
            this.tauFraction = tauFraction;
            return self();
        }

        public T t2Limit(double t2Limit) {
            this.t2Limit = t2Limit;
            return self();
        }

        public T nReplicates(int nReplicates) {
            this.nReplicates = nReplicates;
            return self();
        }

        public abstract FitSpec build();
    }
}
