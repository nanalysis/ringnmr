package org.comdnmr.modelfree;

import java.util.*;

import org.apache.commons.math3.optim.PointValuePair;
import org.comdnmr.modelfree.models.MFModelIso;


// TODO: if tauM is not specified, it must be estimated using the data.
// Including the Map<String, MolDataValues> into the build method will allow
// this to be computed without the risk if an IllegalStateException
// (see getTauM below).
// estimateTau() is currently defined in FitR1R2NOEModel.
// Makes more sense to me to make this a method of this class
public abstract class FitSpec{

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

    public abstract ModelFitResult fit(String key, MolDataValues data);

    protected static Score runFit(RelaxFit relaxFit, MFModelIso model) {
        double[] start = model.getStart();
        double[] lower = model.getLower();
        double[] upper = model.getUpper();

        Optional<PointValuePair> result = relaxFit.fitResidueToModel(start, lower, upper);
        if (result.isEmpty()) {
            throw new RuntimeException("Could not generate fit result.");
        }
        return relaxFit.score(result.get().getPoint(), true);
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

    public double getLocalTauFraction(MolDataValues data) {
        boolean fitTau = (
            t2Limit < 1.0e-6 ||
            data.getData() .stream().anyMatch(value -> value.R2 > t2Limit)
        );
        return fitTau ? tauFraction : 0.0;
    }

    public boolean tauMNeedsComputing() {
        return tauMNeedsComputing;
    }

    public enum BootstrapMode {
        PARAMETRIC,
        NONPARAMETRIC,
        BAYESIAN
    }

    public BootstrapSampler getBootstrapSampler(MolDataValues data) {
        return switch (bootstrapMode) {
            case PARAMETRIC -> new ParametricSampler(data);
            case NONPARAMETRIC -> new NonparametricSampler(data);
            case BAYESIAN -> new BayesianSampler(data);
        };
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
