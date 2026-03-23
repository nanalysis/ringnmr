package org.comdnmr.modelfree;

import java.util.ArrayList;
import java.util.List;

import org.comdnmr.modelfree.models.MFModelIso;

abstract class ModelSelectionFitSpec extends FitSpec {

    private final List<String> modelNames;

    protected ModelSelectionFitSpec(Builder builder) {
        super(builder);
        if (builder.modelNames.isEmpty()) {
            throw new IllegalArgumentException("`modelNames` must not be empty!");
        }
        this.modelNames = builder.modelNames;
    }

    public List<String> getModelNames() {
        return modelNames;
    }

    public List<MFModelIso> getModels() {
        List<MFModelIso> models = new ArrayList<>();
        for (String name : getModelNames()) {
            models.add(MFModelIso.buildModel(name, fitTauM, getTauM(), tauFraction, fitExchange));
        }
        return models;
    }

    // Concrete Builder for ConcreteFitSpec
    public static abstract class Builder extends FitSpec.Builder<Builder> {

        private final List<String> modelNames = new ArrayList<>();

        public Builder modelNames(List<String> modelNames) {
            this.modelNames.addAll(modelNames);
            return this;
        }

        abstract public ModelSelectionFitSpec build();
    }
}
