package org.comdnmr.modelfree;

import java.util.ArrayList;
import java.util.List;

import org.comdnmr.modelfree.models.MFModelIso;

public abstract class ModelSelectionFitSpec extends FitSpec {

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

    public List<MFModelIso> getModels(MolDataValues data) {
        List<MFModelIso> models = new ArrayList<>();
        boolean fitTauM = fitTauM(data);
        for (String name : getModelNames()) {
            MFModelIso model = MFModelIso.buildModel(
                name,
                fitTauM,
                tauM,
                tauFraction,
                fitExchange);
            models.add(model);
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
