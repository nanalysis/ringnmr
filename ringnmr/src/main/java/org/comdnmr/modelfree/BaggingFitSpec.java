package org.comdnmr.modelfree;

public class BaggingFitSpec extends ModelSelectionFitSpec {

    BaggingFitSpec(Builder builder) {
        super(builder);
    }

    public ModelFitResult fit(String key, MolDataValues data) {
        throw new UnsupportedOperationException("Bagging fitting yet to be implemented");
    }

    public static class Builder extends ModelSelectionFitSpec.Builder {

        public BaggingFitSpec build() {
            return new BaggingFitSpec(this);
        }
    }
}
