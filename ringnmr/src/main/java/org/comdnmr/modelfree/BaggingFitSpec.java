package org.comdnmr.modelfree;

class BaggingFitSpec extends ModelSelectionFitSpec {

    BaggingFitSpec(Builder builder) {
        super(builder);
    }

    public static class Builder extends ModelSelectionFitSpec.Builder {

        public BaggingFitSpec build() {
            return new BaggingFitSpec(this);
        }
    }
}
