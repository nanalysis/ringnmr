package org.comdnmr.modelfree;

class ConventionalFitSpec extends ModelSelectionFitSpec {

    ConventionalFitSpec(Builder builder) {
        super(builder);
    }

    public static class Builder extends ModelSelectionFitSpec.Builder {

        public ConventionalFitSpec build() {
            return new ConventionalFitSpec(this);
        }
    }
}
