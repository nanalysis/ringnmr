package org.comdnmr.modelfree;

public class RegularizationFitSpec extends FitSpec {

    private final double lambdaS;
    private final double lambdaTauF;
    private final double lambdaTauS;

    public static class Builder extends FitSpec.Builder<Builder> {
        private double lambdaS = 0.0;
        private double lambdaTauF = 0.0;
        private double lambdaTauS = 0.0;

        private void validateLambda(String name, double value) {
            if (value < 0.0) {
                throw new IllegalArgumentException(String.format("%s must be >= 0.0", name));
            }
        }

        public Builder lambdaS(double lambdaS) {
            validateLambda("lambdaS", lambdaS);
            this.lambdaS = lambdaS;
            return this;
        }

        public Builder lambdaTauF(double lambdaTauF) {
            validateLambda("lambdaTauF", lambdaTauF);
            this.lambdaTauF = lambdaTauF;
            return this;
        }

        public Builder lambdaTauS(double lambdaTauS) {
            validateLambda("lambdaTauS", lambdaTauS);
            this.lambdaTauS = lambdaTauS;
            return this;
        }

        public RegularizationFitSpec build() {
            return new RegularizationFitSpec(this);
        }
    }

    protected RegularizationFitSpec(Builder builder) {
        super(builder);
        this.lambdaS = builder.lambdaS;
        this.lambdaTauF = builder.lambdaTauF;
        this.lambdaTauS = builder.lambdaTauS;
    }

    public ModelFitResult fit(String key, MolDataValues data) {
        throw new UnsupportedOperationException("Regularization fitting yet to be implemented");
    }

    double getLambdaS() { return lambdaS; }

    double getLambdaTauF() { return lambdaTauF; }

    double getLambdaTauS() { return lambdaTauS; }

}
