package org.comdnmr.modelfree;

import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.MFModelIso2sf;

import org.nmrfx.chemistry.relax.OrderPar;
import org.nmrfx.chemistry.relax.OrderParSet;

public class RegularizationFitSpec extends FitSpec {

    private static final String KEY = "REGULARIZATION";

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

    double getLambdaS() { return lambdaS; }

    double getLambdaTauF() { return lambdaTauF; }

    double getLambdaTauS() { return lambdaTauS; }

    @Override
    public String toToml() {
        StringBuilder builder = getBaseTomlBuilder();
        builder.append(String.format("lambdaS = %s%n", lambdaS));
        builder.append(String.format("lambdaTauF = %s%n", lambdaTauF));
        builder.append(String.format("lambdaTauS = %s", lambdaTauS));
        return builder.toString();
    }

    @Override
    protected RelaxFit initRelaxFit(String key, MolDataValues data) {
        RelaxFit relaxFit = super.initRelaxFit(key, data);
        relaxFit.setUseLambda(true);
        relaxFit.setLambdaS(getLambdaS());
        relaxFit.setLambdaTauF(getLambdaTauF());
        relaxFit.setLambdaTauS(getLambdaTauS());
        return relaxFit;
    }

    private MFModelIso2sf getModel2sf(MolDataValues data) {
        boolean fitTauM = fitTauM(data);
        return (MFModelIso2sf) MFModelIso.buildModel("2sf", fitTauM, tauM, tauMFraction, fitExchange);
    }

    @Override
    public ModelFitResult fit(String key, MolDataValues data, Map<String, OrderParSet> orderParSetMap) {
        RelaxFit relaxFit = initRelaxFit(key, data);
        MFModelIso2sf model = getModel2sf(data);
        data.setTestModel(model);

        int nParameters = model.getNPars();
        int nWeights = data.getNValues();
        double[][] parameters = new double[nParameters][nReplicates];
        double[][] weights = new double[nWeights][nReplicates];
        BootstrapSampler sampler = getBootstrapSampler(data);

        Score[] scores = new Score[nReplicates];
        for (int i = 0; i < nReplicates; i++) {
            MolDataValues replicateData = sampler.sample();
            relaxFit.setRelaxData(key, replicateData);
            Score score = runFit(relaxFit, model);
            scores[i] = score;
            double[] replicateParameters = model.getStandardPars(score.getPars());
            double[] replicateWeights = replicateData.getWeights();
            for (int k = 0; k < nParameters; k++) parameters[k][i] = replicateParameters[k];
            for (int j = 0; j < nWeights; j++) weights[j][i] = replicateWeights[j];
        }

        Pair<double[], double[]> parameterEstimates = computeStatistics(parameters, weights);
        double[] fitParameters = parameterEstimates.getLeft();
        double[] fitErrors = parameterEstimates.getRight();

        orderParSetMap.computeIfAbsent(KEY, ky -> new OrderParSet(ky));
        OrderPar orderPar = makeOrderParSet(
            orderParSetMap.get(KEY),
            sampler.getOriginalData(),
            key,
            // FIXME: not sure what Score should be for makeOrderParSet...
            scores[0],
            model,
            fitParameters,
            fitErrors
        );

        return new ModelFitResult(orderPar, parameters, null);
    }

    protected void appendSubclassState(StringBuilder sb) {
        sb.append("lambdaS=").append(Double.doubleToLongBits(lambdaS)).append('|');
        sb.append("lambdaTauF=").append(Double.doubleToLongBits(lambdaTauF)).append('|');
        sb.append("lambdaTauS=").append(Double.doubleToLongBits(lambdaTauS)).append('|');
    }

}
