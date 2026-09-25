package org.comdnmr.modelfree;

import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.MFModelIso2sf;

import org.comdnmr.util.CoMDOptions;
import org.nmrfx.chemistry.relax.OrderPar;
import org.nmrfx.chemistry.relax.OrderParSet;

/**
 * Regularized model-free fitting strategy using the extended
 * ({@code 2sf}) model exclusively.
 *
 * <p>Unlike {@link BaggingFitSpec}, which selects the best model per replicate
 * from a candidate set, this strategy always fits the {@code 2sf} model
 * (parameters: {@code [τm?,] S²f, τf, S²s, τs}) with
 * regularization terms that penalize deviations of the order parameters and
 * correlation times from physically motivated null values. The regularization
 * discourages overfitting without the need for AICc-based model selection.</p>
 *
 * <p>The four strength scale — {@link #lambdaScale}— is passed through to
 * {@link RelaxFit} via {@link RelaxFit#setLambdaS2F(double)},
 * {@link RelaxFit#setLambdaS2S(double)}, {@link RelaxFit#setLambdaTauF(double)},
 * and {@link RelaxFit#setLambdaTauS(double)} respectively. All lower bounds for
 * the model parameters are set to zero (see {@link #getLower(MFModelIso)}).</p>
 *
 * <p>Bootstrap resampling and statistics follow the same procedure as
 * {@link BaggingFitSpec}: for each replicate, the data are resampled via the
 * configured {@link FitSpec.BootstrapMode}, the 2sf model is fit, and the
 * resulting parameters are post-processed via
 * {@link #processParamsAfterFit(double[], boolean)} before being aggregated
 * across replicates by
 * {@link FitSpec#computeStatistics(double[][], double[][])}.</p>
 *
 * <h2>Example usage:</h2>
 * <pre>{@code
 * FitSpec spec = new RegularizationFitSpec.Builder()
 *     .tauM(17.5)
 *     .bootstrapMode(BootstrapMode.NONPARAMETRIC)
 *     .nReplicates(200)
 *     .lambdaS2F(0.5)
 *     .lambdaS2S(0.5)
 *     .lambdaTauF(0.1)
 *     .lambdaTauS(0.2)
 *     .build();
 *
 * ModelFitResult result = spec.fit(key, data, orderParSetMap);
 * }</pre>
 *
 * @see BaggingFitSpec
 * @see ConventionalFitSpec
 * @see FitSpec.BootstrapMode
 */
public class RegularizationFitSpec extends FitSpec {

    /**
     * Key used when registering results in the {@link OrderParSet} map.
     */
    private static final String KEY = "REGULARIZATION";

    /**
     * Order-parameter threshold above which the corresponding motion is
     * considered absent (S² ≈ 1). When both S² values exceed this limit,
     * the result is mapped to "model 0" (no local motions).
     */
    static final double S2_THOLD = 0.99;

    /**
     * Correlation-time threshold (ns) below which the corresponding motion is
     * treated as instantaneous (τ → 0). Motions with τ less than this limit
     * are assigned τ = 0.
     */
    static final double TAU_THOLD = 1.0e-3;  // 1 ps

    /**
     * Regularization strength for the fast order parameter S²f.
     * Higher values push S²f estimates towards 1 (rigid limit).
     */
    private final double lambdaScale;

    /**
     * SNR threshold, in units of the CRLB, below which a component is dropped
     * by {@link MFModelIso2sf#calcThreshold(double[], double)}.  A component
     * survives only if {@code tau / crlb >= stringency}, so this is literally
     * "how many sigma above zero must a component be".
     *
     * <p>Because a component whose true value is zero has an estimate bounded
     * below at zero, {@code tau_hat/sigma} is approximately half-normal and the
     * false-positive rate is about {@code 2*(1 - Phi(stringency))}: roughly 32%
     * at 1.0, 13% at 1.5, 4.6% at 2.0 and 1.2% at 2.5.</p>
     */
    private final double stringency;

    /**
     * Builder for {@link RegularizationFitSpec}.
     *
     * <p>All base configuration (tau_M, bootstrap mode, number of replicates,
     * etc.) is inherited from {@link FitSpec.Builder}. The three
     * regularization strengths have sensible defaults and may be set
     * individually.</p>
     *
     * <h2>Example:</h2>
     * <pre>{@code
     * FitSpec spec = new RegularizationFitSpec.Builder()
     *     .tauM(17.5)
     *     .bootstrapMode(BootstrapMode.NONPARAMETRIC)
     *     .nReplicates(200)
     *     .lambdaS2F(0.5)
     *     .lambdaS2S(0.5)
     *     .lambdaTauF(0.1)
     *     .lambdaTauS(0.2)
     *     .build();
     * }</pre>
     */
    public static class Builder extends FitSpec.Builder<Builder> {

        private static final double DEFAULT_LAMBDA_SCALE = 1.0;
        private static final double DEFAULT_STRINGENCY = 1.0;

        private double lambdaScale = DEFAULT_LAMBDA_SCALE;
        private double stringency = DEFAULT_STRINGENCY;

        /**
         * Returns the default regularization strength for S²f (0.5).
         */
        public static double getDefaultLambdaScale() {
            return DEFAULT_LAMBDA_SCALE;
        }

        private void validateLambda(String name, double value) {
            if (value < 0.0) {
                throw new IllegalArgumentException(
                        String.format("%s must be >= 0.0, got: %s", name, value)
                );
            }
        }

        /**
         * Sets the regularization strength for the fast order parameter S²f.
         *
         * @param lambdaScale regularization weight; must be &ge; 0
         * @return this builder
         * @throws IllegalArgumentException if {@code lambdaS2F} is negative
         */
        public Builder lambdaScale(double lambdaScale) {
            validateLambda("lambdaS2F", lambdaScale);
            this.lambdaScale = lambdaScale;
            return this;
        }

        /**
         * Returns the default component-retention threshold (1.0 sigma).
         */
        public static double getDefaultStringency() {
            return DEFAULT_STRINGENCY;
        }

        /**
         * Sets the SNR threshold for retaining a component, in units of the
         * CRLB.  1.0 keeps anything one sigma above zero; 2.0 cuts the
         * false-positive rate roughly sevenfold at the cost of about half the
         * detections of genuinely marginal (2 sigma) components.
         *
         * @param stringency retention threshold; must be &ge; 0
         * @return this builder
         * @throws IllegalArgumentException if {@code stringency} is negative
         */
        public Builder stringency(double stringency) {
            validateLambda("stringency", stringency);
            this.stringency = stringency;
            return this;
        }

        /**
         * Validates the configuration and constructs a new
         * {@link RegularizationFitSpec}.
         *
         * @return the configured {@code RegularizationFitSpec}
         * @throws IllegalStateException if the configuration is invalid
         */
        @Override
        public RegularizationFitSpec build() {
            validate();
            return new RegularizationFitSpec(this);
        }
    }

    /**
     * Constructs a new {@code RegularizationFitSpec} from the given builder.
     *
     * @param builder the builder containing all configuration values
     */
    protected RegularizationFitSpec(Builder builder) {
        super(builder);
        this.lambdaScale = builder.lambdaScale;
        this.stringency = builder.stringency;
    }

    /**
     * Returns the S²f regularization strength.
     */
    double getLambdaScale() {
        return lambdaScale;
    }

    /**
     * Returns the component-retention SNR threshold.
     */
    double getStringency() {
        return stringency;
    }


    @Override
    protected double[] getLower(MFModelIso model) {
        return new double[model.getNPars()];
    }

    /**
     * {@inheritDoc}
     *
     * <p>Appends {@code lambdaS2F}, {@code lambdaS2S}, {@code lambdaTauF}, and
     * {@code lambdaTauS} to the canonical state string.</p>
     */
    @Override
    protected void appendSubclassState(StringBuilder sb) {
        sb.append("lambdaScale=").append(Double.doubleToLongBits(lambdaScale)).append('|');
        sb.append("stringency=").append(Double.doubleToLongBits(stringency)).append('|');
    }

    @Override
    public String toToml() {
        StringBuilder builder = getBaseTomlBuilder();
        builder.append(String.format("lambdaScale = %s%n", lambdaScale));
        builder.append(String.format("stringency = %s%n", stringency));
        return builder.toString();
    }

    // ── RelaxFit initialization ───────────────────────────────────────────────

    /**
     * {@inheritDoc}
     *
     * <p>In addition to the base initialization, enables lambda regularization
     * on the {@link RelaxFit} instance and sets the three regularization
     * strengths.</p>
     */
    @Override
    public RelaxFit initRelaxFit(String key, MolDataValues<? extends RelaxDataValue> data) {
        RelaxFit relaxFit = super.initRelaxFit(key, data);
        relaxFit.setUseLambda(true);
        relaxFit.setLambdas(getLambdaScale());
        return relaxFit;
    }

    /**
     * Processes output of 2sf fit.
     *
     * <ul>
     *   <li><em>No local motions</em> (both S² above {@code S2_THOLD}):
     *       S²f = S²s = 1, τf = τs = 0.</li>
     *   <li><em>One fast motion</em> (τ below {@code TAU_THOLD}):
     *       S²f = S², τf = 0, S²s = 1, τs = 0.</li>
     *   <li><em>One fast motion with non-zero τ_f</em>
     *       (τ in (TAU_THOLD, SLOW_LIMIT]):
     *       S²f = S², τf = τ, S²s = 1, τs = 0.</li>
     *   <li><em>One slow motion</em> (τ above {@code SLOW_LIMIT}):
     *       S²f = 1, τf = 0, S²s = S², τs = τ.</li>
     *   <li><em>Two instantaneous motions</em> (both τ below {@code TAU_THOLD}):
     *       S²f = S²1 × S²2, τf = τs = 0, S²s = 1.</li>
     *   <li><em>Two motions, one instantaneous</em>:
     *       The instantaneous timescale is mapped to τf = 0.</li>
     *   <li><em>Two independent motions</em>:
     *       Sorted so that τf &lt; τs.</li>
     * </ul>
     *
     * @param params parameter values from the optimizer
     * @return the processed parameter array
     */
    protected double[] processParamsAfterFit(double[] params, boolean fitTau, double slowLimit) {
        int start = fitTau ? 1 : 0;
        double s1 = params[start];      // sf2
        double tau1 = params[start + 1];  // tau_f
        double s2 = params[start + 2];  // ss2
        double tau2 = params[start + 3];  // tau_s

        double sf2, tauf, ss2, taus;

        if (moietyType == MoietyType.AMIDE) {
            // For amide (sN=1): sf2≈1 means no fast-motion contribution.
            // The two-motion spectral density terms are symmetric under a
            // parameter swap, so sorting tau_f ≤ tau_s is valid.
            if (s1 > S2_THOLD && s2 > S2_THOLD) {
                // No local motions
                sf2 = ss2 = 1.0;
                tauf = taus = 0.0;
            } else if (s1 > S2_THOLD || s2 > S2_THOLD) {
                // One local motion
                double s = (s1 > S2_THOLD) ? s2 : s1;
                double tau = (s1 > S2_THOLD) ? tau2 : tau1;
                if (tau < TAU_THOLD) {
                    sf2 = s;
                    tauf = 0.0;
                    ss2 = 1.0;
                    taus = 0.0;
                } else if (tau < slowLimit) {
                    sf2 = s;
                    tauf = tau;
                    ss2 = 1.0;
                    taus = 0.0;
                } else {
                    sf2 = 1.0;
                    tauf = 0.0;
                    ss2 = s;
                    taus = tau;
                }
            } else {
                if (tau1 < TAU_THOLD && tau2 < TAU_THOLD) {
                    // Both motions are instantaneous: collapse to effective order parameter (Model 1)
                    sf2 = s1 * s2;
                    tauf = 0.0;
                    ss2 = 1.0;
                    taus = 0.0;
                } else if (tau1 < TAU_THOLD) {
                    // tau1 is instantaneous: assign it to the fast slot (Model 2s)
                    sf2 = s1;
                    tauf = 0.0;
                    ss2 = s2;
                    taus = tau2;
                } else if (tau2 < TAU_THOLD) {
                    // tau2 is instantaneous: assign it to the fast slot (Model 2s)
                    sf2 = s2;
                    tauf = 0.0;
                    ss2 = s1;
                    taus = tau1;
                } else {
                    // Both timescales are resolvable: sort so that tauf < taus (Model 2sf)
                    if (tau1 < tau2) {
                        sf2 = s1;
                        tauf = tau1;
                        ss2 = s2;
                        taus = tau2;
                    } else {
                        sf2 = s2;
                        tauf = tau2;
                        ss2 = s1;
                        taus = tau1;
                    }
                }
            }
        } else {
            // For deuterium (sN=9): sf2/sN is never near zero for physical sf2
            // values, so sf2≈1 does NOT suppress tau_f.  The spectral density
            // terms are asymmetric in sf2/ss2 due to sN, so sorting is invalid.
            // Only suppress slow motion if ss2≈1, or if tau_s < tau_f (unphysical
            // ordering that cannot be resolved by sorting).
            if (s2 > S2_THOLD || (tau2 > TAU_THOLD && tau2 < tau1)) {
                sf2 = s1;
                tauf = tau1;
                ss2 = 1.0;
                taus = 0.0;
            } else {
                sf2 = s1;
                tauf = tau1;
                ss2 = s2;
                taus = tau2;
            }
        }

        params[start] = sf2;
        params[start + 1] = tauf;
        params[start + 2] = ss2;
        params[start + 3] = taus;

        return params;
    }

    public record FitOnceResult(Score score, double[] crlb) {}

    public FitOnceResult doFit(MFModelIso2sf model, RelaxFit relaxFit, MolDataValues<? extends RelaxDataValue> replicateData, String key, double[] start, int nTry) {
        CoMDOptions options = new CoMDOptions(true);
        model.applyThreshold(null);
        relaxFit.setRelaxData(key, replicateData);


        relaxFit.setLambdas(0.0);
        double[] crlb0 = relaxFit.calcCRLB(replicateData, model);
        model.updateCRLB(crlb0);
        Score unpen = runFit(relaxFit, model, null, nTry);
        double[] up = unpen.getPars();

        relaxFit.setLambdas(lambdaScale);          // the builder's lambdaScale
        model.pars(up);
        model.updateTauWeights();          // w = c'(tau_unpenalized)  ← the whole point
        double[] crlb = relaxFit.calcCRLB(replicateData, model, up);
        model.updateCRLB(crlb);
        Score score = runFit(relaxFit, model, up, 1);


        crlb = relaxFit.calcCRLB(replicateData, model, score.pars);
        model.updateCRLB(crlb);
        model.updateTauWeights();

        score = runFit(relaxFit, model, score.pars, 1);
        crlb = relaxFit.calcCRLB(replicateData, model, score.pars);
        model.updateTauWeights();

        score = runFit(relaxFit, model, score.pars, 1);
        crlb = relaxFit.calcCRLB(replicateData, model, score.pars);
        MFModelIso2sf.ThresholdedPars tPars = model.calcThreshold(crlb, stringency);
        if (tPars.anyChanged()) {
            model.applyThreshold(tPars);
            double[] pars = model.getPars();
            score = runFit(relaxFit, model, pars, 1);
            crlb = relaxFit.calcCRLB(replicateData, model, score.pars);
        }
        return new FitOnceResult(score, crlb);
    }

    /**
     * Performs a regularized model-free fit for a single residue.
     *
     * <p>The algorithm proceeds as follows:</p>
     * <ol>
     *   <li>For each of the {@code nReplicates} bootstrap replicates:
     *     <ol type="a">
     *       <li>Draw a bootstrap sample from {@code data}.</li>
     *       <li>Fit the {@code 2sf} model to the sample with regularization.</li>
     *       <li>Post-process the parameters via
     *           {@link #processParamsAfterFit(double[], boolean)}.</li>
     *       <li>Record the parameters and bootstrap weights.</li>
     *     </ol>
     *   </li>
     *   <li>Aggregate parameters and weights across replicates using
     *       {@link #computeStatistics(double[][], double[][])} to obtain
     *       final estimates and errors.</li>
     *   <li>Construct and register an {@link OrderPar} under the key
     *       {@code "REGULARIZATION"}.</li>
     * </ol>
     *
     * <p><strong>Note on Score:</strong> The {@link Score} stored in the
     * returned {@link ModelFitResult} and passed to {@link #makeOrderPar} is
     * taken from the first replicate. A meaningful aggregate score has not yet
     * been defined.</p>
     *
     * @param key            identifier for the relaxation data set
     * @param data           measured relaxation values for the residue
     * @param orderParSetMap mutable map to which the resulting
     *                       {@link OrderParSet} entry will be added
     * @return the {@link ModelFitResult} with aggregated parameters and errors
     * @throws IllegalStateException if tau_M has not been set
     */
    @Override
    public ModelFitResult fit(String key, MolDataValues<?> data, Map<String, OrderParSet> orderParSetMap) {
        RelaxFit relaxFit = initRelaxFit(key, data);
        MFModelIso2sf model = (MFModelIso2sf) getModel("2sf", data);

        data.setTestModel(model);
        CoMDOptions options = new CoMDOptions(true);
        int nTry = options.getNTries();

        FitOnceResult initialFitResult = doFit(model, relaxFit, data, key, null, nTry);

        int nParameters = model.getNPars();
        int nWeights = data.getNSpectralDensities();
        double[][] parameters = new double[nParameters][nReplicates];
        double[][] weights = new double[nWeights][nReplicates];
        BootstrapSampler<? extends RelaxDataValue> sampler = getBootstrapSampler(data);

        Score[] scores = new Score[nReplicates];
        double[] replicateTimes = new double[nReplicates];
        double[] start = null;
        double[] crossResiduals = new double[nReplicates];
        for (int i = 0; i < nReplicates; i++) {
            long startNs = System.nanoTime();
            MolDataValues<? extends RelaxDataValue> replicateData = sampler.sample();
            FitOnceResult fitOnceResult = doFit(model, relaxFit, replicateData, key, start, nTry);
            scores[i] = fitOnceResult.score;

            start = scores[i].pars.clone();
            double[] replicateParameters = processParamsAfterFit(scores[i].getPars(), model.fitTau(),model.slowLimit());
            double[] replicateWeights = replicateData.getWeights();
            for (int k = 0; k < nParameters; k++) {
                parameters[k][i] = replicateParameters[k];
            }
            for (int j = 0; j < nWeights; j++) {
                weights[j][i] = replicateWeights[j];
            }
            replicateTimes[i] = (System.nanoTime() - startNs) / 1_000_000.0;
            relaxFit.setRelaxData(key, data);
            crossResiduals[i] = relaxFit.maxNormalizedResidual(replicateParameters);
        }

        Pair<double[], double[]> parameterEstimates = computeStatistics(parameters, weights);
        double[] fitParameters = initialFitResult.score.pars;
        double[] fitErrors = parameterEstimates.getRight();

        orderParSetMap.computeIfAbsent(KEY, ky -> new OrderParSet(ky));
        // FIXME: the Score used here (scores[0]) is from the first replicate.
        // For bootstrap fitting, there is no single score over the original
        // data; this should be revisited.
        OrderPar orderPar = makeOrderPar(
                orderParSetMap.get(KEY),
                sampler.getOriginalData(),
                key,
                scores[0],
                model,
                fitParameters,
                fitErrors
        );

        return new ModelFitResult(orderPar, parameters, null, replicateTimes, flagSpuriousReplicates(crossResiduals));
    }
}
