package org.comdnmr.modelfree;

import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.MFModelIso2sf;

/**
 * Abstract intermediate class for model-free fitting strategies that use
 * bootstrap resampling and express results in the extended model-free
 * ({@code 2sf}) form.
 *
 * <p>Currently extended by {@link BaggingFitSpec} and
 * {@link RegularizationFitSpec}. The shared
 * {@link #processParamsAfterFit(MFModelIso, double[])} logic converts the raw
 * optimizer output from any candidate model into the unified
 * {@code [tauM?,] S²f, τf, S²s, τs} layout.</p>
 */
abstract class BootstrappedFitSpec extends FitSpec {

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
    static final double TAU_THOLD = 1.0e-3;

    /**
     * Constructs a new {@code BootstrappedFitSpec} from the given builder.
     *
     * @param builder the builder containing all configuration values
     */
    protected BootstrappedFitSpec(FitSpec.Builder<?> builder) {
        super(builder);
    }

    /**
     * Converts the raw optimizer output from any candidate model into the
     * canonical two-slow-motions ({@code 2sf}) parameter representation:
     * {@code [tauM?,] sf², τ_f, ss², τ_s}.
     *
     * <p>The method first calls {@link MFModelIso#getStandardPars(double[])} to
     * extend the parameters from any simpler model to the full 4-parameter
     * (or 5-parameter when tau_M is free) 2sf layout, then classifies the
     * result into one of the following forms:</p>
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
     *   <li><em>Two motions, one instantaneous</em>:
     *       The instantaneous timescale is mapped to τf = 0.</li>
     *   <li><em>Two independent motions</em>:
     *       Sorted so that τf &lt; τs.</li>
     * </ul>
     *
     * @param model  the model whose parameters are being converted; used to
     *               obtain the standard 2sf-layout parameters via
     *               {@link MFModelIso#getStandardPars(double[])}
     * @param params raw parameter values from the optimizer, in the native
     *               layout of {@code model}
     * @return the (possibly re-ordered) parameter array in canonical 2sf layout
     */
    @Override
    protected double[] processParamsAfterFit(MFModelIso model, double[] params) {
        params = model.getStandardPars(params);
        boolean fitTau = model.fitTau();
        int start = fitTau ? 1 : 0;
        double s1   = params[start];
        double tau1 = params[start + 1];
        double s2   = params[start + 2];
        double tau2 = params[start + 3];

        double sf2, tauf, ss2, taus;

        // No local motions — both order parameters are effectively 1
        if (s1 > S2_THOLD && s2 > S2_THOLD) {
            sf2 = ss2 = 1.0;
            tauf = taus = 0.0;
        }

        // One local motion — the other degree of freedom is frozen out
        else if (s1 > S2_THOLD || s2 > S2_THOLD) {
            double s   = (s1 > S2_THOLD) ? s2 : s1;
            double tau = (s1 > S2_THOLD) ? tau2 : tau1;
            if (tau < TAU_THOLD) {
                // Fast motion with effectively zero correlation time (Model 1)
                sf2 = s;   tauf = 0.0; ss2 = 1.0; taus = 0.0;
            } else if (tau < MFModelIso2sf.SLOW_LIMIT) {
                // Resolvable fast motion (Model 1f)
                sf2 = s;   tauf = tau; ss2 = 1.0; taus = 0.0;
            } else {
                // Slow motion (Model 1s)
                sf2 = 1.0; tauf = 0.0; ss2 = s;   taus = tau;
            }
        }

        // Two local motions
        else {
            if (tau1 < TAU_THOLD) {
                // tau1 is instantaneous: assign it to the fast slot (Model 2s)
                sf2 = s1; tauf = 0.0; ss2 = s2; taus = tau2;
            } else if (tau2 < TAU_THOLD) {
                // tau2 is instantaneous: assign it to the fast slot (Model 2s)
                sf2 = s2; tauf = 0.0; ss2 = s1; taus = tau1;
            } else {
                // Both timescales are resolvable: sort so that tauf < taus (Model 2sf)
                if (tau1 < tau2) {
                    sf2 = s1; tauf = tau1; ss2 = s2; taus = tau2;
                } else {
                    sf2 = s2; tauf = tau2; ss2 = s1; taus = tau1;
                }
            }
        }

        params[start]     = sf2;
        params[start + 1] = tauf;
        params[start + 2] = ss2;
        params[start + 3] = taus;

        return params;
    }
}
