package org.comdnmr.modelfree;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import org.comdnmr.modelfree.models.MFModelIso;

import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.relax.OrderPar;
import org.nmrfx.chemistry.relax.OrderParSet;
import org.nmrfx.chemistry.relax.ResonanceSource;
import org.nmrfx.chemistry.relax.SpectralDensity;


public class ConventionalFitSpec extends ModelSelectionFitSpec {

    private static final String BEST_KEY = "BEST";

    public static class Builder extends ModelSelectionFitSpec.Builder {
        public ConventionalFitSpec build() { return new ConventionalFitSpec(this); }
    }

    ConventionalFitSpec(Builder builder) { super(builder); }

    public ModelFitResult fit(String key, MolDataValues data) {
        double localTauFraction = getLocalTauFraction(data);
        RelaxFit relaxFit = new RelaxFit();
        relaxFit.setRelaxData(key, data);

        // Determine the optimal model using the AICc
        List<MFModelIso> models = getModels();
        MFModelIso bestModel = null;
        Score bestScore = null;
        double aicc = Double.MAX_VALUE;
        for (MFModelIso model : models) {
            model.setTauFraction(localTauFraction);
            data.setTestModel(model);
            Score score = runFit(relaxFit, model);

            // TODO: currently assuming that the number of datapoints is
            // sufficiently large for AICc to be computed (i.e. just using
            // Optional::get without error handling)
            if (score.aicc().get() < aicc) {
                bestModel = model;
                bestScore = score;
                aicc = score.aicc().get();
            }
        }

        if (bestModel == null || bestScore == null) {
            throw new AssertionError("`bestModel` and `bestScore` should not be null here!");
        }

        double[][] replicates = new double[nReplicates][bestModel.getNPars()];
        BootstrapSampler sampler = getBootstrapSampler(data);

        MoleculeBase moleculeBase = MoleculeFactory.getActive();
        Map<String, OrderParSet> orderParSetMap = moleculeBase.orderParSetMap();
        orderParSetMap.computeIfAbsent(BEST_KEY, ky -> new OrderParSet(ky));

        System.out.printf("key: %s%n", key);
        System.out.printf("bestModel: %s%n", bestModel);
        data.setTestModel(bestModel);
        for (int i = 0; i < nReplicates; i++) {
            MolDataValues replicateData = sampler.sample();
            relaxFit.setRelaxData(key, replicateData);
            Score score = runFit(relaxFit, bestModel);
            double[] parameters = score.getPars();
            System.out.printf("parameters:%n%s%n", Arrays.toString(parameters));
            for (int k = 0; k < parameters.length; k++) {
                replicates[i][k] = parameters[k];
            }
        }

        OrderPar orderPar = makeOrderParSet(
            orderParSetMap.get(BEST_KEY),
            sampler.getOriginalData(),
            key,
            bestScore,
            bestModel,
            replicates
        );

        ModelFitResult result = new ModelFitResult(
            orderPar,
            replicates,
            null,
            new Score[]{bestScore}
        );

        return result;
    }

    private OrderPar makeOrderParSet(
        OrderParSet orderParSet,
        MolDataValues data,
        String key,
        Score score,
        MFModelIso model,
        double[][] replicates
    ) {
        ResonanceSource resSource = new ResonanceSource(data.atom);
        Atom atom = data.atom;
        List<String> parameterNames = model.getParNames();
        int nParameters = parameterNames.size();
        OrderPar orderPar = new OrderPar(orderParSet, resSource, score.rss, score.nValues, score.nPars, model.getName());
        double[] parameterBootstraps = new double[nReplicates];
        for (int k = 0; k < nParameters; k++) {
            String parameterName = parameterNames.get(k);
            for (int i = 0; i < nReplicates; i++) {
                parameterBootstraps[i] = replicates[i][k];
            }
            DescriptiveStatistics statistics = new DescriptiveStatistics(parameterBootstraps);
            double parameterMean = statistics.getMean();
            double parameterError = statistics.getStandardDeviation();
            orderPar = orderPar.set(parameterName, parameterMean, parameterError);
        }

        // If tauM is fixed, set it to the fixed value with zero error
        if (!model.fitTau())  orderPar = orderPar.set("Tau_e", model.getTau(), 0.0);

        orderPar = orderPar.set("model", (double) model.getNumber(), null);
        atom.addOrderPar(orderParSet, orderPar);

        SpectralDensity spectralDensity = new SpectralDensity(key, data.getJValues());
        atom.addSpectralDensity(key, spectralDensity);

        return orderPar;
    }
}
