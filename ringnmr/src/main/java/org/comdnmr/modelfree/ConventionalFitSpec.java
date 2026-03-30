package org.comdnmr.modelfree;

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
        RelaxFit relaxFit = new RelaxFit();
        // FIXME: currently weighting doesn't apply to R1/R2/NOE fitting
        relaxFit.setFitJ(true);
        relaxFit.setRelaxData(key, data);

        // Determine the optimal model using the AICc
        List<MFModelIso> models = getModels(data);
        MFModelIso bestModel = null;
        Score bestScore = null;
        double aicc = Double.MAX_VALUE;
        for (MFModelIso model : models) {
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

        data.setTestModel(bestModel);
        for (int i = 0; i < nReplicates; i++) {
            MolDataValues replicateData = sampler.sample();
            relaxFit.setRelaxData(key, replicateData);
            Score score = runFit(relaxFit, bestModel);
            double[] parameters = score.getPars();
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

}
