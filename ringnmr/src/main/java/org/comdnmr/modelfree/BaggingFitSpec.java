package org.comdnmr.modelfree;

import java.util.List;
import java.util.Map;
import java.util.Optional;

import org.apache.commons.lang3.tuple.Pair;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.modelfree.models.MFModelIso2sf;

import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.MoleculeBase;
import org.nmrfx.chemistry.MoleculeFactory;
import org.nmrfx.chemistry.relax.OrderPar;
import org.nmrfx.chemistry.relax.OrderParSet;
import org.nmrfx.chemistry.relax.ResonanceSource;
import org.nmrfx.chemistry.relax.SpectralDensity;

public class BaggingFitSpec extends ModelSelectionFitSpec {

    private static final String KEY = "BAGGING";

    BaggingFitSpec(Builder builder) {
        super(builder);
    }

    public ModelFitResult fit(String key, MolDataValues data) {
        RelaxFit relaxFit = new RelaxFit();
        // FIXME: currently weighting doesn't apply to R1/R2/NOE fitting
        relaxFit.setFitJ(true);

        List<MFModelIso> models = getModels(data);
        MFModelIso2sf model2sf = (MFModelIso2sf) getModel("2sf", data);
        double[][] replicates = new double[nReplicates][5];
        BootstrapSampler sampler = getBootstrapSampler(data);

        MoleculeBase moleculeBase = MoleculeFactory.getActive();
        Map<String, OrderParSet> orderParSetMap = moleculeBase.orderParSetMap();
        orderParSetMap.computeIfAbsent(KEY, ky -> new OrderParSet(ky));

        Score[] bestScores = new Score[nReplicates];
        for (int i = 0; i < nReplicates; i++) {
            MolDataValues replicateData = sampler.sample();
            relaxFit.setRelaxData(key, replicateData);
            Optional<Pair<MFModelIso, Score>> bestModelScore = Optional.empty();
            for (MFModelIso model : models) {
                data.setTestModel(model);
                Score score = runFit(relaxFit, model);
                if (
                    bestModelScore.isEmpty() ||
                    score.aicc().get() < bestModelScore.get().getRight().aicc().get()
                ) {
                    bestModelScore = Optional.of(Pair.of(model, score));
                }
            }

            if (bestModelScore.isEmpty()) {
                throw new AssertionError("`bestModelScore` should not be empty here!");
            }

            MFModelIso bestModel = bestModelScore.get().getLeft();
            Score bestScore = bestModelScore.get().getRight();
            bestScores[i] = bestScore;
            double[] parameters = bestModel.getStandardPars(bestScore.getPars());
            for (int k = 0; k < parameters.length; k++) {
                replicates[i][k] = parameters[k];
            }
        }

        // FIXME: not sure what Score should be for makeOrderParSet...
        OrderPar orderPar = makeOrderParSet(
            orderParSetMap.get(KEY),
            sampler.getOriginalData(),
            key,
            bestScores[0],
            model2sf,
            replicates
        );

        ModelFitResult result = new ModelFitResult(
            orderPar,
            replicates,
            null,
            bestScores
        );

        return result;
    }

    public static class Builder extends ModelSelectionFitSpec.Builder {

        public BaggingFitSpec build() {
            return new BaggingFitSpec(this);
        }
    }
}
