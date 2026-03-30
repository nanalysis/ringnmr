package org.comdnmr.modelfree;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.comdnmr.modelfree.models.MFModelIso;
import org.nmrfx.chemistry.Atom;
import org.nmrfx.chemistry.relax.OrderPar;
import org.nmrfx.chemistry.relax.OrderParSet;
import org.nmrfx.chemistry.relax.ResonanceSource;
import org.nmrfx.chemistry.relax.SpectralDensity;

public abstract class ModelSelectionFitSpec extends FitSpec {

    private final List<String> modelNames;

    protected ModelSelectionFitSpec(Builder builder) {
        super(builder);
        if (builder.modelNames.isEmpty()) {
            throw new IllegalArgumentException("`modelNames` must not be empty!");
        }
        this.modelNames = builder.modelNames;
    }

    public List<String> getModelNames() {
        return modelNames;
    }

    public MFModelIso getModel(String name, MolDataValues data) {
        boolean fitTauM = fitTauM(data);
        return MFModelIso.buildModel(name, fitTauM, tauM, tauFraction, fitExchange);
    }

    public List<MFModelIso> getModels(MolDataValues data) {
        return getModelNames()
            .stream()
            .map(name -> getModel(name, data))
            .collect(Collectors.toList());
    }

    protected OrderPar makeOrderParSet(
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


    // Concrete Builder for ConcreteFitSpec
    public static abstract class Builder extends FitSpec.Builder<Builder> {

        private final List<String> modelNames = new ArrayList<>();

        public Builder modelNames(List<String> modelNames) {
            this.modelNames.addAll(modelNames);
            return this;
        }

        abstract public ModelSelectionFitSpec build();
    }
}
