package org.comdnmr.datasets;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * Generates a simulated relaxation dataset CSV (same form as e.g.
 * datasets/GCN4_NH_600-700-800-900.csv) from a list of per-residue
 * model-free parameter sets and a list of spectrometer fields. For each
 * (residue, field) pair, "true" R1, R2 and NOE values are computed from the
 * model-free parameters, then corrupted with noise and re-derived via
 * {@link RelaxationSimulator}, mimicking a real experimental measurement.
 */
class SimulatedDataset {

    private static final String HEADER = "Residue,ResName,Field,R1,R1_err,R2,R2_err,NOE,NOE_err";
    private static final String DHEADER = "Residue,ResName,Field,R1,R1_err,R2,R2_err,RQ,RQ_err,RAP,RAP_err";

    private SimulatedDataset() {}

    static void generate(
        List<ParameterSet> parameterSets,
        String modelName,
        double[] fields,
        double noiseLevel,
        long seed,
        Path outputFile
    ) throws IOException {
        Random rng = new Random(seed);
        List<ParameterSet> sorted = parameterSets.stream()
            .sorted(Comparator.comparingInt(ParameterSet::residueNumber))
            .toList();

        boolean deuteriumMode = modelName.startsWith("D");
        StringBuilder csv = deuteriumMode ? new StringBuilder(DHEADER) : new StringBuilder(HEADER);
        csv.append("\n");
        for (ParameterSet params : sorted) {
            for (double field : fields) {
                ParameterSet.PredictionInterface truth = params.predictRelaxation(modelName, field, deuteriumMode);
                if (truth instanceof ParameterSet.Prediction predTruth) {
                    RelaxationSimulator.FitResult r1 = RelaxationSimulator.simulateExpDecay(predTruth.r1(), noiseLevel, rng);
                    RelaxationSimulator.FitResult r2 = RelaxationSimulator.simulateExpDecay(predTruth.r2(), noiseLevel, rng);
                    RelaxationSimulator.FitResult noe = RelaxationSimulator.simulateNOE(predTruth.noe(), noiseLevel, rng);
                    csv.append(String.format(
                            Locale.ROOT,
                            "%d,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f%n",
                            params.residueNumber(), params.residueName(), field,
                            r1.value(), r1.error(), r2.value(), r2.error(), noe.value(), noe.error()
                    ));
                } else if (truth instanceof ParameterSet.PredictionDeuterium predictionDeuterium) {
                    RelaxationSimulator.FitResult r1 = RelaxationSimulator.simulateExpDecay(predictionDeuterium.r1(), noiseLevel, rng);
                    RelaxationSimulator.FitResult r2 = RelaxationSimulator.simulateExpDecay(predictionDeuterium.r2(), noiseLevel, rng);
                    RelaxationSimulator.FitResult rQ = RelaxationSimulator.simulateExpDecay(predictionDeuterium.rQ(), noiseLevel, rng);
                    RelaxationSimulator.FitResult rAP = RelaxationSimulator.simulateExpDecay(predictionDeuterium.rAP(), noiseLevel, rng);
                    csv.append(String.format(
                            Locale.ROOT,
                            "%d,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f%n",
                            params.residueNumber(), params.residueName(), field,
                            r1.value(), r1.error(), r2.value(), r2.error(), rQ.value(), rQ.error(), rAP.value(), rAP.error()
                    ));
                }
            }
        }

        Path parent = outputFile.toAbsolutePath().getParent();
        if (parent != null) {
            Files.createDirectories(parent);
        }
        Files.writeString(outputFile, csv.toString());
    }
}
