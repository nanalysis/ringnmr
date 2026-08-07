package org.comdnmr.datasets;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.comdnmr.modelfree.RelaxEquations;
import org.comdnmr.modelfree.models.MFModelIso;

/**
 * A set of model-free parameters for a single residue, in SI units
 * (correlation times in seconds). The set of parameters required beyond the
 * overall correlation time (tauM) depends on the chosen model ("1", "1f",
 * "1s", "2s" or "2sf") and is determined directly from the corresponding
 * MFModelIso subclass (see {@link #requiredColumns(String)}).
 */
public record ParameterSet(int residueNumber, String residueName, double tauM, Map<String, Double> modelParams) {

    public interface PredictionInterface {

    }

    public record Prediction(double r1, double r2, double noe) implements  PredictionInterface {}
    public record PredictionDeuterium(double r1, double r2, double rQ, double rAP) implements  PredictionInterface {}

    /**
     * The CSV columns required to describe a parameter set for a given
     * model, inferred from that model's own parameter names (as reported by
     * {@link MFModelIso#getParNames()}).
     *
     * @param modelName String. One of "1", "1f", "1s", "2s", "2sf".
     * @return List&lt;String&gt;. The required column names, in a sensible order.
     */
    public static List<String> requiredColumns(String modelName) {
        MFModelIso model = MFModelIso.buildModel(modelName, false, 1.0, 0.0, false);
        List<String> columns = new ArrayList<>(List.of("Residue", "ResName", "TauM"));
        columns.addAll(model.getParNames());
        return columns;
    }

    /**
     * Compute the "true" R1, R2 and NOE relaxation values implied by this
     * parameter set at a given 1H spectrometer frequency, using the same
     * model-free spectral density classes as are used elsewhere in this
     * package for fitting.
     *
     * @param modelName String. One of "1", "1f", "1s", "2s", "2sf".
     * @param fieldMHz double. 1H spectrometer frequency, in MHz.
     * @return Prediction. The predicted R1, R2 and NOE values.
     */
    public PredictionInterface predictRelaxation(String modelName, double fieldMHz, boolean deuteriumMode) {
        RelaxEquations relaxObj;
        if (deuteriumMode) {
            relaxObj = RelaxEquations.getRelaxEquations(fieldMHz * 1e6, "D", "C");
        }  else {
            relaxObj = RelaxEquations.getRelaxEquations(fieldMHz * 1.0e6, "H", "N");
        }

        MFModelIso model = MFModelIso.buildModel(modelName, false, tauM * 1.0e9, 0.0, false);

        List<String> parNames = model.getParNames();
        double[] pars = new double[parNames.size()];
        for (int i = 0; i < parNames.size(); i++) {
            String name = parNames.get(i);
            Double value = modelParams.get(name);
            if (value == null) {
                throw new IllegalArgumentException(String.format(
                    "Residue %d: missing value for parameter '%s', required by model '%s'",
                    residueNumber, name, modelName
                ));
            }
            // Correlation times are supplied in seconds; MFModelIso classes expect nanoseconds.
            pars[i] = name.startsWith("Tau") ? value * 1.0e9 : value;
        }

        double[] jValues = model.calc(relaxObj.getOmegas(), pars);
        if (!model.checkParConstraints()) {
            throw new IllegalArgumentException(String.format(
                "Residue %d: parameters are not physically valid for model '%s' " +
                "(internal correlation times must be smaller than TauM)",
                residueNumber, modelName
            ));
        }

        double r1 = relaxObj.R1(jValues);
        double r2 = relaxObj.R2(jValues, 0.0);
        if (deuteriumMode) {
            double rQ = relaxObj.RQ_D(jValues);
            double rAP = relaxObj.Rap_D(jValues);
            return new PredictionDeuterium(r1, r2, rQ, rAP);
        } else {
            double noe = relaxObj.NOE(jValues);
            return new Prediction(r1, r2, noe);
        }

    }

    public static List<ParameterSet> fromFile(File file, String modelName) throws IOException {
        List<String> required = requiredColumns(modelName);
        List<String> paramNames = required.subList(3, required.size());

        List<ParameterSet> result = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
            String header = reader.readLine();
            if (header == null) {
                throw new IllegalArgumentException("Parameter set file is empty: " + file.getName());
            }
            List<String> columns = List.of(header.split(","));
            Map<String, Integer> indices = new HashMap<>();
            for (String name : required) {
                int index = columns.indexOf(name);
                if (index == -1) {
                    throw new IllegalArgumentException(String.format(
                        "Parameter set file must have a header containing the columns: %s (model '%s'). Missing: %s",
                        String.join(",", required), modelName, name
                    ));
                }
                indices.put(name, index);
            }

            String line;
            int lineNumber = 1;
            while ((line = reader.readLine()) != null) {
                lineNumber++;
                if (line.isBlank()) continue;
                String[] fields = line.split(",");
                try {
                    int residueNumber = Integer.parseInt(fields[indices.get("Residue")].trim());
                    String residueName = fields[indices.get("ResName")].trim();
                    double tauM = Double.parseDouble(fields[indices.get("TauM")].trim());
                    Map<String, Double> modelParams = new HashMap<>();
                    for (String name : paramNames) {
                        modelParams.put(name, Double.parseDouble(fields[indices.get(name)].trim()));
                    }
                    result.add(new ParameterSet(residueNumber, residueName, tauM, modelParams));
                } catch (NumberFormatException | ArrayIndexOutOfBoundsException e) {
                    throw new IllegalArgumentException(String.format(
                        "Error parsing parameter set file %s at line %d: %s", file.getName(), lineNumber, e.getMessage()
                    ));
                }
            }
        }
        return result;
    }
}
