package org.comdnmr.datasets;

import org.comdnmr.modelfree.BaggingFitSpec;
import org.comdnmr.modelfree.ConventionalFitSpec;
import org.comdnmr.modelfree.FitSpec;
import org.comdnmr.modelfree.R1R2NOEMolDataValues;
import org.comdnmr.modelfree.RegularizationFitSpec;
import org.comdnmr.modelfree.models.MFModelIso;
import org.comdnmr.util.CoMDPreferences;

import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

public class RingCmd {

    private static final Set<String> ALLOWED_MODELS = new HashSet<>();

    static {
        for (String model : MFModelIso.getAllModelNames()) {
            if (!model.equals("1sf")) ALLOWED_MODELS.add(model);
        }
    }

    public static void main(String[] args) throws IOException, ParseException {
        CoMDPreferences.setOptimizer("CMAES");

        if (args.length > 0 && args[0].equals("simulate")) {
            runSimulate(Arrays.copyOfRange(args, 1, args.length));
        } else if (args.length > 0 && args[0].equals("fit")) {
            runFit(Arrays.copyOfRange(args, 1, args.length));
        }
    }

    public static void printUsage() {
        HelpFormatter hf = new HelpFormatter();
        String syntax = "ringcmd fit|simulate <data-file> <output-dir> [options]";
        hf.printHelp(syntax, buildCliOptions());
    }

    private static void runFit(String[] args) throws IOException {

        try {
            CommandLineParser parser = new DefaultParser();
            Options options = buildCliOptions();
            CommandLine cmd = parser.parse(options, args);
            System.out.println("cmd " + cmd + " args " + cmd.getArgs());

            String[] positional = cmd.getArgs();
            for (var pos : positional) {
                System.out.println("par " + pos);
            }

            if (positional.length != 2) {
                throw new ParseException("Expected exactly two positional arguments: <data-file> <output-dir>");
            }
            File dataFile = new File(positional[0]);
            Path outDir = Path.of(positional[1]);

            if (!dataFile.exists() || !dataFile.isFile()) {
                throw new ParseException(String.format("Data file does not exist or is not a file: %s", dataFile));
            }
            if (Files.exists(outDir) && !Files.isDirectory(outDir)) {
                throw new ParseException(String.format("Output path exists and is not a directory: %s", outDir));
            }

            Dataset dataset = Dataset.fromFile(dataFile);
            List<FitSpec> fitSpec = List.of(getFitSpec(cmd, dataset));
            dataset.saveToToml(outDir);
            dataset.fit(fitSpec, outDir);
        } catch (ParseException e) {
            System.err.println("Error: " + e.getMessage());
            printUsage();
            System.exit(2);
            return;
        }
    }

    private static void runSimulate(String[] args) {
        try {
            CommandLineParser parser = new DefaultParser();
            Options options = buildSimulateCliOptions();
            CommandLine cmd = parser.parse(options, args);

            String[] positional = cmd.getArgs();
            for (var pos : positional) {
                System.out.println("pars " + pos);
            }
            if (positional.length != 2) {
                throw new ParseException("Expected exactly two positional arguments: <params-file> <output-file>");
            }
            File paramsFile = new File(positional[0]);
            Path outputFile = Path.of(positional[1]);

            if (!paramsFile.exists() || !paramsFile.isFile()) {
                throw new ParseException(String.format("Parameter set file does not exist or is not a file: %s", paramsFile));
            }

            if (!cmd.hasOption("model")) {
                throw new ParseException("--model is required");
            }
            String modelName = cmd.getOptionValue("model");
            if (!ALLOWED_MODELS.contains(modelName)) {
                throw new ParseException(String.format(
                        "Invalid --model specified: %s.%nMust be one of: %s",
                        modelName, String.join(", ", ALLOWED_MODELS)
                ));
            }

            if (cmd.hasOption("deuterium")) {
                modelName = "D" + modelName;
            }

            if (!cmd.hasOption("fields")) {
                throw new ParseException("--fields is required");
            }
            double[] fields = parseFields(cmd.getOptionValue("fields"));

            double noise = cmd.hasOption("noise") ? parsePositiveDouble(cmd, "noise", true) : 0.02;
            long seed = cmd.hasOption("seed") ? Long.parseLong(cmd.getOptionValue("seed")) : 1L;

            List<ParameterSet> parameterSets = ParameterSet.fromFile(paramsFile, modelName);
            if (parameterSets.isEmpty()) {
                throw new ParseException("No parameter sets found in: " + paramsFile);
            }

            System.out.println("model " + modelName);
            SimulatedDataset.generate(parameterSets, modelName, fields, noise, seed, outputFile);
        } catch (ParseException e) {
            System.err.println("Error: " + e.getMessage());
            printSimulateUsage();
            System.exit(2);
        } catch (IllegalArgumentException e) {
            System.err.println("Error: " + e.getMessage());
            System.exit(2);
        } catch (IOException e) {
            System.err.println("Error writing output file: " + e.getMessage());
            System.exit(1);
        }
    }

    public static void printSimulateUsage() {
        HelpFormatter hf = new HelpFormatter();
        String syntax = "ringcmd simulate <params-file> <output-file> [options]\"";
        hf.printHelp(syntax, buildSimulateCliOptions());
    }

    private static Options buildSimulateCliOptions() {
        Options opts = new Options();

        opts.addOption(
                Option
                        .builder()
                        .longOpt("model")
                        .hasArg()
                        .argName("model")
                        .desc("The model-free model used to generate 'true' relaxation values from the parameter set file. Must be one of: 1, 1f, 1s, 2s, 2sf")
                        .build()
        );

        opts.addOption(
                Option
                        .builder()
                        .longOpt("fields")
                        .hasArg()
                        .argName("fields")
                        .desc("Comma-separated list of 1H spectrometer frequencies, in MHz (e.g. 500,700,900)")
                        .build()
        );

        opts.addOption(
                Option
                        .builder()
                        .longOpt("noise")
                        .hasArg()
                        .argName("noise")
                        .desc("Relative noise level applied to simulated signal intensities, as a fraction of the reference intensity (default: 0.02)")
                        .build()
        );

        opts.addOption(
                Option
                        .builder()
                        .longOpt("seed")
                        .hasArg()
                        .argName("seed")
                        .desc("Seed for the random number generator used to corrupt simulated measurements with noise (default: 1)")
                        .build()
        );

        opts.addOption(
                Option
                        .builder()
                        .longOpt("deuterium")
                        .argName("d")
                        .desc("Use deuterium models")
                        .build()
        );

        return opts;
    }

    private static double[] parseFields(String fieldsString) throws ParseException {
        String[] parts = fieldsString.split(",\\s*");
        double[] fields = new double[parts.length];
        try {
            for (int i = 0; i < parts.length; i++) {
                fields[i] = Double.parseDouble(parts[i]);
                if (fields[i] <= 0.0) {
                    throw new NumberFormatException("Must be greater than 0.0");
                }
            }
        } catch (NumberFormatException e) {
            throw new ParseException("Error parsing --fields: " + e.getMessage());
        }
        return fields;
    }

    private static Options buildCliOptions() {
        Options opts = new Options();

        opts.addOption(
                Option
                        .builder()
                        .longOpt("method")
                        .hasArg()
                        .argName("method")
                        .desc("Must be one of: conventional, bagging, regularization (default: conventional)")
                        .build()
        );

        opts.addOption(
                Option
                        .builder()
                        .longOpt("tauM")
                        .hasArg()
                        .argName("tauM")
                        .desc("The initial guess of the global correlation time. Must be a positive double (default: estimated based on the R1 and R2 values in the dataset)")
                        .build()
        );

        opts.addOption(
                Option
                        .builder()
                        .longOpt("fitTauM")
                        .desc("Treat tauM as a parameter to be optimized.")
                        .build()
        );

        opts.addOption(
                Option
                        .builder()
                        .longOpt("tauMFraction")
                        .hasArg()
                        .argName("tauMFraction")
                        .desc("If --fitTauM is provided, specifies the bounding of tauM in the optimization. Must be between 0.0 and 1.0 (default: 0.25)")
                        .build()
        );

        opts.addOption(
                Option
                        .builder()
                        .longOpt("r2Limit")
                        .hasArg()
                        .argName("r2Limit")
                        .desc("If --fitTauM is provided, if a residue has R2 values which are below this threshold, tauM will not be fit for this residue. Must be a positive number (default 0.0)")
                        .build()
        );

        opts.addOption(
                Option
                        .builder()
                        .longOpt("bootstrapMode")
                        .hasArg()
                        .argName("bootstrapMode")
                        .desc("Must be one of: parametric, nonparametric, bayesian (default: parametric)")
                        .build()
        );

        opts.addOption(
                Option
                        .builder()
                        .longOpt("nReplicates")
                        .hasArg()
                        .argName("nReplicates")
                        .desc("Number of bootstrap replicates to run. Larger values will lead to slower performnace but more robust statistics. If --bootstrapMode is nonparametric, the following are the largest permitted values: 2 fields: 27, 3 fields: 343, 4 fields: 6859 (default: 25).")
                        .build()
        );

        opts.addOption(
                Option
                        .builder()
                        .longOpt("models")
                        .hasArg()
                        .argName("models")
                        .desc("Only valid when --method is conventional or bagging. Specifies the models to fit to the data (default: 1,1f,1s,2s,2sf).")
                        .build()
        );

        opts.addOption(
                Option
                        .builder()
                        .longOpt("lambdaScale")
                        .hasArg()
                        .argName("lambdaScale")
                        .desc("Only valid with --method regularization. Specifies the scale of the lambda regularization (default: 1.0)")
                        .build()
        );


        opts.addOption(
                Option
                        .builder()
                        .longOpt("useMedian")
                        .desc("If --method is bagging or reglarization, --useMedian will result in the reported parameters being the median of the bootstrap results. If not provided, the mean will be used.")
                        .build()
        );

        opts.addOption(
                Option
                        .builder()
                        .longOpt("j0Mode")
                        .hasArg()
                        .argName("j0Mode")
                        .desc("Only applicable to amide data. Controls how J(0) is treated: independent (one J(0) per field) or averaged_jackknife (a single J(0) averaged across fields with jackknife uncertainty). Must be one of: independent, averaged_jackknife (default: averaged_jackknife)")
                        .build()
        );

        return opts;
    }

    private static FitSpec getFitSpec(CommandLine cmd, Dataset dataset) throws ParseException {
        FitSpec.Builder<?> builder;
        if (cmd.hasOption("method")) {
            String methodString = cmd.getOptionValue("method");
            switch (methodString) {
                case "conventional" -> {
                    builder = new ConventionalFitSpec.Builder()
                            .modelNames(getModels(cmd));
                }
                case "bagging" -> {
                    builder = new BaggingFitSpec.Builder()
                            .useMedian(cmd.hasOption("useMedian"));
                }
                case "regularization" -> {
                    double lambdaScale = (cmd.hasOption("lambdaScale")) ? parsePositiveDouble(cmd, "lambdaScale", true) : 0.0;
                    builder = new RegularizationFitSpec.Builder()
                            .useMedian(cmd.hasOption("useMedian"))
                            .lambdaScale(lambdaScale);
                }
                default -> throw new ParseException(
                        String.format(
                                "Invalid --method specified: %s.%nMust be one of conventional, bagging, regularization",
                                methodString
                        )
                );
            }
        } else {
            builder = new ConventionalFitSpec.Builder().modelNames(getModels(cmd));
        }

        if (dataset.getClass() == AmideDataset.class) {
            builder.moietyType(FitSpec.MoietyType.AMIDE);
            String j0ModeString = cmd.hasOption("j0Mode") ? cmd.getOptionValue("j0Mode") : "averaged_jackknife";
            R1R2NOEMolDataValues.J0Mode j0Mode = switch (j0ModeString) {
                case "independent" -> R1R2NOEMolDataValues.J0Mode.INDEPENDENT;
                case "averaged_jackknife" -> R1R2NOEMolDataValues.J0Mode.AVERAGED_JACKKNIFE;
                default -> throw new ParseException(
                        String.format(
                                "Invalid --j0Mode specified: %s.%nMust be one of independent, averaged_jackknife",
                                j0ModeString
                        )
                );
            };
            builder.j0Mode(j0Mode);
        } else if (dataset.getClass() == DeuteriumDataset.class) {
            builder.moietyType(FitSpec.MoietyType.DEUTERATED_METHYL);
            if (cmd.hasOption("j0Mode")) {
                throw new ParseException("--j0Mode is only applicable to amide data");
            }
        } else {
            throw new AssertionError("Unexpected Dataset class.");
        }

        if (cmd.hasOption("tauM")) {
            builder.tauM(parsePositiveDouble(cmd, "tauM", false));
        } else {
            builder.tauMUnset();
        }

        if (cmd.hasOption("fitTauM")) {
            builder.fitTauM(true);
            if (cmd.hasOption("tauMFraction")) {
                builder.tauMFraction(parsePositiveDouble(cmd, "tauMFraction", true));
            }
            if (cmd.hasOption("r2Limit")) {
                builder.r2Limit(parsePositiveDouble(cmd, "r2Limit", true));
            }
        } else {
            builder.fitTauM(false);
        }

        if (cmd.hasOption("bootstrapMode")) {
            String bootstrapModeString = cmd.getOptionValue("bootstrapMode");
            FitSpec.BootstrapMode bootstrapMode = switch (bootstrapModeString) {
                case "parametric" -> FitSpec.BootstrapMode.PARAMETRIC;
                case "nonparametric" -> FitSpec.BootstrapMode.NONPARAMETRIC;
                case "bayesian" -> FitSpec.BootstrapMode.BAYESIAN;
                default -> throw new ParseException(
                        String.format(
                                "Invalid --bootstrapMode specified: %s.%nMust be one of parametric, nonparametric, bayesian",
                                bootstrapModeString
                        )
                );
            };
            builder.bootstrapMode(bootstrapMode);
        }
        if (cmd.hasOption("nReplicates")) {
            builder.nReplicates(parsePositiveInt(cmd, "nReplicates"));
        }
        builder.fixedSeed(true);

        return builder.build();
    }

    private static double parsePositiveDouble(CommandLine cmd, String name, boolean canBeZero) {
        try {
            double x = Double.parseDouble(cmd.getOptionValue(name));
            boolean valid = (canBeZero) ? x >= 0.0 : x > 0.0;
            if (!valid) {
                String message = String.format("Must be greater than%s 0.0", canBeZero ? " or equal to" : "");
                throw new NumberFormatException(message);
            }
            return x;
        } catch (NumberFormatException e) {
            throw new NumberFormatException(String.format("Error parsing --%s: %s", name, e.getMessage()));
        }
    }

    private static int parsePositiveInt(CommandLine cmd, String name) {
        try {
            int x = Integer.parseInt(cmd.getOptionValue(name));
            if (x <= 0) {
                String message = String.format("Must be greater than 0");
                throw new NumberFormatException(message);
            }
            return x;
        } catch (NumberFormatException e) {
            throw new NumberFormatException(String.format("Error parsing --%s: %s", name, e.getMessage()));
        }
    }

    private static List<String> getModels(CommandLine cmd) throws ParseException {
        List<String> models = new ArrayList<>();
        if (cmd.hasOption("models")) {
            String[] modelStrings = cmd.getOptionValue("models").split(",\s*");
            for (String modelString : modelStrings) {
                if (!ALLOWED_MODELS.contains(modelString)) {
                    throw new ParseException("Invalid model specified: %s.%nValid models are: 1, 1f, 1s, 2s, 2sf");
                }
                models.add(modelString);
            }
        } else {
            for (String modelString : ALLOWED_MODELS) {
                models.add(modelString);
            }
        }
        return models;
    }
}
