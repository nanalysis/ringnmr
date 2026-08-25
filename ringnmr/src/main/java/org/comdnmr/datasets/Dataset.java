package org.comdnmr.datasets;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.function.ToDoubleFunction;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.comdnmr.modelfree.FitSpec;
import org.comdnmr.modelfree.ModelFitResult;
import org.comdnmr.modelfree.MolDataValues;
import org.comdnmr.modelfree.R1R2NOEMolDataValues;
import org.comdnmr.modelfree.RelaxDataValue;
import org.comdnmr.modelfree.StructureValues;


public abstract class Dataset implements Iterable<Map.Entry<String, MolDataValues<? extends RelaxDataValue>>> {

    private static final Pattern KEY_PATTERN = Pattern.compile("\\d+:(\\d+)\\.[A-Za-z0-9]+");
    protected final StructureValues data;
    private Optional<List<Map.Entry<String, MolDataValues<? extends RelaxDataValue>>>> sortedData = Optional.empty();

    Dataset(StructureValues data) { this.data = data; }

    public StructureValues getData() {
        return data;
    }

    public int getNFields() {
        return data.entrySet().iterator().next().getValue().getData().size();
    }

    public Map<String, Double> getMaxR2s() {
        Map<String, Double> result = new TreeMap<>();
        for (var residue : getSortedData()) {
            double maxR2 = residue.getValue().getData()
                .stream()
                .mapToDouble(RelaxDataValue::getR2)
                .max()
                .getAsDouble();
            result.put(residue.getKey(), maxR2);
        }
        return result;
    }

    List<Map.Entry<String, MolDataValues<? extends RelaxDataValue>>> getSortedData() {
        if (sortedData.isEmpty()) {
            sortedData = Optional.of(new ArrayList<>(data.entrySet()));
            sortedData.get().sort(Comparator.comparingInt(Dataset::getResidueNumber));
        }
        return sortedData.get();
    }

    @Override
    public Iterator<Map.Entry<String, MolDataValues<? extends RelaxDataValue>>> iterator() {
        return getSortedData().iterator();
    }

    protected Map<String, double[]> buildRelaxMap(ToDoubleFunction<RelaxDataValue> extractor) {
        Map<String, double[]> map = new TreeMap<>();
        int n = getNFields();
        for (var residue : this) {
            double[] arr = new double[n];
            List<? extends RelaxDataValue> dataValues = residue.getValue().getData();
            dataValues.sort(Comparator.comparingDouble(RelaxDataValue::getB0));
            int index = 0;
            for (RelaxDataValue value : dataValues) {
                arr[index++] = extractor.applyAsDouble(value);
            }
            map.put(residue.getKey(), arr);
        }
        return map;
    }

    public Map<String, double[]> getB0() {
        return buildRelaxMap(RelaxDataValue::getB0);
    }

    public Map<String, double[]> getR1() {
        return buildRelaxMap(RelaxDataValue::getR1);
    }

    public Map<String, double[]> getR1err() {
        return buildRelaxMap(RelaxDataValue::getR1err);
    }

    public Map<String, double[]> getR2() {
        return buildRelaxMap(RelaxDataValue::getR2);
    }

    public Map<String, double[]> getR2err() {
        return buildRelaxMap(RelaxDataValue::getR2err);
    }

    protected Map<String, String> residueToml() {
        Map<String, String> map = new TreeMap<>();
        for (var residue : this) {
            String moietyType = getMoietyType();
            int residueNumber = getResidueNumber(residue);
            String residueName = getResidueName(residue);
            String toml = String.format(
                "[spin_system]%nmoiety_type = \"%s\"%nresidue_number = %d%nresidue_name = \"%s\"",
                moietyType, residueNumber, residueName);
            map.put(residue.getKey(), toml);
        }
        return map;
    }

    abstract Map<String, String> getToml();

    public void saveToToml(Path rootDirectory) throws IOException {
        Map<String, String> tomlMap = getToml();
        for (var residue : this) {
            String directoryName = getDirectoryName(residue);
            Path parent = rootDirectory.resolve(directoryName);
            Path path = parent.resolve("dataset.toml");
            Files.createDirectories(parent);
            Files.writeString(path, tomlMap.get(residue.getKey()) + "\n");
        }
    }

    static int getResidueNumber(Map.Entry<String, MolDataValues<? extends RelaxDataValue>> residue) {
        String key = residue.getKey();
        Matcher matcher = KEY_PATTERN.matcher(key);
        if (matcher.matches()) {
            return Integer.parseInt(matcher.group(1));
        }
        throw new IllegalArgumentException("Key does not match pattern: " + key);
    }

    static String getResidueName(Map.Entry<String, MolDataValues<? extends RelaxDataValue>> residue) {
        return residue.getValue().getAtom().getResidueName();
    }

    static char getAtom(Map.Entry<String, MolDataValues<? extends RelaxDataValue>> residue) {
        return residue.getValue().getAtom().getName().charAt(0);
    }

    protected abstract String getMoietyType();

    protected abstract String getAtomDescriptor(Map.Entry<String, MolDataValues<? extends RelaxDataValue>> residue);

    String getDirectoryName(Map.Entry<String, MolDataValues<? extends RelaxDataValue>> residue) {
        return String.format(getDirectoryFormat() + "-%s",
            getResidueNumber(residue), getResidueName(residue), getAtomDescriptor(residue));
    }

    String getDirectoryFormat() {
        var sortedData = getSortedData();
        int lastResidueValue = getResidueNumber(sortedData.get(sortedData.size() - 1));
        int nDigits = (int) Math.floor(Math.log10(lastResidueValue)) + 1;
        return String.format("%%0%dd-%%s", nDigits);
    }

    public abstract Map<String, ModelFitResult> fit(FitSpec fitSpec);

    public Map<FitSpec, Map<String, ModelFitResult>> fit(List<FitSpec> fitSpecs) {
        Map<FitSpec, Map<String, ModelFitResult>> result = new HashMap<>();
        for (FitSpec fitSpec : fitSpecs) {
            result.put(fitSpec, fit(fitSpec));
        }
        return result;
    }

    public void fit(List<FitSpec> fitSpecs, Path rootDirectory) throws IOException {
        var result = fit(fitSpecs);
        writeFitResults(result, rootDirectory);
    }

    void writeFitResults(Map<FitSpec, Map<String, ModelFitResult>> result, Path rootDirectory)
    throws IOException
    {
        for (var res : result.entrySet()) {
            FitSpec fitSpec = res.getKey();
            String fileName = fitSpec.filenameWithFingerprint("fit", 16, "toml");
            for (var residue : this) {
                ModelFitResult modelFitResult = res.getValue().get(residue.getKey());
                Path parent = rootDirectory.resolve(getDirectoryName(residue)).resolve("fits");
                Files.createDirectories(parent);
                Path path = parent.resolve(fileName);
                String fitSpecToml = fitSpec.toToml();
                String fitResultToml = modelFitResult.toToml(true);
                String toml = String.join("\n\n", fitSpecToml, fitResultToml);
                Files.writeString(path, toml);
            }
        }
    }

    public static Dataset fromFile(File file) throws IOException {
        try (var reader = new BufferedReader(new FileReader(file))) {
            String header = reader.readLine();
            if (header != null) {
                if (header.contains("RQ") && header.contains("RQ_err")
                        && header.contains("RAP") && header.contains("RAP_err")) {
                    return DeuteriumDataset.fromFile(file);
                }
                if (header.contains("NOE") && header.contains("NOE_err")) {
                    return AmideDataset.fromFile(file);
                }
            }
        }
        throw new IllegalArgumentException("Unrecognised dataset format in file: " + file.getName());
    }
}
