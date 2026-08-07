package org.comdnmr.datasets;

import java.io.File;
import java.io.IOException;
import java.util.*;

import org.comdnmr.data.DataIO;
import org.comdnmr.modelfree.DeuteriumDataValue;
import org.comdnmr.modelfree.DeuteriumMolDataValues;
import org.comdnmr.modelfree.DeuteriumStructureValues;
import org.comdnmr.modelfree.FitDeuteriumModel;
import org.comdnmr.modelfree.FitSpec;
import org.comdnmr.modelfree.ModelFitResult;
import org.comdnmr.modelfree.MolDataValues;
import org.comdnmr.modelfree.RelaxDataValue;
import org.comdnmr.modelfree.SpectralDensityCalculator;

public class DeuteriumDataset extends Dataset {

    private Optional<Map<String, double[][]>> spectralDensityData = Optional.empty();

    private DeuteriumDataset(DeuteriumStructureValues data) { super(data); }

    public static DeuteriumDataset fromFile(File file) throws IOException {
        DataIO.loadRelaxationTextFile(file);
        DeuteriumStructureValues data = FitDeuteriumModel.getData(false);
        return new DeuteriumDataset(data);
    }

    protected String getMoietyType() { return "CDH2"; }

    // Carbon position in amino acid, where 1 represents the α-carbon
    static int getCarbonPosition(Map.Entry<String, MolDataValues<? extends RelaxDataValue>> residue) {
        char letter = residue.getValue().getAtom().getName().charAt(1);
        return switch (letter) {
            // 'G' -> γ
            case 'G' -> 3;
            // Others 'B' -> β, 'D' -> δ, 'E' -> ε can be determined based on arithmetic
            default -> letter - 'A' + 1;
        };
    }

    // If there are two methyl groups in the amino acid (i.e. leucine, with two
    // δ-methyl groups), this gives each chemically equivalent methyl a unique
    // index.
    static Optional<Integer> getCarbonIndex(Map.Entry<String, MolDataValues<? extends RelaxDataValue>> residue) {
        String atomDescriptor = residue.getValue().getAtom().getName();
        if (atomDescriptor.length() != 3) return Optional.empty();
        return Optional.of(atomDescriptor.charAt(2) - '0');
    }

    static String buildAtomDescriptor(Map.Entry<String, MolDataValues<? extends RelaxDataValue>> residue) {
        StringBuilder descriptor = new StringBuilder("C" + (char) ('A' + getCarbonPosition(residue) - 1));
        getCarbonIndex(residue).ifPresent(i -> descriptor.append(i));
        return descriptor.toString();
    }

    @Override
    protected String getAtomDescriptor(Map.Entry<String, MolDataValues<? extends RelaxDataValue>> residue) {
        return buildAtomDescriptor(residue);
    }

    public Map<String, double[]> getRQ() {
        return buildRelaxMap(v -> ((DeuteriumDataValue) v).getObservables()[2]);
    }

    public Map<String, double[]> getRQerr() {
        return buildRelaxMap(v -> ((DeuteriumDataValue) v).getObservableErrors()[2]);
    }

    public Map<String, double[]> getRAP() {
        return buildRelaxMap(v -> ((DeuteriumDataValue) v).getObservables()[3]);
    }

    public Map<String, double[]> getRAPErr() {
        return buildRelaxMap(v -> ((DeuteriumDataValue) v).getObservableErrors()[3]);
    }

    private Map<String, double[][]> getSpectralDensityData() {
        if (spectralDensityData.isEmpty()) {
            Map<String, double[][]> map = new TreeMap<>();
            for (var residue : this) {
                DeuteriumMolDataValues molData = (DeuteriumMolDataValues) residue.getValue();
                List<DeuteriumDataValue> dataValues = molData.getData();
                dataValues.sort(Comparator.comparingDouble(RelaxDataValue::getB0));
                double[][] jData = SpectralDensityCalculator.calcJDeuterium(dataValues);
                // rows 0/1/2 are omega/J/Jerr at each unique spectral density frequency,
                // for both the single-field (independentMapping) and multi-field (jointMapping) cases
                map.put(residue.getKey(), new double[][]{jData[0], jData[1], jData[2]});
            }
            spectralDensityData = Optional.of(map);
        }
        return spectralDensityData.get();
    }

    private Map<String, double[]> getSpectralDensityVector(int index) {
        Map<String, double[][]> specDenData = getSpectralDensityData();
        Map<String, double[]> map = new TreeMap<>();
        for (var residue : this) {
            String key = residue.getKey();
            map.put(key, specDenData.get(key)[index]);
        }
        return map;
    }

    public Map<String, double[]> getOmega() {
        return getSpectralDensityVector(0);
    }

    public Map<String, double[]> getSpectralDensities() {
        return getSpectralDensityVector(1);
    }

    public Map<String, double[]> getSpectralDensityErrors() {
        return getSpectralDensityVector(2);
    }

    protected Map<String, String> residueToml() {
        Map<String, String> map = super.residueToml();
        for (var residue : this) {
            String key = residue.getKey();
            StringBuilder toAddBuilder = new StringBuilder();
            toAddBuilder.append(String.format("%ncarbon_position = %d", getCarbonPosition(residue)));
            Optional<Integer> carbonIndex = getCarbonIndex(residue);
            if (carbonIndex.isPresent()) {
                toAddBuilder.append(String.format("%ncarbon_index = %d", carbonIndex.get()));
            }
            String toAdd = toAddBuilder.toString();
            map.put(key, map.get(key) + toAdd);
        }
        return map;
    }

    protected Map<String, String> relaxationRateToml() {
        Map<String, double[]> fieldMap = getB0();
        Map<String, double[]> r1Map = getR1();
        Map<String, double[]> r1ErrMap = getR1err();
        Map<String, double[]> r2Map = getR2();
        Map<String, double[]> r2ErrMap = getR2err();
        Map<String, double[]> rqMap = getRQ();
        Map<String, double[]> rqErrMap = getRQerr();
        Map<String, double[]> rapMap = getRAP();
        Map<String, double[]> rapErrMap = getRAPErr();

        Map<String, String> map = new TreeMap<>();
        for (var residue : this) {
            String key = residue.getKey();
            StringBuilder builder = new StringBuilder("[relaxation_rates]\n");
            builder.append(String.format("field = %s%n", Arrays.toString(fieldMap.get(key))));
            builder.append(String.format("R1 = %s%n", Arrays.toString(r1Map.get(key))));
            builder.append(String.format("R1_err = %s%n", Arrays.toString(r1ErrMap.get(key))));
            builder.append(String.format("R2 = %s%n", Arrays.toString(r2Map.get(key))));
            builder.append(String.format("R2_err = %s%n", Arrays.toString(r2ErrMap.get(key))));
            builder.append(String.format("RQ = %s%n", Arrays.toString(rqMap.get(key))));
            builder.append(String.format("RQ_err = %s%n", Arrays.toString(rqErrMap.get(key))));
            builder.append(String.format("RAP = %s%n", Arrays.toString(rapMap.get(key))));
            builder.append(String.format("RAP_err = %s", Arrays.toString(rapErrMap.get(key))));
            map.put(key, builder.toString());
        }
        return map;
    }

    protected Map<String, String> spectralDensityToml() {
        Map<String, double[]> omegaMap = getOmega();
        Map<String, double[]> jMap = getSpectralDensities();
        Map<String, double[]> jErrMap = getSpectralDensityErrors();

        Map<String, String> map = new TreeMap<>();
        for (var residue : this) {
            String key = residue.getKey();
            StringBuilder builder = new StringBuilder("[spectral_densities]\n");
            builder.append(String.format("omega = %s%n", Arrays.toString(omegaMap.get(key))));
            builder.append(String.format("J = %s%n", Arrays.toString(jMap.get(key))));
            builder.append(String.format("J_err = %s", Arrays.toString(jErrMap.get(key))));
            map.put(key, builder.toString());
        }
        return map;
    }

    @Override
    public Map<String, ModelFitResult> fit(FitSpec fitSpec) {
        FitDeuteriumModel fitModel = new FitDeuteriumModel();
        fitModel.setData(data);
        fitModel.setFitSpec(fitSpec);
        return fitModel.testIsoModel();
    }

    @Override
    Map<String, String> getToml() {
        Map<String, String> residueTomlMap = residueToml();
        Map<String, String> relaxationRateTomlMap = relaxationRateToml();
        Map<String, String> spectralDensityTomlMap = spectralDensityToml();
        Map<String, String> map = new TreeMap<>();
        for (var residue : this) {
            String key = residue.getKey();
            String toml = String.join(
                "\n\n",
                residueTomlMap.get(key),
                relaxationRateTomlMap.get(key),
                spectralDensityTomlMap.get(key)
            );
            map.put(key, toml);
        }
        return map;
    }
}
