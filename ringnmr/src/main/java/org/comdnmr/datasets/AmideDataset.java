package org.comdnmr.datasets;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.ToDoubleFunction;

import org.comdnmr.data.DataIO;
import org.comdnmr.modelfree.FitR1R2NOEModel;
import org.comdnmr.modelfree.FitSpec;
import org.comdnmr.modelfree.ModelFitResult;
import org.comdnmr.modelfree.R1R2NOEDataValue;
import org.comdnmr.modelfree.MolDataValues;
import org.comdnmr.modelfree.R1R2NOEMolDataValues;
import org.comdnmr.modelfree.R1R2NOEStructureValues;
import org.comdnmr.modelfree.RelaxDataValue;
import org.comdnmr.modelfree.SpectralDensityCalculator;

public class AmideDataset extends Dataset {

    private Optional<Map<String, double[][]>> spectralDensityData = Optional.empty();

    private AmideDataset(R1R2NOEStructureValues data) { super(data); }

    public static AmideDataset fromFile(File file) throws IOException {
        DataIO.loadRelaxationTextFile(file);
        FitR1R2NOEModel fitModel = new FitR1R2NOEModel();
        R1R2NOEStructureValues data = fitModel.getData(false);
        return new AmideDataset(data);
    }

    protected String getMoietyType() { return "NH"; }

    @Override
    protected String getAtomDescriptor(Map.Entry<String, MolDataValues<? extends RelaxDataValue>> residue) {
        return "N";
    }

    private Map<String, double[][]> getSpectralDensityData() {
        if (spectralDensityData.isEmpty()) {
            Map<String, double[][]> map = new TreeMap<>();
            int nFields = getNFields();
            int nOmegas = 3 * nFields;
            for (var residue : this) {
                double[][] arr = new double[3][nOmegas];
                R1R2NOEMolDataValues molData = (R1R2NOEMolDataValues) residue.getValue();
                List<R1R2NOEDataValue> dataValues = molData.getData();
                dataValues.sort(Comparator.comparingDouble(RelaxDataValue::getB0));
                double[][] jData = SpectralDensityCalculator.calcJR1R2NOE(dataValues);
                for (int typeIndex = 0; typeIndex < 3; typeIndex++) {
                    for (int fieldIndex = 0; fieldIndex < nFields; fieldIndex++) {
                        // ω = 0
                        arr[typeIndex][fieldIndex] = jData[typeIndex][3 * fieldIndex];
                        // ωN
                        arr[typeIndex][nFields + fieldIndex] = jData[typeIndex][3 * fieldIndex + 2];
                        // 0.87ωH
                        arr[typeIndex][2 * nFields + fieldIndex] = jData[typeIndex][3 * fieldIndex + 1];
                    }
                }
                map.put(residue.getKey(), arr);
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

    private Map<String, double[]> buildR1R2NOEMap(ToDoubleFunction<R1R2NOEDataValue> extractor) {
        Map<String, double[]> map = new TreeMap<>();
        int n = getNFields();
        for (var residue : this) {
            double[] arr = new double[n];
            R1R2NOEMolDataValues molData = (R1R2NOEMolDataValues) residue.getValue();
            List<R1R2NOEDataValue> dataValues = molData.getData();
            dataValues.sort(Comparator.comparingDouble(RelaxDataValue::getB0));
            int index = 0;
            for (R1R2NOEDataValue value : dataValues) {
                arr[index++] = extractor.applyAsDouble(value);
            }
            map.put(residue.getKey(), arr);
        }
        return map;
    }

    public Map<String, double[]> getNOE() {
        return buildR1R2NOEMap(R1R2NOEDataValue::getNOE);
    }

    public Map<String, double[]> getNOEerr() {
        return buildR1R2NOEMap(R1R2NOEDataValue::getNOEerr);
    }

    protected Map<String, String> relaxationRateToml() {
        Map<String, double[]> fieldMap = getB0();
        Map<String, double[]> r1Map = getR1();
        Map<String, double[]> r1ErrMap = getR1err();
        Map<String, double[]> r2Map = getR2();
        Map<String, double[]> r2ErrMap = getR2err();
        Map<String, double[]> noeMap = getNOE();
        Map<String, double[]> noeErrMap = getNOEerr();

        Map<String, String> map = new TreeMap<>();
        for (var residue : this) {
            String key = residue.getKey();
            StringBuilder builder = new StringBuilder("[relaxation_rates]\n");
            builder.append(String.format("field = %s%n", Arrays.toString(fieldMap.get(key))));
            builder.append(String.format("R1 = %s%n", Arrays.toString(r1Map.get(key))));
            builder.append(String.format("R1_err = %s%n", Arrays.toString(r1ErrMap.get(key))));
            builder.append(String.format("R2 = %s%n", Arrays.toString(r2Map.get(key))));
            builder.append(String.format("R2_err = %s%n", Arrays.toString(r2ErrMap.get(key))));
            builder.append(String.format("NOE = %s%n", Arrays.toString(noeMap.get(key))));
            builder.append(String.format("NOE_err = %s", Arrays.toString(noeErrMap.get(key))));
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
        FitR1R2NOEModel fitModel = new FitR1R2NOEModel();
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
