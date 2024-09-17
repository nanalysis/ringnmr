package org.comdnmr.util.traindata;

import org.comdnmr.util.DataUtil;
import org.comdnmr.util.traindata.samplers.*;
import org.comdnmr.util.traindata.dependants.*;

import org.comdnmr.eqnfit.CESTEquation;
import org.comdnmr.eqnfit.CPMGEquation;
import org.comdnmr.eqnfit.EquationType;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.text.StringSubstitutor;

import java.io.File;
import java.io.IOException;

import java.util.*;

public class DataGenerator {

    static final ObjectMapper objectMapper = new ObjectMapper();
    final JsonNode config;
    final int nDatasets;

    Iterator<String[]> iterator;
    public String currentExperimentName;
    public String currentEnumName;

    String rootDirectoryTemplate;
    String dataDirectoryTemplate;
    String dataPathTemplate;

    public DataGenerator(String configFile) throws IOException, ReflectiveOperationException {
        config = objectMapper.readTree(new File(configFile));

        nDatasets = config.get("n_datasets").asInt();
        rootDirectoryTemplate = config.get("root_dir").asText();
        dataDirectoryTemplate = config.get("data_dir").asText();
        dataPathTemplate = config.get("ring_output_path").asText();

        ArrayList<String[]> dataIter = new ArrayList<>();
        Iterator<Map.Entry<String, JsonNode>> dataTypes = config.get("data_types").fields();
        while (dataTypes.hasNext()) {
            Map.Entry<String, JsonNode> info = dataTypes.next();
            String experimentName = info.getKey();
            Iterator<Map.Entry<String, JsonNode>> enums = info.getValue().get("enums").fields();
            while (enums.hasNext()) {
                String enumName = enums.next().getKey();
                String[] dataInfo = new String[2];
                dataInfo[0] = experimentName;
                dataInfo[1] = enumName;
                dataIter.add(dataInfo);
            }
        }

        this.iterator = dataIter.iterator();
    }

    public boolean hasAnotherDataset() {
        return iterator.hasNext();
    }

    public void nextDataset() {
        if (iterator.hasNext()) {
            String[] info = iterator.next();
            currentExperimentName = info[0];
            currentEnumName = info[1];
        } else {
            currentExperimentName = null;
            currentEnumName = null;
        }
    }

    public HashMap<String, Object> generateDataset() {
        if (currentExperimentName == null) {
            return null;
        }

        JsonNode experimentInfo = config.get("data_types").get(currentExperimentName);
        JsonNode parameterInfo = experimentInfo.get("enums").get(currentEnumName);

        MultiSampler<Double> xValuesSampler = fetchSampler(experimentInfo.get("x_values"));

        JsonNode variableInfo = experimentInfo.get("variable");
        MultiSampler<Double> variableSampler = fetchSampler(experimentInfo.get("variable"));

        List<DependantGenerator> dependantGenerators = null;
        JsonNode dependantInfo = variableInfo.get("dependants");
        if (dependantInfo != null) {
            dependantGenerators = getDependantGenerators(dependantInfo);
        }

        List<JsonNode> constantInfoList = objectMapper.convertValue(
            experimentInfo.get("constants"),
            new TypeReference<ArrayList<JsonNode>>() { }
        );
        ArrayList<MultiSampler<Double>> constantSamplers = new ArrayList<MultiSampler<Double>>();
        for (JsonNode constantInfo : constantInfoList) {
            constantSamplers.add(fetchSampler(constantInfo));
        }

        List<RandomDoubleMultiSampler> parameterSamplers = getParameterSamplers(parameterInfo);

        // Stores a list of 1D double lists
        DataList<ArrayList<Double>> variableList = new DataList<>(variableSampler.name);

        // Stores a list of 2D double lists
        DataList<ArrayList<ArrayList<Double>>> profileList = new DataList<>("profile");

        ArrayList<DataList<Double>> parameterLists = new ArrayList<>();
        for (RandomDoubleMultiSampler parameterSampler : parameterSamplers) {
            DataList<Double> parameterList = new DataList<>(parameterSampler.name);
            parameterLists.add(parameterList);
        }

        ArrayList<DataList<Double>> constantLists = new ArrayList<>();
        for (MultiSampler<Double> constantSampler : constantSamplers) {
            DataList<Double> constantList = new DataList<>(constantSampler.name);
            constantLists.add(constantList);
        }

        EquationType cls = fetchEquationClass();

        int xSize = 2 + constantSamplers.size();
        int nDeps = 0;
        if (dependantGenerators != null) {
            nDeps = dependantGenerators.size();
            xSize += nDeps;
        }

        double[] x = new double[xSize];

        int parSize = parameterSamplers.size();
        int[] map = new int[parSize];
        for (int i = 0; i < parSize; i++) {
            map[i] = i;
        }
        double[] par = new double[parSize];

        // >>> interpolation information >>>
        ArrayList<Double> xValuesInterpolation = fetchSampler(experimentInfo.get("x_values_interpolation")).sample();
        DataList<ArrayList<ArrayList<Double>>> profileInterpolationList = new DataList<>("profile_interpolation");
        // <<< interpolation information <<<

        for (int n = 0; n < this.nDatasets; n++) {
            // >>> Make `par` >>>
            int parIdx = 0;
            for (RandomDoubleMultiSampler parameterSampler : parameterSamplers) {
                par[parIdx] = getDoubleFromSampler(parameterSampler);
                parameterLists.get(parIdx).add(par[parIdx]);
                parIdx += 1;
            }

            // >>> Add constants to `x` >>>
            int constantIdx = 2 + nDeps;
            for (MultiSampler<Double> constantSampler : constantSamplers) {
                x[constantIdx] = getDoubleFromSampler(constantSampler);
                constantLists.get(constantIdx - 2 - nDeps).add(x[constantIdx]);
                constantIdx += 1;
            }
            // <<< Add constants to `x` <<<

            ArrayList<Double> xValues = xValuesSampler.sample();
            // >>> Add constants to `x` >>>
            ArrayList<Double> variables = variableSampler.sample();

            variableList.add(variables);

            ArrayList<ArrayList<Double>> profile2D = new ArrayList<>();
            for (int varIdx = 0; varIdx < variables.size(); varIdx++) {
                ArrayList<Double> profile1D = new ArrayList<>();
                x[1 + nDeps] = variables.get(varIdx);
                for (int depIdx = 0; depIdx < nDeps; depIdx++) {
                    x[1 + depIdx] = dependantGenerators.get(depIdx).fetch(x[1 + nDeps]);
                }
                for (int xValueIdx = 0; xValueIdx < xValues.size(); xValueIdx++) {
                    x[0] = xValues.get(xValueIdx);

                    profile1D.add(cls.calculate(par, map, x, 0));
                }
                profile2D.add(profile1D);
            }

            profileInterpolationList.add(
                getInterpolatedProfile(
                    profile2D,
                    xValues,
                    xValuesInterpolation
                )
            );

            profileList.add(profile2D);

            if ((n + 1) % 1000 == 0 && (n + 1) != nDatasets) {
                System.out.println(
                    String.format(
                        "[%s - %s]: %d/%d datasets created",
                        currentExperimentName,
                        currentEnumName,
                        n + 1,
                        nDatasets
                    )
                );
            }
        }

        HashMap<String, Object> data = makeData(
            parameterLists,
            constantLists,
            variableList,
            profileList,
            profileInterpolationList
        );

        System.out.println(
            String.format(
                "[%s - %s]: Finished making %d datasets",
                currentExperimentName,
                currentEnumName,
                nDatasets
            )
        );

        return data;
    }

    public String getRootDirectory() {
        Map<String, String> sub = new HashMap<>();
        sub.put("name", String.format("%s-%s", currentExperimentName, currentEnumName));
        return StringSubstitutor.replace(rootDirectoryTemplate, sub, "{", "}");
    }

    public String getDataDirectory() {
        Map<String, String> sub = new HashMap<>();
        sub.put("root_dir", getRootDirectory());
        return StringSubstitutor.replace(dataDirectoryTemplate, sub, "{", "}");
    }

    public String getDataPath() {
        Map<String, String> sub = new HashMap<>();
        sub.put("data_dir", getDataDirectory());
        return StringSubstitutor.replace(dataPathTemplate, sub, "{", "}");
    }

    ArrayList<ArrayList<Double>> getInterpolatedProfile(
        ArrayList<ArrayList<Double>> profile,
        ArrayList<Double> xValuesProfile,
        ArrayList<Double> xValuesInterplation
    ) {
        ArrayList<ArrayList<Double>> profileInterpolated = new ArrayList<>();
        ArrayList<Double> profInterp;
        ArrayList<Double> prof;
        double[] xValues = ArrayUtils.toPrimitive(xValuesProfile.toArray(new Double[xValuesProfile.size()]));
        for (int i = 0; i < profile.size(); i++) {
            prof = profile.get(i);
            double[] yValues = ArrayUtils.toPrimitive(prof.toArray(new Double[prof.size()]));
            SplineInterpolator splineInterpolator = new SplineInterpolator();
            PolynomialSplineFunction spline = splineInterpolator.interpolate(xValues, yValues);

            profInterp = new ArrayList<Double>();

            for (Double xInterp : xValuesInterplation) {
                profInterp.add(spline.value(xInterp));
            }

            profileInterpolated.add(profInterp);
        }


        return profileInterpolated;
    }

    HashMap<String, Object> makeData(
        ArrayList<DataList<Double>> parameterLists,
        ArrayList<DataList<Double>> constantLists,
        DataList<ArrayList<Double>> variableList,
        DataList<ArrayList<ArrayList<Double>>> profileList,
        DataList<ArrayList<ArrayList<Double>>> profileInterpolationList
    ) {
        HashMap<String, Object> data = new HashMap<>();

        ArrayList<HashMap<String, Object>> parameterMapList = new ArrayList<>();
        for (DataList<Double> parameterList : parameterLists) {
            parameterMapList.add(parameterList.asHashMap());
        }

        ArrayList<HashMap<String, Object>> constantMapList = new ArrayList<>();
        for (DataList<Double> constantList : constantLists) {
            constantMapList.add(constantList.asHashMap());
        }

        data.put("parameters", parameterMapList);
        data.put("constants", constantMapList);
        data.put("variable", variableList.asHashMap());
        data.put("profile", profileList.asHashMap());
        data.put("profile_interpolation", profileInterpolationList.asHashMap());

        return data;
    }


    MultiSampler<Double> fetchSampler(JsonNode info) {
        // Bit hacky, need to initialize `sampler` to something to avoid the
        // compiler complaining.
        MultiSampler<Double> sampler = new ExplicitSampler<Double>("DUMMY", 0.0);
        String name = info.get("name").asText();
        JsonNode samplerInfo = info.get("sampling");

        switch (samplerInfo.get("type").asText()) {
            case "RandomDouble":
                sampler = new RandomDoubleMultiSampler(
                    name,
                    samplerInfo.get("minValue").asDouble(),
                    samplerInfo.get("maxValue").asDouble(),
                    samplerInfo.get("minSamples").asInt(),
                    samplerInfo.get("maxSamples").asInt()
                );
                break;

            case "Choose":
                sampler = new ChooseSampler<Double>(
                    name,
                    objectMapper.convertValue(samplerInfo.get("options"), ArrayList.class),
                    samplerInfo.get("minSamples").asInt(),
                    samplerInfo.get("maxSamples").asInt()
                );
                break;

            case "RandomRange":
                sampler = new RangeSampler(
                    name,
                    samplerInfo.get("minStart").asDouble(),
                    samplerInfo.get("maxStart").asDouble(),
                    samplerInfo.get("minStop").asDouble(),
                    samplerInfo.get("maxStop").asDouble(),
                    objectMapper.convertValue(samplerInfo.get("stepOptions"), ArrayList.class)
                );
                break;

            case "RandomLinspace":
                sampler = new LinspaceSampler(
                    name,
                    samplerInfo.get("minFirst").asDouble(),
                    samplerInfo.get("maxFirst").asDouble(),
                    samplerInfo.get("minLast").asDouble(),
                    samplerInfo.get("maxLast").asDouble(),
                    samplerInfo.get("minSamples").asInt(),
                    samplerInfo.get("maxSamples").asInt()
                );
                break;

            case "RandomSegmentedLinspace":
                sampler = new SegmentedLinspaceSampler(
                    name,
                    objectMapper.convertValue(samplerInfo.get("minFirsts"), ArrayList.class),
                    objectMapper.convertValue(samplerInfo.get("maxFirsts"), ArrayList.class),
                    objectMapper.convertValue(samplerInfo.get("minLasts"), ArrayList.class),
                    objectMapper.convertValue(samplerInfo.get("maxLasts"), ArrayList.class),
                    objectMapper.convertValue(samplerInfo.get("minSampless"), ArrayList.class),
                    objectMapper.convertValue(samplerInfo.get("maxSampless"), ArrayList.class)
                );
                break;

            case "RandomLogspace":
                sampler = new LogspaceSampler(
                    name,
                    samplerInfo.get("minFirst").asDouble(),
                    samplerInfo.get("maxFirst").asDouble(),
                    samplerInfo.get("scale").asDouble(),
                    samplerInfo.get("minSamples").asInt(),
                    samplerInfo.get("maxSamples").asInt()
                );
                break;

            case "ExplicitValue":
                sampler = new ExplicitSampler<Double>(
                    name,
                    samplerInfo.get("value").asDouble()
                );
                break;

            case "ExplicitList":
                sampler = new ExplicitSamplerList<Double>(
                    name,
                    objectMapper.convertValue(samplerInfo.get("value"), ArrayList.class)
                );
                break;

            case "ExplicitLinspace":
                sampler = new ExplicitSamplerLinspace(
                    name,
                    samplerInfo.get("first").asDouble(),
                    samplerInfo.get("last").asDouble(),
                    samplerInfo.get("nSamples").asInt()
                );
                break;

            default: throw new RuntimeException("Unknown sampler specification");
        }

        return sampler;
    }

    List<RandomDoubleMultiSampler> getParameterSamplers(JsonNode parametersInfo) {
        Iterator<JsonNode> parameterIterator = parametersInfo.elements();
        List<RandomDoubleMultiSampler> parameterSamplers = new ArrayList<RandomDoubleMultiSampler>();

        while (parameterIterator.hasNext()) {
            JsonNode parameterInfo = parameterIterator.next();
            String parameterName = parameterInfo.get(0).asText();
            double parameterMin = parameterInfo.get(1).asDouble();
            double parameterMax = parameterInfo.get(2).asDouble();
            parameterSamplers.add(new RandomDoubleMultiSampler(parameterName, parameterMin, parameterMax, 1, 1));
        }

        return parameterSamplers;
    }

    EquationType fetchEquationClass() {
        // Assign `cls` to any old class initially, so the compiler does not complain with:
        // "variable might not have been initialized"
        EquationType cls = CPMGEquation.CPMGMQ;

        switch (currentExperimentName) {
            case "CPMGEquation":
                switch (currentEnumName) {
                    case "CPMGMQ": cls = CPMGEquation.CPMGMQ; break;
                    case "CPMGFAST": cls = CPMGEquation.CPMGFAST; break;
                    case "CPMGSLOW": cls = CPMGEquation.CPMGSLOW; break;
                    default: throw new RuntimeException();
                }
                break;

            case "CESTEquation":
                switch (currentEnumName) {
                    case "EXACT0": cls = CESTEquation.EXACT0; break;
                    default: throw new RuntimeException();
                }
                break;

            default: throw new RuntimeException();
        }

        return cls;
    }

    double getDoubleFromSampler(MultiSampler<Double> sampler) {
        ArrayList<Double> valueList = sampler.sample();
        double value = valueList.get(0);
        return value;
    }

    List<DependantGenerator> getDependantGenerators(JsonNode info) {
        Iterator<JsonNode> dependantIterator = info.elements();
        List<DependantGenerator> dependantGenerators = new ArrayList<DependantGenerator>();

        while (dependantIterator.hasNext()) {
            JsonNode dependantInfo = dependantIterator.next();
            dependantGenerators.add(fetchDependantGenerator(dependantInfo));
        }

        return dependantGenerators;
    }

    DependantGenerator fetchDependantGenerator(JsonNode info) {
        // Bit hacky, need to initialize `generator` to something to avoid the
        // compiler complaining.
        DependantGenerator generator = new Multiplier(1.0);

        String name = info.get("name").asText();
        String generatorInfo = info.get("relation").asText();
        String[] arr = generatorInfo.split(":");
        String generatorType = arr[0];
        String generatorParams = arr[1];

        switch (generatorType) {
            case "Mul":
                double mul = Double.parseDouble(generatorParams);
                generator = new Multiplier(mul);
                break;

            default: throw new RuntimeException("Unknown dependant specification");
        }

        return generator;
    }
}
