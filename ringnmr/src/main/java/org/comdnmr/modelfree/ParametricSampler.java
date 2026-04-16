package org.comdnmr.modelfree;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.ZigguratSampler;

/**
 * Parametric bootstrap sampler that generates synthetic observations by adding
 * Gaussian noise to the original measured values.
 *
 * <p>On each call to {@link #sample()}, every relaxation observable (R1, R2,
 * and either NOE for {@link R1R2NOEDataValue} or rQ/rAP for
 * {@link DeuteriumDataValue}) is independently perturbed by a zero-mean Gaussian
 * with standard deviation equal to the corresponding reported experimental error.
 * The original values are preserved internally and restored by
 * {@link #getOriginalData()}.
 *
 * <p>Normal deviates are drawn using the Ziggurat algorithm via Apache Commons RNG.
 *
 * <p><strong>Note:</strong> the underlying {@link MolDataValues} object is mutated
 * in-place. After calling {@link #sample()}, it holds perturbed values until the
 * next call to {@link #sample()} or {@link #getOriginalData()}.
 */
public class ParametricSampler extends BootstrapSampler {

    private final NormalizedGaussianSampler gaussian;

    private Class<? extends RelaxDataValue> dataType;
    private Map<String, List<Double>> originalDataMap;

    /**
     * Constructs a {@code ParametricSampler} and snapshots the original observed
     * values from {@code data} for later restoration.
     *
     * @param data the relaxation data to be resampled; must contain at least one entry
     * @throws IllegalStateException if {@code data} has no data values
     * @throws AssertionError if the data values are an unrecognised subtype
     */
    public ParametricSampler(MolDataValues data) {
        super(data);
        gaussian = ZigguratSampler.NormalizedGaussian.of(rng);
        setOriginalData();
    }

    private void setOriginalData() {
        try {
            dataType = data.dataValues.get(0).getClass();
        } catch (IndexOutOfBoundsException e) {
            throw new IllegalStateException("`MolDataValues` object has no dataValues!");
        }

        originalDataMap = new HashMap<>();
        originalDataMap.put("R1", new ArrayList<>());
        originalDataMap.put("R2", new ArrayList<>());

        if (dataType.equals(R1R2NOEDataValue.class)) {
            originalDataMap.put("NOE", new ArrayList<>());
            for (RelaxDataValue dataValue : data.dataValues) {
                originalDataMap.get("R1").add(dataValue.R1);
                originalDataMap.get("R2").add(dataValue.R2);
                originalDataMap.get("NOE").add(((R1R2NOEDataValue) dataValue).NOE);
            }
        } else if (dataType.equals(DeuteriumDataValue.class)) {
            originalDataMap.put("rQ", new ArrayList<>());
            originalDataMap.put("rAP", new ArrayList<>());
            for (RelaxDataValue dataValue : data.dataValues) {
                originalDataMap.get("R1").add(dataValue.R1);
                originalDataMap.get("R2").add(dataValue.R2);
                DeuteriumDataValue value = (DeuteriumDataValue) dataValue;
                originalDataMap.get("rQ").add(value.rQ);
                originalDataMap.get("rAP").add(value.rAP);
            }
        } else {
            throw new AssertionError("Unexpected subclass of `RelaxDataValue` detected!");
        }
    }

    /**
     * {@inheritDoc}
     *
     * <p>Perturbs each observable (R1, R2, and NOE or rQ/rAP depending on data type)
     * by independent Gaussian noise scaled by the per-field experimental error.
     * Clears the cached J-values so they are recomputed on next access.
     */
    public MolDataValues sample() {
        List<RelaxDataValue> dataValues = data.getData();
        for (int i = 0; i < getNFields(); i++) {
            RelaxDataValue dataValue = dataValues.get(i);
            dataValue.R1 = originalDataMap.get("R1").get(i) + dataValue.R1err * gaussian.sample();
            dataValue.R2 = originalDataMap.get("R2").get(i) + dataValue.R2err * gaussian.sample();

            if (dataType.equals(R1R2NOEDataValue.class)) {
                R1R2NOEDataValue value = (R1R2NOEDataValue) dataValue;
                value.NOE = originalDataMap.get("NOE").get(i) + value.NOEerr * gaussian.sample();
            } else if (dataType.equals(DeuteriumDataValue.class)) {
                DeuteriumDataValue value = (DeuteriumDataValue) dataValue;
                value.rQ  = originalDataMap.get("rQ").get(i)  + value.rQError  * gaussian.sample();
                value.rAP = originalDataMap.get("rAP").get(i) + value.rAPError * gaussian.sample();
            } else {
                throw new AssertionError("Unexpected subclass of `RelaxDataValue` detected!");
            }
        }
        // Force J-values to be recomputed on next access
        data.jValues = null;
        return data;
    }

    /**
     * {@inheritDoc}
     *
     * <p>Restores every observable to the value snapshotted at construction and
     * clears the cached J-values.
     */
    public MolDataValues getOriginalData() {
        for (int i = 0; i < getNFields(); i++) {
            RelaxDataValue dataValue = data.dataValues.get(i);
            dataValue.R1 = originalDataMap.get("R1").get(i);
            dataValue.R2 = originalDataMap.get("R2").get(i);

            if (dataType.equals(R1R2NOEDataValue.class)) {
                R1R2NOEDataValue value = (R1R2NOEDataValue) dataValue;
                value.NOE = originalDataMap.get("NOE").get(i);
            } else if (dataType.equals(DeuteriumDataValue.class)) {
                DeuteriumDataValue value = (DeuteriumDataValue) dataValue;
                value.rQ  = originalDataMap.get("rQ").get(i);
                value.rAP = originalDataMap.get("rAP").get(i);
            } else {
                throw new AssertionError("Unexpected subclass of `RelaxDataValue` detected!");
            }
        }
        // Force J-values to be recomputed on next access
        data.jValues = null;
        return data;
    }
}
