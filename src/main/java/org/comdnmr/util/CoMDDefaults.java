/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.util;

/**
 *
 * @author brucejohnson
 */
public class CoMDDefaults {

    static final Double REF_FIELD = 500.0;
    static final Double CPMG_MAX_FREQ = 3000.0;
    static final Double REX_RATIO = 3.0;
    static final Double START_RADIUS = 20.0;
    static final Double FINAL_RADIUS = -5.0;
    static final Double TOLERANCE = -5.0;
    static final Double DELTA_AB_DIFF = 0.1;
    static final Boolean WEIGHT_FIT = true;
    static final Boolean ABS_VALUE_FIT = false;// false
    static final Boolean NON_PARAMETRIC_BOOTSTRAP = true;
    static final Boolean NEURAL_NETWORK_GUESS = true;
    static final Integer SAMPLE_SIZE = 50; // 
    static final String OPTIMIZER = "CMA-ES";
    static final String BOOTSTRAP_OPTIMIZER = "CMA-ES";

    /**
     * @return the refField
     */
    public static Double getRefField() {
        return REF_FIELD;
    }

    /**
     * @return the cpmgMaxFreq
     */
    public static Double getCpmgMaxFreq() {
        return CPMG_MAX_FREQ;
    }

    /**
     * @return the rexRatio
     */
    public static Double getRexRatio() {
        return REX_RATIO;
    }

      /**
     * @return the AB Diff limit
     */
    public static Double getDeltaABDiff() {
        return DELTA_AB_DIFF;
    }


    /**
     * @return the startRadius
     */
    public static Double getStartRadius() {
        return START_RADIUS;
    }

    /**
     * @return the finalRadius
     */
    public static Double getFinalRadius() {
        return FINAL_RADIUS;
    }

    /**
     * @return the tolerance
     */
    public static Double getTolerance() {
        return TOLERANCE;
    }

    /**
     * @return the weightFit
     */
    public static Boolean getWeightFit() {
        return WEIGHT_FIT;
    }

    /**
     * @return the absValueFit
     */
    public static Boolean getAbsValueFit() {
        return ABS_VALUE_FIT;
    }

    /**
     * @return the nonParametricBootstrap
     */
    public static Boolean getNonParametricBootstrap() {
        return NON_PARAMETRIC_BOOTSTRAP;
    }

    /**
     * @return the neuralNetworkGuess
     */
    public static Boolean getNeuralNetworkGuess() {
        return NEURAL_NETWORK_GUESS;
    }

    /**
     * @return the sampleSize
     */
    public static Integer getSampleSize() {
        return SAMPLE_SIZE;
    }

    /**
     * @return the optimizer
     */
    public static String getOptimizer() {
        return OPTIMIZER;
    }

    /**
     * @return the bootStrapOptimizer
     */
    public static String getBootStrapOptimizer() {
        return BOOTSTRAP_OPTIMIZER;
    }

}
