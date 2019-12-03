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
public class CoMDOptions {

    final boolean usePrefs;

    public CoMDOptions(boolean usePrefs) {
        this.usePrefs = usePrefs;
    }

    /**
     * @return the refField
     */
    public Double getRefField() {
        return usePrefs ? CoMDPreferences.getRefField() : CoMDDefaults.getRefField();
    }

    /**
     * @return the cpmgMaxFreq
     */
    public Double getCpmgMaxFreq() {
        return usePrefs ? CoMDPreferences.getCPMGMaxFreq() : CoMDDefaults.getCpmgMaxFreq();
    }

    /**
     * @return the rexRatio
     */
    public Double getRexRatio() {
        return usePrefs ? CoMDPreferences.getRexRatio() : CoMDDefaults.getRexRatio();
    }

    public Double getDeltaABDiff() {
        return usePrefs ? CoMDPreferences.getDeltaABDiff() : CoMDDefaults.getDeltaABDiff();

    }

    /**
     * @return the startRadius
     */
    public Double getStartRadius() {
        return usePrefs ? CoMDPreferences.getStartRadius() : CoMDDefaults.getStartRadius();
    }

    /**
     * @return the finalRadius
     */
    public Double getFinalRadius() {
        return usePrefs ? CoMDPreferences.getFinalRadius() : CoMDDefaults.getFinalRadius();
    }

    /**
     * @return the tolerance
     */
    public Double getTolerance() {
        return usePrefs ? CoMDPreferences.getTolerance() : CoMDDefaults.getTolerance();
    }

    /**
     * @return the weightFit
     */
    public Boolean getWeightFit() {
        return usePrefs ? CoMDPreferences.getWeightFit() : CoMDDefaults.getWeightFit();
    }

    /**
     * @return the absValueFit
     */
    public Boolean getAbsValueFit() {
        return usePrefs ? CoMDPreferences.getAbsValueFit() : CoMDDefaults.getAbsValueFit();
    }

    /**
     * @return the nonParametricBootstrap
     */
    public Boolean getNonParametricBootstrap() {
        return usePrefs ? CoMDPreferences.getNonParametricBootstrap() : CoMDDefaults.getNonParametricBootstrap();
    }

    /**
     * @return the neuralNetworkGuess
     */
    public Boolean getNeuralNetworkGuess() {
        return usePrefs ? CoMDPreferences.getNeuralNetworkGuess() : CoMDDefaults.getNeuralNetworkGuess();
    }

    /**
     * @return the sampleSize
     */
    public Integer getSampleSize() {
        return usePrefs ? CoMDPreferences.getSampleSize() : CoMDDefaults.getSampleSize();
    }

    /**
     * @return the optimizer
     */
    public String getOptimizer() {
        return usePrefs ? CoMDPreferences.getOptimizer() : CoMDDefaults.getOptimizer();
    }

    /**
     * @return the bootStrapOptimizer
     */
    public String getBootStrapOptimizer() {
        return usePrefs ? CoMDPreferences.getBootStrapOptimizer() : CoMDDefaults.getBootStrapOptimizer();
    }

}
