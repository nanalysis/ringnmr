/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.gui;

import java.util.ArrayList;
import java.util.List;
import org.comdnmr.data.ExperimentResult;
import org.comdnmr.data.ExperimentSet;
import org.comdnmr.data.ValueSet;

/**
 *
 * @author brucejohnson
 */
public class ChartInfo {

    String mapName;
    String state;
    String equationName;
    String[] currentResidues;
    ExperimentResult experimentalResult;
    ValueSet valueSet;
    List<int[]> currentStates = new ArrayList<>();

    public boolean hasExperiments() {
        return valueSet != null;
    }

    public ValueSet getExperiments() {
        return valueSet;
    }

    public boolean hasResult() {
        return experimentalResult != null;
    }

    public boolean hasResidue() {
        return experimentalResult != null;
    }

    public boolean hasResidues() {
        return currentResidues != null;
    }
    
    public int getResNum() {
        return experimentalResult.getResNum();
    }

    public String getResNumStr() {
        return String.valueOf(experimentalResult.getResNum());
    }

    public String[] getResidues() {
        return currentResidues;
    }

    public void setResidues(String[] residues) {
        this.currentResidues = residues;
    }

    public ExperimentResult getResult() {
        return experimentalResult;
    }

    public String getEquationName() {
        return equationName;
    }
    
    public void setStates(List<int[]> states) {
        currentStates.clear();
        currentStates.addAll(states);        
    }

    public void clear() {
        currentResidues = null;
        experimentalResult = null;
        valueSet = null;
        currentStates = new ArrayList<>();
    }

}
