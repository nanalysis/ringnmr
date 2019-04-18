/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.fit.gui;

import java.util.List;
import javafx.scene.control.ChoiceBox;
import javafx.scene.layout.VBox;
import org.comdnmr.fit.calc.ExperimentData;
import org.comdnmr.fit.calc.ExperimentData.Nuclei;
import org.comdnmr.fit.calc.ParValueInterface;

/**
 *
 * @author Bruce Johnson
 */
public abstract class EquationControls {

    ChoiceBox<ExperimentData.Nuclei> nucleiSelector = new ChoiceBox<>();
    ChoiceBox<String> equationSelector = new ChoiceBox<>();
    ChoiceBox<String> stateSelector = new ChoiceBox<>();

    public ChoiceBox<ExperimentData.Nuclei> getNucChoiceBox() {
        return nucleiSelector;
    }

    abstract String getEquation();

    abstract List<String> getParNames();

    abstract public double[] getExtras();

    public VBox init() {
        VBox vBox = new VBox();
        nucleiSelector.getItems().clear();
        for (ExperimentData.Nuclei nuc : ExperimentData.Nuclei.values()) {
            nucleiSelector.getItems().add(nuc);
        }
        nucleiSelector.getSelectionModel().select(ExperimentData.Nuclei.N15);
        return vBox;
    }

    public Nuclei getNucleus() {
        return nucleiSelector.getSelectionModel().getSelectedItem();
    }
    
    public void setNucleus(String nucName) {
        nucleiSelector.getSelectionModel().select(ExperimentData.Nuclei.get(nucName));
    }

    abstract public VBox makeControls(PyController controller);

    abstract void updateSliders(List<ParValueInterface> parValues, String equationName);

    abstract void updateStates(List<int[]> allStates);

    abstract public void simSliderAction(String label);

    abstract public double[] sliderGuess(String equationName, int[][] map);

    public void updateEquations(ChoiceBox<String> equationChoice, List<String> equationNames) {
        equationChoice.getItems().clear();
        equationChoice.getItems().add("+");
        equationChoice.getItems().addAll(equationNames);
        equationChoice.setValue(equationNames.get(0));
    }

}
