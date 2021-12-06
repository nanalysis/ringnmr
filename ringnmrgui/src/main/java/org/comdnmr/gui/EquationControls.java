/*
 * CoMD/NMR Software : A Program for Analyzing NMR Dynamics Data
 * Copyright (C) 2018-2019 Bruce A Johnson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
package org.comdnmr.gui;

import java.util.List;
import javafx.scene.control.ChoiceBox;
import javafx.scene.layout.VBox;
import org.comdnmr.eqnfit.ParValueInterface;
import org.nmrfx.datasets.Nuclei;

/**
 *
 * @author Bruce Johnson
 */
public abstract class EquationControls {

    ChoiceBox<Nuclei> nucleiSelector = new ChoiceBox<>();
    ChoiceBox<String> equationSelector = new ChoiceBox<>();
    ChoiceBox<String> stateSelector = new ChoiceBox<>();

    public ChoiceBox<Nuclei> getNucChoiceBox() {
        return nucleiSelector;
    }

    abstract String getEquation();

    abstract List<String> getParNames();

    abstract public double[] getExtras();

    public VBox init() {
        VBox vBox = new VBox();
        nucleiSelector.getItems().clear();
        for (Nuclei nuc : Nuclei.values()) {
            nucleiSelector.getItems().add(nuc);
        }
        nucleiSelector.getSelectionModel().select(Nuclei.N15);
        return vBox;
    }

    public Nuclei getNucleus() {
        return nucleiSelector.getSelectionModel().getSelectedItem();
    }
    
    public void setNucleus(String nucName) {
        nucleiSelector.getSelectionModel().select(Nuclei.findNuclei(nucName));
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
